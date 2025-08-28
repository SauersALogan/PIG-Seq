#!/usr/bin/env python3
"""
sample_expression_matrix.py - Merge normalized samples into expression matrix
==========================================

Function to merge multiple normalized and annotated sample files into a 
comparative expression matrix, handling the challenge of matching genes 
across different assemblies with different coordinates and identifiers.

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import pandas as pd
import re
import shutil
import sys
import numpy as np
import glob

# =============================================================================
# Import utility functions
# =============================================================================
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_pairing import extract_file_identifiers, pair_files_by_sample

# =============================================================================
# Actual functions to test
# =============================================================================
def create_gene_identifier(row):
    """
    Create unique identifier for cross-sample gene matching
    """
    if row['product'] == 'hypothetical protein':
        return f"{row['binning']}|hypothetical|{row['Length']}"
    else:
        return f"{row['binning']}|{row['product'].strip()}"

def create_expression_matrix(normalized_files, output_file="expression_matrix.txt", discard_unique_hypotheticals=True):
    """
    Create expression matrix from multiple normalized sample files using exact size matching
    """
    all_sample_data = []
    sample_names = []
    for file_path in normalized_files:
        df = pd.read_csv(file_path, sep='\t')
        sample_name = os.path.basename(file_path).split('_')[0]
        sample_names.append(sample_name)
        df['gene_id'] = df.apply(create_gene_identifier, axis=1)
        sample_data = df[['gene_id', 'binning', 'product', 'Length', 'TPM']].copy()
        sample_data['sample'] = sample_name
        all_sample_data.append(sample_data)
        print(f"Loaded {len(df)} genes from {sample_name}")
    combined_data = pd.concat(all_sample_data, ignore_index=True)
    if discard_unique_hypotheticals:
        hyp_genes = combined_data[combined_data['product'] == 'hypothetical protein']
        gene_counts = hyp_genes['gene_id'].value_counts()
        universal_hyp_genes = gene_counts[gene_counts == len(normalized_files)].index
        annotated_genes = combined_data[combined_data['product'] != 'hypothetical protein']
        universal_hyp_data = combined_data[combined_data['gene_id'].isin(universal_hyp_genes)]
        filtered_data = pd.concat([annotated_genes, universal_hyp_data], ignore_index=True)
        print(f"Kept {len(annotated_genes['gene_id'].unique())} unique annotated genes")
        print(f"Kept {len(universal_hyp_genes)} universally matched hypothetical proteins")
        print(f"Discarded {len(gene_counts) - len(universal_hyp_genes)} unique hypothetical proteins")
    else:
        filtered_data = combined_data
    matrix = filtered_data.pivot_table(
        index=['binning', 'product', 'gene_id'],
        columns='sample',
        values='TPM',
        fill_value=0,
        aggfunc='sum'
    ).reset_index()
    sample_columns = [col for col in matrix.columns if col in sample_names]
    matrix['total'] = matrix[sample_columns].sum(axis=1)
    matrix = matrix.sort_values(['binning', 'total'], ascending=[True, False])
    matrix = matrix.reset_index(drop=True)
    matrix.to_csv(output_file, sep='\t', index=False)
    print(f"\nExpression matrix created: {output_file}")
    print(f"Matrix shape: {matrix.shape}")
    print(f"Unique gene functions: {len(matrix)}")
    print(f"Bins represented: {matrix['binning'].nunique()}")
    return matrix

def run_sample_merging(input_files, output_directory="./"):
    """
    Run the sample merging pipeline on multiple input formats
    """
    normalized_files = []
    for pattern in input_files:
        if os.path.isdir(pattern):
            normalized_files.extend(glob.glob(os.path.join(pattern, "*_normalized.txt")))
            normalized_files.extend(glob.glob(os.path.join(pattern, "*normalized*.txt")))
        else:
            expanded = glob.glob(pattern)
            if expanded:
                normalized_files.extend(expanded)
            else:
                print(f"Warning: No files found for pattern: {pattern}")
    normalized_files = sorted(list(set(normalized_files))) # Remove duplicates just in case glob acts funny
    if not normalized_files:
        print("ERROR: No normalized sample files found")
        return None
    print(f"Found {len(normalized_files)} normalized sample files:")
    for file in normalized_files:
        print(f" - {file}")
    os.makedirs(output_directory, exist_ok=True)
    output_file = os.path.join(output_directory, "expression_matrix.txt")
    print(f"\nCreating expression matrix...")
    matrix = create_expression_matrix(
        normalized_files,
        output_file,
        discard_unique_hypotheticals=discard_unique_hypotheticals)
    if matrix is None:
        print("ERROR: Failed to create expression matrix")
        return None
    print(f"\n=== EXPRESSION MATRIX SUMMARY ===")
    print(f"Output file: {output_file}")
    print(f"Matrix dimensions: {matrix.shape}")
    sample_cols = [col for col in matrix.columns if col not in ['binning', 'product', 'gene_id', 'total']]
    print(f"Samples included: {len(sample_cols)}")
    print(f"Sample names: {', '.join(sample_cols)}")
    total_genes = len(matrix)
    annotated_genes = len(matrix[matrix['product'] != 'hypothetical protein'])
    hypothetical_genes = total_genes - annotated_genes
    print(f"Total gene functions: {total_genes}")
    print(f"Annotated functions: {annotated_genes}")
    print(f"Hypothetical proteins: {hypothetical_genes}")
    bins = matrix['binning'].unique()
    print(f"Bins represented: {len(bins)}")
    total_tpm = matrix['total'].sum()
    print(f"Total TPM across all genes: {total_tpm:.2f}")
    return output_file

# =============================================================================
# Command line interfacing
# =============================================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Merge normalized sample files into expression matrix")
    parser.add_argument("--input_files", nargs="+", required=True,
        help="Normalized sample files or directories containing them")
    parser.add_argument("--output", default="./",
        help="Output directory (default: current directory)")
    parser.add_argument("--keep_unique_hypotheticals", action="store_true",
        help="Keep hypothetical proteins found in only some samples (default: discard)")

    args = parser.parse_args()

    result = run_sample_merging(args.input_files,
        args.output,
        discard_unique_hypotheticals=not args.keep_unique_hypotheticals)
    if result:
        print(f"\nExpression matrix successfully created: {result}")
    else:
        print("\nFailed to create expression matrix")
        exit(1)

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def temp_test_dir():
    """Create a temp directory for test files"""
    temp_dir = tempfile.mkdtemp(prefix='Sample_merging_test')
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def normalized_sample_files(temp_test_dir):
    """Create multiple normalized sample files for testing"""

    # Sample 1 data
    sample1_content = """Geneid	Chr	Start	End	Strand	Length	Sample1_mapped.sam	binning	CPM	RPK	TPM	product
PROKKA_00001	Contig00001	149	1174	+	1026	201	Bin1.fa		5.07	195.90	4.79	tRNA-specific adenosine deaminase
PROKKA_00002	Contig00001	1327	2295	+	969	190	Bin1.fa		4.79	196.08	4.80	hypothetical protein
PROKKA_00003	Contig00001	2653	3135	-	483	75	Bin2.fa		1.89	155.28	3.80	hypothetical protein
PROKKA_00004	Contig00002	100	800	+	701	120	unbinned	3.02	171.18	4.19	DNA polymerase
"""

    # Sample 2 data - same functions, different PROKKA IDs and coordinates
    sample2_content = """Geneid	Chr	Start	End	Strand	Length	Sample2_mapped.sam	binning	CPM	RPK	TPM	product
PROKKA_00008	Contig00001	200	1300	+	1101	180	Bin1.fa		4.50	163.49	4.20	tRNA-specific adenosine deaminase
PROKKA_00009	Contig00001	1500	2400	+	901	165	Bin2.fa		4.12	183.13	4.70	hypothetical protein
PROKKA_00010	Contig00001	3000	3400	-	401	95	Bin2.fa		2.37	236.91	6.08	hypothetical protein
PROKKA_00011	Contig00002	500	1200	+	701	140	unbinned	3.50	199.72	5.13	DNA polymerase
"""

    # Sample 3 data
    sample3_content = """Geneid	Chr	Start	End	Strand	Length	Sample3_mapped.sam	binning	CPM	RPK	TPM	product
PROKKA_00015	Contig00001	300	1400	+	1101	220	Bin1.fa		5.50	199.82	5.10	tRNA-specific adenosine deaminase
PROKKA_00016	Contig00001	1600	2100	+	501	85	Bin2.fa		2.12	169.66	4.33	hypothetical protein
PROKKA_00017	Contig00002	200	900	+	701	160	unbinned	4.00	228.25	5.83	DNA polymerase
PROKKA_00018	Contig00003	100	2000	-	1901	280	Bin3.fa		7.00	147.29	3.76	Ferrienterobactin receptor
"""

    sample_files = []
    for i, content in enumerate([sample1_content, sample2_content, sample3_content], 1):
        sample_file = os.path.join(temp_test_dir, f'S2{i}_Pos_mapped_binned_normalized.txt')
        with open(sample_file, 'w') as f:
            f.write(content)
        sample_files.append(sample_file)
    return sample_files

@pytest.fixture
def expected_matrix_output(temp_test_dir):
    """Expected expression matrix output"""
    expected_content = """binning	product	gene_id	Sample1	Sample2	Sample3	total
Bin1.fa	tRNA-specific adenosine deaminase	Bin3.fa|tRNA-specific adenosine deaminase	4.79	4.20	5.10	14.09
Bin2.fa	hypothetical protein	Bin2.fa|hypothetical|969	4.80	4.70	0.00	9.50
Bin2.fa	hypothetical protein	Bin2.fa|hypothetical|401	3.80	6.08	0.00	9.88
Bin3.fa	Ferrienterobactin receptor	Bin3.fa|Ferrienterobactin receptor	0.00	0.00	3.76	3.76
unbinned	DNA polymerase	unbinned|DNA polymerase	4.19	5.13	5.83	15.15
"""

    expected_file = os.path.join(temp_test_dir, 'expected_matrix.txt')
    with open(expected_file, 'w') as f:
        f.write(expected_content)
    return expected_file

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_gene_categorization():
    """Test that gene identifiers are created correctly"""
    row_annotated = pd.Series({
        'binning': 'Bin1.fa',
        'product': 'DNA polymerase',
        'Length': 1000
    })
    assert create_gene_identifier(row_annotated) == "Bin1.fa|DNA polymerase"
    print(f"Test categorizer produced the following annotations {row_annotated}")
    row_hypothetical = pd.Series({
        'binning': 'Bin1.fa',
        'product': 'hypothetical protein',
        'Length': 969
    })
    assert create_gene_identifier(row_hypothetical) == "Bin1.fa|hypothetical|969"
    print(f"Test categorizer produced the following hypotheticals {row_hypothetical}")

def test_run_sample_merging(normalized_sample_files):
    """Test the main sample merging function"""
    

@pytest.fixture(autouse=True)
def cleanup_files(request):
    def cleanup():
        import glob
        for file in glob.glob("*.txt"):
            if os.path.exists(file):
                os.unlink(file)
    request.addfinalizer(cleanup)
