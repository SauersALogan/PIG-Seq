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
    # TODO: Handle different input formats (files vs directory)
    # TODO: Call create_expression_matrix
    # TODO: Generate summary statistics
    # TODO: Return output path

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
    row_hypothetical = pd.Series({
        'binning': 'Bin1.fa',
        'product': 'hypothetical protein',
        'Length': 969
    })
    assert create_gene_identifier(row_hypothetical) == "Bin1.fa|hypothetical|969"

def test_single_sample_processing(normalized_sample_files):
    """Test processing of a single normalized sample file"""
    # TODO: Test loading and categorizing genes from one file
    pass

def test_expression_matrix_creation(normalized_sample_files, expected_matrix_output):
    """Test creation of expression matrix from multiple samples"""
    # TODO: Test full matrix creation pipeline
    # TODO: Compare with expected output
    # TODO: Verify TPM summation properties
    pass

def test_run_sample_merging(normalized_sample_files):
    """Test the main sample merging function"""
    # TODO: Test the main pipeline function
    # TODO: Verify output file creation
    # TODO: Test summary statistics
    pass

@pytest.fixture(autouse=True)
def cleanup_files(request):
    def cleanup():
        import glob
        for file in glob.glob("*.txt"):
            if os.path.exists(file):
                os.unlink(file)
    request.addfinalizer(cleanup)
