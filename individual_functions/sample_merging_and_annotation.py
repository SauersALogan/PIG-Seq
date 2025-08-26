#!/usr/bin/env python3
"""
feature_parsing.py - Use the output parsed PAFs to filter sam files
==========================================

Function to filter feature counts by assembly specific contig
to bin map files generated earlier in this pipeline

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import pandas as pd
import re
import shutil
import sys

# =============================================================================
# Import utility functions
# =============================================================================
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_pairing import extract_file_identifiers, pair_files_by_sample

# =============================================================================
# Actual functions to test
# =============================================================================
def normalize_and_annotate(count_file, gff_file):
    count_file = pd.read_csv(count_file, delimiter='\t', header=True)
    count_col = count_file.columns[-2]
    total_counts = count_file[count_col].sum()
    count_file['CPM'] = (count_file[count_col] / total_counts) * 1_000_000
    count_file['RPK'] = count_file[count_col] / (count_file['Length'] / 1000)
    total_RPK = count_file['RPK'].sum()
    count_file['TPM'] = (count_file['RPK'] / total_RPK) * 1_000_000

    gff = pd.read_csv(gff_file, delimiter='\t', comment='#', header=None)
    cds_df = gff_df[gff_df[2] == 'CDS'].copy()
    cds_df['locus_tag'] = cds_df[8].str.extract(r'locus_tag=([^;]+)')
    cds_df['product'] = cds_df[8].str.extract(r'product=([^;]+)')
    cds_df['product'] = cds_df['product'].fillna('hypothetical protein')
    annotations = dict(zip(cds_df['locus_tag'], cds_df['product']))
    count_file['product'] = count_file['Geneid'].map(annotations)

    count_base = os.path.basename(count_file_path)
    count_name = os.path.splitext(count_base)[0]
    output_name = count_name + "_normalized.txt"
    output_path = output_name

    count_file.to_csv(output_path, sep='\t', index=False)
    return output_path

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def temp_test_dir():
    """Create a temp directory for test files"""
    temp_dir = tempfile.mkdtemp(prefix='Feature_parse_unit_test')
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def single_sample(temp_test_dir):
    """Creating a single sample fixture for testing"""
    content = """Geneid	Chr	Start	End	Strand	Length	TIP-Seq_mapping/S20_Native_mapping_renamed_contigs/S20_Native_mapped.sam	binning
PROKKA_00003	Contig00001	626	835	-	210	65	S29_Nativebin.10.fa
PROKKA_00004	Contig00001	1099	2196	-	1098	179	S29_Nativebin.10.fa
PROKKA_00005	Contig00001	2302	2595	-	294	27	S29_Nativebin.10.fa
PROKKA_00006	Contig00002	200	800	+	601	45	unbinned
PROKKA_00007	Contig00002	1000	1500	-	501	23	unbinned
"""

    sample_file = os.path.join(temp_test_dir, 'Test_sample1.txt')
    with open(sample_file, 'w') as f:
        f.write(content)
    return sample_file

@pytest.fixture
def multiple_sample(temp_test_dir):
    """Create multiple sample files for unit testing."""
    sample1_content = """Geneid	Chr	Start	End	Strand	Length	TIP-Seq_mapping/S20_Native_mapping_renamed_contigs/S20_Native_mapped.sam	binning
PROKKA_00003	Contig00001	626	835	-	210	65	S29_Nativebin.10.fa
PROKKA_00004	Contig00001	1099	2196	-	1098	179	S29_Nativebin.10.fa
PROKKA_00005	Contig00001	2302	2595	-	294	27	S29_Nativebin.10.fa
PROKKA_00006	Contig00002	200	800	+	601	45	unbinned
"""

    sample2_content = """Geneid	Chr	Start	End	Strand	Length	TIP-Seq_mapping/S21_Native_mapping_renamed_contigs/S21_Native_mapped.sam	binning
PROKKA_00008	Contig00001	500	900	+	401	145	S29_Nativebin.10.fa
PROKKA_00009	Contig00001	1200	2000	-	801	203	S29_Nativebin.10.fa
PROKKA_00010	Contig00002	300	700	+	401	67	unbinned
PROKKA_00011	Contig00003	100	500	-	401	23	S45_Nativebin.05.fa
"""

    sample3_content = """Geneid	Chr	Start	End	Strand	Length	TIP-Seq_mapping/S22_Native_mapping_renamed_contigs/S22_Native_mapped.sam	binning
PROKKA_00012	Contig00001	700	1100	+	401	87	S29_Nativebin.10.fa
PROKKA_00013	Contig00001	1500	2300	-	801	156	S29_Nativebin.10.fa
PROKKA_00014	Contig00002	400	900	+	501	12	unbinned
PROKKA_00015	Contig00003	200	600	-	401	34	S45_Nativebin.05.fa
"""

    sample_files = []
    for i, content in enumerate([sample1_content, sample2_content, sample3_content], 1):
        sample_file = os.path.join(temp_test_dir, f'Test_sample{i}.txt')
        with open(sample_file, 'w') as f:
            f.write(content)
        sample_files.append(sample_file)
    return sample_files

@pytest.fixture
def single_gff(temp_test_dir):
    """Create a mock gff file for annotation."""
    gff_content = """##gff-version 3
##sequence-region Contig00001 1 50000
##sequence-region Contig00002 1 30000
Contig00001	prokka	CDS	626	835	.	-	0	ID=PROKKA_00003;locus_tag=PROKKA_00003;product=hypothetical protein;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	1099	2196	.	-	0	ID=PROKKA_00004;locus_tag=PROKKA_00004;product=DNA polymerase;gene=polA;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	2302	2595	.	-	0	ID=PROKKA_00005;locus_tag=PROKKA_00005;product=ribosomal protein L1;gene=rplA;inference=ab initio prediction:Prodigal:2.6
Contig00002	prokka	CDS	200	800	.	+	0	ID=PROKKA_00006;locus_tag=PROKKA_00006;product=tRNA synthetase;gene=trpS;inference=ab initio prediction:Prodigal:2.6
Contig00002	prokka	CDS	1000	1500	.	-	0	ID=PROKKA_00007;locus_tag=PROKKA_00007;product=membrane protein;inference=ab initio prediction:Prodigal:2.6
"""

    gff_file = os.path.join(temp_test_dir, 'Sample1.gff')
    with open(gff_file, 'w') as f:
        f.write(gff_content)
    return gff_file

@pytest.fixture
def multiple_gffs(temp_test_dir):
    """Create a multiple gff files for testing."""
    gff1_content = """##gff-version 3
##sequence-region Contig00001 1 50000
##sequence-region Contig00002 1 30000
Contig00001	prokka	CDS	626	835	.	-	0	ID=PROKKA_00003;locus_tag=PROKKA_00003;product=hypothetical protein;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	1099	2196	.	-	0	ID=PROKKA_00004;locus_tag=PROKKA_00004;product=DNA polymerase;gene=polA;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	2302	2595	.	-	0	ID=PROKKA_00005;locus_tag=PROKKA_00005;product=ribosomal protein L1;gene=rplA;inference=ab initio prediction:Prodigal:2.6
Contig00002	prokka	CDS	200	800	.	+	0	ID=PROKKA_00006;locus_tag=PROKKA_00006;product=tRNA synthetase;gene=trpS;inference=ab initio prediction:Prodigal:2.6
"""

    gff2_content = """##gff-version 3
##sequence-region Contig00001 1 50000
##sequence-region Contig00002 1 30000
##sequence-region Contig00003 1 20000
Contig00001	prokka	CDS	500	900	.	+	0	ID=PROKKA_00008;locus_tag=PROKKA_00008;product=DNA polymerase;gene=polA;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	1200	2000	.	-	0	ID=PROKKA_00009;locus_tag=PROKKA_00009;product=ribosomal protein L1;gene=rplA;inference=ab initio prediction:Prodigal:2.6
Contig00002	prokka	CDS	300	700	.	+	0	ID=PROKKA_00010;locus_tag=PROKKA_00010;product=hypothetical protein;inference=ab initio prediction:Prodigal:2.6
Contig00003	prokka	CDS	100	500	.	-	0	ID=PROKKA_00011;locus_tag=PROKKA_00011;product=ATP synthase;gene=atpA;inference=ab initio prediction:Prodigal:2.6
"""

    gff3_content = """##gff-version 3
##sequence-region Contig00001 1 50000
##sequence-region Contig00002 1 30000
##sequence-region Contig00003 1 20000
Contig00001	prokka	CDS	700	1100	.	+	0	ID=PROKKA_00012;locus_tag=PROKKA_00012;product=ribosomal protein L1;gene=rplA;inference=ab initio prediction:Prodigal:2.6
Contig00001	prokka	CDS	1500	2300	.	-	0	ID=PROKKA_00013;locus_tag=PROKKA_00013;product=DNA polymerase;gene=polA;inference=ab initio prediction:Prodigal:2.6
Contig00002	prokka	CDS	400	900	.	+	0	ID=PROKKA_00014;locus_tag=PROKKA_00014;product=tRNA synthetase;gene=trpS;inference=ab initio prediction:Prodigal:2.6
Contig00003	prokka	CDS	200	600	.	-	0	ID=PROKKA_00015;locus_tag=PROKKA_00015;product=ATP synthase;gene=atpA;inference=ab initio prediction:Prodigal:2.6
"""

    gff_files = []
    for i, content in enumerate([gff1_content, gff2_content, gff3_content], 1):
        gff_file = os.path.join(temp_test_dir, f'Sample{i}.gff')
        with open(gff_file, 'w') as f:
            f.write(content)
        gff_files.append(gff_file)
    return gff_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_assembly(single_sample, single_gff):
    """Test that the parsing function works on a single sample/assembly"""

def test_multiple_assemblies(multiple_samples, multiple_gffs):
    """Test that parsing function works on multiple samples/assemblies"""

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for sam_file in glob.glob("*.txt"):
            if os.path.exists(sam_file):
                os.unlink(sam_file)
    request.addfinalizer(cleanup)
