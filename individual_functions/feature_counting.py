#!/usr/bin/env python3
"""
SAM_parsing.py - Use the output parsed PAFs to filter sam files
==========================================

Function to filter and separate SAM files mapped to assemblies
using PAF files created from earlier in this pipeline

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import subprocess

# =============================================================================
# Actual functions to test
# =============================================================================
def build_featureCounts (gtf_files, output_path, sam_files, threads, gene, gene_id):
    """Count SAM file features using GTF annotation"""
    featureCounts = [
        "featureCounts",
        "-T", str(threads),
        "-t", gene, # The prokka classification of the read, usually default to gene or CDS
        "-g", gene_id, # The assigned gene_id, used to group features with similar gene_ids 
        "-a", gtf_files,
        "-o", output_path,
        sam_files
    ]
    return featureCounts

def execute_feature_counting(sam_file, gtf_file, threads = 2, gene = "CDS", gene_id = "gene_id")
    count_tables = []
    sam_base=os.path.basename(sam_file)
    sam_name=os.path.splitext(sam_base)[0]
    output_name=sam_name+".txt"
    output_path=(output_name)
    Counter = build_featureCounts(gtf_file, output_path, sam_file, threads, gene, gene_id)
    try:
        counts = subprocess.run(Counter, check=True)
        count_tables.append(output_path)
    except subprocess.CalledProcessError:
        print("featureCounts failed to process {sam_file}")
    return(count_tables)

def run_counter(sam_files, gtf_files, threads = 2, gene = "CDS", gene_id = "gene_id")
    """Runs the counting function imported above on the actual data"""
    if isinstance(sam_files, str):
        sam_file = sam_files
        gtf_file = gtf_files
        execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
    elif isinstance(sam_files, list):
        if isinstance(gtf_files, str):
        gtf_files = gtf_file
        for sam_file in sam_files:
            execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
        if isinstance(gtf_files, list):
            for gtf_file in gtf_files:
                for sam_file in sam_files:
                    execture_feature_counting(sam_file, gtf_file, threads, gene, gene_id)

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_sam():
    """Create a single test SAM file for unit testing."""
    sam_content = """@SQ     SN:Contig1   LN:341341
@SQ     SN:Contig2   LN:325562
@SQ     SN:Contig3   LN:317386
@SQ     SN:Contig4   LN:301055
@SQ     SN:Contig5   LN:298171
@SQ     SN:Contig6   LN:282893
read1:AA 147     Contig1     3554    60      150M    =       3376 -328     GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACATCCCATGGACAGTATATGTTTGTATTACCAAACAGATAGCCGTCATACCACTCCTTTACAGCGCCGTACTGGTCCGTAAATCCATAATACGCCAGCATCT        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII        NM:i:0  MD:Z:150        MC:Z:150M       AS:i:150     XS:i:21
read2:LL 83      Contig2    72558   60      150M    =       72369-339     CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCAGTTCTCCTGTATATTCATCGGTAAAACGAATCTTTAAATCCTTTTCCACATAGACAAAGCTGACATTGACCTTGCCTGATGCCGTCTCCGAAGTCTC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII        NM:i:1  MD:Z:103G46     MC:Z:150M       AS:i:145     XS:i:20
read3:ZZ 163     Contig3    72369   60      150M    =       72558339      CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACCTTATCCTTTACTGTAACCTGCTGCGGAACCGACGGGAACGTATATTTCTCTGTGGAGGTTATGATCGCATCAAACACGCCCGCGTTCACATTCTGCGC        IIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII9IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIII        NM:i:0  MD:Z:150        MC:Z:150M       AS:i:150     XS:i:0
read4:YY 99      Contig4    77025   60      150M    =       77222326      CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACACAGTTTGAGGAGATGGGGCTATTTAACCGCGGTATCAAAACCAGAATCGGCAGCAATGGAAACTCAGATTATTATGGCATCATACGTGAGTGTTCGA        IIII-IIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9III        NM:i:0  MD:Z:150        MC:Z:129M       AS:i:150     XS:i:74
read5:XX 147     Contig5    77222   60      129M    =       77025-326     GTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTGTTCGTGACGCACAGGCAGTCGCACGGTACTTTGGGCTTGTGTCTAAAGATAAAACAAAGGATTACAGCAATTATGCT     99-I9IIIIIIIIIIIIIIIIIIIIIIII--IIII9IIIIIIIIIIIIIIIIIIIIIIII-IIIII9IIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIII     NM:i:2  MD:Z:0A29G98    MC:Z:150M       AS:i:123        XS:i:0
read6:TT 163     Contig6    17929   60      150M    =       18129350      CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAATAAAATGCGATGATATCGCCGGTCTTGTATTTATTTCCGTTCAGCGCCACTACGATGTCTTCCTCCTGCAGCGACTCCGTCACGGAAGTTCCGCTGATC        III9IIIIIIII--IIIIIIIIIIIII9II9IIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIII9I9I9IIIII-IIII9IIIIIIIIII9II-999I99I9I--9IIIII99-IIII-II        NM:i:2  MD:Z:113A19T16  MC:Z:150M       AS:i:140     XS:i:0
"""

@pytest.fixture
def multiple_sams():
    """Create multiple SAM files for unit testing."""
    sam_1_content ="""@SQ     SN:Contig1   LN:341341
@SQ     SN:Contig2   LN:325562
@SQ     SN:Contig3   LN:317386
read1:AA 147     Contig1     3554    60      150M    =       3376 -328     GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACA>
read2:LL 83      Contig2    72558   60      150M    =       72369-339     CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGAT>
read3:ZZ 163     Contig3    72369   60      150M    =       72558339      CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCAC>
"""

    sam_2_content ="""@SQ     SN:Contig4   LN:301055
@SQ     SN:Contig5   LN:298171
@SQ     SN:Contig6   LN:282893
read4:YY 99      Contig4    77025   60      150M    =       77222326      CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGA>
read5:XX 147     Contig5    77222   60      129M    =       77025-326     GTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGT>
read6:TT 163     Contig6    17929   60      150M    =       18129350      CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAA>
"""

    sam_files = []
    for i, content in enumerate([sam_1_content, sam_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.paf') as tmp:
            tmp.write(content)
            sam_files.append(tmp.name)
    return sam_files

@pytest.fixture
def mock_gtf():
"""Create a mock gtf file for FeatureCounts."""
    gtf_content = gtf_content = """Contig1	Prokka	gene	1000	2500	.	+	.	gene_id "gene_001"; gene_name "hypothetical_protein_001";
Contig1	Prokka	CDS	1200	2300	.	+	0	gene_id "gene_001"; gene_name "hypothetical_protein_001";
Contig1	Prokka	gene	3000	4500	.	-	.	gene_id "gene_002"; gene_name "ABC_transporter";
Contig1	Prokka	CDS	3200	4300	.	-	0	gene_id "gene_002"; gene_name "ABC_transporter";
Contig2	Prokka	gene	50000	52000	.	+	.	gene_id "gene_003"; gene_name "ribosomal_protein";
Contig2	Prokka	CDS	50200	51800	.	+	0	gene_id "gene_003"; gene_name "ribosomal_protein";
Contig2	Prokka	gene	70000	75000	.	+	.	gene_id "gene_004"; gene_name "DNA_polymerase";
Contig2	Prokka	CDS	70200	74800	.	+	0	gene_id "gene_004"; gene_name "DNA_polymerase";
Contig3	Prokka	gene	60000	65000	.	-	.	gene_id "gene_005"; gene_name "membrane_protein";
Contig3	Prokka	CDS	60200	64800	.	-	0	gene_id "gene_005"; gene_name "membrane_protein";
Contig3	Prokka	gene	70000	75000	.	+	.	gene_id "gene_006"; gene_name "kinase_domain";
Contig3	Prokka	CDS	70200	74800	.	+	0	gene_id "gene_006"; gene_name "kinase_domain";
Contig4	Prokka	gene	75000	80000	.	+	.	gene_id "gene_007"; gene_name "transcriptional_regulator";
Contig4	Prokka	CDS	75200	79800	.	+	0	gene_id "gene_007"; gene_name "transcriptional_regulator";
Contig5	Prokka	gene	76000	78000	.	-	.	gene_id "gene_008"; gene_name "oxidoreductase";
Contig5	Prokka	CDS	76200	77800	.	-	0	gene_id "gene_008"; gene_name "oxidoreductase";
Contig6	Prokka	gene	15000	20000	.	+	.	gene_id "gene_009"; gene_name "hydrolase";
Contig6	Prokka	CDS	15200	19800	.	+	0	gene_id "gene_009"; gene_name "hydrolase";
"""
    return gtf_content

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_sam(single_sam):
    """Test that the parsing function works on a single sam file"""

def test_multiple_pafs(multiple_pafs):
    """Test that parsing function works on multiple sam files"""

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for sam_file in glob.glob("*.sam"):
            if os.path.exists(sam_file):
                os.unlink(sam_file)
    request.addfinalizer(cleanup)
