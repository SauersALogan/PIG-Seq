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
        "-p", "--countReadPairs", "-B", "-P",
        "-T", str(threads),
        "-t", gene, # The prokka classification of the read, usually default to gene or CDS
        "-g", gene_id, # The assigned gene_id, used to group features with similar gene_ids 
        "-a", gtf_files,
        "-o", output_path,
        sam_files
    ]
    return featureCounts

def execute_feature_counting(sam_file, gtf_file, threads = 2, gene = "CDS", gene_id = "gene_id"):
    sam_base=os.path.basename(sam_file)
    sam_name=os.path.splitext(sam_base)[0]
    output_name=sam_name+".txt"
    output_path=(output_name)
    print(f"DEBUG: Output path will be: {output_path}")
    Counter = build_featureCounts(gtf_file, output_path, sam_file, threads, gene, gene_id)
    try:
        subprocess.run(Counter, check=True)
    except subprocess.CalledProcessError:
        print("featureCounts failed to process {sam_file}")
    return(output_path)

def run_counter(sam_files, gtf_files, threads = 2, gene = "CDS", gene_id = "gene_id"):
    """Runs the counting function imported above on the actual data"""
    count_tables = []
    print(f"DEBUG: run_counter called with sam_files={sam_files}, type={type(sam_files)}")
    print(f"DEBUG: gtf_files={gtf_files}, type={type(gtf_files)}")
    if isinstance(sam_files, str):
        sam_file = sam_files
        gtf_file = gtf_files
        counts = execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
        print(f"DEBUG: execute_feature_counting returned: {counts}")
        count_tables.append(counts)
    elif isinstance(sam_files, list) and isinstance(gtf_files, str):
        gtf_file = gtf_files
        for sam_file in sam_files:
            counts = execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
            print(f"DEBUG: execute_feature_counting returned: {counts}")
            count_tables.append(counts)
    elif isinstance(sam_files, list) and isinstance(gtf_files, list):
        if len(gtf_files) == 1 and len(sam_files) > 1:
            print("It seems you have a single GTF file being read as a list, we can fix this")
            gtf_file = gtf_files[0]
            for sam_file in sam_files:
                counts = execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
                print(f"DEBUG: execute_feature_counting returned: {counts}")
                count_tables.append(counts)
        elif len(sam_files) == len(gtf_files):
            print("I will sort the files for pairing")
            sam_files.sort()
            gtf_files.sort()
            for sam_file, gtf_file in zip(sam_files, gtf_files):
                counts = execute_feature_counting(sam_file, gtf_file, threads, gene, gene_id)
                print(f"DEBUG: execute_feature_counting returned: {counts}")
                count_tables.append(counts)
    return(count_tables)

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_sam():
    """Create a single test SAM file for unit testing."""
    sam_content = """@SQ	SN:Contig1	LN:341341
@SQ	SN:Contig2	LN:325562
@SQ	SN:Contig3	LN:317386
@SQ	SN:Contig4	LN:301055
@SQ	SN:Contig5	LN:298171
@SQ	SN:Contig6	LN:282893
read1	99	Contig1	1500	60	150M	=	1700	350	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read1	147	Contig1	1700	60	150M	=	1500	-350	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	99	Contig1	3200	60	150M	=	3450	400	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	Contig1	3450	60	150M	=	3200	-400	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	99	Contig2	51000	60	150M	=	51250	400	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	147	Contig2	51250	60	150M	=	51000	-400	CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	99	Contig2	72000	60	150M	=	72250	400	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	147	Contig2	72250	60	150M	=	72000	-400	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	99	Contig3	62000	60	150M	=	62300	450	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	147	Contig3	62300	60	150M	=	62000	-450	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read6	99	Contig3	72100	60	150M	=	72350	400	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read6	147	Contig3	72350	60	150M	=	72100	-400	CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	99	Contig4	77000	60	150M	=	77300	450	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	147	Contig4	77300	60	150M	=	77000	-450	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	99	Contig5	76800	60	150M	=	77100	450	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	147	Contig5	77100	60	150M	=	76800	-450	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	99	Contig6	17000	60	150M	=	17300	450	CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	147	Contig6	17300	60	150M	=	17000	-450	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_sample1.sam') as tmp:
        tmp.write(sam_content)
        sam_file = tmp.name
    return [sam_file]

@pytest.fixture
def multiple_sams():
    """Create multiple SAM files for unit testing."""
    sam_1_content ="""@SQ     SN:Contig1   LN:341341
@SQ     SN:Contig2   LN:325562
@SQ     SN:Contig3   LN:317386
read1	99	Contig1	1500	60	150M	=	1700	350	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read1	147	Contig1	1700	60	150M	=	1500	-350	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	99	Contig2	51000	60	150M	=	51250	400	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	Contig2	51250	60	150M	=	51000	-400	CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam_2_content ="""@SQ     SN:Contig4   LN:301055
@SQ     SN:Contig5   LN:298171
@SQ     SN:Contig6   LN:282893
read3	99	Contig4	77000	60	150M	=	77300	450	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	147	Contig4	77300	60	150M	=	77000	-450	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAAGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	99	Contig5	76800	60	150M	=	77100	450	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	147	Contig5	77100	60	150M	=	76800	-450	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACGTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTCTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam_files = []
    for i, content in enumerate([sam_1_content, sam_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_sample{i}.paf') as tmp:
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
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1.gtf') as tmp:
        tmp.write(gtf_content)
        gtf_file = tmp.name
    return [gtf_file]

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_sam(single_sam, mock_gtf):
    """Test that the parsing function works on a single sam file"""
    results = run_counter(single_sam, mock_gtf, threads=2, gene="CDS", gene_id="gene_id")
    for file in results:
        if os.path.exists(file):
            print(f"DEBUG: Output file {file} was created successfully")
            file_size = os.path.getsize(file)
            assert file_size > 0, f"Output file {file} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
        else:
            print(f"DEBUG: WARNING - Output file {file} was NOT created")

def test_multiple_sams(multiple_sams, mock_gtf):
    """Test that parsing function works on multiple sam files"""
    results = run_counter(multiple_sams, mock_gtf, threads=2, gene="CDS", gene_id="gene_id")
    assert len(results) == len(multiple_sams), f"Expected {len(multiple_sams)} output files, got {len(results)}"
    for file in results:
        if os.path.exists(file):
            print(f"DEBUG: Output file {file} was created successfully")
            file_size = os.path.getsize(file)
            assert file_size > 0, f"Output file {file} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
        else:
            print(f"DEBUG: WARNING - Output file {file} was NOT created")

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for file in glob.glob("*.txt"):
            if os.path.exists(file):
                os.unlink(file)
        for file in glob.glob("*.txt.summary"):
            if os.path.exists(file):
                os.unlink(file)
    request.addfinalizer(cleanup)

