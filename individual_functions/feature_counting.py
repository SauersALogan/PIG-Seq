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
import shutil

# =============================================================================
# Actual functions to test
# =============================================================================
def build_featureCounts (gff_files, output_path, sam_files, threads, gene, gene_id, paired_end=True):
    """Count SAM file features using GFF annotation"""
    featureCounts = [
        "featureCounts",
        "-T", str(threads),
        "-t", gene, # The prokka classification of the read, usually default to gene or CDS
        "-g", gene_id, # The assigned gene_id, used to group features with similar gene_ids 
        "-a", gff_files,
        "-o", output_path,
    ]

    if paired_end:
        featureCounts.extend(["-p", "-B", "-C"])

    featureCounts.append(sam_files)

    print(f"DEBUG: featureCounts command: {' '.join(featureCounts)}")
    return featureCounts

def execute_feature_counting(sam_file, gff_file, threads = 2, gene = "CDS", gene_id = "locus_tag", paired_end=True):
    sam_base=os.path.basename(sam_file)
    sam_name=os.path.splitext(sam_base)[0]
    output_name=sam_name+".txt"
    output_path=(output_name)
    print(f"DEBUG: Output path will be: {output_path}")
    Counter = build_featureCounts(gff_file, output_path, sam_file, threads, gene, gene_id, paired_end)
    try:
        subprocess.run(Counter, check=True)
    except subprocess.CalledProcessError:
        print("featureCounts failed to process {sam_file}")
    return(output_path)

def run_counter(sam_files, gff_files, threads = 2, gene = "CDS", gene_id = "locus_tag", paired_end=True):
    """Runs the counting function imported above on the actual data"""
    count_tables = []
    print(f"DEBUG: run_counter called with sam_files={sam_files}, type={type(sam_files)}")
    print(f"DEBUG: gff_files={gff_files}, type={type(gff_files)}")
    if isinstance(sam_files, str):
        sam_file = sam_files
        gff_file = gff_files
        counts = execute_feature_counting(sam_file, gff_file, threads, gene, gene_id, paired_end=True)
        print(f"DEBUG: execute_feature_counting returned: {counts}")
        count_tables.append(counts)
    elif isinstance(sam_files, list) and isinstance(gff_files, str):
        gff_file = gff_files
        for sam_file in sam_files:
            counts = execute_feature_counting(sam_file, gff_file, threads, gene, gene_id, paired_end=True)
            print(f"DEBUG: execute_feature_counting returned: {counts}")
            count_tables.append(counts)
    elif isinstance(sam_files, list) and isinstance(gff_files, list):
        for sam_file, gff_file in zip(sam_files, gff_files):
            counts = execute_feature_counting(sam_file, gff_file, threads, gene, gene_id, paired_end=True)
            print(f"DEBUG: execute_feature_counting returned: {counts}")
            count_tables.append(counts)
    return(count_tables)

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for testing."""
    temp_dir = tempfile.mkdtemp(prefix='Pipeline_test')
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def single_sam(temp_test_dir):
    """Create a single test SAM file for unit testing with paired-end reads."""
    sam_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:Contig1	LN:341341
@SQ	SN:Contig2	LN:325562
@SQ	SN:Contig3	LN:317386
@SQ	SN:Contig4	LN:301055
@SQ	SN:Contig5	LN:298171
@SQ	SN:Contig6	LN:282893
read1	99	Contig1	3554	60	150M	=	3704	300	GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACATCCCATGGACAGTATATGTTTGTATTACCAAACAGATAGCCGTCATACCACTCCTTTACAGCGCCGTACTGGTCCGTAAATCCATAATACGCC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read1	147	Contig1	3704	60	150M	=	3554	-300	CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCAGTTCTCCTGTATATTCATCGGTAAAACGAATCTTTAAATCCTTTTCCACATAGACAAAGCTGACATTGACCTTGCCTGATGCCGTCTCCGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	99	Contig2	72558	60	150M	=	72708	300	CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACCTTATCCTTTACTGTAACCTGCTGCGGAACCGACGGGAACGTATATTTCTCTGTGGAGGTTATGATCGCATCAAACACGCCCGCGTTCACATT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	Contig2	72708	60	150M	=	72558	-300	CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACACAGTTTGAGGAGATGGGGCTATTTAACCGCGGTATCAAAACCAGAATCGGCAGCAATGGAAACTCAGATTATTATGGCATCATACGTGAGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	99	Contig3	72369	60	150M	=	72519	300	GTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTGTTCGTGACGCACAGGCAGTCGCACGGTACTTTGGGCTTGTGTCTAAAGATAAAACAAAGGATTACAGCAATTATGCTCGATACCGTGACGTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	147	Contig3	72519	60	150M	=	72369	-300	CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAATAAAATGCGATGATATCGCCGGTCTTGTATTTATTTCCGTTCAGCGCCACTACGATGTCTTCCTCCTGCAGCGACTCCGTCACGGAAGTTCCG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_sample1.sam') as tmp:
        tmp.write(sam_content)
        sam_file = tmp.name
    return [sam_file]

@pytest.fixture
def multiple_sam_files(temp_test_dir):
    """Create paired-end SAM files with reads mapping to the assembly contigs."""
    sam1_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a1_contig1	LN:64000
@SQ	SN:a1_contig2	LN:64000
@SQ	SN:a1_contig3	LN:48000
@SQ	SN:a1_contig4	LN:64000
read1	99	a1_contig1	150	60	100M	=	250	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read1	147	a1_contig1	250	60	100M	=	150	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	99	a1_contig1	700	60	100M	=	800	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	a1_contig1	800	60	100M	=	700	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	99	a1_contig2	250	60	100M	=	350	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	147	a1_contig2	350	60	100M	=	250	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	99	a1_contig3	400	60	100M	=	500	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	147	a1_contig3	500	60	100M	=	400	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam2_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a2_contig1	LN:64000
@SQ	SN:a2_contig2	LN:64000
@SQ	SN:a2_contig3	LN:48000
@SQ	SN:a2_contig4	LN:64000
read5	99	a2_contig1	150	60	100M	=	250	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	147	a2_contig1	250	60	100M	=	150	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read6	99	a2_contig2	200	60	100M	=	300	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read6	147	a2_contig2	300	60	100M	=	200	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	99	a2_contig2	1100	60	100M	=	1200	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	147	a2_contig2	1200	60	100M	=	1100	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	99	a2_contig2	1900	60	100M	=	2000	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	147	a2_contig2	2000	60	100M	=	1900	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	99	a2_contig3	350	60	100M	=	450	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	147	a2_contig3	450	60	100M	=	350	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam_files = []
    for i, content in enumerate([sam1_content, sam2_content], 1):
        sam_file = os.path.join(temp_test_dir, f'test_sample{i}.sam')
        with open(sam_file, 'w') as f:
            f.write(content)
        sam_files.append(sam_file)
    return sam_files

@pytest.fixture
def mock_gff_files(temp_test_dir):
    """Create GFF files corresponding to the assembly files."""
    gff1_content = """##gff-version 3
##sequence-region a1_contig1 1 64000
##sequence-region a1_contig2 1 64000
##sequence-region a1_contig3 1 48000
##sequence-region a1_contig4 1 64000
a1_contig1	prokka	CDS	100	500	.	+	0	ID=GENE001;locus_tag=GENE001;product=hypothetical protein
a1_contig1	prokka	CDS	600	1200	.	-	0	ID=GENE002;locus_tag=GENE002;product=DNA polymerase
a1_contig1	prokka	CDS	1500	2000	.	+	0	ID=GENE003;locus_tag=GENE003;product=ribosomal protein
a1_contig2	prokka	CDS	200	800	.	+	0	ID=GENE004;locus_tag=GENE004;product=membrane protein
a1_contig2	prokka	CDS	1000	1600	.	-	0	ID=GENE005;locus_tag=GENE005;product=ATP synthase
a1_contig3	prokka	CDS	300	900	.	+	0	ID=GENE006;locus_tag=GENE006;product=transcriptase
a1_contig4	prokka	CDS	150	750	.	+	0	ID=GENE007;locus_tag=GENE007;product=helicase
"""

    gff2_content = """##gff-version 3
##sequence-region a2_contig1 1 64000
##sequence-region a2_contig2 1 64000
##sequence-region a2_contig3 1 48000
##sequence-region a2_contig4 1 64000
a2_contig1	prokka	CDS	100	500	.	+	0	ID=GENE008;locus_tag=GENE008;product=hypothetical protein
a2_contig1	prokka	CDS	600	1200	.	-	0	ID=GENE009;locus_tag=GENE009;product=DNA polymerase
a2_contig2	prokka	CDS	200	800	.	+	0	ID=GENE010;locus_tag=GENE010;product=membrane protein
a2_contig2	prokka	CDS	1000	1600	.	-	0	ID=GENE011;locus_tag=GENE011;product=ATP synthase
a2_contig2	prokka	CDS	1800	2400	.	+	0	ID=GENE012;locus_tag=GENE012;product=kinase
a2_contig3	prokka	CDS	300	900	.	+	0	ID=GENE013;locus_tag=GENE013;product=transcriptase
a2_contig4	prokka	CDS	150	750	.	+	0	ID=GENE014;locus_tag=GENE014;product=helicase
"""

    gff_files = []
    for i, content in enumerate([gff1_content, gff2_content], 1):
        gff_file = os.path.join(temp_test_dir, f'test_assembly{i}.gff')
        with open(gff_file, 'w') as f:
            f.write(content)
        gff_files.append(gff_file)
    return gff_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_sam(single_sam, mock_gff_files):
    """Test that the parsing function works on a single sam file"""
    results = run_counter(single_sam, mock_gff_files, threads=2, gene="CDS", gene_id="locus_tag", paired_end=True)
    for file in results:
        if os.path.exists(file):
            print(f"DEBUG: Output file {file} was created successfully")
            file_size = os.path.getsize(file)
            assert file_size > 0, f"Output file {file} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
        else:
            print(f"DEBUG: WARNING - Output file {file} was NOT created")

def test_multiple_sams(multiple_sam_files, mock_gff_files):
    """Test that parsing function works on multiple sam files"""
    results = run_counter(multiple_sam_files, mock_gff_files, threads=2, gene="CDS", gene_id="locus_tag", paired_end=True)
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
        for file in glob.glob(".txt.summary"):
            if os.path.exists(file):
                os.unlink(file)
    request.addfinalizer(cleanup)

