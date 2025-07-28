#!/usr/bin/env python3
"""
PAF_parsing.py - Parse the output PAF
==========================================

Function to parse the output PAF from contig_mapping.py selecting bins mapped to contigs with the specified identity and coverage
This file contains the function and its tests
"""

import tempfile
import os
import pytest
import subprocess # Needed to run external command sin python

# =============================================================================
# Actual functions to test
# =============================================================================
def PAF_parsing(paf_files, identity_threshold = 0.95, coverage_threshold = 0.8, quality_threshold = 40):
    """Parse PAF files to select bins meeting user requirement specifications"""
    if isinstance(paf_files, str):
        paf = pd.read_csv(paf_files, delimiter='\t')
        paf_query_end = paf.iloc[:,3]
        paf_query_start = paf.iloc[:,2]
        paf_alignment_length = paf.iloc[:,10]
        paf_matches = paf.iloc[:,9]
        paf_query_length = paf.iloc[:,1]
        paf_alignment_quality = paf.iloc[:,8]
        paf_base=os.path.basename(paf_files)
        paf_name=os.path.splitext(paf_files)[0]
        output_name=paf_name
        output_path=(output_name)
        paf['Identity']=paf_matches/paf_alignment_length
        paf['Coverage']=(paf_query_end - paf_query_start)/paf_query_length
        good_alignments = paf[(paf['Identity'] > identity_threshold) & 
        (paf['Coverage'] > coverage_threshold) & 
        (paf_alignment_quality > quality_threshold)]  
        good_alignments.to_csv(f"{output_path}.tsv",
            sep='\t',header=False,index=False) 
    elif isinstance(paf_files, list):
        for file in paf_files:
            paf = pd.read_csv(file, delimiter='\t')
            paf_query_end = paf.iloc[:,3]
            paf_query_start = paf.iloc[:,2]
            paf_alignment_length = paf.iloc[:,10]
            paf_matches = paf.iloc[:,9]
            paf_query_length = paf.iloc[:,1]
            paf_alignment_quality = paf.iloc[:,8]
            paf_base=os.path.basename(file)
            paf_name=os.path.splitext(file)[0]
            output_name=paf_name
            output_path=(output_name)
            paf['Identity']=paf_matches/paf_alignment_length
            paf['Coverage']=(paf_query_end - paf_query_start)/paf_query_length
            good_alignments = paf[(paf['Identity'] > identity_threshold) & 
            (paf['Coverage'] > coverage_threshold) & 
            (paf_alignment_quality > quality_threshold)]  
            good_alignments.to_csv(f"{output_path}.tsv",
                sep='\t',header=False,index=False) 
    else:
        print("There is no valid PAF file found")

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_paf():
    """Create a single test paf file for unit testing."""
    paf_content = """
	read_001\t100000\t0\t99995\t+\tcontig_1\t750000\t500000
		\t599995\t99995\t99995\t55
	read_002\t1000\t0\t695\t+\tcontig_2\t375000\t375
		\t1070\t675\t695\t60
	read_003\t10000\t0\t8995\t+\tcontig_3\t175000\t10000
		\t18995\t2995\t8995\t54
	read_004\t7000\t3000\t6750\t+\tcontig_4\t55000\t7500
		\t11250\t3713\t3750\t21
	"""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1.paf') as tmp:
        tmp.write(assembly_content)
        assembly_file = tmp.name
    return [paf_file]

@pytest.fixture
def multiple_pafs():
    """Create multiple paf files for integration testing."""
    paf_1_content = """
	read_001\t100000\t0\t99995\t+\tcontig_1\t750000\t500000
		\t599995\t99995\t99995\t55
	read_002\t1000\t0\t695\t+\tcontig_2\t375000\t375
		\t1070\t675\t695\t60
	read_003\t10000\t0\t8995\t+\tcontig_3\t175000\t10000
		\t18995\t2995\t8995\t54
	read_004\t7000\t3000\t6750\t+\tcontig_4\t55000\t7500
		\t11250\t3713\t3750\t21
	"""
    paf_2_content = """
    read_001\t1000000\t0\t909995\t+\tcontig_4\t2500000\t550000
		\t1459995\t891795\t909995\t55
	read_002\t4000\t1000\t1695\t+\tcontig_5\t375000\t375
		\t1070\t675\t695\t60
	read_003\t10000\t0\t8995\t+\tcontig_6\t175000\t20000
		\t28995\t2995\t8995\t54
	read_004\t7000\t3000\t6750\t+\tcontig_7\t55000\t7500
		\t11250\t3713\t3750\t21
	"""

    paf_files = []
    for i, content in enumerate([paf_1_content, paf_2_content, paf_3_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.paf') as tmp:
            tmp.write(content)
            paf_files.append(tmp.name)
    return paf_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_paf(single_paf):
    """Test that the parsing function works on a single paf file"""
    good_read = "read_001"
    poor_reads = ["read_002", "read_003", "read_004"]
    PAF_parsing(single_paf)
    for input_file in single_paf:
        expected_output = os.path.splitext(input_file)[0] + ".tsv"
        output_file = pd.read_csv(expected_output, delimiter='\t', header=None)
        read_names = output_file.iloc[:, 0].tolist()
        assert good_read in read_names, "read_001 is in the output" 
        for poor in poor_reads:
            assert poor not in read_names, f"'{poor}' was found in output!"

def test_multiple_pafs(multiple_pafs):
    """Test that parsing function works on multiple pafs"""
    good_read = "read_001"
    poor_reads = ["read_002", "read_003", "read_004"]
    PAF_parsing(multiple_pafs)
    for input_file in multiple_pafs:
        expected_output = os.path.splitext(input_file)[0] + ".tsv"
        output_file = pd.read_csv(expected_output, delimiter='\t', header=None)
        read_names = output_file.iloc[:, 0].tolist()
        assert good_read in read_names, "read_001 is in the output" 
        for poor in poor_reads:
            assert poor not in read_names, f"'{poor}' was found in output!"

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for tsv_file in glob.glob("*.tsv"):
            if os.path.exists(tsv_file):
                os.unlink(tsv_file)
    request.addfinalizer(cleanup)
