#!/usr/bin/env python3
"""
PAF_parsing.py - Parse the output PAF
==========================================

Function to parse the output PAF from contig_mapping.py selecting bins mapped to
contigs with the specified identity and coverage

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import pandas as pd

# =============================================================================
# Actual functions to test
# =============================================================================
def PAF_parsing(paf_file, assembly_file, identity_threshold = 0.95, coverage_threshold = 0.8, quality_threshold = 40):
    paf = pd.read_csv(paf_file, delimiter='\t', header=None, dtype={
        1: 'int64', # query_length
        2: 'int64', # query_start
        3: 'int64', # query_end
        8: 'int64', # target_end
        9: 'int64', # matches
        10: 'int64', # alignment_length
        11: 'int64' # mapping_quality
    })
    paf_query_end = paf.iloc[:,3]
    paf_query_start = paf.iloc[:,2]
    paf_alignment_length = paf.iloc[:,10]
    paf_matches = paf.iloc[:,9]
    paf_query_length = paf.iloc[:,1]
    paf_alignment_quality = paf.iloc[:,11]
    paf_base=os.path.basename(paf_file)
    paf_name=os.path.splitext(paf_file)[0]
    output_name=paf_name
    output_path=(output_name)
    paf['Identity']=paf_matches/paf_alignment_length
    paf['Coverage']=(paf_query_end - paf_query_start)/paf_query_length
    good_alignments = paf[(paf['Identity'] > identity_threshold) &
    (paf['Coverage'] > coverage_threshold) &
    (paf_alignment_quality > quality_threshold)]
    print(f"DEBUG: Content in good alignments is:")
    print(f"{good_alignments}")
    mapping = good_alignments[[0, 5]].drop_duplicates()
    mapping.columns = ['Contig', 'Bin']
    print(f"DEBUG: Content in mapping is:")
    print(f"DEBUG: {mapping}")
    passed_contigs = mapping['Contig'].unique()
    paf_contigs = paf[0].unique()
    assembly_contigs = set()
    with open(assembly_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_name = line[1:].split()[0]  # Take first word after '>'
                assembly_contigs.add(contig_name)
    all_contigs = paf_contigs | assembly_contigs
    failed_contigs = set(all_contigs) - set(passed_contigs)
    unbinned_contigs = pd.DataFrame({
        'Contig': list(failed_contigs),
        'Bin': ['unbinned'] * len(failed_contigs)
    })
    final_mapping = pd.concat([mapping, unbinned_contigs], ignore_index = True)
    mapping_output = f"{output_path}_contigs_to_bin_mapping.txt"
    final_mapping.to_csv(mapping_output, sep="\t", header=True, index=False)
    print(f"Contig-to-bin mapping file written to: {mapping_output}")
    print(f"The map has the following content:")
    print(f"{final_mapping}")

def run_paf_parsing(paf_files, identity_threshold = 0.95, coverage_threshold = 0.8, quality_threshold = 40):
    """Parse PAF files to select bins meeting user requirement specifications"""
    if isinstance(paf_files, str):
        paf_file = paf_files
        PAF_parsing(paf_file)
    elif isinstance(paf_files, list):
        for paf_file in paf_files:
            PAF_parsing(paf_file)
    else:
        print("There is no valid PAF file found")

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_paf():
    """Create a single test paf file for unit testing."""
    paf_content = """read_001\t100000\t0\t99995\t+\tcontig_1\t750000\t500000\t599995\t99995\t99995\t55
read_002\t1000\t0\t695\t+\tcontig_2\t375000\t375\t1070\t675\t695\t60
read_003\t10000\t0\t8995\t+\tcontig_3\t175000\t10000\t18995\t2995\t8995\t54
read_004\t7000\t3000\t6750\t+\tcontig_4\t55000\t7500\t11250\t3713\t3750\t21"""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1.paf') as tmp:
        tmp.write(paf_content)
        paf_file = tmp.name
    return [paf_file]

@pytest.fixture
def multiple_pafs():
    """Create multiple paf files for integration testing."""
    paf_1_content = """read_001\t100000\t0\t99995\t+\tcontig_1\t750000\t500000\t599995\t99995\t99995\t55
read_002\t1000\t0\t695\t+\tcontig_2\t375000\t375\t1070\t675\t695\t60
read_003\t10000\t0\t8995\t+\tcontig_3\t175000\t10000\t18995\t2995\t8995\t54
read_004\t7000\t3000\t6750\t+\tcontig_4\t55000\t7500\t11250\t3713\t3750\t21"""

    paf_2_content = """read_001\t1000000\t0\t909995\t+\tcontig_4\t2500000\t550000\t1459995\t891795\t909995\t55
read_002\t4000\t1000\t1695\t+\tcontig_5\t375000\t375\t1070\t675\t695\t60
read_003\t10000\t0\t8995\t+\tcontig_6\t175000\t20000\t28995\t2995\t8995\t54
read_004\t7000\t3000\t6750\t+\tcontig_7\t55000\t7500\t11250\t3713\t3750\t21"""

    paf_files = []
    for i, content in enumerate([paf_1_content, paf_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.paf') as tmp:
            tmp.write(content)
            paf_files.append(tmp.name)
    return paf_files



# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_paf(single_paf):
    """Test that the parsing function works on a single paf file"""
    binned_read = "read_001"
    unbinned_reads = ["read_002", "read_003", "read_004"]
    run_paf_parsing(single_paf)
    for input_file in single_paf:
        expected_output = os.path.splitext(input_file)[0] + "_contigs_to_bin_mapping.txt"
        output_file = pd.read_csv(expected_output, delimiter='\t', header=0)
        read_names = output_file.iloc[:, 0].tolist()
        print(f"DEBUG: Content in output file is:")
        print(f"{output_file}")
        binned = output_file.loc[output_file['Contig'] == binned_read]
        unbinned = output_file.loc[output_file['Contig'].isin(unbinned_reads)]
        assert binned['Bin'].iloc[0] in ["contig_1", "contig_4"], f"read_001 incorrectly assigned to {binned['Bin'].iloc[0]}"
        for _, row in unbinned.iterrows():
            assert row['Bin'] == "unbinned", f"{row['Contig']} should be unbinned but assigned to {row['Bin']}!"

def test_multiple_pafs(multiple_pafs):
    """Test that parsing function works on multiple pafs"""
    binned_read = "read_001"
    unbinned_reads = ["read_002", "read_003", "read_004"]
    run_paf_parsing(multiple_pafs)
    for input_file in multiple_pafs:
        expected_output = os.path.splitext(input_file)[0] + "_contigs_to_bin_mapping.txt"
        output_file = pd.read_csv(expected_output, delimiter='\t', header=0)
        read_names = output_file.iloc[:, 0].tolist()
        binned = output_file.loc[output_file['Contig'] == binned_read]
        unbinned = output_file.loc[output_file['Contig'].isin(unbinned_reads)]
        assert binned['Bin'].iloc[0] in ["contig_1", "contig_4"], f"read_001 incorrectly assigned to {binned['Bin'].iloc[0]}"
        for _, row in unbinned.iterrows():
            assert row['Bin'] == "unbinned", f"{row['Contig']} should be unbinned but assigned to {row['Bin']}!"

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for tsv_file in glob.glob("*.tsv"):
            if os.path.exists(tsv_file):
                os.unlink(tsv_file)
    request.addfinalizer(cleanup)
