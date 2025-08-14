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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import re
import shutil

# =============================================================================
# Actual functions to test
# =============================================================================
def extract_file_identifiers(filename):
    """Extracts the similar file identifier for all input files to ensure proper pairing.
    First it looks for common patterns such as Sample_XXX, then numerical suffixes such as sample1. If these fail
    it will finish by looking at the last part of the file after controlling for common delimiters"""
    basename = os.path.basename(filename)
    basename = re.sub(r'\.(txt|sam|bam|gtf|gff|paf|fasta|fa|fastq|fq)$', '', basename, flags=re.IGNORECASE)
    match = re.search(r'(\d+)', basename)
    if match:
        return match.group(1)
    return basename

def pair_files_by_sample(files_list1, files_list2, list1_name="files1", list2_name="files2"):
    """This function will pair files based upon the layouts controlled for in the function above"""
    list1_files_by_id = {}
    for file in files_list1:
        identifier = extract_file_identifiers(file)
        list1_files_by_id[identifier] = file
    list2_files_by_id = {}
    for file in files_list2:
        identifier = extract_file_identifiers(file)
        list2_files_by_id[identifier] = file
    paired_files = []
    unmatched_list1 = []
    unmatched_list2 = []
    for identifier in list1_files_by_id.keys():
        if identifier in list2_files_by_id:
            file1 = list1_files_by_id[identifier]
            file2 = list2_files_by_id[identifier]
            paired_files.append((file1, file2))
            print(f"DEBUG: Paired {os.path.basename(file1)} with {os.path.basename(file2)} (ID: {identifier})")
        else:
            unmatched_list1.append(list1_files_by_id[identifier])
    for identifier in list2_files_by_id.keys():
        if identifier not in list1_files_by_id:
            unmatched_list2.append(list2_files_by_id[identifier])
    if unmatched_list1:
        print(f"WARNING: Unmatched {list1_name}: {[os.path.basename(f) for f in unmatched_list1]}")
    if unmatched_list2:
        print(f"WARNING: Unmatched {list2_name}: {[os.path.basename(f) for f in unmatched_list2]}")
    return paired_files

# =============================================================================
# Create the mock data for testing
# =============================================================================

# =============================================================================
# Setup the actual tests
# =============================================================================

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for tsv_file in glob.glob("*.tsv"):
            if os.path.exists(tsv_file):
                os.unlink(tsv_file)
    request.addfinalizer(cleanup)
