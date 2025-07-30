#!/usr/bin/env python3
"""
Paired immunoglobulin sequencing - Test-Driven Development 
========================================================

Starting with basic concepts in a TDD approach to learn the basics of dev-ops. 
Building simple bioinformatics functions step by step.

ULTIMATE GOALS:
1. Write simple tests first
2. Implement basic functions
3. Understand the Red-Green-Refactor cycle
4. Build confidence with testing
5. Impliment more complex error handling

As functions are built they will replace the TODO sections 
"""

import pytest
import tempfile
import os
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

# =============================================================================
# Integration test fixtures - Shared test data
# =============================================================================
def generate_random_sequence(length, gc_content=0.5):
    bases=['A', 'T', 'G', 'C']
    weights = [(1-gc_content)/2, (1-gc_content)/2, gc_content/2, gc_content/2]
    sequence = ''.join(random.choices(bases, weights=weights, k=length))
    return sequence

# Generate the sequences with a consistent seed
random.seed(42)
perfect_match = generate_random_sequence(64000) # 64kb
bad_match = generate_random_sequence(64000, gc_content=0.52)
random.seed(101010)
bad_coverage_start = generate_random_sequence(8000)
bad_coverage_middle = perfect_match[16000:48000]
random.seed(20547391)
bad_coverage_end = generate_random_sequence(8000)
bad_coverage = bad_coverage_start + bad_coverage_middle + bad_coverage_end
bad_quality_bases = list(perfect_match)
for i in range(0, len(bad_quality_bases), 50):
    if i < len(bad_quality_bases):
        bad_quality_bases[i] = 'N'
        if i + 1 < len(bad_quality_bases):
            bad_quality_bases[i + 1] = 'N'
bad_quality = "".join(bad_quality_bases)

@pytest.fixture
def multiple_assemblies():
    """Create assembly files where contig1 will match bin sequences well, others poorly."""
    assembly_1_content = f">a1_contig1\n{perfect_match}\n>a1_contig2\n{bad_match}\n>a1_contig3\n{bad_coverage}\n>a1_contig4\n{bad_quality}"
    assembly_2_content = f">a2_contig1\n{bad_match}\n>a2_contig2\n{perfect_match}\n>a2_contig3\n{bad_coverage}\n>a2_contig4\n{bad_quality}"

    assembly_files = []
    for i, content in enumerate([assembly_1_content, assembly_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.fasta') as tmp:
            tmp.write(content)
            assembly_files.append(tmp.name)
    return assembly_files

@pytest.fixture
def sample_bins():
    """Create bin files where one sequence matches contig1 perfectly, others don't match well."""
    bin1_content = f">bin1_scaffold1\n{perfect_match}\n>bin1_scaffold2\n{generate_random_sequence(3000)}"
    bin2_content = f">bin2_scaffold1\n{generate_random_sequence(5000)}\n>bin2_scaffold2\n{perfect_match[:20000]}"

    bin_files = []
    for i, content in enumerate([bin1_content, bin2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_bin{i}.fasta') as tmp:
            tmp.write(content)
            bin_files.append(tmp.name)
    return bin_files

# =============================================================================
# Import the functions from the individual_functions folder
# =============================================================================
from individual_functions.contig_mapping import build_minimap2_command, run_alignment
from individual_functions.PAF_parsing import PAF_parsing

# =============================================================================
# Integration test - Test functions working together
# =============================================================================
class TestFullPipeline:
    def test_complete_pipeline(self, sample_bins, multiple_assemblies):
        """
        Test to check if we're ready for integration testing.
        This test checks if individual functions can be imported.
        """
        # TODO: Uncomment as individual functions are completed

        # Step 1: Align
        aligned_paf_files = run_alignment(sample_bins, multiple_assemblies)
        for paf_file in aligned_paf_files:
            paf_out = pd.read_csv(paf_file, delimiter='\t', header=None)
            query_names = paf_out.iloc[:, 0].tolist()
            mapped_reads = ["a1_contig1", "a2_contig2"]
            found_reads = set(mapped_reads) & set(query_names)
            assert found_reads, f"None of {mapped_reads} found in PAF file {paf_file} queries: {set(query_names)}"

        # Step 2: Parse the PAF output
        good_reads = ["a1_contig1", "a2_contig2"]
        poor_reads = ["a1_contig2", "a1_contig3", "a1_contig4", "a2_contig1", "a2_contig3", "a2_contig4"]
        for paf_file in aligned_paf_files:
            PAF_parsing(aligned_paf_files)
        all_read_names = []
        for input_file in aligned_paf_files:
            print(input_file)
            expected_output = os.path.splitext(input_file)[0] + ".tsv"
            output_file = pd.read_csv(expected_output, delimiter='\t', header=None)
            read_names = output_file.iloc[:, 0].tolist()
            all_read_names.extend(read_names)
            found_good = set(good_reads) & set(all_read_names)
            assert found_good, "a1_contig1 and a2_contig2 are NOT in the output"
            for poor in poor_reads:
                assert poor not in read_names, f"'{poor}' was found in output!"

        print(f"✅ Pipeline completed and tests all passed!")

@pytest.fixture(autouse=True)
def cleanup_test_files(request):
    def cleanup():
        import glob
        for paf_file in glob.glob("*.paf"):
            if os.path.exists(paf_file):
                os.unlink(paf_file)
        for tsv_file in glob.glob("*.tsv"):
            if os.path.exists(tsv_file):
                os.unlink(tsv_file)
    request.addfinalizer(cleanup)

