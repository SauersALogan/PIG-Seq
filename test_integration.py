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

# =============================================================================
# Integration test fixtures - Shared test data
# =============================================================================
@pytest.fixture
def multiple_assemblies():
    """Create assembly files where contig1 will match bin sequences well, others poorly."""
    assembly_1_content = """
    >contig1 length=1000
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >contig2 length=800
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    >contig3 length=1200
    TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
    >contig4 length=600
    AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT
    """

    assembly_2_content = """
    >contig1 length=1000
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCAATCGATCGATCGATCGATCGATCG
    >contig2 length=800
    GCTAGCTAGCTAGCTAGCtAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    >contig3 length=1200
    TTTTAAAACCCCGGGGTTTTAAACCCCGGGGTTTTAAAACCCCGGGG
    >contig4 length=600
    AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT
    """

    assembly_files = []
    for i, content in enumerate([assembly_1_content, assembly_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.fasta') as tmp:
            tmp.write(content)
            assembly_files.append(tmp.name)
    return assembly_files

@pytest.fixture
def sample_bins():
    """Create bin files where one sequence matches contig1 perfectly, others don't match well."""
    bin1_content = """
    >bin1_scaffold1 length=1000
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >bin1_scaffold2 length=500
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    """

    bin2_content = """
    >bin2_scaffold1 length=800
    TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
    """

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
            assert "contig1" in query_names, f"contig1 not found in PAF file {paf_file} queries: {set(query_names)}"

        # Step 2: Parse the PAF output
        good_read = "contig1"
        poor_reads = ["contig2", "contig3", "contig4"]
        for paf_file in aligned_paf_files:
            PAF_parsing(paf_file)
        for input_file in aligned_paf_files:
            print(input_file)
            expected_output = os.path.splitext(input_file)[0] + ".tsv"
            output_file = pd.read_csv(expected_output, delimiter='\t', header=None)
            read_names = output_file.iloc[:, 0].tolist()
            assert good_read in read_names, "contig1 is in the output"
            for poor in poor_reads:
                assert poor not in read_names, f"'{poor}' was found in output!"

        print(f"✅ Pipeline processed {len(pipeline_results)} combinations")

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
