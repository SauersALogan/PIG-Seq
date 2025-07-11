#!/usr/bin/env python3
"""
contig_mapping.py - Map assembly contigs to bins
==========================================

Function to map assembly contigs to bins
This file contains the function and its tests.
"""

import tempfile
import os
import pytest
import subprocess # Needed to run external command sin python

# =============================================================================
# Actual functions to test
# =============================================================================

def build_minimap2_command(assembly_file, bin_files):
    """"Build the minimap2 command as a list - Need run_alignment to execute."""
    output_paf = "alignment_output.paf"

    cmd = [ 
        "minimap2",
        "-x", "asm5", 
        "-o", output_paf, 
    ]

    if isinstance(bin_files, list): 
        cmd.extend(bin_files) 
    else: 
        cmd.append(bin_files)

    return cmd, output_paf

def run_alignment(assembly_file, bin_files):
    """Run minimap2 alignment and return PAF files"""
    try:
        cmd, output_paf = build_minimap2_command(assembly_file, bin_files)
        results = subprocess.run(cmd, check=True)
        return output_paf

    except subprocess.CalledProcessError:
        print("Minimap2 alignment failed")
        return None

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_assembly():
    """Create a sample assembly file for unit testing."""
    assembly_content = """
    >contig1 length=1000
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >contig2 length=800
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    >contig3 length=1200
    TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
    >contig4 length=600
    AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT
    """

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
        tmp.write(assembly_content)
        return tmp.name

@pytest.fixture
def multiple_assemblies():
    """Create multiple assembly files for integration testing."""
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

    assembly_3_content = """
    >contig1 length=1000
    ATCGATCGATCGATCGATCGATGATCGATCGATCAATCGATCGATCGATCGATCGATCG
    >contig2 length=800
    GCTAGCTAGCTAGCTAGCtAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    >contig3 length=1200
    TTTTAAAACCCCGGGGTTACCCCGGGGTTTTAAAACCCCGGGG
    >contig4 length=600
    AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT
    """

    assembly_files = []
    for i, content in enumerate([assembly_1_content, assembly_2_content, assembly_3_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.fasta') as tmp:
            tmp.write(content)
            assembly_files.append(tmp.name)
    return assembly_files

@pytest.fixture
def sample_bins():
    """Create sample bin files for testing."""
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
# Setup the actual tests
# =============================================================================

def test_single_assembly_alignment(single_assembly, sample_bins):
    """Test that minimap2 actually runs and produces output on a single assembly"""

    paf_file = run_alignment(single_assembly, sample_bins)
    assert os.path.exists(paf_file), "PAF file should be created"
    assert os.path.getsize(paf_file) > 0, "PAF file should not be empty"

def test_multiple_assemblies(multiple_assemblies, sample_bins):
    results = []
    for assembly in multiple_assemblies:
        paf_file = run_alignment(assembly, sample_bins)
        assert os.path.exists(paf_file), "PAF file should be created"
        assert os.path.getsize(paf_file) > 0, "PAF file should not be empty"

if __name__ == "__main__":
    # Run tests when file is executed directly
    success = run_all_tests()
    
    if success:
        print("\n" + "=" * 50)
        print("FUNCTION READY FOR INTEGRATION!")
        print("=" * 50)
        print("You can now:")
        print("1. Import this function in integration_test.py")
        print("2. Uncomment the relevant integration tests")
        print("3. Run: pytest integration_test.py -v")
    else:
        exit(1)
