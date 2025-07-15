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
def build_minimap2_command(assembly, bin, output_path):
    """"Build the minimap2 command as a list - Need run_alignment to execute."""
    aligner = [ 
        "minimap2",
        "-x", "asm5",
        assembly, bin,
        "-o", output_path
    ]
    return aligner

def run_alignment(assembly, bin, output_path):
    """Run minimap2 alignment and return PAF files"""
    try:
        aligner = build_minimap2_command(assembly, bin, output_path)
        result = subprocess.run(aligner, check=True)
        return output_path
    except subprocess.CalledProcessError:
        print("Minimap2 alignment failed for {assembly}")
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

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1.fasta') as tmp:
        tmp.write(assembly_content)
        assembly_file = tmp.name
    return [assembly_file]

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

def test_single_assembly(single_assembly, sample_bins):
    """Test that minimap2 actually runs and produces output on a single assembly"""
    for assembly in single_assembly:
        for bin in sample_bins:
            assembly_base=os.path.basename(assembly)
            assembly_name=os.path.splitext(assembly_base)[0]
            bin_base=os.path.basename(bin)
            bin_name=os.path.splitext(bin_base)[0]
            output_name=assembly_name+"_"+bin_name+".paf"
            output_path=(output_name)
            results = run_alignment(assembly, bin, output_path)
            assert os.path.exists(results), "PAF file should be created"
            assert os.path.getsize(results) > 0, "PAF file should not be empty"

def test_multiple_assemblies(multiple_assemblies, sample_bins):
    """Test that minimap2 runs and produces proper output with multiple assemblies"""
    for assembly in multiple_assemblies:
        for bin in sample_bins:
            assembly_base=os.path.basename(assembly)
            assembly_name=os.path.splitext(assembly_base)[0]
            bin_base=os.path.basename(bin)
            bin_name=os.path.splitext(bin_base)[0]
            output_name=assembly_name+"_"+bin_name+".paf"
            output_path=(output_name)
            results = run_alignment(assembly, bin, output_path)
            assert os.path.exists(results), "PAF file should be created"
            assert os.path.getsize(results) > 0, "PAF file should not be empty"

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for paf_file in glob.glob("*.paf"):
            if os.path.exists(paf_file):
                os.unlink(paf_file)
    request.addfinalizer(cleanup)
