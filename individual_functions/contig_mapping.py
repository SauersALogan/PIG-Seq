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
import pandas as pd

# =============================================================================
# Actual functions to test
# =============================================================================
def build_minimap2_command(bins, assemblies, output_path):
    """"Build the minimap2 command as a list - Need run_alignment to execute."""
    minimap = [
        "minimap2",
        "-x", "asm5",
        bins, assemblies,
        "-o", output_path
    ]
    return minimap

def run_alignment(bins, assemblies):
    """Run minimap2 alignment and return PAF files"""
    if isinstance(assemblies, str):
        assembly = assemblies
        if isinstance(bins, str):
            bin = bins
            assembly_base=os.path.basename(assembly)
            assembly_name=os.path.splitext(assembly_base)[0]
            output_name=assembly_name+".paf"
            output_path=(output_name)
            aligner = build_minimap2_command(bin, assembly, output_path)
            try:
                alignment = subprocess.run(aligner, check=True)
                return output_path
            except subprocess.CalledProcessError:
                print("Minimap2 alignment failed for {assembly}")
            return None
        elif isinstance(bins, list):
            with open("merged_bins.fa", "w") as merged:
                for bin in bins:
                    with open(bin, "r") as input_bin:
                        bin_content = input_bin.read()
                        merged.write(bin_content)
                        merged.write("\n")
            bins = "merged_bins.fa"
            assembly_base=os.path.basename(assembly)
            assembly_name=os.path.splitext(assembly_base)[0]
            output_name=assembly_name+".paf"
            output_path=(output_name)
            aligner = build_minimap2_command(bins, assembly, output_path)
            try:
                alignment = subprocess.run(aligner, check=True)
                return output_path
            except subprocess.CalledProcessError:
                print("Minimap2 alignment failed for {assembly}")
            return None
    if isinstance(assemblies, list):
        if isinstance(bins, str):
            bin = bins
            for assembly in assemblies:
                assembly_base=os.path.basename(assembly)
                assembly_name=os.path.splitext(assembly_base)[0]
                output_name=assembly_name+".paf"
                output_path=(output_name)
                aligner = build_minimap2_command(bin, assembly, output_path)
                try:
                    alignment = subprocess.run(aligner, check=True)
                    return output_path
                except subprocess.CalledProcessError:
                    print("Minimap2 alignment failed for {assembly}")
                return None
        elif isinstance(bins, list):
            with open("merged_bins.fa", "w") as merged:
                for bin in bins:
                    with open(bin, "r") as input_bin:
                        bin_content = input_bin.read()
                        merged.write(bin_content)
                        merged.write("\n")
            bins = "merged_bins.fa"
            for assembly in assemblies:
                assembly_base=os.path.basename(assembly)
                assembly_name=os.path.splitext(assembly_base)[0]
                output_name=assembly_name+".paf"
                output_path=(output_name)
                aligner = build_minimap2_command(bins, assembly, output_path)
                try:
                    alignment = subprocess.run(aligner, check=True)
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
    results = run_alignment(single_assembly, sample_bins)
    paf_out = pd.read_csv(results, delimiter='\t', header=None)
    query_names = paf_out.iloc[:, 0].tolist()
    assert "contig1" in query_names, f"contig1 not found in PAF file queries: {set(query_names)}"

def test_multiple_assemblies(multiple_assemblies, sample_bins):
    """Test that minimap2 runs and produces proper output with multiple assemblies"""
    results = []
    result = run_alignment(multiple_assemblies, sample_bins)
    results.append(result)
    for paf_file in results:
        paf_out = pd.read_csv(paf_file, delimiter='\t', header=None)
        query_names = paf_out.iloc[:, 0].tolist()
        assert "contig1" in query_names, f"contig1 not found in PAF file {paf_file} queries: {set(query_names)}"

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for paf_file in glob.glob("*.paf"):
            if os.path.exists(paf_file):
                os.unlink(paf_file)
    request.addfinalizer(cleanup)
