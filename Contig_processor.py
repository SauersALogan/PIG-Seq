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
from pathlib import Path
from unittest.mock import patch, MagicMock
import subprocess

# =============================================================================
# FUNCTION 1: Count Lines Starting with '>'
# =============================================================================

def count_fasta_sequences(filename):
    """
    Count how many lines start with '>' in a file.
    
    Args:
        filename (str): Path to file
        
    Returns:
        int: Number of lines starting with '>'
    """
    # TODO: Write your implementation here
    pass

# SIMPLE TESTS FOR FUNCTION 1
class TestCountSequences:
    """Simple tests - just the happy path first!"""
    
    def test_count_three_sequences(self):
        """Test counting 3 sequences."""
        # Create a simple test file
        content = """
>seq1
ATCG
>seq2  
GCTA
>seq3
TTTT
"""
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            # Test our function
            result = count_fasta_sequences(tmp_path)
            assert result == 3
        finally:
            os.unlink(tmp_path)  # Clean up
    
    def test_count_zero_sequences(self):
        """Test with no sequences (no > lines)."""
        content = """
Just some text
No sequences here
ATCGATCG"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            result = count_fasta_sequences(tmp_path)
            assert result == 0
        finally:
            os.unlink(tmp_path)

# =============================================================================
# FUNCTION 2: Extract Sequence Names
# =============================================================================

def get_sequence_names(filename):
    """
    Get all sequence names from a FASTA file.
    Sequence names are the parts after '>' but before any spaces.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        list: List of sequence names
    """
    # TODO: Write your implementation here
    pass

class TestGetSequenceNames:
    """Test extracting sequence names."""
    
    def test_get_simple_names(self):
        """Test getting simple sequence names."""
        content = """>seq1
ATCG
>seq2
GCTA
>seq3
TTTT"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            result = get_sequence_names(tmp_path)
            expected = ['seq1', 'seq2', 'seq3']
            assert result == expected
        finally:
            os.unlink(tmp_path)
    
    def test_get_names_with_descriptions(self):
        """Test names with descriptions (ignore part after space)."""
        content = """>contig1 length=1000 coverage=10.5
ATCG
>contig2 length=500
GCTA"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            result = get_sequence_names(tmp_path)
            expected = ['contig1', 'contig2']  # Just the name part
            assert result == expected
        finally:
            os.unlink(tmp_path)

# =============================================================================
# FUNCTION 3: Simple PAF Parser
# =============================================================================

def parse_simple_paf(filename):
    """
    Parse a PAF file and return basic alignment info.
    PAF format has these columns (we only care about first few):
    0: query_name
    1: query_length  
    9: matches
    10: alignment_length
    
    Args:
        filename (str): Path to PAF file
        
    Returns:
        list: List of dictionaries with alignment info
    """
    # TODO: Write your implementation here
    pass

class TestParsePAF:
    """Test PAF parsing - keep it simple!"""
    
    def test_parse_one_alignment(self):
        """Test parsing one alignment."""
        # Simple PAF line (tab-separated)
        content = "contig1\t1000\t0\t950\t+\tbin1\t2000\t100\t1050\t900\t950\t60"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.paf') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            result = parse_simple_paf(tmp_path)
            expected = [{
                'query_name': 'contig1',
                'query_length': 1000,
                'matches': 900,
                'alignment_length': 950
            }]
            assert result == expected
        finally:
            os.unlink(tmp_path)
    
    def test_parse_multiple_alignments(self):
        """Test parsing multiple alignments."""
        content = """contig1\t1000\t0\t950\t+\tbin1\t2000\t100\t1050\t900\t950\t60
contig2\t800\t0\t400\t+\tbin2\t1500\t200\t600\t350\t400\t40"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.paf') as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            result = parse_simple_paf(tmp_path)
            assert len(result) == 2
            assert result[0]['query_name'] == 'contig1'
            assert result[1]['query_name'] == 'contig2'
        finally:
            os.unlink(tmp_path)

# =============================================================================
# FUNCTION 4: Build Command Lists
# =============================================================================

def build_minimap_command(reference, query, output):
    """
    Build a minimap2 command as a list (for subprocess).
    
    Args:
        reference (str): Reference file path
        query (str): Query file path  
        output (str): Output file path
        
    Returns:
        list: Command as list of strings
    """
    # TODO: Write your implementation here
    pass

class TestBuildCommands:
    """Test command building - no actual execution!"""
    
    def test_build_minimap_command(self):
        """Test building a minimap2 command."""
        result = build_minimap_command("ref.fa", "query.fa", "output.paf")
        expected = ["minimap2", "-x", "asm5", "-c", "ref.fa", "query.fa"]
        assert result == expected
    
    def test_build_minimap_with_threads(self):
        """Test building command with thread option."""
        # We'll add this later - start simple!
        pass

