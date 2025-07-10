#!/usr/bin/env python3
"""
count_sequences.py - Count FASTA Sequences
==========================================

Function to count sequences in FASTA files.
This file contains the function and its tests.
"""

import tempfile
import os

def count_fasta_sequences(filename):
    """
    Count how many lines start with '>' in a file.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        int: Number of lines starting with '>'
        
    Example:
        >>> with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
        ...     tmp.write(">seq1\\nATCG\\n>seq2\\nGCTA\\n")
        ...     tmp_path = tmp.name
        >>> count_fasta_sequences(tmp_path)
        2
        >>> os.unlink(tmp_path)
    """
    count = 0
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count

# =============================================================================
# TESTS FOR THIS FUNCTION ONLY
# =============================================================================

def test_count_three_sequences():
    """Test counting 3 sequences."""
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
        result = count_fasta_sequences(tmp_path)
        assert result == 3, f"Expected 3 sequences, got {result}"
        print("✅ test_count_three_sequences passed")
    finally:
        os.unlink(tmp_path)

def test_count_zero_sequences():
    """Test with no sequences (no > lines)."""
    content = """Just some text
No sequences here
ATCGATCG"""
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    
    try:
        result = count_fasta_sequences(tmp_path)
        assert result == 0, f"Expected 0 sequences, got {result}"
        print("✅ test_count_zero_sequences passed")
    finally:
        os.unlink(tmp_path)

def test_count_one_sequence():
    """Test with single sequence."""
    content = """>single_sequence
ATCGATCGATCGATCG"""
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    
    try:
        result = count_fasta_sequences(tmp_path)
        assert result == 1, f"Expected 1 sequence, got {result}"
        print("✅ test_count_one_sequence passed")
    finally:
        os.unlink(tmp_path)

def test_count_multiline_sequences():
    """Test sequences that span multiple lines."""
    content = """>seq1
ATCGATCGATCG
GCTAGCTAGCTA
>seq2
TTTTAAAAAAA
CCCCGGGGGGGG
AAAATTTTCCCC"""
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    
    try:
        result = count_fasta_sequences(tmp_path)
        assert result == 2, f"Expected 2 sequences, got {result}"
        print("✅ test_count_multiline_sequences passed")
    finally:
        os.unlink(tmp_path)

def test_count_empty_file():
    """Test with completely empty file."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
        tmp.write("")  # Empty file
        tmp_path = tmp.name
    
    try:
        result = count_fasta_sequences(tmp_path)
        assert result == 0, f"Expected 0 sequences in empty file, got {result}"
        print("✅ test_count_empty_file passed")
    finally:
        os.unlink(tmp_path)

def run_all_tests():
    """Run all tests for this function."""
    print("Running tests for count_fasta_sequences...")
    print("=" * 50)
    
    tests = [
        test_count_three_sequences,
        test_count_zero_sequences,
        test_count_one_sequence,
        test_count_multiline_sequences,
        test_count_empty_file
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"❌ {test_func.__name__} failed: {e}")
            failed += 1
    
    print("=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All tests passed! Function is ready for integration.")
        return True
    else:
        print("💥 Some tests failed. Fix the function before integration.")
        return False

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
