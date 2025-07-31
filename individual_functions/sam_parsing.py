#!/usr/bin/env python3
"""
SAM_parsing.py - Use the output parsed PAFs to filter sam files
==========================================

Function to filter and separate SAM files mapped to assemblies
using PAF files created from earlier in this pipeline

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import pandas as pd

# =============================================================================
# Actual functions to test
# =============================================================================
def SAM_parsing():
    """Filter SAM files"""

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_sam():
    """Create a single test SAM file for unit testing."""
    sam_content = """@SQ     SN:Contig1   LN:341341
@SQ     SN:Contig2   LN:325562
@SQ     SN:Contig3   LN:317386
@SQ     SN:Contig4   LN:301055
@SQ     SN:Contig5   LN:298171
@SQ     SN:Contig6   LN:282893
read1:AA 147     Contig1     3554    60      150M    =       3376 -328     GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACATCCCATGGACAGTATATGTTTGTATTACCAAACAGATAGCCGTCATACCACTCCTTTACAGCGCCGTACTGGTCCGTAAATCCATAATACGCCAGCATCT        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII        NM:i:0  MD:Z:150        MC:Z:150M       AS:i:150     XS:i:21
read2:LL 83      Contig2    72558   60      150M    =       72369-339     CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGATCAGTTCTCCTGTATATTCATCGGTAAAACGAATCTTTAAATCCTTTTCCACATAGACAAAGCTGACATTGACCTTGCCTGATGCCGTCTCCGAAGTCTC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII        NM:i:1  MD:Z:103G46     MC:Z:150M       AS:i:145     XS:i:20
read3:ZZ 163     Contig3    72369   60      150M    =       72558339      CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCACCTTATCCTTTACTGTAACCTGCTGCGGAACCGACGGGAACGTATATTTCTCTGTGGAGGTTATGATCGCATCAAACACGCCCGCGTTCACATTCTGCGC        IIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII9IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIII        NM:i:0  MD:Z:150        MC:Z:150M       AS:i:150     XS:i:0
read4:YY 99      Contig4    77025   60      150M    =       77222326      CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGACACAGTTTGAGGAGATGGGGCTATTTAACCGCGGTATCAAAACCAGAATCGGCAGCAATGGAAACTCAGATTATTATGGCATCATACGTGAGTGTTCGA        IIII-IIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9III        NM:i:0  MD:Z:150        MC:Z:129M       AS:i:150     XS:i:74
read5:XX 147     Contig5    77222   60      129M    =       77025-326     GTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGTGTTCGTGACGCACAGGCAGTCGCACGGTACTTTGGGCTTGTGTCTAAAGATAAAACAAAGGATTACAGCAATTATGCT     99-I9IIIIIIIIIIIIIIIIIIIIIIII--IIII9IIIIIIIIIIIIIIIIIIIIIIII-IIIII9IIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIII     NM:i:2  MD:Z:0A29G98    MC:Z:150M       AS:i:123        XS:i:0
read6:TT 163     Contig6    17929   60      150M    =       18129350      CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAATAAAATGCGATGATATCGCCGGTCTTGTATTTATTTCCGTTCAGCGCCACTACGATGTCTTCCTCCTGCAGCGACTCCGTCACGGAAGTTCCGCTGATC        III9IIIIIIII--IIIIIIIIIIIII9II9IIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIII9I9I9IIIII-IIII9IIIIIIIIII9II-999I99I9I--9IIIII99-IIII-II        NM:i:2  MD:Z:113A19T16  MC:Z:150M       AS:i:140     XS:i:0
"""

def multiple_pafs():
    """Create multiple SAM files for unit testing."""
    sam_1_content ="""@SQ     SN:Contig1   LN:341341
@SQ     SN:Contig2   LN:325562
@SQ     SN:Contig3   LN:317386
read1:AA 147     Contig1     3554    60      150M    =       3376 -328     GCGGTTCTGACACACTTCCGTCACGTAAATCACCGCAGTAATTTATTACA>
read2:LL 83      Contig2    72558   60      150M    =       72369-339     CTTCTTTCCTTTGGCATTCGTCAAAGTCACCTCAAAAACAGTTCCTGTGAT>
read3:ZZ 163     Contig3    72369   60      150M    =       72558339      CTGCTTCTCTGTCTTGACCTCATCTTGCACATTGATAACGACGTACTCCAC>
"""

    sam_2_content ="""@SQ     SN:Contig4   LN:301055
@SQ     SN:Contig5   LN:298171
@SQ     SN:Contig6   LN:282893
read4:YY 99      Contig4    77025   60      150M    =       77222326      CAACGGGGCAGTTGTACAGTCAGGGATATTCTGCTGCCAATGAATTTCTGA>
read5:XX 147     Contig5    77222   60      129M    =       77025-326     GTCAATGACACGCCATATATTCAGTCCCCGTAGAAGCTCAGAGAGTTTGGT>
read6:TT 163     Contig6    17929   60      150M    =       18129350      CAGTCTCCCGCATATGCAATAACCCGTTTGACCAGGATATTGTTGTTATAA>
"""

    sam_files = []
    for i, content in enumerate([sam_1_content, sam_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.paf') as tmp:
            tmp.write(content)
            sam_files.append(tmp.name)
    return sam_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_sam(single_sam):
    """Test that the parsing function works on a single sam file"""

def test_multiple_pafs(multiple_pafs):
    """Test that parsing function works on multiple sam files"""

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for sam_file in glob.glob("*.sam"):
            if os.path.exists(sam_file):
                os.unlink(sam_file)
    request.addfinalizer(cleanup)
