#!/usr/bin/env python3
"""
feature_parsing.py - Use the output parsed PAFs to filter sam files
==========================================

Function to filter feature counts by assembly specific contig
to bin map files generated earlier in this pipeline

This file contains the function and its tests
"""

import tempfile
import os
import pytest
import pandas as pd
import re
import shutil
import sys

# =============================================================================
# Import utility functions
# =============================================================================
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_pairing import extract_file_identifiers, pair_files_by_sample

# =============================================================================
# Actual functions to test
# =============================================================================

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def temp_test_dir():
    """Create a temp directory for test files"""
    temp_dir = tempfile.mkdtemp(prefix='Feature_parse_unit_test')
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def single_sample(temp_test_dir):
    """Creating a single sample fixture for testing"""

@pytest.fixture
def multiple_sample(temp_test_dir):
    """Create multiple sample files for unit testing."""


feature_files = []
    for i, content in enumerate([feature_1_content, feature_2_content], 1):
        feature_file = os.path.join(temp_test_dir, f'test_assembly{i}_counts.txt')
        with open(feature_file, 'w') as f:
            f.write(content)
        feature_files.append(feature_file)
    return feature_files

@pytest.fixture
def single_gff(temp_test_dir):
    """Create a mock gff file for annotation."""

@pytest.fixture
def multiple_gffs(temp_test_dir):
    """Create a multiple gff files for testing."""


    map_files = []
    for i, content in enumerate([map_1_content, map_2_content], 1):
        map_file = os.path.join(temp_test_dir, f'test_assembly{i}_map.txt')
        with open(map_file, 'w') as f:
            f.write(content)
        map_files.append(map_file)
    return map_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_assembly(single_sample, single_gff):
    """Test that the parsing function works on a single sample/assembly"""

def test_multiple_assemblies(multiple_samples, multiple_gffs):
    """Test that parsing function works on multiple samples/assemblies"""

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for sam_file in glob.glob("*.txt"):
            if os.path.exists(sam_file):
                os.unlink(sam_file)
    request.addfinalizer(cleanup)
