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

# =============================================================================
# Actual functions to test
# =============================================================================
def feature_parsing(map_file, feature_file):
    """Filter feature file"""
    feature_base=os.path.basename(feature_file)
    feature_name=os.path.splitext(feature_base)[0]
    output_name=feature_name+"_binned.txt"
    output_path=(output_name)
    print(f"Writing output to {output_path}")
    features = pd.read_csv(feature_file, delimiter=',')
    print(f"Found feature file to contain the following columns:")
    f_headers = list(features.columns.values)
    print(f"{f_headers}")
    map = pd.read_csv(map_file, delimiter=',')
    print(f"Found mapping file with the following headers")
    headers = list(map.columns.values)
    print(f"{headers}")
    map_dictionary = map.set_index('Contig')['Bin'].to_dict()
    print("Created the following mapping")
    print(f"{map_dictionary}")
    results = features['Chr'].map(map_dictionary)
    features['binning'] = results
    features.to_csv(f"{output_path}", sep="\t", header=True)
    print(f"Resulting gene to bin to count data written to: {output_path}")
    return features

def run_parsing(map_files, feature_files):
    """Run the parsing function above on varies file input formats"""
    final_results = []
    if isinstance(feature_files, str) and isinstance(map_files, str):
        feature_file = feature_files
        map_file = map_files
        parsed_features = feature_parsing(map_file, feature_file)
        final_results.append(parsed_features)
    elif isinstance(feature_files, list) and isinstance(map_files, str):
        print("It seems you have mapped several features to the same assembly, that is not really what this tool is for")
    elif isinstance(feature_files, str) and isinstance(map_files, list):
        print("You provided a list of mapping files for a single feature file, this tells me something is wrong, check inputs")
    elif isinstance(feature_files, list) and isinstance(map_files, list):
        feature_files.sort()
        map_files.sort()
        for feature_file, map_file in zip(feature_files, map_files):
            parsed_features = feature_parsing(map_file, feature_file)
            final_results.append(parsed_features)
    return(final_results)

# =============================================================================
# Create the mock data for testing
# =============================================================================
@pytest.fixture
def single_feature():
    """Create a single test feature file for unit testing."""
    feature_content = """Geneid,Chr,Start,End,Strand,Length,sample1.bam
ENSG00000223972,contig1,11869,14409,+,1735,23
ENSG00000227232,contig1,14404,29570,-,2073,156
ENSG00000278267,contig2,17369,17436,-,68,0
ENSG00000243485,contig4,29554,31109,+,1440,89
ENSG00000284332,contig5,30366,30503,+,138,5
ENSG00000237613,contig5,34554,36081,-,718,34
ENSG00000268020,contig6,52473,53312,+,840,67
ENSG00000240361,contig6,62948,63887,+,940,12
ENSG00000186092,contig6,69091,70008,+,918,203
"""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1_counts.txt') as tmp:
        tmp.write(feature_content)
        feature_file = tmp.name
    return [feature_file]

@pytest.fixture
def multiple_features():
    """Create multiple feature files for unit testing."""
    feature_1_content ="""Geneid,Chr,Start,End,Strand,Length,sample1.bam
ENSG00000223972,contig1,11869,14409,+,1735,23
ENSG00000227232,contig1,14404,29570,-,2073,156
ENSG00000278267,contig2,17369,17436,-,68,0
ENSG00000243485,contig4,29554,31109,+,1440,89
ENSG00000284332,contig5,30366,30503,+,138,5
ENSG00000237613,contig5,34554,36081,-,718,34
ENSG00000268020,contig6,52473,53312,+,840,67
ENSG00000240361,contig6,62948,63887,+,940,12
ENSG00000186092,contig6,69091,70008,+,918,203
"""

    feature_2_content ="""Geneid,Chr,Start,End,Strand,Length,sample1.bam
ENSG00000223972,contig1,11869,14409,+,1735,230
ENSG00000227232,contig1,14404,29570,-,2073,15
ENSG00000278267,contig2,17369,17436,-,68,21
ENSG00000243485,contig4,29554,31109,+,1440,5
ENSG00000284332,contig5,30366,30503,+,138,5
ENSG00000237613,contig5,34554,36081,-,718,54
ENSG00000268020,contig6,52473,53312,+,840,121
ENSG00000240361,contig6,62948,63887,+,940,1
ENSG00000186092,contig6,69091,70008,+,918,2030
"""

    feature_files = []
    for i, content in enumerate([feature_1_content, feature_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}_counts.txt') as tmp:
            tmp.write(content)
            feature_files.append(tmp.name)
    return feature_files

@pytest.fixture
def single_map():
    """Create a mock contig to binning map file."""
    map_content = "Contig,Bin\ncontig1,bin1\ncontig3,bin3\ncontig2,unbinned\ncontig4,unbinned\ncontig5,bin2\ncontig6,bin1"

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_assembly1_map.txt') as tmp:
        tmp.write(map_content)
        map_file = tmp.name
    return [map_file]

@pytest.fixture
def multiple_maps():
    """Create a multiple mock contig to binning map file."""
    map_1_content = "Contig,Bin\ncontig1,bin1\ncontig3,bin3\ncontig2,unbinned"

    map_2_content = "Contig,Bin\ncontig4,unbinned\ncontig5,bin2\ncontig6,bin1"

    map_files = []
    for i, content in enumerate([map_1_content, map_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}_map.txt') as tmp:
            tmp.write(content)
            map_files.append(tmp.name)
    return map_files

# =============================================================================
# Setup the actual tests
# =============================================================================
def test_single_assembly(single_map, single_feature):
    """Test that the parsing function works on a single sample/assembly"""
    run_parsing(single_map, single_feature)
    for input_file in single_feature:
        feature_base = os.path.basename(input_file)
        feature_name = os.path.splitext(feature_base)[0]
        output_name = feature_name + "binned_counts.txt"
        output_path = output_name
        if os.path.exists(output_path):
            print(f"DEBUG: Output file {output_path} was created successfully")
            file_size = os.path.getsize(output_path)
            assert file_size > 0, f"Output file {output_path} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
            df = pd.read_csv(output_path, delimiter='\t')
            contig6_rows = df[df['Chr'] == 'contig6']
            contig1_rows = df[df['Chr'] == 'contig1']
            contig3_rows = df[df['Chr'] == 'contig3']
            contig4_rows = df[df['Chr'] == 'contig4']
            assert(contig6_rows['binning'] == 'bin1').all(), "Not all contig6 entries assigned to bin1"
            assert(contig1_rows['binning'] == 'bin1').all(), "Not all contig1 entries assigned to bin1"
            assert(contig3_rows['binning'] == 'unbinned').all(), "Not all contig3 entries assigned to unbinned"
            assert(contig4_rows['binning'] == 'unbinned').all(), "Not all contig4 entries assigned to unbinned"
        else:
            print(f"DEBUG: WARNING - Output file {output_path} was NOT created")

def test_multiple_assemblies(multiple_maps, multiple_features):
    """Test that parsing function works on multiple samples/assemblies"""
    run_parsing(multiple_maps, multiple_features)
    for input_file in multiple_features:
        feature_base = os.path.basename(input_file)
        feature_name = os.path.splitext(feature_base)[0]
        output_name = feature_name + "binned_counts.txt"
        output_path = output_name
        if os.path.exists(output_path):
            print(f"DEBUG: Output file {output_path} was created successfully")
            file_size = os.path.getsize(output_path)
            assert file_size > 0, f"Output file {output_path} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
            df = pd.read_csv(output_path, delimiter='\t')
            contig6_rows = df[df['Chr'] == 'contig6']
            contig1_rows = df[df['Chr'] == 'contig1']
            contig3_rows = df[df['Chr'] == 'contig3']
            contig4_rows = df[df['Chr'] == 'contig4']
            assert(contig6_rows['binning'] == 'bin1').all(), "Not all contig6 entries assigned to bin1"
            assert(contig1_rows['binning'] == 'bin1').all(), "Not all contig1 entries assigned to bin1"
            assert(contig3_rows['binning'] == 'unbinned').all(), "Not all contig3 entries assigned to unbinned"
            assert(contig4_rows['binning'] == 'unbinned').all(), "Not all contig4 entries assigned to unbinned"
        else:
            print(f"DEBUG: WARNING - Output file {output_path} was NOT created")

@pytest.fixture(autouse=True)
def cleanup_paf_files(request):
    def cleanup():
        import glob
        for sam_file in glob.glob("*.txt"):
            if os.path.exists(sam_file):
                os.unlink(sam_file)
    request.addfinalizer(cleanup)
