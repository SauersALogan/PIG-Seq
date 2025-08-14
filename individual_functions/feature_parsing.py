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

# =============================================================================
# Actual functions to test
# =============================================================================
def extract_file_identifiers(filename):
    """Extracts the similar file identifier for all input files to ensure proper pairing.
    First it looks for common patterns such as Sample_XXX, then numerical suffixes such as sample1. If these fail
    it will finish by looking at the last part of the file after controlling for common delimiters"""
    basename = os.path.basename(filename)
    basename = re.sub(r'\.(txt|sam|bam|gtf|gff|paf|fasta|fa|fastq|fq)$', '', basename, flags=re.IGNORECASE)
    match = re.search(r'(\d+)', basename)
    if match:
        return match.group(1)
    return basename

def pair_files_by_sample(feature_files, map_files):
    """This function will pair files based upon the layouts controlled for in the function above"""
    feature_dict = {}
    for feature_file in feature_files:
        identifier = extract_file_identifiers(feature_file)
        feature_dict[identifier] = feature_file
    map_dict = {}
    for map_file in map_files:
        identifier = extract_file_identifiers(map_file)
        map_dict[identifier] = map_file
    paired_files = []
    unmatched_features = []
    unmatched_maps = []
    for identifier in feature_dict.keys():
        if identifier in map_dict:
            feature_file = feature_dict[identifier]
            map_file = map_dict[identifier]
            paired_files.append((feature_file, map_file))
            print(f"DEBUG: Paired {os.path.basename(feature_file)} with {os.path.basename(map_file)} (ID: {identifier})")
        else:
            unmatched_features.append(feature_dict[identifier])
    for identifier in map_dict.keys():
        if identifier not in feature_dict:
            unmatched_maps.append(map_dict[identifier])
    if unmatched_features:
        print(f"WARNING: Unmatched feature files: {[os.path.basename(f) for f in unmatched_features]}")
    if unmatched_maps:
        print(f"WARNING: Unmatched mapping files: {[os.path.basename(f) for f in unmatched_maps]}")
    return paired_files

def feature_parsing(map_file, feature_file):
    """Filter feature file"""
    feature_base=os.path.basename(feature_file)
    feature_name=os.path.splitext(feature_base)[0]
    output_name=feature_name+"_binned.txt"
    output_path=(output_name)
    print(f"Writing output to {output_path}")
    features = pd.read_csv(feature_file, delimiter='\t', comment='#')
    print(f"Found feature file to contain the following columns:")
    f_headers = list(features.columns.values)
    print(f"{f_headers}")
    map = pd.read_csv(map_file, delimiter='\t', header=0)
    print(f"Found mapping file with the following headers")
    headers = list(map.columns.values)
    print(f"{headers}")
    map_dictionary = map.set_index('Contig')['Bin'].to_dict()
    print("Created the following mapping")
    print(f"{map_dictionary}")
    results = features['Chr'].map(map_dictionary)
    features['binning'] = results
    features.to_csv(f"{output_path}", sep="\t", header=True, index = False)
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
        paired_files = pair_files_by_sample(feature_files, map_files)
        if not paired_files:
            print(f"DEBUG: No matching file pairs found, please check your file naming conventions or check the documentation")
            return final_results
        for feature_file, map_file in paired_files:
            parsed_features = feature_parsing(map_file, feature_file)
            final_results.append(parsed_features)
    return(final_results)

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
def single_feature(temp_test_dir):
    """Create a single test feature file for unit testing."""
    feature_content = """Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample1.bam
ENSG00000223972\tcontig1\t11869\t14409\t+\t1735\t23
ENSG00000227232\tcontig1\t14404\t29570\t-\t2073\t156
ENSG00000278267\tcontig2\t17369\t17436\t-\t68\t0
ENSG00000243485\tcontig4\t29554\t31109\t+\t1440\t89
ENSG00000284332\tcontig5\t30366\t30503\t+\t138\t5
ENSG00000237613\tcontig5\t34554\t36081\t-\t718\t34
ENSG00000268020\tcontig6\t52473\t53312\t+\t840\t67
ENSG00000240361\tcontig6\t62948\t63887\t+\t940\t12
ENSG00000186092\tcontig6\t69091\t70008\t+\t918\t203
"""

    feature_file = os.path.join(temp_test_dir, 'test_assembly1_counts.txt')
    with open(feature_file, 'w') as f:
        f.write(feature_content)
    return [feature_file]

@pytest.fixture
def multiple_features(temp_test_dir):
    """Create multiple feature files for unit testing."""
    feature_1_content ="""Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample1.bam
ENSG00000223972\tcontig1\t11869\t14409\t+\t1735\t23
ENSG00000227232\tcontig1\t14404\t29570\t-\t2073\t156
ENSG00000278267\tcontig2\t17369\t17436\t-\t68\t0
ENSG00000243485\tcontig4\t29554\t31109\t+\t1440\t89
ENSG00000284332\tcontig5\t30366\t30503\t+\t138\t5
ENSG00000237613\tcontig5\t34554\t36081\t-\t718\t34
ENSG00000268020\tcontig6\t52473\t53312\t+\t840\t67
ENSG00000240361\tcontig6\t62948\t63887\t+\t940\t12
ENSG00000186092\tcontig6\t69091\t70008\t+\t918\t203
"""

    feature_2_content ="""Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample1.bam
ENSG00000223972\tcontig1\t11869\t14409\t+\t1735\t230
ENSG00000227232\tcontig1\t14404\t29570\t-\t2073\t15
ENSG00000278267\tcontig2\t17369\t17436\t-\t68\t21
ENSG00000243485\tcontig4\t29554\t31109\t+\t1440\t5
ENSG00000284332\tcontig5\t30366\t30503\t+\t138\t5
ENSG00000237613\tcontig5\t34554\t36081\t-\t718\t54
ENSG00000268020\tcontig6\t52473\t53312\t+\t840\t121
ENSG00000240361\tcontig6\t62948\t63887\t+\t940\t1
ENSG00000186092\tcontig6\t69091\t70008\t+\t918\t2030
"""

    feature_files = []
    for i, content in enumerate([feature_1_content, feature_2_content], 1):
        feature_file = os.path.join(temp_test_dir, f'test_assembly{i}_counts.txt')
        with open(feature_file, 'w') as f:
            f.write(content)
        feature_files.append(feature_file)
    return feature_files

@pytest.fixture
def single_map(temp_test_dir):
    """Create a mock contig to binning map file."""
    map_content = "Contig\tBin\ncontig1\tbin1\ncontig3\tbin3\ncontig2\tunbinned\ncontig4\tunbinned\ncontig5\tbin2\ncontig6\tbin1"
    map_file = os.path.join(temp_test_dir, 'test_assembly1_map.txt')
    with open(map_file, 'w') as f:
        f.write(map_content)
    return [map_file]

@pytest.fixture
def multiple_maps(temp_test_dir):
    """Create a multiple mock contig to binning map file."""
    map_1_content = "Contig\tBin\ncontig1\tbin1\ncontig3\tbin3\ncontig2\tunbinned\ncontig4\tunbinned\ncontig5\tbin2\ncontig6\tbin1"

    map_2_content = "Contig\tBin\ncontig1\tbin1\ncontig3\tbin3\ncontig2\tunbinned\ncontig4\tunbinned\ncontig5\tbin2\ncontig6\tbin1"

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
def test_single_assembly(single_map, single_feature):
    """Test that the parsing function works on a single sample/assembly"""
    run_parsing(single_map, single_feature)
    for input_file in single_feature:
        feature_base = os.path.basename(input_file)
        feature_name = os.path.splitext(feature_base)[0]
        output_name = feature_name + "_binned.txt"
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
        output_name = feature_name + "_binned.txt"
        output_path = output_name
        if os.path.exists(output_path):
            print(f"DEBUG: Output file {output_path} was created successfully")
            file_size = os.path.getsize(output_path)
            assert file_size > 0, f"Output file {output_path} is empty!"
            print(f"DEBUG: File size: {file_size} bytes")
            df = pd.read_csv(output_path, delimiter='\t')
            print(f"DEBUG: Content in output file is:")
            print(f"{df}")
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
