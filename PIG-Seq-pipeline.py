#!/usr/bin/env python3

##################################################################################
# Imports
##################################################################################

import pytest
import tempfile
import os
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import re
import shutil
import argparse
import glob

##################################################################################
# Constants
##################################################################################

##################################################################################
# Core Logic
##################################################################################
from individual_functions.contig_mapping import build_minimap2_command, run_alignment
from individual_functions.PAF_parsing import PAF_parsing, run_paf_parsing
from individual_functions.feature_counting import build_featureCounts, execute_feature_counting, run_counter
from individual_functions.feature_parsing import extract_file_identifiers, pair_files_by_sample, feature_parsing, run_parsing
from utils.file_pairing import extract_file_identifiers, pair_files_by_sample

##################################################################################
# Command line interfacing
##################################################################################
parser = argparse.ArgumentParser(description="Script for processing contigs from assemblies and splitting into binned and unbinned sets")
parser.add_argument("--assemblies", nargs="+", required=True, help="The assembly files")
parser.add_argument("--bins", nargs="+", required=True, help="The bin files")
parser.add_argument("--output", required=True, help="Output directory")
parser.add_argument("--sam_files", nargs="+", required=True, help="The directory containing SAM files")
parser.add_argument("--gtf_files", nargs="+", required=True, help="The directory containing GTF files")
parser.add_argument("--pattern_source", required=False, help="File or string for how identifiers should be handled")

##################################################################################
# Main workflow
##################################################################################
if __name__ == "__main__":
    args = parser.parse_args()

    # Expand the wildcard
    assembly_files = []
    for pattern in args.assemblies:
        assembly_files.extend(glob.glob(pattern))

    bin_files = []
    for pattern in args.bins:
        if os.path.isdir(pattern):
            bin_files.extend(glob.glob(os.path.join(pattern, "*.fasta")))
            bin_files.extend(glob.glob(os.path.join(pattern, "*.fa")))
            bin_files.extend(glob.glob(os.path.join(pattern, "*.fas")))
        else:
            expanded = glob.glob(pattern)
            if expanded:
                bin_files.extend(expanded)
            else:
                print(f"Warning: No files found for pattern: {pattern}")

    sam_files = []
    for pattern in args.sam_files:
        sam_files.extend(glob.glob(pattern))
    sam_files.sort()

    gtf_files = []
    for pattern in args.gtf_files:
        gtf_files.extend(glob.glob(pattern))
    gtf_files.sort()

    # Load the pattern_source
    pattern_source = args.pattern_source

    # Display the files
    if assembly_files:
        print("Assembly files:")
        for file in assembly_files:
            print(f" - {file}")
    else:
        print("No assembly files found")

    if bin_files:
        print("Bin files:")
        for file in bin_files:
            print(f" - {file}")
    else:
        print("No bin files found")

    if sam_files:
        print("Sam files:")
        for file in sam_files:
            print(f" - {file}")
    else:
        print("No sam files found")

    if gtf_files:
        print("GTF files:")
        for file in gtf_files:
            print(f" - {file}")
    else:
        print("No GTF files found")

    if pattern_source is None:
        print(f"You have not selected a method for finding file identifiers, proceeding with default logic")
    elif os.path.isfile(pattern_source):
        print(f"You have selected to use an input file for indentifiers, using file: {pattern_source}")
    else:
        print(f"You've input a custom string for identifiers, I will search using: {pattern_source}")

    os.makedirs(args.output, exist_ok=True)

    print("\n=== Running alignment ===")
    aligned_paf_files = run_alignment(bin_files, assembly_files)
    for paf_file in aligned_paf_files:
        print(f" - {paf_file:3}")
        print(f" - ... and {len(paf_file) - 3} additional lines")

    print("\n=== Running PAF parsing ===")
    mapping_files = run_paf_parsing(aligned_paf_files, assembly_files, bin_files)
    print(f"Created mapping files: {mapping_files}")
    for mapping_file in mapping_files:
        print(f"\n=== First 5 lines of {mapping_file} ===")
        try:
            with open(mapping_file, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 5:  # Stop after 5 lines
                        break
                    print(line.strip())
        except FileNotFoundError:
            print(f"Warning: File {mapping_file} not found")
        except Exception as e:
            print(f"Error reading {mapping_file}: {e}")

    print("\n=== Running feature counting ===")
    results = run_counter(sam_files, gtf_files, threads=2, gene="CDS", gene_id="locus_tag")
    print(f"Feature counting - Type: {type(results)}, Count: {len(results) if results else 0}")
    print(f" - Files contain: {results[:2]} and ... (+{len(results)-2} additional lines)")

    print("\n=== Preparing expected outputs ===")
    expected_feature_outputs = []
    for count_file in results:
        count_base = os.path.basename(count_file)
        count_name = os.path.splitext(count_base)[0]
        feature_output = count_name + "_binned.txt"
        expected_feature_outputs.append(feature_parsing_results)

    print("\n=== Running feature parsing ===")
    feature_parsing_results = run_parsing(mapping_files, results)
    print(f"Feature parsing - Type: {type(feature_parsing_results)}, Result: {feature_parsing_results}")

