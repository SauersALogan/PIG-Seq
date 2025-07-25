#!/usr/bin/env python3

##################################################################################
# Imports
##################################################################################

import tempfile
import os
import pytest
import subprocess
import argparse
import glob

##################################################################################
# Constants
##################################################################################



##################################################################################
# Core Logic
##################################################################################

from individual_functions.contig_mapping import build_minimap2_command, run_alignment

##################################################################################
# Command line interfacing
##################################################################################

parser = argparse.ArgumentParser(description="Script for processing contigs from assemblies and splitting into binned and unbinned sets")
parser.add_argument("--assemblies", nargs="+", required=True, help="The assembly files")
parser.add_argument("--bins", nargs="+", required=True, help="The bin files")
parser.add_argument("--output", required=True, help="Output directory")

##################################################################################
# Main workflow
##################################################################################

if __name__ == "__main__":
    args = parser.parse_args()

    # Expand the wildcard
    assembly_files=[]
    for pattern in args.assemblies:
        assembly_files.extend(glob.glob(pattern))

    bin_files=[]
    for pattern in args.bins:
        bin_files.extend(glob.glob(pattern))

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

    os.makedirs(args.output, exists_ok=True)

    failed_alignments=[]

    for assembly in assembly_files:
        for bin in bin_files:
            assembly_base=os.path.basename(assembly)
            assembly_name=os.path.splitext(assembly_base)[0]
            bin_base=os.path.basename(bin)
            bin_name=os.path.splitext(bin_base)[0]
            output_name=assembly_name+"_"+bin_name+".paf"
            output_path=(output_name)
            results = run_alignment(assembly, bin_files, output_path)

            if result:
                print(f"Successful output {output}")
            else:
                print(f"Failed for {assembly}")
                failed_alignments.append(assembly)
