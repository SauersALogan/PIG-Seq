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

def build_minimap2_command(assembly_files, bin_files):
    """"Build the minimap2 command as a list - Need run_alignment to execute."""

    cmd = [ 
        "minimap2",
        "-x", "asm5", 
        "-o", output_paf, 
    ]

    if isinstance(bin_files, list): 
        cmd.extend(bin_files) 
    else: 
        cmd.append(bin_files)

    return cmd, output_paf

def run_alignment(assembly_files, bin_files):
    """Run minimap2 alignment and return PAF files"""
    try:
        cmd, output_paf = build_minimap2_command(assembly_files, bin_files)
        results = subprocess.run(cmd, check=True)
        return output_paf
    except subprocess.CalledProcessError:
        print("Minimap2 alignment failed")
        return None

##################################################################################
# Command line interfacing
##################################################################################

parser = argparse.ArgumentParser(description="Script for processing contigs from assemblies and splitting into binned and unbinned sets")
parser.add_argument("--assemblies", nargs="+", required=True, help="The assembly files")
parser.add_argument("--bins", nargs="+", required=True, help="The bin files")
parser.add_argument("--output", required=True, help="Output path file")

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

    results = run_alignment(assembly_files, bin_files, args.output)

    if results:
        print(f"Processing completed successfuly. Output: {results}@
    else:
        print("Processing failed")
    exit(1)
