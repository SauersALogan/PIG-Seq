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
    """Build the minimap2 command as a list - Need run_alignment to execute."""
    aligner = [
        "minimap2",
        "-x", "asm5",
        "-o", output_paf,
    ]
    if isinstance(bin_files, list):
        aligner.extend(bin_files)
    else:
        aligner.append(bin_files)
    return aligner

def run_alignment(assembly_files, bin_files):
    """Run minimap2 alignment and return PAF files"""
    try:
        aligner, output_paf = build_minimap2_command(assembly_files, bin_files)
        output_paf = subprocess.run(aligner, check=True)
        return output_path
    except subprocess.CalledProcessError:
        print(f"Minimap2 alignment failed for {assembly_file}")
        return None

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
        basename=os.path.basename(assembly)
        name=os.path.splitext(basename)[0]
        output_name=name+".paf"
        output_path=os.path.join(args.output, output_name)
        
        results = run_alignment(assembly, bin_files, output_path)

        if result:
            print(f"Successful output {output}")
        else:
            print(f"Failed for {assembly}")
            failed_alignments.append(assembly)

    if failed_alignments:
        print("Not all assemblies produced alignments")
        exit(1)
    else:
        exit(0)
