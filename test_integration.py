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
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

# =============================================================================
# Integration test fixtures - Shared test data
# =============================================================================
def generate_random_sequence(length, gc_content=0.5):
    bases=['A', 'T', 'G', 'C']
    weights = [(1-gc_content)/2, (1-gc_content)/2, gc_content/2, gc_content/2]
    sequence = ''.join(random.choices(bases, weights=weights, k=length))
    return sequence

# Generate the sequences with a consistent seed
random.seed(42)
perfect_match = generate_random_sequence(64000) # 64kb
bad_match = generate_random_sequence(64000, gc_content=0.52)
random.seed(101010)
bad_coverage_start = generate_random_sequence(8000)
bad_coverage_middle = perfect_match[16000:48000]
random.seed(20547391)
bad_coverage_end = generate_random_sequence(8000)
bad_coverage = bad_coverage_start + bad_coverage_middle + bad_coverage_end
bad_quality_bases = list(perfect_match)
for i in range(0, len(bad_quality_bases), 50):
    if i < len(bad_quality_bases):
        bad_quality_bases[i] = 'N'
        if i + 1 < len(bad_quality_bases):
            bad_quality_bases[i + 1] = 'N'
bad_quality = "".join(bad_quality_bases)

@pytest.fixture
def multiple_assemblies():
    """Create assembly files where contig1 will match bin sequences well, others poorly."""
    assembly_1_content = f">a1_contig1\n{perfect_match}\n>a1_contig2\n{bad_match}\n>a1_contig3\n{bad_coverage}\n>a1_contig4\n{bad_quality}"
    assembly_2_content = f">a2_contig1\n{bad_match}\n>a2_contig2\n{perfect_match}\n>a2_contig3\n{bad_coverage}\n>a2_contig4\n{bad_quality}"

    assembly_files = []
    for i, content in enumerate([assembly_1_content, assembly_2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.fasta') as tmp:
            tmp.write(content)
            assembly_files.append(tmp.name)
    return assembly_files

@pytest.fixture
def sample_bins():
    """Create bin files where one sequence matches contig1 perfectly, others don't match well."""
    bin1_content = f">bin1_scaffold1\n{perfect_match}\n>bin1_scaffold2\n{generate_random_sequence(3000)}"
    bin2_content = f">bin2_scaffold1\n{generate_random_sequence(5000)}\n>bin2_scaffold2\n{perfect_match[:20000]}"

    bin_files = []
    for i, content in enumerate([bin1_content, bin2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_bin{i}.fasta') as tmp:
            tmp.write(content)
            bin_files.append(tmp.name)
    return bin_files

@pytest.fixture
def mock_gtf_files(multiple_assemblies):
    """Create GTF files corresponding to the assembly files."""
    gtf_files = []

    # GTF content for assembly 1 (a1_contig1, a1_contig2, a1_contig3, a1_contig4)
    gtf1_content = """a1_contig1	prokka	CDS	100	500	.	+	0	gene_id "GENE001"; product "hypothetical protein"
a1_contig1	prokka	CDS	600	1200	.	-	0	gene_id "GENE002"; product "DNA polymerase"
a1_contig1	prokka	CDS	1500	2000	.	+	0	gene_id "GENE003"; product "ribosomal protein"
a1_contig2	prokka	CDS	200	800	.	+	0	gene_id "GENE004"; product "membrane protein"
a1_contig2	prokka	CDS	1000	1600	.	-	0	gene_id "GENE005"; product "ATP synthase"
a1_contig3	prokka	CDS	300	900	.	+	0	gene_id "GENE006"; product "transcriptase"
a1_contig4	prokka	CDS	150	750	.	+	0	gene_id "GENE007"; product "helicase"
"""

    # GTF content for assembly 2 (a2_contig1, a2_contig2, a2_contig3, a2_contig4)
    gtf2_content = """a2_contig1	prokka	CDS	100	500	.	+	0	gene_id "GENE008"; product "hypothetical protein"
a2_contig1	prokka	CDS	600	1200	.	-	0	gene_id "GENE009"; product "DNA polymerase"
a2_contig2	prokka	CDS	200	800	.	+	0	gene_id "GENE010"; product "membrane protein"
a2_contig2	prokka	CDS	1000	1600	.	-	0	gene_id "GENE011"; product "ATP synthase"
a2_contig2	prokka	CDS	1800	2400	.	+	0	gene_id "GENE012"; product "kinase"
a2_contig3	prokka	CDS	300	900	.	+	0	gene_id "GENE013"; product "transcriptase"
a2_contig4	prokka	CDS	150	750	.	+	0	gene_id "GENE014"; product "helicase"
"""

    for i, content in enumerate([gtf1_content, gtf2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_assembly{i}.gtf') as tmp:
            tmp.write(content)
            gtf_files.append(tmp.name)
    return gtf_files

@pytest.fixture
def mock_sam_files(multiple_assemblies):
    """Create SAM files with reads mapping to the assembly contigs."""
    sam_files = []

    # SAM header and content for assembly 1
    sam1_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a1_contig1	LN:64000
@SQ	SN:a1_contig2	LN:64000
@SQ	SN:a1_contig3	LN:48000
@SQ	SN:a1_contig4	LN:64000
read1	0	a1_contig1	150	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	0	a1_contig1	300	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	0	a1_contig1	700	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	0	a1_contig2	250	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	0	a1_contig3	400	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam2_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a2_contig1	LN:64000
@SQ	SN:a2_contig2	LN:64000
@SQ	SN:a2_contig3	LN:48000
@SQ	SN:a2_contig4	LN:64000
read6	0	a2_contig1	150	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	0	a2_contig2	200	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	0	a2_contig2	1100	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	0	a2_contig2	1900	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read10	0	a2_contig3	350	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    for i, content in enumerate([sam1_content, sam2_content], 1):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_sample{i}.sam') as tmp:
            tmp.write(content)
            sam_files.append(tmp.name)
    return sam_files

# =============================================================================
# Import the functions from the individual_functions folder
# =============================================================================
from individual_functions.contig_mapping import build_minimap2_command, run_alignment
from individual_functions.PAF_parsing import PAF_parsing, run_paf_parsing
from individual_functions.feature_counting.py build_featureCounts, execute_feature_counting, run_counter
from individual_functions.feature_parsing.py feature_parsing, run_parsing

# =============================================================================
# Integration test - Test functions working together
# =============================================================================
class TestFullPipeline:
    def test_complete_pipeline(self, sample_bins, multiple_assemblies):
        """
        Test to check if we're ready for integration testing.
        This test checks if individual functions can be imported.
        """
        # TODO: Uncomment as individual functions are completed

        # Step 1: Align
        aligned_paf_files = run_alignment(sample_bins, multiple_assemblies)
        for paf_file in aligned_paf_files:
            paf_out = pd.read_csv(paf_file, delimiter='\t', header=None)
            query_names = paf_out.iloc[:, 0].tolist()
            mapped_reads = ["a1_contig1", "a2_contig2"]
            found_reads = set(mapped_reads) & set(query_names)
            assert found_reads, f"None of {mapped_reads} found in PAF file {paf_file} queries: {set(query_names)}"

        # Step 2: Parse the PAF output
        binned_reads = ["a1_contig1", "a2_contig2"]
        unbinned_reads = ["a1_contig2", "a1_contig3", "a1_contig4", "a2_contig1", "a2_contig3", "a2_contig4"]
        for paf_file in aligned_paf_files:
            run_paf_parsing(aligned_paf_files)
        for input_file in aligned_paf_files:
            expected_output = os.path.splitext(input_file)[0] + "_contigs_to_bin_mapping.txt"
            output_file = pd.read_csv(expected_output, delimiter='\t', header=0)
            read_names = output_file.iloc[:, 0].tolist()
            print(f"DEBUG: Content in output file is:")
            print(f"{output_file}")
            binned = output_file.loc[output_file['Contig'].isin(binned_reads)]
            unbinned = output_file.loc[output_file['Contig'].isin(unbinned_reads)]
            for _, row in binned.iterrows():
                assert row['Bin'] == "bin1_scaffold1", f"{row['Contig']} incorrectly assigned to {binned['Bin'].iloc[0]}"
            for _, row in unbinned.iterrows():
                assert row['Bin'] == "unbinned", f"{row['Contig']} should be unbinned but assigned to {row['Bin']}!"

        # Step 3: Count the features using a SAM and GTF file
        results = run_counter(mock_sam_files, mock_gtf_files, threads=2, gene="CDS", gene_id="gene_id")
        for file in results:
            if os.path.exists(file):
                print(f"DEBUG: Output file {file} was created successfully")
                file_size = os.path.getsize(file)
                assert file_size > 0, f"Output file {file} is empty!"
                print(f"DEBUG: File size: {file_size} bytes")
            else:
                print(f"DEBUG: WARNING - Output file {file} was NOT created")

        # Step 4: Parse the output featureCounts to assign the gene ids to bins contigs align well to, if no alignment than
        # It is returned as unbinned as above
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

        print(f"✅ Pipeline completed and tests all passed!")

@pytest.fixture(autouse=True)
def cleanup_test_files(request):
    def cleanup():
        import glob
        for paf_file in glob.glob("*.paf"):
            if os.path.exists(paf_file):
                os.unlink(paf_file)
        for tsv_file in glob.glob("*.tsv"):
            if os.path.exists(tsv_file):
                os.unlink(tsv_file)
        for fasta in glob.glob("*.fa"):
            if os.path.exists(fasta):
                os.unlink(fasta)
        for txt_file in glob.glob("*.txt"):
            if os.path.exists(txt_file):
                os.unlink(txt_file)
        for summary_file in glob.glob("*.txt.summary"):
            if os.path.exists(summary_file):
                os.unlink(summary_file)
    request.addfinalizer(cleanup)

