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
import re
import shutil

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
def temp_test_dir():
    temp_dir = tempfile.mkdtemp(prefix='Pipeline_test')
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def multiple_assemblies(temp_test_dir):
    """Create assembly files where contig1 will match bin sequences well, others poorly."""
    assembly_1_content = f">a1_contig1\n{perfect_match}\n>a1_contig2\n{bad_match}\n>a1_contig3\n{bad_coverage}\n>a1_contig4\n{bad_quality}"
    assembly_2_content = f">a2_contig1\n{bad_match}\n>a2_contig2\n{perfect_match}\n>a2_contig3\n{bad_coverage}\n>a2_contig4\n{bad_quality}"

    assembly_files = []
    for i, content in enumerate([assembly_1_content, assembly_2_content], 1):
        assembly_file = os.path.join(temp_test_dir, f'test_assembly{i}.fasta')
        with open(assembly_file, 'w') as f:
            f.write(content)
        assembly_files.append(assembly_file)
    return assembly_files

@pytest.fixture
def sample_bins(temp_test_dir):
    """Create bin files where one sequence matches contig1 perfectly, others don't match well."""
    bin1_content = f">bin1_scaffold1\n{perfect_match}\n>bin1_scaffold2\n{generate_random_sequence(3000)}"
    bin2_content = f">bin2_scaffold1\n{generate_random_sequence(5000)}\n>bin2_scaffold2\n{perfect_match[:20000]}"

    bin_files = []
    for i, content in enumerate([bin1_content, bin2_content], 1):
        bin_file = os.path.join(temp_test_dir, f'test_bin{i}.fasta')
        with open(bin_file, 'w') as f:
            f.write(content)
        bin_files.append(bin_file)
    return bin_files

@pytest.fixture
def mock_gff_files(temp_test_dir):
    """Create GFF files corresponding to the assembly files."""
    gff1_content = """##gff-version 3
##sequence-region a1_contig1 1 64000
##sequence-region a1_contig2 1 64000
##sequence-region a1_contig3 1 48000
##sequence-region a1_contig4 1 64000
a1_contig1	prokka	CDS	100	500	.	+	0	ID=GENE001;locus_tag=GENE001;product=hypothetical protein
a1_contig1	prokka	CDS	600	1200	.	-	0	ID=GENE002;locus_tag=GENE002;product=DNA polymerase
a1_contig1	prokka	CDS	1500	2000	.	+	0	ID=GENE003;locus_tag=GENE003;product=ribosomal protein
a1_contig2	prokka	CDS	200	800	.	+	0	ID=GENE004;locus_tag=GENE004;product=membrane protein
a1_contig2	prokka	CDS	1000	1600	.	-	0	ID=GENE005;locus_tag=GENE005;product=ATP synthase
a1_contig3	prokka	CDS	300	900	.	+	0	ID=GENE006;locus_tag=GENE006;product=transcriptase
a1_contig4	prokka	CDS	150	750	.	+	0	ID=GENE007;locus_tag=GENE007;product=helicase
"""

    gff2_content = """##gff-version 3
##sequence-region a2_contig1 1 64000
##sequence-region a2_contig2 1 64000
##sequence-region a2_contig3 1 48000
##sequence-region a2_contig4 1 64000
a2_contig1	prokka	CDS	100	500	.	+	0	ID=GENE008;locus_tag=GENE008;product=hypothetical protein
a2_contig1	prokka	CDS	600	1200	.	-	0	ID=GENE009;locus_tag=GENE009;product=DNA polymerase
a2_contig2	prokka	CDS	200	800	.	+	0	ID=GENE010;locus_tag=GENE010;product=membrane protein
a2_contig2	prokka	CDS	1000	1600	.	-	0	ID=GENE011;locus_tag=GENE011;product=ATP synthase
a2_contig2	prokka	CDS	1800	2400	.	+	0	ID=GENE012;locus_tag=GENE012;product=kinase
a2_contig3	prokka	CDS	300	900	.	+	0	ID=GENE013;locus_tag=GENE013;product=transcriptase
a2_contig4	prokka	CDS	150	750	.	+	0	ID=GENE014;locus_tag=GENE014;product=helicase
"""

    gff_files = []
    for i, content in enumerate([gff1_content, gff2_content], 1):
        gff_file = os.path.join(temp_test_dir, f'test_assembly{i}.gff')
        with open(gff_file, 'w') as f:
            f.write(content)
        gff_files.append(gff_file)
    return gff_files

@pytest.fixture
def mock_sam_files(temp_test_dir):
    """Create SAM files with reads mapping to the assembly contigs."""
    sam1_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a1_contig1	LN:64000
@SQ	SN:a1_contig2	LN:64000
@SQ	SN:a1_contig3	LN:48000
@SQ	SN:a1_contig4	LN:64000
read1	99	a1_contig1	150	60	100M	=	250	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read1	147	a1_contig1	250	60	100M	=	150	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	99	a1_contig1	300	60	100M	=	400	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	a1_contig1	400	60	100M	=	300	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	99	a1_contig1	700	60	100M	=	800	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	147	a1_contig1	800	60	100M	=	700	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	99	a1_contig2	250	60	100M	=	350	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	147	a1_contig2	350	60	100M	=	250	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	99	a1_contig3	400	60	100M	=	500	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read5	147	a1_contig3	500	60	100M	=	400	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam2_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:a2_contig1	LN:64000
@SQ	SN:a2_contig2	LN:64000
@SQ	SN:a2_contig3	LN:48000
@SQ	SN:a2_contig4	LN:64000
read6	99	a2_contig1	150	60	100M	=	250	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read6	147	a2_contig1	250	60	100M	=	150	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	99	a2_contig2	200	60	100M	=	300	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read7	147	a2_contig2	300	60	100M	=	200	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	99	a2_contig2	1100	60	100M	=	1200	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read8	147	a2_contig2	1200	60	100M	=	1100	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	99	a2_contig2	1900	60	100M	=	2000	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read9	147	a2_contig2	2000	60	100M	=	1900	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read10	99	a2_contig3	350	60	100M	=	450	200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read10	147	a2_contig3	450	60	100M	=	350	-200	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

    sam_files = []
    for i, content in enumerate([sam1_content, sam2_content], 1):
        sam_file = os.path.join(temp_test_dir, f'test_sample{i}.sam')
        with open(sam_file, 'w') as f:
            f.write(content)
        sam_files.append(sam_file)
    return sam_files

# =============================================================================
# Import the functions from the individual_functions folder
# =============================================================================
from individual_functions.contig_mapping import build_minimap2_command, run_alignment
from individual_functions.PAF_parsing import PAF_parsing, run_paf_parsing
from individual_functions.feature_counting import build_featureCounts, execute_feature_counting, run_counter
from individual_functions.feature_parsing import extract_file_identifiers, pair_files_by_sample, feature_parsing, run_parsing
from utils.file_pairing import extract_file_identifiers, pair_files_by_sample

# =============================================================================
# Integration test - Test functions working together
# =============================================================================
class TestFullPipeline:
    def test_complete_pipeline(self, sample_bins, multiple_assemblies, mock_sam_files, mock_gff_files):
        """
        Test to check if we're ready for integration testing.
        This test checks if individual functions can be imported.
        """
        pattern_source = None
        # Step 1: Align
        aligned_paf_files = run_alignment(sample_bins, multiple_assemblies)
        if aligned_paf_files and all(os.path.exists(paf_file) for paf_file in aligned_paf_files):
            for paf_file in aligned_paf_files:
                paf_out = pd.read_csv(paf_file, delimiter='\t', header=None)
                query_names = paf_out.iloc[:, 0].tolist()
                mapped_reads = ["a1_contig1", "a2_contig2"]
                found_reads = set(mapped_reads) & set(query_names)
                assert found_reads, f"None of {mapped_reads} found in PAF file {paf_file} queries: {set(query_names)}"
            print(f"✅ Alignment completed and tests all passed!")
        else:
            print(f"{chr(0x2757)} ERROR: This step did not generate any output")

        # Step 2: Parse the PAF output
        binned_reads = ["a1_contig1", "a2_contig2"]
        unbinned_reads = ["a1_contig2", "a1_contig3", "a1_contig4", "a2_contig1", "a2_contig3", "a2_contig4"]
        run_paf_parsing(aligned_paf_files, multiple_assemblies, sample_bins)
        expected_outputs = [os.path.splitext(input_file)[0] + "_contigs_to_bin_mapping.txt" for input_file in aligned_paf_files]
        if all(os.path.exists(output_file) for output_file in expected_outputs):
            for input_file in aligned_paf_files:
                expected_output = os.path.splitext(input_file)[0] + "_contigs_to_bin_mapping.txt"
                output_file = pd.read_csv(expected_output, delimiter='\t', header=0)
                read_names = output_file.iloc[:, 0].tolist()
                print(f"DEBUG: Content in output file is:")
                print(f"{output_file}")
                binned = output_file.loc[output_file['Contig'].isin(binned_reads)]
                unbinned = output_file.loc[output_file['Contig'].isin(unbinned_reads)]
                for _, row in binned.iterrows():
                    assert row['Bin'] == "test_bin1.fasta", f"{row['Contig']} incorrectly assigned to {binned['Bin'].iloc[0]}"
                for _, row in unbinned.iterrows():
                    assert row['Bin'] == "unbinned", f"{row['Contig']} should be unbinned but assigned to {row['Bin']}!"
            print(f"✅ PAF parsing completed and tests all passed!")
        else:
            print(f"{chr(0x2757)} ERROR: PAF parsing step did not generate expected output files")

        # Step 3: Count the features using a SAM and GTF file
        results = run_counter(mock_sam_files, mock_gff_files, threads=2, gene="CDS", gene_id="locus_tag", paired_end=False)
        if all(os.path.exists(file) for file in results):
            assembly1_genes = ['GENE001', 'GENE002', 'GENE003', 'GENE004', 'GENE005', 'GENE006', 'GENE007']
            assembly2_genes = ['GENE008', 'GENE009', 'GENE010', 'GENE011', 'GENE012', 'GENE013', 'GENE014']
            for i, file in enumerate(results):
                print(f"Output file {file} was created successfully")
                file_size = os.path.getsize(file)
                assert file_size > 0, f"Output file {file} is empty!"
                print(f"File size: {file_size} bytes")

                df = pd.read_csv(file, delimiter='\t', comment='#')
                print(f"The content of the feature count file is:")
                print(f"{df}")
                if i == 0:
                    expected_genes = assembly1_genes
                else:
                    expected_genes = assembly2_genes

                found_genes = df['Geneid'].tolist()
                missing_genes = set(expected_genes) - set(found_genes)
                assert len(missing_genes) == 0, f"Missing genes in output: {missing_genes}"
                print(f"All expected genes found in {file}")

                count_column = df.columns[-1]
                total_counts = df[count_column].sum()
                assert total_counts > 0, f"No reads counted in {file}"
                print(f"Total read counts in {file}: {total_counts}")
            print(f"✅ Feature counting completed and tests all passed!")
        else:
            print(f"{chr(0x2757)} ERROR: Feature counting step did not generate expected output files")

        # Step 4: Parse the output featureCounts to assign the gene ids to bins contigs align well to, if no alignment than
        # It is returned as unbinned as above
        feature_parsing_results = run_parsing(expected_outputs, results)
        print(f"{feature_parsing_results}")
        expected_feature_outputs = []
        for count_file in results:
            count_base = os.path.basename(count_file)
            count_name = os.path.splitext(count_base)[0]
            feature_output = count_name + "_binned.txt"
            expected_feature_outputs.append(feature_output)
        if all(os.path.exists(output_file) for output_file in expected_feature_outputs):
            for feature_file in expected_feature_outputs:
                print(f"DEBUG: Output file {feature_file} was created successfully")
                file_size = os.path.getsize(feature_file)
                assert file_size > 0, f"Output file {feature_file} is empty!"
                print(f"DEBUG: File size: {file_size} bytes")
                binned_counts = pd.read_csv(feature_file, delimiter='\t', header=0)
                binned_genes = binned_counts[binned_counts['Chr'].isin(['a1_contig1', 'a2_contig2'])]
                unbinned_genes = binned_counts[binned_counts['Chr'].isin(['a1_contig2', 'a1_contig3', 'a1_contig4', 'a2_contig1', 'a2_contig3', 'a2_contig4'])]
                for _, row in binned_genes.iterrows():
                    assert row['binning'] in ["test_bin1.fasta"], f"Gene {row['Geneid']} incorrectly assigned: {row['binning']}"
                for _, row in unbinned_genes.iterrows():
                    assert row['binning'] == "unbinned", f"Gene {row['Geneid']} should be unbinned but assigned to {row['binning']}"
            print(f"✅ Feature parsing completed and tests all passed!")
        else:
            print(f"{chr(0x2757)} ERROR: Feature parsing step did not generate expected output files")
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

