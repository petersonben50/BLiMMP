#################################
# ABA_GID.py
# Benjamin D. Peterson
#################################


###########################
# Set up environment
###########################
import io
import os
import sys
import argparse
import glob
from Bio import AlignIO
from Bio import SearchIO
from Bio import SeqIO
from multiprocessing import Pool
import pandas as pd
import subprocess



######################################################
######################################################
# Parse commands
######################################################
######################################################

###########################
# Set up an argument parser
###########################
parser = argparse.ArgumentParser()

# Input data
parser.add_argument('--orf_file')
parser.add_argument('--orf_folder')

# Read in other inputs
parser.add_argument('--hmm')
parser.add_argument('--cluster_cutoff', default='0.97')
parser.add_argument('--n_value_cdhit', default='5')
parser.add_argument('--metagenome_list', default='Do_not_run')
parser.add_argument('--metagenomes_location', default='Do_not_run')
parser.add_argument('--metatranscriptome_location', default='Do_not_run')
parser.add_argument('--read_depth_cutoff', default=150)
parser.add_argument('--number_threads', default=2)
parser.add_argument('--reference_aa_dataset', default='none_provided')
parser.add_argument('--do_not_use_super5_for_alignment', default='super5', action='store_const', const='align')

# Output info
parser.add_argument('--output_location')
parser.add_argument('--output_prefix')

# Other
parser.add_argument('--testing', action='store_true')

# Skip commands
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--skip_orf_concat', action='store_true')
parser.add_argument('--skip_hmm_search', action='store_true')
parser.add_argument('--skip_pull_out_aa', action='store_true')
parser.add_argument('--skip_aa_alignment', action='store_true')
parser.add_argument('--skip_clustering_seqs', action='store_true')
parser.add_argument('--skip_generate_g2A', action='store_true')


###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
ORF_FILE = inputs.orf_file
ORF_FOLDER = inputs.orf_folder

# Read in other inputs
HMM = inputs.hmm
CLUSTER_CUTOFF = inputs.cluster_cutoff
N_VALUE_CDHIT = inputs.n_value_cdhit
METAGENOME_LIST = inputs.metagenome_list
METAGENOMES_LOCATION = inputs.metagenomes_location
METATRANSCRIPTOMES_LOCATION = inputs.metatranscriptome_location
READ_DEPTH_CUTOFF = inputs.read_depth_cutoff
NUMBER_THREADS = inputs.number_threads
REFERENCE_AA_DATASET = inputs.reference_aa_dataset
SUPER5_INPUT = inputs.do_not_use_super5_for_alignment

# Output info
OUTPUT_LOCATION = inputs.output_location
OUTPUT_LOCATION = OUTPUT_LOCATION + '/'
OUTPUT_PREFIX = inputs.output_prefix


# Other
TESTING = inputs.testing

# Skip commands
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
SKIP_ORF_CONCAT = inputs.skip_orf_concat
SKIP_HMM_SEARCH = inputs.skip_hmm_search
SKIP_PULL_OUT_AA = inputs.skip_pull_out_aa
SKIP_AA_ALIGNMENT = inputs.skip_aa_alignment
SKIP_CLUSTERING_SEQS = inputs.skip_clustering_seqs
SKIP_GENERATE_G2A = inputs.skip_generate_g2A

######################################################
######################################################
# Verify commands and establish modes
######################################################
######################################################

###########################
# Make sure we have only one input, set the type
###########################
if ORF_FILE is not None:
    if ORF_FOLDER is not None:
        print("You've supplied both an concatenated ORF file and a folder with the ORFs for the bins. Only supply one.")
        sys.exit()
    else:
        print("Input is single ORFs file")
        input_type = "single_orf"
elif ORF_FOLDER is not None:
    print("Input is a folder of ORF files")
    input_type = "orf_folder"
else:
    print("You haven't supplied an input for the bins")
    sys.exit()


###########################
# Check if output directory already exists
###########################
# Set location for temporary working directory
working_directory = OUTPUT_LOCATION + '/working_directory/'
if os.path.isdir(OUTPUT_LOCATION) == True:
    if SKIP_NEW_DIRECTORY:
        print("Directory exists, but you said that's okay")
        if os.path.isdir(working_directory) == False:
            os.mkdir(working_directory)
    else:
        print("Hey dummy, " + OUTPUT_LOCATION + " is already a directory. Please give me an empty directory")
        sys.exit()
else:
    os.mkdir(OUTPUT_LOCATION)
    os.mkdir(working_directory)


######################################################
######################################################
# Prepare ORFs
######################################################
######################################################

# Variable for ORFs to use:
concat_orf_to_use = working_directory + 'all_ORFs_concat.faa'
g2a_file = working_directory + OUTPUT_PREFIX + '_ORFs_G2A.tsv'


###########################
# Input set-up: using folder of ORFs
###########################
# If we supplied folder of ORFs, concatenate them to a folder
if SKIP_ORF_CONCAT:
    print("Skipping ORF concatenation")
else:
    if input_type == "orf_folder":
        print("Concatenating ORFs and generating G2A file from all assemblies")
        concat_cmd = "cat "
        g2akey = dict()
        genome_files = ORF_FOLDER + "/**"
        genomes = glob.glob(genome_files)
        for genome in genomes:
            if genome.endswith('.faa'):
                # Add assembly to list for concatenation
                concat_cmd = concat_cmd + " " + genome
        concat_cmd = concat_cmd + " > " + concat_orf_to_use
        os.system(concat_cmd)
        # Generate gene-to-assembly file
        g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + ORF_FOLDER + " > " + g2a_file
        os.system(g2a_cmd)


###########################
# Input set-up: using single ORF file
###########################
# If we supplied concatenated ORF, move it over:
if input_type == "single_orf":
    print("Copying " + ORF_FILE + " to " + concat_orf_to_use)
    cp_command = "cp " + ORF_FILE + " " + concat_orf_to_use
    os.system(cp_command)
    # Still going to make the gene-to-assembly file
    g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + working_directory + " > " + g2a_file
    os.system(g2a_cmd)


###########################
# Set other variables
###########################
hmmer_results_file_name = working_directory + OUTPUT_PREFIX + '_HMM.out'
g2a_for_gene = OUTPUT_LOCATION + OUTPUT_PREFIX + '_G2A.tsv'
fasta_output_for_hits = OUTPUT_LOCATION + '/' + OUTPUT_PREFIX + '.faa'
PA_files = [i for i in glob.glob(METATRANSCRIPTOMES_LOCATION + "/*.tsv")]
MT_output = OUTPUT_LOCATION + '/' + OUTPUT_PREFIX + '_MT_coverage.tsv'


######################################################
######################################################
# Run HMM search
######################################################
######################################################

###########################
# Run all HMMs on ORF file
###########################
if SKIP_HMM_SEARCH:
    print("Simon says skip the HMM run")
else:
    print("Running HMM-based search for " + OUTPUT_PREFIX)
    hmmer_log_file_name = working_directory + OUTPUT_PREFIX + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu ' + NUMBER_THREADS + ' -E 0.000001 ' + HMM + " " + concat_orf_to_use + " > " + hmmer_log_file_name 
    os.system(hmm_cmd)


###########################
# Save out tsv file with the gene-to-assembly information
###########################
if SKIP_GENERATE_G2A:
    print("Simon say skip the G2A file generation")
else:
    hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
    for sampleID in hmmer_output:
        g2b_for_gene_cmd = "awk '$1 == \"" + sampleID.id + "\" { print $0 }' " + g2a_file + " >> " + g2a_for_gene
        os.system(g2b_for_gene_cmd)
g2a_data = pd.read_csv(g2a_for_gene, delimiter="\t", names=['gene', 'assembly'])


###########################
# Pull out amino acid sequences
###########################
if SKIP_PULL_OUT_AA:
    print("Simon says skip the HMM run")
else:
    hmmer_results_file_length = subprocess.check_output('wc -l < ' + hmmer_results_file_name, shell=True)
    if int(hmmer_results_file_length) > 13:
        print("Extracting AA sequences for " + OUTPUT_PREFIX)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for seq_record in SeqIO.parse(concat_orf_to_use, "fasta"):
                for sampleID in hmmer_output:
                    if sampleID.id == seq_record.id:
                        resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(seq_record.seq).replace("*","") + '\n')
    else:
        print('No hits for ' + OUTPUT_PREFIX + '. Ending the script now.')
        sys.exit()


###########################
# Align amino acid sequences to HMM
###########################
if SKIP_AA_ALIGNMENT:
    print("Simon says skip the AA alignment")
else:
    if os.path.isfile(fasta_output_for_hits):
        sto_output = working_directory + OUTPUT_PREFIX + '.sto'
        afa_output = OUTPUT_LOCATION + OUTPUT_PREFIX + '.afa'
        print("Aligning sequences of " + OUTPUT_PREFIX + " to HMM")
        hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
        os.system(hmmalign_cmd)
        # Write out afa alignment
        with open(afa_output, 'w') as alignOut:
            for align_record in AlignIO.read(sto_output, "stockholm"):
                alignOut.write('>' + align_record.id + '\n')
                alignOut.write(str(align_record.seq) + '\n')
        os.remove(sto_output)
    else:
        print("No fasta file with the HMM hits in it.")
        sys.exit()


######################################################
######################################################
# Cluster sequences
######################################################
######################################################
if SKIP_CLUSTERING_SEQS:
    print("Simon says skip clustering seqs")
else:
    derep_fasta = working_directory + OUTPUT_PREFIX + '_derep.faa'
    clustering_info_output = OUTPUT_LOCATION + OUTPUT_PREFIX + '_cluster_data.tsv'
    clustering_info_output_log = working_directory + OUTPUT_PREFIX + '_cluster_data_log.txt'
    cdhit_cmd = "cd-hit -g 0 -i " + fasta_output_for_hits
    cdhit_cmd = cdhit_cmd + " -o " + derep_fasta
    cdhit_cmd = cdhit_cmd + " -c " + CLUSTER_CUTOFF
    cdhit_cmd = cdhit_cmd + " -n " + N_VALUE_CDHIT
    cdhit_cmd = cdhit_cmd + " -d 0 "
    cdhit_cmd = cdhit_cmd + " > " + clustering_info_output_log
    os.system(cdhit_cmd)
    cdhit_parsing_cmd = "clstr2txt.pl " + derep_fasta + ".clstr > " + clustering_info_output
    os.system(cdhit_parsing_cmd)



######################################################
######################################################
# Pull out MG depth information
######################################################
######################################################
if METAGENOME_LIST == "Do_not_run":
    print("List of metagenomes not provided")
if METAGENOMES_LOCATION == "Do_not_run":
    print("Folder of metagenomes not provided")
if METAGENOME_LIST != "Do_not_run" and METAGENOMES_LOCATION != "Do_not_run":
    print("Pulling out mapping information for " + OUTPUT_PREFIX)
    # Set up G2A key
    for index, row in g2a_data.iterrows():
        scaffold_of_interest = row.gene.rsplit("_", 1)[0]
        #print("Mapping data for " + scaffold_of_interest)
        with open(METAGENOME_LIST, 'r') as mg_list:
            for metagenome_nl in mg_list.readlines():
                metagenome = metagenome_nl.strip()
                mg_cov_out_raw = working_directory + metagenome + "_" + OUTPUT_PREFIX + "_coverage_raw.tsv"
                mapping_file = METAGENOMES_LOCATION + "/" + metagenome + "_to_" + row.assembly + ".bam"
                if os.path.isfile(mapping_file):
                    sam_cmd = "samtools depth -a -r " + scaffold_of_interest + " " + mapping_file + " >> " + mg_cov_out_raw
                    os.system(sam_cmd)
                else:
                    print("   " + metagenome + " not mapped to " + row.assembly)
    # Aggregate the coverage within metagenomes
    with open(METAGENOME_LIST, 'r') as mg_list:
        for metagenome_nl in mg_list.readlines():
            metagenome = metagenome_nl.strip()
            mg_cov_out_raw = working_directory + metagenome + "_" + OUTPUT_PREFIX + "_coverage_raw.tsv"
            mg_cov_out = working_directory + metagenome + "_" + OUTPUT_PREFIX + "_coverage.tsv"
            # Open up coverage table
            mg_cov_data_raw = pd.read_table(mg_cov_out_raw, names = ['contigs', 'locus', 'depth'])
            mg_cov_data_raw = mg_cov_data_raw[mg_cov_data_raw['locus'] >= READ_DEPTH_CUTOFF]
            # Filter out residues at end of contig
            lengthOfContig = mg_cov_data_raw[['contigs', 'locus']].groupby('contigs').max()
            lengthOfContig.rename(columns = {'locus':'lengthOfContig'}, inplace = True)
            # Save out original length
            lengthOfContigOriginal = lengthOfContig
            # Substract the length to filter out
            lengthOfContig['maxLengthToInclude'] = lengthOfContig['lengthOfContig'] - READ_DEPTH_CUTOFF
            # Then, join max contig length DF with depth table
            mg_cov_data_raw = pd.merge(mg_cov_data_raw, lengthOfContig, on='contigs', how='outer')
            # Filter out end of contig
            mg_cov_data_raw = mg_cov_data_raw[mg_cov_data_raw['locus'] <= mg_cov_data_raw['maxLengthToInclude']]
            mg_cov_data_raw = mg_cov_data_raw[['contigs', 'depth', 'lengthOfContig']]
            # Aggregate depth by contig
            mg_cov_data = mg_cov_data_raw.groupby('contigs').mean()
            # Read out data
            mg_cov_data.to_csv(mg_cov_out, sep='\t', header = False)
            # Add column with metagenome name
            add_name_column = 'sed -i "s/$/\t' + metagenome + '/" ' + mg_cov_out
            os.system(add_name_column)
    all_mg_cov = OUTPUT_LOCATION + OUTPUT_PREFIX + "_MG_coverage.tsv"
    concat_cov_cmd = "cat " + working_directory + "*" + OUTPUT_PREFIX + "_coverage.tsv > " + all_mg_cov
    print(concat_cov_cmd)
    os.system(concat_cov_cmd)



######################################################
######################################################
# Pull out MT coverage information
######################################################
######################################################

def retrieve_RNA_pseudoalignment_counts(MT_2_A_PA_file):
    MT_ID = MT_2_A_PA_file.rsplit("/", 1)[1].split("_to_")[0]
    kallisto_data = pd.DataFrame(columns=['geneID', 'length_gene', 'effective_length', 'counts', 'tpm'])
    print(MT_2_A_PA_file)
    with open(MT_2_A_PA_file, 'r') as PA_data:
        for kallisto_entry_NL in PA_data.readlines():
            kallisto_entry = kallisto_entry_NL.rstrip()
            kallisto_geneID = kallisto_entry.split('\t')[0]
            if kallisto_geneID in list(g2a_data.gene):
                kallisto_entry_df = pd.read_csv(io.StringIO(kallisto_entry), delimiter = '\t', header = None, names=['geneID', 'length_gene', 'effective_length', 'counts', 'tpm'])
                kallisto_data = pd.concat([kallisto_data, kallisto_entry_df], ignore_index = True)
        kallisto_data['mtID'] = MT_ID
        return kallisto_data

def combine_RNA_pseudoalignment_counts(PA_file_list, MT_output_file):
    MT_results_list = list()
    print("Pulling out MT read counts")
    with Pool(int(NUMBER_THREADS)) as pool:
        for result in pool.map(retrieve_RNA_pseudoalignment_counts, PA_file_list):
            MT_results_list.append(result)
    MT_results_df = pd.concat(MT_results_list, ignore_index = True)
    MT_results_df.to_csv(MT_output_file, sep = '\t', index = False, header = True)
if METATRANSCRIPTOMES_LOCATION == "Do_not_run":
    print("No MT location provided")
elif os.path.isdir(METATRANSCRIPTOMES_LOCATION) == True:
    combine_RNA_pseudoalignment_counts(PA_files, MT_output)
else:
    print("Provided MT location is not a folder.")

######################################################
######################################################
# Generate rough tree with reference amino acid dataset
######################################################
######################################################
if REFERENCE_AA_DATASET == 'none_provided':
    print("No reference dataset provided")
else:
    if os.path.isfile(REFERENCE_AA_DATASET):
        print(" Generating a rough tree for " + OUTPUT_PREFIX)
        tree_orfs_to_use = working_directory + OUTPUT_PREFIX + "_orfs_for_tree.faa"
        tree_align_to_use = working_directory + OUTPUT_PREFIX + "_orfs_for_tree.afa"
        tree_align_to_use_cleaned = working_directory + OUTPUT_PREFIX + "_orfs_for_tree_cleaned.afa"
        tree_output = OUTPUT_LOCATION + OUTPUT_PREFIX + ".tree"
        cat_cmd = "cat " + fasta_output_for_hits + " " + REFERENCE_AA_DATASET + " > " + tree_orfs_to_use
        align_cmd = "muscle -" + SUPER5_INPUT + " " + tree_orfs_to_use + " -output " + tree_align_to_use
        clean_cmd = "trimal -in " + tree_align_to_use + " -out " + tree_align_to_use_cleaned + " -gt 0.5"
        tree_cmd = "FastTree " + tree_align_to_use_cleaned + " > " + tree_output
        print(" Concatenating " + OUTPUT_PREFIX + " sequences.")
        os.system(cat_cmd)
        print(" Aligning " + OUTPUT_PREFIX + " sequences.")
        os.system(align_cmd)
        os.system(clean_cmd)
        print(" Running FastTree on " + OUTPUT_PREFIX + " sequences.")
        os.system(tree_cmd)
    else:
        print("Provided file doesn't exist. No tree will be generated.")



######################################################
######################################################
# Pull out gene neighborhood
######################################################
######################################################


######################################################
######################################################
# Lines to break script while testing, if needed
######################################################
######################################################
if TESTING:
    print(str(TESTING) + ", I'm testing")
    sys.exit()
