import glob
import pandas as pd

PA_folder='/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/output'
b2a='/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/bin_to_assembly.tsv'
g2b='/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/hgcA_bins_G2B.tsv'
output_file='/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/bin_MT_data.tsv'
# Read in bin to assembly file
bin_to_assembly = pd.read_csv(b2a, sep='\t', header=None, names=['bin_id', 'assembly_id'])
# Get list of all assembly IDs
assemblyIDs = list(set(bin_to_assembly['assembly_id']))
assemblyIDs.sort()
# Read in gene to bin file
gene_to_bin = pd.read_csv(g2b, sep='\t', header=None, names=['gene_id', 'bin_id'])
# Add column with anvio or Das Tool
gene_to_bin['tool'] = gene_to_bin['bin_id'].str.split('_').str[2]

# Set up a dataframe to store all data
all_data = pd.DataFrame()

# Aggregate data
# Loop over all assemblies
for assembly in assemblyIDs:
    # Set up a temporary dataframe to store all data for this assembly
    temp_anvio_df = pd.DataFrame()
    # Look for all PA files associated with this assembly in anvio
    anvio_PA_files = glob.glob(f"{PA_folder}/*_to_{assembly}_kallisto_anvio.tsv")
    # If no anvio PA files, skip
    # Sort anvio PA files
    anvio_PA_files.sort()
    # Loop over all anvio PA files
    for anvio_PA_file in anvio_PA_files:
        # Get metatranscriptome id
        metatranscriptome_id = anvio_PA_file.split('/')[-1].split('_to_')[0]
        print(f"Working on {metatranscriptome_id} pseudoalignment file for {assembly} (anvio)")
        # Read in anvio PA file
        anvio_PA = pd.read_csv(anvio_PA_file, sep='\t', header=0, names=['gene_id', 'length', 'eff_length', 'est_counts', 'tpm'])
        # Calculate mapped reads per base pair, name with MT id
        anvio_PA[metatranscriptome_id] = anvio_PA['est_counts'] / anvio_PA['eff_length']
        anvio_PA = anvio_PA[['gene_id', metatranscriptome_id]]
        # Add tool info
        anvio_PA.loc[:, ('tool')] = 'anvio'
        # Merge anvio PA file with gene to bin file
        anvio_PA = pd.merge(anvio_PA, gene_to_bin, on=['gene_id', 'tool'], how='left')
        # Drop tool column
        anvio_PA = anvio_PA.drop(columns=['tool'])
        anvio_PA = anvio_PA[['gene_id', 'bin_id', metatranscriptome_id]]
        # Add to temporary dataframe
        if temp_anvio_df.empty:
            temp_anvio_df = anvio_PA
        else:
            temp_anvio_df = pd.merge(temp_anvio_df, anvio_PA, on=['gene_id', 'bin_id'], how='outer')
    # Set up a temporary dataframe to store all data for this assembly
    temp_dasTool_df = pd.DataFrame()
    # Look for all PA files associated with this assembly in Das Tool
    dasTool_PA_files = glob.glob(f"{PA_folder}/*_to_{assembly}_kallisto_dasTool.tsv")
    # If no dasTool PA files, skip
    # Sort dasTool PA files
    dasTool_PA_files.sort()
    # Loop over all dasTool PA files
    for dasTool_PA_file in dasTool_PA_files:
        # Get metatranscriptome id
        metatranscriptome_id = dasTool_PA_file.split('/')[-1].split('_to_')[0]
        print(f"Working on {metatranscriptome_id} pseudoalignment file for {assembly} (dasTool)")
        # Read in dasTool PA file
        dasTool_PA = pd.read_csv(dasTool_PA_file, sep='\t', header=0, names=['gene_id', 'length', 'eff_length', 'est_counts', 'tpm'])
        # Calculate mapped reads per base pair, name with MT id
        dasTool_PA[metatranscriptome_id] = dasTool_PA['est_counts'] / dasTool_PA['eff_length']
        dasTool_PA = dasTool_PA[['gene_id', metatranscriptome_id]]
        # Add tool info
        dasTool_PA.loc[:, ('tool')] = 'dasTool'
        # Merge dasTool PA file with gene to bin file
        dasTool_PA = pd.merge(dasTool_PA, gene_to_bin, on=['gene_id', 'tool'], how='left')
        # Drop tool column
        dasTool_PA = dasTool_PA.drop(columns=['tool'])
        dasTool_PA = dasTool_PA[['gene_id', 'bin_id', metatranscriptome_id]]
        # Add to temporary dataframe
        if temp_dasTool_df.empty:
            temp_dasTool_df = dasTool_PA
        else:
            temp_dasTool_df = pd.merge(temp_dasTool_df, dasTool_PA, on=['gene_id', 'bin_id'], how='outer')
    # Merge anvio and dasTool dataframes if neither is empty
    if temp_anvio_df.empty and temp_dasTool_df.empty:
        continue
    elif temp_anvio_df.empty:
        temp_assembly_df = temp_dasTool_df
    elif temp_dasTool_df.empty:
        temp_assembly_df = temp_anvio_df
    else:
        temp_assembly_df = pd.concat([temp_anvio_df, temp_dasTool_df], axis=0, ignore_index=True)
    # Add to all data
    if all_data.empty:
        all_data = temp_assembly_df
    else:
        all_data = pd.concat([all_data, temp_assembly_df], axis=0, ignore_index=True)

# Write out data
all_data.to_csv(output_file, sep='\t', header=True, index=False)
