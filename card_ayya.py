import os
import pandas as pd


# Directory containing the gene mapping data files
input_dir = '/net/mraid20/export/jasmine/card/NICU'
output_file = 'sampleXgenes.csv'
metadata_file = 'gene_metadata.csv'


# Initialize an empty DataFrame for storing the samples X genes data
samples_genes_df = pd.DataFrame()


# Initialize an empty DataFrame for storing the metadata
metadata_df = pd.DataFrame()


# Loop through each file in the directory
for filename in os.listdir(input_dir):
    if filename.endswith('.gene_mapping_data.txt'):
        sample_id = filename.split('.')[0]  # Extract the sample ID from the filename (including "_v2_fullrun")


        # Read the file into a DataFrame
        file_path = os.path.join(input_dir, filename)
        df = pd.read_csv(file_path, sep='\t')


        # Calculate RPK (Reads Per Kilobase)
        df['Reference Length (kb)'] = df['Reference Length'] / 1000
        df['RPK'] = df['All Mapped Reads'] / df['Reference Length (kb)']


        # Create a Series with ARO Term as the index and RPK as the value
        rpk_series = df.set_index('ARO Term')['RPK']


        # Add the Series to the samples_genes_df DataFrame
        samples_genes_df[sample_id] = rpk_series


        # Extract metadata information and add to metadata_df if not already present
        metadata_subset = df[['ARO Term', 'ARO Accession', 'AMR Gene Family', 'Drug Class', 'Resistance Mechanism']]
        metadata_subset.set_index('ARO Term', inplace=True)
        metadata_df = pd.concat([metadata_df, metadata_subset]).drop_duplicates()


# Transpose the samples_genes_df so that samples are rows and genes are columns
samples_genes_df = samples_genes_df.transpose()


# Save the resulting DataFrames as CSV files
samples_genes_df.to_csv(output_file)
metadata_df.to_csv(metadata_file)
