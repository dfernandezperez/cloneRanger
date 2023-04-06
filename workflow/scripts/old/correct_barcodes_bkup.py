import gzip
import re
import pandas as pd
from Bio import SeqIO

# input_file  = "data/bc_filt/collapsed/CTR-rep1_FB_S1_L001_R2_001_GFP_collapsed-hd4.fastq.gz"
# output_file = "data/bc_filt/collapsed/corrected.fastq.gz"
# feature_ref   = "data/bc_filt/collapsed/feature_ref.csv"
# color = "GFP"
# read = "R2"
input_file    = snakemake.input[0]
output_fastq  = snakemake.output[0]
feature_ref   = snakemake.output[1]
color         = snakemake.wildcards.larry_color
read          = snakemake.wildcards.read

cluster_dict_id        = dict()
record_list            = []
consensus_read_pattern = re.compile(r'cluster_size=')

#------------------------------------------------------------------------------------------------------------------------------------
# Declare functions
#------------------------------------------------------------------------------------------------------------------------------------

def create_feature_ref(bc_dict, read):
    """Get dictionary with barcode sequences and their corresponding larry color to 
    generate the feature reference csv file required by cellranger."""
    
    bc_df                 = pd.DataFrame.from_dict(bc_dict, orient='index', columns = ['id'])
    bc_df                 = bc_df.rename_axis("sequence").reset_index()
    bc_df['name']         = bc_df.groupby('id').cumcount() + 1
    bc_df['name']         = bc_df['id'] + "_" + bc_df['name'].astype(str)
    bc_df['id']           = bc_df['name']
    bc_df['read']         = read
    bc_df['pattern']      = "(BC)"
    bc_df['feature_type'] = "Custom"
    return(bc_df)

#------------------------------------------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------------------------------------------

# Open the compressed fastq file and store in a dictionary the reference barcode for every cluster of similar barcodes
with gzip.open(input_file, 'rt') as input_handle:
    for record in SeqIO.parse(input_handle, 'fastq'):
        # Store records in a list to access them again later instead of parsing the fastq 2 times
        record_list += record,
        # Look for reads containing consensus barcode
        match = consensus_read_pattern.search(record.description)
        if match:
            # Create a dictionary with the cluster id and the corresponding consensus sequence
            cluster_id = re.sub('^.* cluster_id=([0-9]+).*', r'\1', record.description)
            cluster_dict_id[cluster_id] = record.seq

# Loop through the fastq entries again and substite every sequence by the corresponding reference sequence of every cluster of reads
with gzip.open(output_fastq, 'wt') as output_handle:
    for record in record_list:
    # Get cluster id for read, use that cluster id to access cluster_dict_id dict
    # which contains the reference sequence. This way I generate a second dict in which
    # Every read id is a associated to a specific sequence
        cluster_id = re.sub('^.* cluster_id=([0-9]+).*', r'\1', record.description)
        record.seq = cluster_dict_id[cluster_id]
        count = SeqIO.write(record, output_handle, 'fastq')

# Generate a feature referece csv file with the format specified by 10X cellranger
bc_seqs_dict = {str(key): color for key in cluster_dict_id.values()}
create_feature_ref(bc_seqs_dict, read).to_csv(feature_ref, index = False)