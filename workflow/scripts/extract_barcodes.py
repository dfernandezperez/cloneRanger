import gzip
import re
import pandas as pd
import numpy as np
from Bio import SeqIO

#------------------------------------------------------------------------------------------------------------------------------------
# Declare functions
#------------------------------------------------------------------------------------------------------------------------------------

def create_feature_ref(bc_dict):
    """Get dictionary with barcode sequences and their corresponding larry color to 
    generate the feature reference csv file required by cellranger."""
    
    bc_df = pd.DataFrame.from_dict(bc_dict, orient='index', columns = ['id'])
    bc_df = bc_df.rename_axis("sequence").reset_index()
    bc_df['name'] = bc_df.groupby('id').cumcount() + 1
    bc_df['name'] = bc_df['id'] + "_" + bc_df['name'].astype(str)
    bc_df['id']    = bc_df['name']
    bc_df['read']    = "R2"
    bc_df['pattern'] = "(BC)"
    bc_df['feature_type'] = "Custom"
    return(bc_df)

#------------------------------------------------------------------------------------------------------------------------------------
# Declare input/output files
#------------------------------------------------------------------------------------------------------------------------------------

# Define the input and output file names
input_file    = snakemake.input[0]
output_fastq  = snakemake.output[0]
feature_ref   = snakemake.output[1]
barcodes_dict = snakemake.params[0]

#------------------------------------------------------------------------------------------------------------------------------------
# Generate regex patterns to filer fastq
#------------------------------------------------------------------------------------------------------------------------------------

# Create a dictionary of barcode_patterns : color. This will be used to know to which color
# corresponds each pattern during the regex matching parsing the fastq file
patterns      = [re.compile(r'{}'.format(re.escape(barcode))) for barcode in barcodes_dict.keys()]
patterns_dict = {key: value for key, value in zip(patterns, barcodes_dict.values())}

# Store matched sequences with their corresponding barcode
matched_seqs = {}

#------------------------------------------------------------------------------------------------------------------------------------
# Parse fastq file and subset sequences containing valid barcodes. Generate feature reference csv file
#------------------------------------------------------------------------------------------------------------------------------------

# Open the compressed fastq file & output file
with gzip.open(input_file, 'rt') as input_handle, gzip.open(output_fastq, 'wt') as output_handle:
    # Iterate over each record in the fastq file
    for record in SeqIO.parse(input_handle, 'fastq'):
        # Look for barcode patterns
        for pattern in patterns:
            match = pattern.search(str(record.seq))
            if match:
                # Update record to contain just matched sequence 
                new_record = record[match.start() : match.end()]
                # Save matched seq into a dictionary with the associated barcode 
                if str(record.seq) not in matched_seqs: matched_seqs[str(record.seq)] = patterns_dict[pattern]
                # Write output to fastq. Just 1 match can happen (different barcodes). Once a match is found pass to next fastq entry
                count = SeqIO.write(new_record, output_handle, 'fastq')
                break 

create_feature_ref(matched_seqs).to_csv(feature_ref, index = False)