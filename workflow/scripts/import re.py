import gzip
import re
import pandas as pd
from Bio import SeqIO

# Define the input and output file names
input_file  = '/stemcell/scratch/dfernandezp/IndranilSingh/aml_dnmt3_larry/data/scRNAseq_LARRY/config_files/tmp_fq/FB_Day5_exp1_1-SCI7T030_H2CNKDSX5_S4_L004_R2_001.fastq.gz'
output_file = 'output.fastq'

# Create a dictionary of barcode_patterns : color. This will be used to know to which color
# corresponds each pattern during the regex matching parsing the fastq file
barcodes_dict = {
    "GAGAGACCATATGTGTGTGAATCCAGCCTACCGCTGA": "GFP",
    "GAGTAATGAGCTGCATGTCTGGATCCGA" : "Sapphire"
}
patterns      = [re.compile(r'{}'.format(re.escape(barcode))) for barcode in barcodes_dict.keys()]
patterns_dict = {key: value for key, value in zip(patterns, barcodes_dict.values())}

# Store matched sequences with their corresponding barcode
matched_seqs = {}

# Open the compressed fastq file & output file
with gzip.open(input_file, 'rt') as input_handle, gzip.open(output_file, 'wt') as output_handle:
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



def create_feature_ref(bc_dict):
    """Get dataframe with barcode counts and generate the feature reference csv file required by cellranger."""
    bc_df = pd.DataFrame.from_dict(matched_seqs, orient='index', columns = 'id')

    bc_df = bc_df.rename_axis("sequence").reset_index()
    bc_df['id']      = [ id + '_' + str(i) for i in np.arange(len(bc_df)) ]
    bc_df['name']    = 'str' + df['col'].astype(str)
    bc_df['read']    = "R2"
    bc_df['pattern'] = "(BC)"
    bc_df['feature_type'] = "Custom"

    return(bc_df)
