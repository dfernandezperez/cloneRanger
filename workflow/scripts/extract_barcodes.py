import gzip
import re
from Bio import SeqIO

#------------------------------------------------------------------------------------------------------------------------------------
# Declare input/output files
#------------------------------------------------------------------------------------------------------------------------------------

# Define the input and output file names
input_file    = snakemake.input[0]
output_fastqs = snakemake.output
barcodes_dict = snakemake.params[0]
read_feat_bc  = snakemake.params[1]

#------------------------------------------------------------------------------------------------------------------------------------
# Generate regex patterns to filer fastq
#------------------------------------------------------------------------------------------------------------------------------------

# Create a dictionary of barcode_patterns : color. This will be used to know to which color
# corresponds each pattern during the regex matching parsing the fastq file
patterns      = [re.compile(r'{}'.format(barcode)) for barcode in barcodes_dict.keys()]
patterns_dict = {key: value for key, value in zip(patterns, barcodes_dict.values())}

# Store matched sequences with their corresponding barcode
matched_seqs = {key: [] for key in barcodes_dict.values()}

#------------------------------------------------------------------------------------------------------------------------------------
# Parse fastq file and subset sequences containing valid barcodes. Generate feature reference csv file
#------------------------------------------------------------------------------------------------------------------------------------

# Open the compressed fastq file & output file
with gzip.open(input_file, 'rt') as input_handle:
    # Iterate over each record in the fastq file
    for record in SeqIO.parse(input_handle, 'fastq'):
        # Look for barcode patterns
        for pattern in patterns:
            match = pattern.search(str(record.seq))
            if match:
                # Update record to contain just matched sequence 
                new_record = record[match.start() : match.end()]
                # Save updated records to dictionary. Every key is a different barcode color, the content are all the records
                # corresponding to that color.
                matched_seqs[patterns_dict[pattern]] += new_record,
                break

i = 0 # Order of output files is the same as the keys of the dictionary since both of them are taken from the same config variable
for key in matched_seqs.keys():
    with gzip.open(output_fastqs[i], 'wt') as output_handle:
        count = SeqIO.write(matched_seqs[key], output_handle, 'fastq')
    i += 1