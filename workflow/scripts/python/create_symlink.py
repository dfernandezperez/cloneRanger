import os

# output fastq files
output_fw = snakemake.output[0]
output_rv = snakemake.output[1]

# Transform input paths to absolute in case they are relative
input_fw = os.path.abspath(snakemake.input[0])
input_rv = os.path.abspath(snakemake.input[1])

# Create symbolic links
os.symlink(input_fw, output_fw)
os.symlink(input_rv, output_rv)