import os

# Transform input paths to absolute in case they are relative
input = os.path.abspath(snakemake.input)

# Create symbolic links
os.symlink(input, snakemake.output)