
import os

"""Create the content of library.csv for the feature barcoding and multiome pipeline
from cellranger.

Based on what is written in the column lib_type from the units file,
write the following: "fastq folder,fastq name,library type".
"""

abs_path        = os.getcwd()
sample          = snakemake.wildcards["sample"]
is_feature_bc   = snakemake.params["is_feature_bc"]
is_cell_hashing = snakemake.params["is_cell_hashing"]
lib_type        = snakemake.params["library_type"]

lib_csv = dict()


########################################################################################################
# Create library dict
########################################################################################################
if is_feature_bc:
    lib_csv['FB'] = f'{abs_path}/data/clean,{sample}_FB,Custom'

if is_cell_hashing:
    lib_csv['CH'] = f'{abs_path}/data/clean,{sample}_CH,Antibody Capture'

if lib_type == 'GEX':
    lib_csv['GEX'] = f'{abs_path}/data/clean,{sample}_GEX,Gene Expression'

if lib_type == 'ATAC':
    lib_csv['ATAC'] = f'{abs_path}/data/clean,{sample}_ATAC,Chromatin Accessibility'

if lib_type == 'ARC' and not is_feature_bc:
    lib_csv['ATAC'] = f'{abs_path}/data/clean,{sample}_ATAC,Chromatin Accessibility'
    lib_csv['GEX']  = f'{abs_path}/data/clean,{sample}_GEX,Gene Expression'
    
if lib_type == 'ARC' and is_feature_bc:
    lib_csv['GEX']  = f'{abs_path}/data/clean,{sample}_GEX,Gene Expression'

# Save library dict to csv
out_csv = '\n'.join(lib_csv.values())
with open(snakemake.output["library"], 'wt') as output_handle:
    output_handle.write("fastqs,sample,library_type\n")
    output_handle.write(out_csv)


########################################################################################################
# Create library dict for arc
########################################################################################################
# If the library is multiome and there is feature barcode data, the pipeline must be run 2 times
# One for rna + larry and the other for rna+atac. Due to this, we need 2 different library csv files.
if lib_type == 'ARC' and is_feature_bc:

    lib_csv2 = dict()
    lib_csv2['ATAC'] = f'{abs_path}/data/clean,{sample}_ATAC,Chromatin Accessibility'
    lib_csv2['GEX']  = f'{abs_path}/data/clean,{sample}_GEX,Gene Expression'

    # Save library (arc) dict to csv
    out_csv2 = '\n'.join(lib_csv2.values())
    with open(snakemake.output["library_arc"], 'wt') as output_handle:
        output_handle.write("fastqs,sample,library_type\n")
        output_handle.write(out_csv2)

else:
    lib_csv2 = "If library type is not ARC, ignore this file."
    with open(snakemake.output["library_arc"], 'wt') as output_handle:
        output_handle.write(lib_csv2)
