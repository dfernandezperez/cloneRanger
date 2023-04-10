from snakemake.utils import validate
from pathlib import Path
from yaml import safe_load
import os 
import sys

# # Validate structures
# with open(config["samples"], "r") as file:
#     metadata = safe_load(file)
# validate(metadata, schema="../schema/samples.schema.yaml")
# validate(config, schema="../schema/config.schema.yaml")

# # Some useful variables
# LANES = list(metadata.keys())
# SAMPLES = set()
# for _, v in metadata.items():
#     SAMPLES = SAMPLES.union(v.keys())

# wildcard constraints to be sure that output names are defined properly
wildcard_constraints:
    # hamming distance can just be a value
    hd="\d+" ,
    # reads corresponding to cellular barcode and feature barcode can only be R1 or R2
    read_fb="R[12]",
    read_cb="R[12]",
    # lib type can be just ATAC, GEX or FB
    lib_type='|'.join([x for x in ["ATAC", "GEX", "FB"]])


# Input functions
def get_fastqs(wildcards):
    """Return all FASTQS specified in sample metadata."""
    return units.loc[(wildcards.sample, wildcards.lib_type, wildcards.lane), ["R1", "R2"]].dropna()

def convert_introns():
    """Specify whether introns should be counted

    For ease of use, the user only specifies True or False
    This function handles the conversion.
    """
    if config["10x_pipeline"] == "GEX":
        if not config["cellranger_count"]["introns"]:
            return "--include-introns False"
        else:
            return ""

    elif config["10x_pipeline"] == "GEX_ATAC":
        if not config["cellranger_count"]["introns"]:
            return "--gex-exclude-introns"
        else:
            return ""
    
    else: 
        sys.exit("Intronic use can be only specified if 10x_pipeline == GEX or GEX_ATAC")


def is_feature_bc():
    """Specify whether feature barcoding has been performed
    or not
    """
    if config["feature_bc_config"]["feature_bc"]:
        return True
    else:
        return False

def get_library_type(wildcards):
    """Create the content of library.csv for the feature barcoding and multiome pipeline
    from cellranger.
    
    Based on what is written in the column lib_type from the units file,
    write the following: "fastq folder,fastq name,library type".
    """
    abs_path            = os.getcwd()
    feature_barcode_col = set( units.loc[(wildcards.sample), "lib_type"] )
    lib_types           = dict()

    for fb in feature_barcode_col:
        if fb == 'GEX':
            lib_types[fb] = f'{abs_path}/data/clean,{wildcards.sample}_{fb},Gene Expression'
        if fb == 'FB':
            lib_types[fb] = f'{abs_path}/data/clean,{wildcards.sample}_{fb},Custom'
        if fb == 'ATAC':
            lib_types[fb] = f'{abs_path}/data/clean,{wildcards.sample}_{fb},Chromatin Accessibility'

    fb_names = [lib_types.get(fb) for fb in feature_barcode_col]
    return '\n'.join(fb_names)


def agg_fastqc():
    """Aggregate input from FASTQC for MultiQC.

    A simple snakemake expand is not used to avoid the cross product,
    as not all samples will be in all lanes.

    A checkpoint could be used. However, wildcard constraints and the sample sheet
    guarantee we can know this a priori.
    """
    fastq_out = []
    for lane, samples in metadata.items():
        for sample in samples.keys():
            fastq_out = fastq_out + [
                f"results/fastqc/{lane}_{sample}_{read}_fastqc.zip"
                for read in ["R1", "R2"]
            ]
    return dict(
        fastqc=fastq_out,
    )


def get_qc_data(wildcards):
    """Get all QC'd data, grouped by lane, for SOLO."""
    keys_only = {k: list(v.keys()) for k, v in metadata.items()}
    keys_only = {
        k: [f"results/qc/{val}/adata.h5ad" for val in v] for k, v in keys_only.items()
    }
    return dict(data=keys_only[wildcards.lane])
