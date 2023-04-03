from snakemake.utils import validate
from pathlib import Path
from yaml import safe_load

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


# # wildcard restraints - as reads can only be R1 or R2
# wildcard_constraints:
#     read="R[12]",


# Input functions
def get_fastqs(wildcards):
    """Return all FASTQS specified in sample metadata."""
    return units.loc[(wildcards.sample, wildcards.feature_bc, wildcards.lane), ["R1", "R2"]].dropna()

def get_library_type(wildcards):
    """Create the content of library.csv for the feature barcoding pipeline
    from cellranger.
    
    Based on what is written in the column feature_bc from the units file,
    write the following: "fastq folder,fastq name,library type".
    """
    feature_barcode_col = set( units.loc[(wildcards.sample), "feature_bc"] )
    lib_types = {
        'GEX': 'data,' + wildcards.sample + '_GEX,Gene Expression', 
        'FB' : 'data,' + wildcards.sample + '_FB,Custom'
    }

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


def convert_introns():
    """Specify whether introns should be counted

    For ease of use, the user only specifies True or False
    This function handles the conversion.
    """
    if not config["cellranger_count"]["introns"]:
        return "--include-introns False"
    else:
        return ""

def is_feature_bc(wildcards):
    """Specify whether feature barcoding has been performed
    or not
    """
    if not config["feature_bc_config"]["feature_bc"]:
        return "works"
        # return "--libraries=${library} --feature-ref=${feature_ref}"
    else:
        return ""

def get_qc_data(wildcards):
    """Get all QC'd data, grouped by lane, for SOLO."""
    keys_only = {k: list(v.keys()) for k, v in metadata.items()}
    keys_only = {
        k: [f"results/qc/{val}/adata.h5ad" for val in v] for k, v in keys_only.items()
    }
    return dict(data=keys_only[wildcards.lane])
