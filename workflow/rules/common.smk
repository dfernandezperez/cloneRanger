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
# def get_fastqs(wildcards):
#     """Return all FASTQS specified in sample metadata."""
#     return (metadata[wildcards.lane][wildcards.sample][wildcards.read],)

def get_fastqs(wildcards):
    """Return all FASTQS specified in sample metadata."""
    return samples.loc[(wildcards.sample, wildcards.lane), ["R1", "R2"]].dropna()

def get_bcl_results():
    """Retrieve the bcl results files."""
    bcls = list(Path("data").glob("*.bcl2fastq.zip"))
    return dict(
        bcl_zip=bcls,
    )


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
    if config["counts"]["introns"]:
        return "--include-introns"
    else:
        return ""


def get_sample_reads(wildcards):
    """Get all reads for Cellranger.

    Cellranger doesn't actually need the file names,
    but we need to mark them as dependencies and introduce the sample wildcard.
    Given we know the lane (metadata.key()) and the read (R1 or R2),
    we use the following expand statement to do so.
    """
    files = expand(
        "data/{{sample}}_S1_L00{lane}_{read}_001.fastq.gz",
        lane=metadata.keys(),
        read=["R1", "R2"],
    )
    return dict(
        fastqs=files,
    )


def get_qc_data(wildcards):
    """Get all QC'd data, grouped by lane, for SOLO."""
    keys_only = {k: list(v.keys()) for k, v in metadata.items()}
    keys_only = {
        k: [f"results/qc/{val}/adata.h5ad" for val in v] for k, v in keys_only.items()
    }
    return dict(data=keys_only[wildcards.lane])
