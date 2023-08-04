from pathlib import Path
import os 
import sys

# wildcard constraints to be sure that output names are defined properly
wildcard_constraints:
    # hamming distance can just be a value
    hd="\d+" ,
    # reads corresponding to cellular barcode and feature barcode can only be R1 or R2
    read_fb="R[12]",
    read_cb="R[12]",
    # lib type can be just ATAC, GEX or FB
    lib_type='|'.join([x for x in ["ATAC", "GEX", "FB", "CH"]])


# Input functions
def get_fastqs(wildcards):
    """Return all FASTQS specified in sample metadata."""
    return units.loc[(wildcards.sample, wildcards.lib_type, wildcards.lane), ["R1", "R2"]].dropna()


def get_atac_fastqs(wildcards):
    """Return all FASTQS specified in sample metadata."""
    return units.loc[(wildcards.sample, wildcards.lib_type, wildcards.lane), ["R3"]].dropna()


def convert_introns():
    """Specify whether introns should be counted

    For ease of use, the user only specifies True or False
    This function handles the conversion.
    """
    if config["cellranger_count"]["10x_pipeline"] == "GEX":
        if not config["cellranger_count"]["introns"]:
            return "--include-introns False"
        else:
            return ""

    elif config["cellranger_count"]["10x_pipeline"] == "ARC":
        if not config["cellranger_count"]["introns"]:
            return "--gex-exclude-introns"
        else:
            return ""
    
    else: 
        sys.exit("Intronic use can be only specified if 10x_pipeline == GEX or GEX_ATAC")


def get_expected_cells(wildcards):
    """Function to set the number of expected cells if it was set by the user,
    otherwise estimate them automatically by cellranger.
    """
    if sample_config["UMI_cutoff"][wildcards.sample] == "auto":
        return ""
    else:
        return f'--expect-cells {sample_config["UMI_cutoff"][wildcards.sample]}'


def get_feature_ref(wildcards):
    if is_feature_bc() or is_cell_hashing(wildcards.sample):
        return f"--feature-ref data/feature_reference/Feature-reference_{wildcards.sample}.csv"        
    else:
         return ""


def is_feature_bc():
    """Specify whether feature barcoding has been performed
    or not
    """
    if LARRY["feature_bc"]:
        return True
    else:
        return False


def is_cell_hashing(sample):
    """Specify if a sample has been processed with cell hashing (like totalseq)
    and should be processed with the specific cellranger multi pipeline
    """
    if CELL_HASHING["assignments"] is None:
        return False
    elif sample in CELL_HASHING["assignments"]:
        return True
    else:
        return False


def get_seurat_rds(wildcards):
    if is_feature_bc():
        return expand("results/02_createSeurat/seurat_{sample}_noDoublets-larry-filt.rds", sample = SAMPLES)
    else:
        return expand("results/02_createSeurat/seurat_{sample}_noDoublets.rds", sample = SAMPLES)


def get_cellranger_input(wildcards):
    if is_feature_bc() or is_cell_hashing(wildcards.sample):
        return {
            "fq" : expand(
                        "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                        sample   = wildcards.sample,
                        lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                        read     = ["R1", "R2"]
                        ),
            "libraries"   : "data/feature_bc_libraries/{sample}_library.csv",
            "feature_ref" : "data/feature_reference/Feature-reference_{sample}.csv",            
        }
    else:
         return {
            "fq" : expand(
                        "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                        sample   = wildcards.sample,
                        lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                        read     = ["R1", "R2"]
                        ),
            "libraries"   : "data/feature_bc_libraries/{sample}_library.csv"
        }

def get_cellranger_output(wildcards):
    if is_feature_bc() and config["cellranger_count"]["10x_pipeline"] == "ARC":
        return {
            "arc"   : "results/01_arc/{sample}/outs/filtered_feature_bc_matrix",
            "counts" : "results/01_counts/{sample}/outs/filtered_feature_bc_matrix",
            "fragments" : "results/01_arc/{sample}/outs/atac_fragments.tsv.gz",
        }

    elif not is_feature_bc() and config["cellranger_count"]["10x_pipeline"] == "ARC":
        return {
            "arc" : "results/01_arc/{sample}/outs/filtered_feature_bc_matrix",
            "counts" : "results/01_arc/{sample}/outs/filtered_feature_bc_matrix", # repeat arc input just to have the count variable
            "fragments" : "results/01_arc/{sample}/outs/atac_fragments.tsv.gz",

        }

    else:
        return { 
            "counts" : "results/01_counts/{sample}/outs/filtered_feature_bc_matrix" 
        }
        

def get_feature_ref_input(wildcards):
    if is_feature_bc():
    
        if is_cell_hashing(wildcards.sample):

            return {
                "cell_hash_ref" : "data/cellhashing/cellhashing-reference_{sample}.csv",
                "larry_ref"     : "data/feature_reference/Feature_reference_larry.csv"
            }     

        else:

            return {
                "larry_ref" : "data/feature_reference/Feature_reference_larry.csv"
            }

    else:
        # If no larry/cellhashing is performed, this file has
        return {
            "cell_hash_ref" : "data/cellhashing/cellhashing-reference_{sample}.csv"
        }