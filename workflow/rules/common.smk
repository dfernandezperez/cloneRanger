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
    if config["10x_pipeline"] == "GEX":
        if not config["cellranger_count"]["introns"]:
            return "--include-introns False"
        else:
            return ""

    elif config["10x_pipeline"] == "ARC":
        if not config["cellranger_count"]["introns"]:
            return "--gex-exclude-introns"
        else:
            return ""
    
    else: 
        sys.exit("Intronic use can be only specified if 10x_pipeline == GEX or GEX_ATAC")


def get_expected_cells():
    """Function to set the number of expected cells if it was set by the user,
    otherwise estimate them automatically by cellranger.
    """
    if config["cellranger_count"]["n_cells"] == "auto":
        return ""
    else:
        return f'--expect-cells {config["cellranger_count"]["n_cells"]}'


def is_feature_bc():
    """Specify whether feature barcoding has been performed
    or not
    """
    if config["feature_bc_config"]["feature_bc"]:
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


def get_library_type(wildcards):
    """Create the content of library.csv for the feature barcoding and multiome pipeline
    from cellranger.
    
    Based on what is written in the column lib_type from the units file,
    write the following: "fastq folder,fastq name,library type".
    """
    abs_path     = os.getcwd()
    lib_type_col = set( units.loc[(wildcards.sample), "lib_type"] )
    lib_types    = dict()

    for lib_type in lib_type_col:
        if lib_type == 'GEX':
            lib_types[lib_type] = f'{wildcards.sample}_{lib_type},{abs_path}/data/clean,Gene Expression'
        elif lib_type == 'ATAC':
            lib_types[lib_type] = f'{wildcards.sample}_{lib_type},{abs_path}/data/clean,Chromatin Accessibility'
        # Use and if statement for larry and cellhashing data, in case we want to analzye the libs
        # without some of them despite having the larry/ch fastq files in units.tsv.  
        elif lib_type == 'FB':
            if not is_feature_bc():
                next
            else:
                lib_types[lib_type] = f'{wildcards.sample}_{lib_type},{abs_path}/data/clean,Custom'
        elif lib_type == 'CH':
            if not is_cell_hashing(wildcards.sample):
                next
            else:
                lib_types[lib_type] = f'{wildcards.sample}_{lib_type},{abs_path}/data/clean,Multiplexing Capture'
        else:
            sys.exit("Library type must be GEX, FB, ATAC or CH.")

    return '\n'.join(lib_types.values())


def get_library_input(wildcards):
    if is_feature_bc() and is_cell_hashing(wildcards.sample):
        return {
            "fq"      : expand(
                                "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                                sample   = wildcards.sample,
                                lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                                read     = ["R1", "R2"]
                            ),
            "larry"   : "data/feature_reference/Feature_reference.csv",
            "cmo_set" : "data/feature_bc_libraries/cmo-set.csv"
        }
    elif is_feature_bc() and not is_cell_hashing(wildcards.sample):
        return {
            "fq"      : expand(
                                "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                                sample   = wildcards.sample,
                                lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                                read     = ["R1", "R2"]
                            ),
            "larry"   : "data/feature_reference/Feature_reference.csv"
        }
    elif not is_feature_bc() and not is_cell_hashing(wildcards.sample):
        return {
            "fq"      : expand(
                                "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                                sample   = wildcards.sample,
                                lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                                read     = ["R1", "R2"]
                            ),
        }
    elif not is_feature_bc() and is_cell_hashing(wildcards.sample):
        return {
            "fq"      : expand(
                                "data/clean/{sample}_{lib_type}_S1_L001_{read}_001.fastq.gz", 
                                sample   = wildcards.sample,
                                lib_type = SAMPLE_LIB_DICT[wildcards.sample],
                                read     = ["R1", "R2"]
                            ),
            "cmo_set" : "data/feature_bc_libraries/cmo-set.csv"
        }


def get_sample_matrices():

    sample_paths = expand(
                "results/01_counts/{sample}/outs/per_sample_outs/{subsample}/count/sample_filtered_feature_bc_matrix/",
                zip,
                sample    = SAMPLES_LONG,
                subsample = SUBSAMPLES
                )
    
    return(sample_paths)