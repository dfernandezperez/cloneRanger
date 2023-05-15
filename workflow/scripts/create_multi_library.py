import os

#-----------------------------------------------------------------------------------------------------------------------
# Read parameters
#-----------------------------------------------------------------------------------------------------------------------
genome            = snakemake.params["genome"]
library_types     = snakemake.params["library_types"]
cell_hashing      = snakemake.params["cell_hashing"]
params_cellranger = snakemake.params["cellranger_params"]
sample            = snakemake.wildcards["sample"]
abs_path          = os.getcwd()

if hasattr(snakemake.input, 'larry'):
    larry_ref = snakemake.input["larry"]
else:
    larry_ref = False

if hasattr(snakemake.input, 'cmo_set'):
    cmo_set = snakemake.input["cmo_set"]
else:
    cmo_set = False

#-----------------------------------------------------------------------------------------------------------------------
# Declare functions
#-----------------------------------------------------------------------------------------------------------------------
def samples_cell_hashing(cell_hashing_dict, sample, cmo_set):
    """Adapt hashing-id : sample pairs to csv format"""
    if not cmo_set:
        return ""
    else:
        samples_cell_hashing_csv = ""
        sample_assignments = cell_hashing_dict[sample]

        for key, value in sample_assignments.items():
            samples_cell_hashing_csv += f"{value},{key}\n"

        params_samples_cell_hashing_csv = f"""
[samples] # for Cell Multiplexing libraries only
sample_id,cmo_ids
{samples_cell_hashing_csv}
        """

        return params_samples_cell_hashing_csv


def samples_larry(larry_ref, abs_path):
    """Prepare larry barcode fastq paths to cellranger multi csv format"""

    if not larry_ref:
        return ""
    else:
        params_samples_larry = f"""
[feature] # For Feature Barcode libraries only
reference,{abs_path}/{larry_ref}
        """

        return params_samples_larry


def samples_gene_expression(genome, abs_path, cmo_set, params_cellranger_csv):

    if not cmo_set:
        params_samples_gex = f"""
[gene-expression]
reference,{genome}
{params_cellranger_csv}
    """
    else:
        params_samples_gex = f"""
[gene-expression]
reference,{genome}
cmo-set,{abs_path}/{cmo_set}
{params_cellranger_csv}
    """
        
    return params_samples_gex


def samples_libraries(library_types):
    params_samples_libraries = f"""
[libraries]
fastq_id,fastqs,feature_types
{library_types}
    """

    return params_samples_libraries


def params_cellranger_to_csv(params_cellranger):
    """Adapt cellranger params to csv format"""
    params_cellranger_csv = ""
    for param in params_cellranger:
        params_cellranger_csv += f"{param}\n"

    return params_cellranger_csv


def generate_multiconfig_csv(params_gene_expression, params_larry, params_libraries, params_cellhasg):

    out_csv = f"""
{params_gene_expression}
{params_larry}
{params_libraries}
{params_cellhasg}
    """

    return out_csv


#-----------------------------------------------------------------------------------------------------------------------
# Main
#-----------------------------------------------------------------------------------------------------------------------

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    params_cellranger_csv    = params_cellranger_to_csv(params_cellranger)

    gex_config = samples_gene_expression(
        genome                = genome,
        abs_path              = abs_path,
        params_cellranger_csv = params_cellranger_csv,
        cmo_set               = cmo_set
    )

    cellhashing_config = samples_cell_hashing(
        cell_hashing_dict = cell_hashing,
        sample            = sample,
        cmo_set           = cmo_set
    )

    larry_config = samples_larry(
        larry_ref=larry_ref,
        abs_path=abs_path
    )

    samples_config = samples_libraries(
        library_types = library_types
    )

    output_csv               = generate_multiconfig_csv(params_gene_expression = gex_config, 
                                                        params_larry     = larry_config,
                                                        params_libraries = samples_config,
                                                        params_cellhasg  = cellhashing_config
                                                        )

    with open(snakemake.output[0], 'wt') as output_handle:
        output_handle.write(output_csv)

