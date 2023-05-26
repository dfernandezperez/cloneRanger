def generate_cellhash_ref_csv(cell_hashing_dict):
    """Generate cmo-set csv file indicating the totalseq names and barcode sequences + positions"""
    cmo_csv = "id,name,read,pattern,sequence,feature_type\n"
    for key, value in cell_hashing_dict.items():
        cmo_csv += f"{key},{key},{value[0]},{value[1]},{value[2]},Antibody Capture\n"
        
    return cmo_csv

cell_hashing_Abs  = snakemake.params["cellhash_ab_names"]
cell_hashing_dict = snakemake.params["cell_hashing"]
cell_hashing_dict = dict((Ab, cell_hashing_dict[Ab]) for Ab in cell_hashing_Abs) # Use just cellhash Abs corresponding to this sample

cellhash_csv = generate_cellhash_ref_csv(cell_hashing_dict)

with open(snakemake.output[0], 'wt') as output_handle:
    output_handle.write(cellhash_csv)