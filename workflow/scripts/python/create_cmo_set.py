def generate_cmo_set_csv(cell_hashing_dict):
    """Generate cmo-set csv file indicating the totalseq names and barcode sequences + positions"""
    cmo_csv = "id,name,read,pattern,sequence,feature_type\n"
    for key, value in cell_hashing_dict.items():
        cmo_csv += f"{key},{key},{value[0]},{value[1]},{value[2]},Multiplexing Capture\n"
        
    return cmo_csv


cmo_csv = generate_cmo_set_csv(snakemake.params["cell_hashing"])

with open(snakemake.output[0], 'wt') as output_handle:
    output_handle.write(cmo_csv)