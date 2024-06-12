rule create_library:
    output:
        library     = "data/feature_bc_libraries/{sample}_library.csv",
        library_arc = "data/feature_bc_libraries/{sample}_library_arc.csv",
    params:
        library_type    = config["cellranger_count"]["10x_pipeline"],
        is_feature_bc   = is_feature_bc(),
        is_cell_hashing = lambda w: is_cell_hashing(w.sample)
    script:
        "../scripts/python/create_library_csv.py"


rule create_cellhashing_ref:
    output:
        "data/cellhashing/cellhashing-reference_{sample}.csv"
    params:
        cell_hashing       = CELL_HASHING["barcodes"],
        cellhash_ab_names  = lambda w: CELL_HASHING["assignments"][w.sample]
    log:
        "results/00_logs/create_cellhashing_ref/{sample}.log"
    script:
        "../scripts/python/generate_cellhashing_ref.py"


rule generate_feature_ref:
    input:
        unpack(get_feature_ref_input)
    output:
        "data/feature_reference/Feature-reference_{sample}.csv"
    log:
        "results/00_logs/generate_feature_ref/{sample}.log"
    conda:
        "../envs/Seurat.yaml"
    script:
        "../scripts/R/generate_feature_ref.R"