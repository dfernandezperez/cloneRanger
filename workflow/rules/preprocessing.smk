rule create_seurat:
    input:
        expand("results/01_counts/{sample}/outs/filtered_feature_bc_matrix/", sample = SAMPLES)
    output:
        no_doublets = "results/02_createSeurat/seurat_noDoublets.rds",
        raw         = "results/02_createSeurat/seurat_allCells.rds",
    params:
        sample_names = SAMPLES,
        min_cells    = config["preprocessing"]["min_cells"]
    conda:
        "R_singlecell"
    threads:
        RESOURCES["create_seurat"]["cpu"]
    resources:
        mem_mb = RESOURCES["create_seurat"]["MaxMem"]
    log:
        "results/00_logs/create_seurat/log"
    benchmark:
        "results/benchmarks/create_seurat/benchmark.txt"
    script:
        "../scripts/create_seurat.R"