rule create_seurat:
    input:
        expand("results/01_counts/{sample}/outs/filtered_feature_bc_matrix/", sample = SAMPLES)
    output:
        no_doublets = "results/02_createSeurat/seurat_noDoublets.rds",
        raw         = "results/02_createSeurat/seurat_allCells.rds"
    params:
        sample_names = SAMPLES,
        min_cells    = config["preprocessing"]["min_cells"],
        mito_pattern = config["preprocessing"]["mito_pattern"],
        ribo_pattern = config["preprocessing"]["ribo_pattern"],
    conda:
         "../envs/Seurat.yaml"
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

rule qc:
    input:
        "results/02_createSeurat/seurat_noDoublets.rds"
    output:
        html          = "results/03_qc/seurat_sneakPeak.html",
        cluster_degs  = "results/03_qc/cluster_degs.tsv",
        sample_degs   = "results/03_qc/sample_degs.tsv"
    params:
        marker_genes = config["preprocessing"]["marker_genes"],
        species      = config["species"]
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["create_seurat"]["cpu"]
    resources:
        mem_mb = RESOURCES["create_seurat"]["MaxMem"]
    log:
        "results/00_logs/qc/log"
    benchmark:
        "results/benchmarks/qc/benchmark.txt"
    script:
        "../scripts/mini_report.Rmd"