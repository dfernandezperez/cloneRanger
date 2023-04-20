rule create_seurat:
    input:
        expand("results/01_counts/{sample}/outs/filtered_feature_bc_matrix/", sample = SAMPLES)
    output:
        no_doublets = "results/02_createSeurat/seurat_noDoublets.rds",
        raw         = "results/02_createSeurat/seurat_allCells.rds"
    params:
        sample_names    = SAMPLES,
        min_cells_gene  = config["preprocessing"]["min_cells_gene"],
        min_cells_larry = config["preprocessing"]["min_cells_larry"],
        mito_pattern    = config["preprocessing"]["mito_pattern"],
        ribo_pattern    = config["preprocessing"]["ribo_pattern"],
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

rule RNA_exploration:
    input:
        "results/02_createSeurat/seurat_noDoublets.rds"
    output:
        html = "results/03_RNA-exploration/RNA_exploration.html"
    params:
        marker_genes  = config["preprocessing"]["marker_genes"],
        species       = config["species"],
        cluster_degs  = "results/03_RNA-exploration/cluster_degs.tsv",
        sample_degs   = "results/03_RNA-exploration/sample_degs.tsv"
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["create_seurat"]["cpu"]
    resources:
        mem_mb = RESOURCES["create_seurat"]["MaxMem"]
    log:
        "results/00_logs/RNA_exploration/log"
    benchmark:
        "results/benchmarks/RNA_exploration/benchmark.txt"
    script:
        "../scripts/RNA_exploration.Rmd"

rule barcode_exploration:
    input:
        "results/02_createSeurat/seurat_noDoublets.rds"
    output:
        html = "results/04_barcode-exploration/barcode_exploration.html"
    params:
        barcodes = LARRY_COLORS
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
        "../scripts/barcode_exploration.Rmd"