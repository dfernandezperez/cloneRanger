rule create_seurat:
    input:
        unpack(get_cellranger_output)
    output:
        temp("results/02_createSeurat/tmp_seurat_{sample}.rds")
    params:
        is_cell_hashing = lambda w: "TRUE" if is_cell_hashing(w.sample) else "FALSE",
        cellhash_names  = lambda w: CELL_HASHING["assignments"][w.sample] if is_cell_hashing(w.sample) else "FALSE",
        is_larry        = "TRUE" if is_feature_bc() else "FALSE",
        min_cells_gene  = config["preprocessing"]["min_cells_gene"],
        umi_cutoff      = lambda w: sample_config["UMI_cutoff"][w.sample],
        mito_pattern    = config["preprocessing"]["mito_pattern"],
        ribo_pattern    = config["preprocessing"]["ribo_pattern"],
        library_type    = config["cellranger_count"]["10x_pipeline"]
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["create_seurat"]["cpu"]
    resources:
        mem_mb = RESOURCES["create_seurat"]["MaxMem"]
    log:
        "results/00_logs/create_seurat/{sample}.log"
    benchmark:
        "results/benchmarks/create_seurat/{sample}_benchmark.txt"
    script:
        "../scripts/R/create_seurat.R"


# If there is larry filter the matrix, otherwise just update .rds filenames.
if is_feature_bc():
    rule barcode_filtering:
        input:
            "results/02_createSeurat/tmp_seurat_{sample}.rds"
        output:
            "results/02_createSeurat/seurat_{sample}.rds",
        params:
            reads_cutoff  = LARRY["reads_cutoff"],
            umi_cutoff    = LARRY["umi_cutoff"],
            molecule_info = lambda w: f"results/01_counts/{w.sample}/outs/molecule_info.h5",
        conda:
            "../envs/Seurat.yaml"
        threads:
            RESOURCES["barcode_filtering"]["cpu"]
        resources:
            mem_mb = RESOURCES["barcode_filtering"]["MaxMem"]
        log:
            "results/00_logs/barcode_filtering/{sample}.log"
        benchmark:
            "results/benchmarks/barcode_filtering/{sample}_benchmark.txt"
        script:
            "../scripts/R/barcode_filtering.R"
else:
    rule update_rds_name:
        input:
            "results/02_createSeurat/tmp_seurat_{sample}.rds"
        output:
            "results/02_createSeurat/seurat_{sample}.rds",
        log:
            "results/00_logs/update_rds_name/{sample}.log"
        benchmark:
            "results/benchmarks/update_rds_name/{sample}_benchmark.txt"
        shell:
            """
            mv {input} {output}
            """


rule merge_seurat:
    input:
        expand("results/02_createSeurat/seurat_{sample}.rds", sample = SAMPLES)
    output:
        seurat = "results/03_mergeSeurat/seurat_merged.rds",
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["merge_seurat"]["cpu"]
    resources:
        mem_mb = RESOURCES["merge_seurat"]["MaxMem"]
    log:
        "results/00_logs/merge_seurat/log"
    benchmark:
        "results/benchmarks/merge_seurat/benchmark.txt"
    script:
        "../scripts/R/merge_seurat.R"


rule RNA_exploration:
    input:
        "results/03_mergeSeurat/seurat_merged.rds"
    output:
        html = "results/04_RNA-exploration/RNA_exploration.html"
    params:
        marker_genes  = config["preprocessing"]["marker_genes"],
        species       = config["preprocessing"]["species"],
        cluster_degs  = "results/04_RNA-exploration/cluster_degs.tsv",
        sample_degs   = "results/04_RNA-exploration/sample_degs.tsv"
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["RNA_exploration"]["cpu"]
    resources:
        mem_mb = RESOURCES["RNA_exploration"]["MaxMem"]
    log:
        "results/00_logs/RNA_exploration/log"
    benchmark:
        "results/benchmarks/RNA_exploration/benchmark.txt"
    script:
        "../scripts/R/RNA_exploration.Rmd"


rule barcode_summary:
    input:
       "results/02_createSeurat/seurat_{sample}.rds"
    output:
        html = "results/05_barcode-exploration/{sample}_barcode-summary.html"
    params:
        molecule_info = lambda w: f"results/01_counts/{w.sample}/outs/molecule_info.h5",
    conda:
         "../envs/Seurat.yaml"
    threads:
        RESOURCES["barcode_filtering"]["cpu"]
    resources:
        mem_mb = RESOURCES["barcode_filtering"]["MaxMem"]
    log:
        "results/00_logs/barcode_filtering/{sample}.log"
    benchmark:
        "results/benchmarks/barcode_filtering/{sample}_benchmark.txt"
    script:
        "../scripts/R/barcode_summary.Rmd"