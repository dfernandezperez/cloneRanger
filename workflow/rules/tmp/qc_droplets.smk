rule filter_empty:
    input:
        mtx=rules.counts.output.mtx,
    output:
        knee=report(
            "results/filter_empty/{sample}/knee_plot.png",
            caption="../reports/knee_plot.rst",
            category="Filter Empty Droplets",
            subcategory="{sample}",
        ),
        tail=report(
            "results/filter_empty/{sample}/tail_end_plot.png",
            caption="../reports/tail_end_plot.rst",
            category="Filter Empty Droplets",
            subcategory="{sample}",
        ),
        empty="results/filter_empty/{sample}/empty.csv",
    params:
        niters=config["filter_empty"]["niters"],
    log:
        "results/logs/filter_empty/{sample}.log",
    benchmark:
        "results/benchmarks/filter_empty/{sample}.txt"
    threads: 1
    conda:
        "../envs/filter_empty.yaml"
    script:
        "../scripts/filter_empty.R"


rule qc:
    input:
        dir=rules.counts.output.results,
        empty=rules.filter_empty.output.empty,
    params:
        introns=config["counts"]["introns"],
        pct_counts_mt=config["qc"]["pct_counts_mt"],
        total_counts=config["qc"]["total_counts"],
        n_genes_by_counts=config["qc"]["n_genes_by_counts"],
    output:
        qc_violin=report(
            "results/qc/{sample}/qc_violin.png",
            caption="../reports/qc_violin.rst",
            category="Cell Quality Control",
            subcategory="{sample}",
        ),
        mt_scatter=report(
            "results/qc/{sample}/mt_scatter.png",
            caption="../reports/qc_scatter_mt.rst",
            category="Cell Quality Control",
            subcategory="{sample}",
        ),
        n_genes=report(
            "results/qc/{sample}/n_gene_scatter.png",
            caption="../reports/qc_scatter_ngene.rst",
            category="Cell Quality Control",
            subcategory="{sample}",
        ),
        n_top=report(
            "results/qc/{sample}/n_top_genes.png",
            caption="../reports/qc_top_genes.rst",
            category="Cell Quality Control",
            subcategory="{sample}",
        ),
        data="results/qc/{sample}/adata.h5ad",
    log:
        "results/logs/qc/{sample}.log",
    benchmark:
        "results/benchmarks/qc/{sample}.txt"
    threads: 1
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/qc.py"


rule doublets:
    input:
        unpack(get_qc_data),
    output:
        table=report(
            "results/doublets/{lane}/table.png",
            caption="../reports/doublets_table.rst",
            category="Find Doublets",
            subcategory="{lane}",
        ),
        data="results/doublets/{lane}/adata.h5ad",
    log:
        "results/logs/doublets/{lane}.log",
    benchmark:
        "results/benchmarks/doublets/{lane}.txt"
    threads: 1
    resources:
        ndivia_gpu=1,
    conda:
        "../envs/doublets.yaml"
    script:
        "../scripts/doublets.py"
