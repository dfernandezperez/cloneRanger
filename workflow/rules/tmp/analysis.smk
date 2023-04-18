rule dim_reduc:
    input:
        data=expand("results/doublets/{lane}/adata.h5ad", lane=LANES),
    params:
        nHVG=config.get("dim_reduc", {}).get("nHVG"),
        var_thresh=config.get("dim_reduc", {}).get("var_thresh", 0.85),
    output:
        data="results/dim_reduc/adata.h5ad",
        hvg=report(
            "results/dim_reduc/hvgs.png",
            caption="../reports/dr_hvg.rst",
            category="Dimensional Reduction",
        ),
        elbow=report(
            "results/dim_reduc/elbow.png",
            caption="../reports/dr_elbow.rst",
            category="Dimensional Reduction",
        ),
    log:
        "results/logs/dim_reduc/dim_reduc.log",
    benchmark:
        "results/benchmarks/dim_reduc/dim_reduc.txt"
    threads: 16
    conda:
        "../envs/dim_reduc.yaml"
    script:
        "../scripts/dim_reduc.py"


rule cluster:
    input:
        data=rules.dim_reduc.output.data,
    params:
        res=config["cluster"]["res"],
        markers=config["cluster"].get("markers", []),
    output:
        data="results/cluster/adata.h5ad",
        umap=report(
            "results/cluster/umap.png",
            caption="../reports/umap.rst",
            category="Clustering",
        ),
        markers=report(
            "results/cluster/markers.pdf",
            caption="../reports/markers.rst",
            category="Clustering",
        ),
        dot=report(
            "results/cluster/dotplot.png",
            caption="../reports/dotplot.rst",
            category="Differential Expression",
        ),
        de=report(
            "results/cluster/de_genes.xlsx",
            caption="../reports/de_genes.rst",
            category="Differential Expression",
        ),
    log:
        "results/logs/cluster/cluster.log",
    benchmark:
        "results/benchmarks/cluster/cluster.txt"
    threads: 16
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/cluster.py"
