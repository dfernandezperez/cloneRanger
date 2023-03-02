rule densities:
    input:
        data=rules.cluster.output.data,
    params:
        features=config.get("densities", {}).get("features", ["sample"]),
    output:
        data="results/densities/adata.h5ad",
        densities=report(
            "results/densities/densities.pdf",
            caption="../reports/densities.rst",
            category="Embedding Densities",
        ),
    log:
        "results/logs/densities/densities.log",
    benchmark:
        "results/benchmarks/densities/densities.txt"
    threads: 1
    conda:
        "../envs/densities.yaml"
    script:
        "../scripts/densities.py"
