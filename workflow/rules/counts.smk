if config["10x_pipeline"] == "GEX":

    if is_feature_bc():
    # No conda env used as cellranger cannot be installed in one
        rule cellranger_count:
            input:
                fq = lambda w: expand(
                                "data/clean/{col.sample_id}_{col.lib_type}_S1_L001_{read}_001.fastq.gz", 
                                col=units.loc[w.sample].itertuples(), 
                                read=["R1", "R2"]
                            ),
                libraries   = "data/feature_bc_libraries/{sample}_library.csv",
                feature_ref = "data/feature_reference/Feature_reference.csv"
            output:
                mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
                html = report(
                    "results/01_counts/{sample}/outs/web_summary.html",
                    caption="../reports/counts.rst",
                    category="Cellranger Counts",
                    subcategory="{sample}",
                ),
            params:
                introns     = convert_introns(),
                n_cells     = config["cellranger_count"]["n_cells"],
                genome      = config["genome_reference_gex"],
                mem_gb      = config["cellranger_count"]["mem"]
            log:
                "results/00_logs/counts/{sample}.log",
            benchmark:
                "results/benchmarks/counts/{sample}.txt"
            threads:
                RESOURCES["cellranger_count"]["cpu"]
            resources:
                mem_mb = RESOURCES["cellranger_count"]["MaxMem"]
            container: 
                config["cellranger_rna_sif"]
            shell:
                """
                cellranger count \
                --nosecondary \
                {params.introns} \
                --id {wildcards.sample} \
                --transcriptome {params.genome} \
                --libraries {input.libraries} \
                --feature-ref {input.feature_ref} \
                --expect-cells {params.n_cells} \
                --localcores {threads} \
                --localmem {params.mem_gb} \
                &> {log} && \
                # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
                # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
                rm -r results/01_counts/{wildcards.sample} && \
                mv {wildcards.sample} results/01_counts/{wildcards.sample}
            """

    else:
        rule cellranger_count:
            input:
                lambda w: expand("data/clean/{{sample}}_GEX_S1_L001_{read}_001.fastq.gz", read=["R1", "R2"])
            output:
                mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
                html = report(
                    "results/01_counts/{sample}/outs/web_summary.html",
                    caption="../reports/counts.rst",
                    category="Cellranger Counts",
                    subcategory="{sample}",
                ),
            params:
                introns     = convert_introns(),
                n_cells     = config["cellranger_count"]["n_cells"],
                genome      = config["genome_reference_gex"],
                mem_gb      = config["cellranger_count"]["mem"]
            log:
                "results/00_logs/counts/{sample}.log",
            benchmark:
                "results/benchmarks/counts/{sample}.txt"
            threads:
                RESOURCES["cellranger_count"]["cpu"]
            resources:
                mem_mb = RESOURCES["cellranger_count"]["MaxMem"]
            container: 
                config["cellranger_rna_sif"]
            shell:
                """
                cellranger count \
                --nosecondary \
                {params.introns} \
                --id {wildcards.sample} \
                --transcriptome {params.genome} \
                --fastqs data/clean \
                --sample {wildcards.sample}_GEX \
                --expect-cells {params.n_cells} \
                --localcores {threads} \
                --localmem {params.mem_gb} \
                &> {log} && \
                # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
                # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
                rm -r results/01_counts/{wildcards.sample} && \
                mv {wildcards.sample} results/01_counts/{wildcards.sample}
                """

elif config["10x_pipeline"] == "ATAC":
   rule cellranger_count:
        input:
            lambda w: expand("data/clean/{{sample}}_ATAC_S1_L001_{read}_001.fastq.gz", read=["R1", "R2"])
        output:
            mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix.h5",
            html = report(
                "results/01_counts/{sample}/outs/web_summary.html",
                caption="../reports/counts.rst",
                category="Cellranger Counts",
                subcategory="{sample}",
            ),
        params:
            genome      = config["genome_refrerence_arc"],
            mem_gb      = config["cellranger_count"]["mem"]
        log:
            "results/00_logs/counts/{sample}.log",
        benchmark:
            "results/benchmarks/counts/{sample}.txt"
        threads:
            RESOURCES["cellranger_count"]["cpu"]
        resources:
            mem_mb = RESOURCES["cellranger_count"]["MaxMem"]
        container: 
            config["cellranger_atac_sif"]
        shell:
            """
            cellranger-atac count \
            --id {wildcards.sample} \
            --reference {params.genome} \
            --fastqs data/clean \
            --localcores {threads} \
            --localmem {params.mem_gb} \
            &> {log} && \
            # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
            # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
            rm -r results/01_counts/{wildcards.sample} && \
            mv {wildcards.sample} results/01_counts/{wildcards.sample}
        """

elif config["10x_pipeline"] == "ARC":
    rule cellranger_count:
        input:
            fq = lambda w: expand(
                            "data/clean/{lane.sample_id}_{lane.lib_type}_S1_L001_{read}_001.fastq.gz", 
                            lane=units.loc[w.sample].itertuples(), 
                            read=["R1", "R2"]
                        ),
            libraries   = "data/feature_bc_libraries/{sample}_library.csv"
        output:
            mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix.h5",
            html = report(
                "results/01_counts/{sample}/outs/web_summary.html",
                caption="../reports/counts.rst",
                category="Cellranger Counts",
                subcategory="{sample}",
            ),
        params:
            introns     = convert_introns(),
            genome      = config["genome_refrerence_arc"],
            mem_gb      = config["cellranger_count"]["mem"]
        log:
            "results/00_logs/counts/{sample}.log",
        benchmark:
            "results/benchmarks/counts/{sample}.txt"
        threads:
            RESOURCES["cellranger_count"]["cpu"]
        resources:
            mem_mb = RESOURCES["cellranger_count"]["MaxMem"]
        container: 
            config["cellranger_multiome_sif"]
        shell:
            """
            cellranger-arc count \
            {params.introns} \
            --id {wildcards.sample} \
            --reference {params.genome} \
            --libraries {input.libraries} \
            --localcores {threads} \
            --localmem {params.mem_gb} \
            &> {log} && \
            # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
            # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
            rm -r results/01_counts/{wildcards.sample} && \
            mv {wildcards.sample} results/01_counts/{wildcards.sample}
            """