if config["cellranger_count"]["10x_pipeline"] == "GEX":

    rule cellranger_count:
        input:
            unpack(get_cellranger_input)
        output:
            mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
            outs = directory("results/01_counts/{sample}/outs/filtered_feature_bc_matrix"),
            html = report(
                "results/01_counts/{sample}/outs/web_summary.html",
                caption  = "../reports/counts.rst",
                category = "Cellranger Counts",
                subcategory = "{sample}",
            ),
        params:
            introns     = convert_introns(),
            bam         = create_bam(),
            n_cells     = get_expected_cells,
            genome      = config["cellranger_count"]["genome_reference_gex"],
            extra_p     = config["cellranger_count"]["extra_parameters_rna"],
            mem_gb      = round(int(RESOURCES["cellranger_count"]["mem_mb"])/1000),
            feature_ref = get_feature_ref
        log:
            "results/00_logs/counts/{sample}.log",
        benchmark:
            "results/benchmarks/counts/{sample}.txt"
        threads:
            RESOURCES["cellranger_count"]["cpu"]
        resources:
            mem_mb  = RESOURCES["cellranger_count"]["mem_mb"],
            runtime = RESOURCES["cellranger_count"]["runtime"]
        container: 
            config["singularity"]["cellranger_rna_sif"]
        shell:
            """
            cellranger count \
            {params.feature_ref} \
            {params.n_cells} \
            {params.extra_p}  \
            {params.introns} \
            {params.bam} \
            --id {wildcards.sample} \
            --transcriptome {params.genome} \
            --libraries {input.libraries} \
            --localcores {threads} \
            --localmem {params.mem_gb} \
            &> {log} && \
            # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
            # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
            rm -r results/01_counts/{wildcards.sample} && \
            mv {wildcards.sample} results/01_counts/{wildcards.sample}
            """


elif config["cellranger_count"]["10x_pipeline"] == "ATAC":
   rule cellranger_count:
        input:
            lambda w: expand("data/clean/{{sample}}_ATAC_S1_L001_{read}_001.fastq.gz", read=["R1", "R2"])
        output:
            mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix.h5",
            atac = "results/01_arc/{sample}/outs/atac_fragments.tsv.gz",
            html = report(
                "results/01_counts/{sample}/outs/web_summary.html",
                caption     = "../reports/counts.rst",
                category    = "Cellranger Counts",
                subcategory = "{sample}",
            ),
        params:
            genome  = config["cellranger_count"]["genome_refrerence_arc"],
            mem_gb  = round(int(RESOURCES["cellranger_count"]["mem_mb"])/1000),
            extra_p = config["cellranger_count"]["extra_parameters_atac"]
        log:
            "results/00_logs/counts/{sample}.log",
        benchmark:
            "results/benchmarks/counts/{sample}.txt"
        threads:
            RESOURCES["cellranger_count"]["cpu"]
        resources:
            mem_mb  = RESOURCES["cellranger_count"]["mem_mb"],
            runtime = RESOURCES["cellranger_count"]["runtime"]
        container: 
            config["singularity"]["cellranger_atac_sif"]
        shell:
            """
            cellranger-atac count \
            {params.extra_p} \
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

elif config["cellranger_count"]["10x_pipeline"] == "ARC":
    rule cellranger_arc_count:
        input:
            fq_gex  = expand(
                            "data/clean/{sample}_GEX_S1_L001_{read}_001.fastq.gz", 
                            sample = SAMPLES,
                            read   = ["R1", "R2"]
                        ),
            fq_atac = expand(
                            "data/clean/{sample}_ATAC_S1_L001_{read}_001.fastq.gz", 
                            sample = SAMPLES,
                            read   = ["R1", "R2", "R3"]
                        ),
            libraries   = "data/feature_bc_libraries/{sample}_library_arc.csv"
        output:
            mtx  = "results/01_arc/{sample}/outs/filtered_feature_bc_matrix.h5",
            atac = "results/01_arc/{sample}/outs/atac_fragments.tsv.gz",
            outs = directory("results/01_arc/{sample}/outs/filtered_feature_bc_matrix"),
            html = report(
                "results/01_arc/{sample}/outs/web_summary.html",
                caption     = "../reports/counts.rst",
                category    = "Cellranger Counts",
                subcategory = "{sample}",
            ),
        params:
            introns = convert_introns(),
            genome  = config["cellranger_count"]["genome_refrerence_arc"],
            mem_gb  = round(int(RESOURCES["cellranger_arc_count"]["mem_mb"])/1000),
            extra_p = config["cellranger_count"]["extra_parameters_arc"]
        log:
            "results/00_logs/cellranger_arc_count/{sample}.log",
        benchmark:
            "results/benchmarks/cellranger_arc_count/{sample}.txt"
        threads:
            RESOURCES["cellranger_arc_count"]["cpu"]
        resources:
            mem_mb  = RESOURCES["cellranger_arc_count"]["mem_mb"],
            runtime = RESOURCES["cellranger_arc_count"]["runtime"]
        container: 
            config["singularity"]["cellranger_multiome_sif"]
        shell:
            """
            cellranger-arc count \
            {params.introns} \
            {params.extra_p} \
            --id {wildcards.sample} \
            --reference {params.genome} \
            --libraries {input.libraries} \
            --localcores {threads} \
            --localmem {params.mem_gb} \
            &> {log} && \
            # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
            # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
            rm -r results/01_arc/{wildcards.sample} && \
            mv {wildcards.sample} results/01_arc/{wildcards.sample}
            """
            
    if is_feature_bc():
        rule cellranger_count:
            input:
                unpack(get_cellranger_input)
            output:
                mtx  = "results/01_counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
                outs = directory("results/01_counts/{sample}/outs/filtered_feature_bc_matrix"),
                html = report(
                    "results/01_counts/{sample}/outs/web_summary.html",
                    caption  = "../reports/counts.rst",
                    category = "Cellranger Counts",
                    subcategory = "{sample}",
                ),
            params:
                introns     = convert_introns(),
                bam         = create_bam(),
                n_cells     = get_expected_cells,
                genome      = config["cellranger_count"]["genome_reference_gex"],
                extra_p     = config["cellranger_count"]["extra_parameters_rna"],
                mem_gb      = round(int(RESOURCES["cellranger_count"]["mem_mb"])/1000),
                feature_ref = get_feature_ref
            log:
                "results/00_logs/counts/{sample}.log",
            benchmark:
                "results/benchmarks/counts/{sample}.txt"
            threads:
                RESOURCES["cellranger_count"]["cpu"]
            resources:
                mem_mb  = RESOURCES["cellranger_count"]["mem_mb"],
                runtime = RESOURCES["cellranger_count"]["runtime"]
            container: 
                config["singularity"]["cellranger_rna_sif"]
            shell:
                """
                cellranger count \
                {params.feature_ref} \
                {params.n_cells} \
                {params.extra_p}  \
                {params.introns} \
                {params.bam} \
                --id {wildcards.sample} \
                --chemistry=ARC-v1 \
                --transcriptome {params.genome} \
                --libraries {input.libraries} \
                --localcores {threads} \
                --localmem {params.mem_gb} \
                &> {log} && \
                # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
                # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
                rm -r results/01_counts/{wildcards.sample} && \
                mv {wildcards.sample} results/01_counts/{wildcards.sample}
                """