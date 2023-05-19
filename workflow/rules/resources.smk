# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = temp("data/symlink/{sample}_{lib_type}_S1_L00{lane}_R1_001.fastq.gz"),
        rv = temp("data/symlink/{sample}_{lib_type}_S1_L00{lane}_R2_001.fastq.gz"),
    container:
        None
    shell:
        """
        ln -s {input[0]} {output.fw}
        ln -s {input[1]} {output.rv}
        """

rule clean_names_atac:
    input:
        get_atac_fastqs,
    output:
        temp("data/symlink/{sample}_{lib_type}_S1_L00{lane}_R3_001.fastq.gz"),
    container:
        None
    shell:
        """
        ln -s {input} {output}
        """

rule merge_lanes:
    input:
        fw = lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_S1_L00{sample.lane}_R1_001.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            ),
        rv = lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_S1_L00{sample.lane}_R2_001.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            )
    output:
        fw = "data/lane_merged/{sample}_{lib_type}_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_{lib_type}_S1_L001_R2_001.fastq.gz"
    log:
        "results/00_logs/merge_lanes/{sample}_{lib_type}.log"
    container:
        None
    shell:
        """
        cat {input.fw} > {output.fw} 2> {log}
        cat {input.rv} > {output.rv} 2>> {log}
        """

rule merge_lanes_atac:
    input:
        lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_S1_L00{sample.lane}_R3_001.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            )
    output:
        "data/lane_merged/{sample}_{lib_type}_S1_L001_R3_001.fastq.gz"
    log:
        "results/00_logs/merge_lanes_atac/{sample}_{lib_type}.log"
    container:
        None
    shell:
        """
        cat {input} > {output} 2> {log}
        """
    
rule move_gex_fq:
    # Dummy rule to move the gex fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_GEX_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_GEX_S1_L001_R2_001.fastq.gz"
    output:
        fw = "data/clean/{sample}_GEX_S1_L001_R1_001.fastq.gz",
        rv = "data/clean/{sample}_GEX_S1_L001_R2_001.fastq.gz"
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        """

rule move_cellhash_fq:
    # Dummy rule to move the gex fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_CH_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_CH_S1_L001_R2_001.fastq.gz"
    output:
        fw = "data/clean/{sample}_CH_S1_L001_R1_001.fastq.gz",
        rv = "data/clean/{sample}_CH_S1_L001_R2_001.fastq.gz"
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        """

rule move_atac_fq:
    # Dummy rule to move the atac fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_ATAC_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_ATAC_S1_L001_R2_001.fastq.gz",
        r3 = "data/lane_merged/{sample}_ATAC_S1_L001_R3_001.fastq.gz"
    output:
        fw = "data/clean/{sample}_ATAC_S1_L001_R1_001.fastq.gz",
        rv = "data/clean/{sample}_ATAC_S1_L001_R2_001.fastq.gz",
        r3 = "data/clean/{sample}_ATAC_S1_L001_R3_001.fastq.gz"
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        mv {input.r3} {output.r3}
        """

rule create_library:
    input:
        expand("data/clean/{{sample}}_{lib_type}_S1_L001_{read}_001.fastq.gz", lib_type = LIB_TYPES, read = ["R1", "R2"]),
    output:
        "data/feature_bc_libraries/{sample}_library.csv"
    params:
        library_types = get_library_type
    shell:
        """
        echo "fastqs,sample,library_type" > {output}
        echo "{params.library_types}" >> {output}
        """

rule create_cmo_set:
    output:
        "data/feature_bc_libraries/cmo-set.csv"
    params:
        cell_hashing = CELL_HASHING["barcodes"],
    script:
        "../scripts/python/create_cmo_set.py"


rule create_library_multi:
    input:
        unpack(get_library_input)
    output:
        "data/feature_bc_libraries/{sample}_library_multi.csv"
    params:
        library_types     = get_library_type,
        cell_hashing      = CELL_HASHING["assignments"],
        cellranger_params = config["cellranger_count"]["extra_parameters_rna"],
        introns           = convert_introns(),
        n_cells           = get_expected_cells(),
        genome            = config["genome_reference_gex"],
        mem_gb            = config["cellranger_count"]["mem"],
    log:
        "results/00_logs/create_library_multi/{sample}.log"
    script:
        "../scripts/python/create_multi_library.py"