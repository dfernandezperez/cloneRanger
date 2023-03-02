# Downloads and extracts cell ranger
# The naming and rules ensure that any verion can be used
# rule get_cellranger:
#     output:
#         cr=directory("resources/cellranger"),
#         bin="resources/cellranger/bin/cellranger",
#     params:
#         url=config["get_cellranger"]["url"],
#     log:
#         "results/logs/get_cellranger/get_cellranger.log",
#     benchmark:
#         "results/benchmarks/get_cellranger/get_cellranger.txt"
#     shell:
#         """
#             wget -O resources/cellranger.tar.gz "{params.url}" &> {log}
#         tar -xzf resources/cellranger.tar.gz -C resources &> {log} && \
#         rm -rf resources/cellranger.tar.gz
#         mv resources/cellranger-* resources/cellranger
#         """


# rule get_reference:
#     output:
#         dir=directory("resources/genome"),
#     params:
#         url=config["get_reference"]["url"],
#     log:
#         "results/logs/get_reference/get_reference.log",
#     benchmark:
#         "results/benchmarks/get_reference/get_reference.txt"
#     shell:
#         """
#             wget -O resources/genome.tar.gz "{params.url}" &> {log}
#         tar -xzf resources/genome.tar.gz -C resources &> {log} && \
#         rm -rf resources/genome.tar.gz
#         mv resources/refdata-* resources/genome
#         """

# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = "data/{sample}_S1_L00{lane}_R1_001.fastq.gz",
        rv = "data/{sample}_S1_L00{lane}_R2_001.fastq.gz",
    log:
        "results/logs/clean_names/{lane}_{sample}.log",
    benchmark:
        "results/benchmarks/clean_names/{lane}_{sample}.txt"
    shell:
        """
        mv {input[0]} {output.fw} &> {log}
        mv {input[1]} {output.rv} &>> {log}
        """
