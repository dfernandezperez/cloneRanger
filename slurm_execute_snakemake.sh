# Usage --> ./execute_pipeline.sh  <rules> (can be nothing or a rule defined in the snakefile without wildcards in the input files)
snakemake -j 1 --unlock
snakemake --profile config/snakemake_profile/slurm "$@" > snakemake.log 2> snakemake.err < /dev/null