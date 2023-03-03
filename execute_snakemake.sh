# Usage --> ./execute_pipeline.sh  <rules> (can be nothing or a rule defined in the snakefile without wildcards in the input files)
# Save PID of snakemake process to kill it in case it is required
snakemake -j 1 --unlock
nohup snakemake --profile config/snakemake_profile "$@" > results/snakemake.log 2> results/snakemake.err < /dev/null &
echo $! > pid.txt