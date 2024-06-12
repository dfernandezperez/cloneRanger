# Snakemake profiles

In Snakemake 4.1 [snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) were introduced. They are supposed to substitute the classic cluster.json file and make the execution of snakemake more simple. The parameters that will be passed to snakemake (i.e: --cluster, --use-singularity...) now are inside a yaml file (`config.yaml`) inside the profile folder (in the case of this repository is `config/snakemake_profile`). There is one config.yaml for local execution (`config/snakemake_profile/local`) and another for execution in HPCs using the slurm scheduler (`config/snakemake_profile/slurm`). The `config.yaml` contains the parameters passed to snakemake. 

If you are running the pipeline locally (personal computer/workstation), some important parameters that you may have to adapt in `config/snakemake_profile/local` are:

```yaml
cores            : 120                # Define total amount of cores that the pipeline can use
resources        : mem_mb=128000      # Define total amount of ram that pipeline can use
singularity-args : "--bind /stemcell" # Volumes to mount for singularity (important to mount the volume where the pipeline will be executed)
```

Adapt `cores`, `resources` and `default-resources` to your computer hardware. Then adapt the path to bind in singularity from `singularity-args` to the folder in which you are running the pipeline form.

In case of running the pipeline using slurm:

```yaml
jobs             : unlimited # Define total amount of parallel jobs that the pipeline can execute
singularity-args : "--bind /scratch --bind /data --bind /home" # Volumes to mount for singularity (important to mount the volume where the pipeline will be executed)
```

Particularly important to adapt the volumes to bind for singularity to your own HPC.