# Snakemake workflow: `10X single-cell + LARRY`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


A Snakemake workflow to process single-cell libraries generated with 10XGenomics platform (RNA, ATAC and RNA+ATAC) together with [LARRY barcoding](https://www.nature.com/articles/s41586-020-2503-6).


## Setup

The following files are located inside the folder `config`. In this folder you will find the files with raw data paths (`units.tsv`), resources configuration (`resources.yaml`) and pipeline parameters (`config.yaml`).

### Raw data

Paths to raw data are located in the file `units.tsv`. The file has the following structure:

| sample_id | lane | lib_type | R1 | R2 | R3 |
|-----------|------|----------|----|----|----|
| name_of_sample | name_of_lane_or_resequencing | library type | path/to/forward.fastq.gz | path/to/reverse.fastq.gz | path/to/ATAC-R3.fastq.gz

* The first field correspond to the sample name. This field has to be identical for all the fastqs corresponding to the same sample, idependently of the library type of the fastq. If a sample is split in 2 different library types (such as ATAC + RNA or RNA and LARRY), both of them must have the same sample_id.

* The second field corresponds to `lane`. The idea of this field is to group fastq files corresponding to the same sample (or to samples that have to be merged). For example, if 1 sample arrived in 2 different lanes from a PE experiment, in total there will be 4 fastqs (2 forward and 2 reverse). In this case, one should enter the same sample 2 times, putting in the `lane` field the corresponding lanes (lane1 and lane2, for example). Actually one can write any word in this field, the idea is to group fastqs from the same sample. All the entries with the same name in the `sample` field with different `lane` will be merged in the same fastq. Here an example of how it would be with 1 sample that arrived in 2 lanes:

| sample_id | lane | lib_type | R1 | R2 | R3 |
|-----------|------|----------|----|----|----|
| foo | lane1 | GEX | path/to/forward_lane1.fastq.gz | path/to/reverse_lane1.fastq.gz | |
| foo | lane2 | GEX | path/to/forward_lane2.fastq.gz | path/to/reverse_lane2.fastq.gz | |

Here I am using lane1 and lane2 for consistency and making things more clear, but the following would also work:

| sample_id | lane | lib_type | R1 | R2 | R3 |
|-----------|------|----------|----|----|----|
| foo | potato | ATAC | path/to/forward_lane1.fastq.gz | path/to/reverse_lane1.fastq.gz | path/to/R3_lane1.fastq.gz
| foo | whatever | ATAC | path/to/forward_lane2.fastq.gz | path/to/reverse_lane2.fastq.gz | path/to/R3_lane2.fastq.gz

* The third field `lib_type` corresponds to the type of library. Basically it can be one of 3: GEX (Gene expression), ATAC (Chromatin accessibility) or FB (Feature barcoding, which has to be set for the LARRY fastqs).

* Finally the last 3 fields `R1`, `R2` and `R3`correspond to the paths to the fastq files. For RNA (GEX) `R1` is the FORWARD read and `R2` the REVERSE. Usually the cellular barcode is present in the R1 and the transcript in the R2, same for LARRY. In the case of the ATAC, `R2` corresponds to the dual illumina indexing and `R3` corresopnds to the reverse read (where the tn5 fragment is present).

### Configuration of pipeline parameters

In this configuration folder (I say this because there's in another folder a file with the same name) there is the file `config.yaml`. This file contains the configuration of the software and parameters used in the pipeline. Modify them as you wish. Check always that you are using the correct genome files corresponding to the version that you want to use. 

Inside the file, and also in the file `workflow/schema/config.schema.yaml` you can find what is controled by each tunable parameter.

### Resources configuration

`config/resources.yaml` contains the per rule resource parameters (ncpus, ram, walltime...). It can be modified as desired.

### Snakemake profiles

In Snakemake 4.1 [snakemake profiles](https://github.com/Snakemake-Profiles) were introduced. They are supposed to substitute the classic cluster.json file and make the execution of snakemake more simple. The parameters that will be passed to snakemake (i.e: --cluster, --use-singularity...) now are inside a yaml file (`config.yaml`) inside the profile folder (in the case of this repository is `snakemake_profile`). The `config.yaml` inside `snakemake_profile` contains the parameters passed to snakemake. So if you were executing snakemake as `snakemake --cluster qsub --use-singularity` the new `config.yaml` would be like this:

```yaml
cluster: qsub
use-singularity: true
```

## Execution of the pipeline

Once you have all the configuration files as desired, it's time to execute the pipeline. For that you have to execute the `execute_pipeline.sh` script, followed by the name of the rule that you want to execute. If no rule is given it will automatically execute the rule `all` (which would execute the standard pipeline). Examples:

```bash
./execute_pipeline.sh all
```

is equivalent to 

```bash
./execute_pipeline.sh
```
