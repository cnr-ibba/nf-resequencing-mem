# Resequencing Nextflow Pipeline

<!-- markdownlint-disable MD014 -->

## Setting up

In order to execute this pipeline, you will need `nextflow` installed and one of this
different executors: `conda`, `singularity` and `docker`. You can choose to clone
this repository if you plan to change this pipeline according your needs:

```text
$ git clone https://github.com/cnr-ibba/nf-resequencing-mem
```

The other way to running this pipeline is described in
[pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html#pipeline-sharing)
nextflow manual. You will need also to define your credentials for private
repositories. See [SCM configuration file](https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file)
for more details. After that, call this pipeline with:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  --input <samplesheet.csv> --genome_fasta <genome_fasta> --outdir <results dir>
```

where the following are `nextflow` specific parameters:

- `-resume`: recover previous attempt
- `-profile`: specify one of `docker`, `singularity` and `conda` profiles. `singularity`
  is the recommended profile in a HPC environment

these are pipeline parameters which are _mandatory_:

- `--input`: (required) specify the samplesheet CSV/TSV file where `sample,fastq_1,fastq_2`
  columns are described (see `assets/samplesheet.csv` for an example). In the
  `fastq_1` and `fastq_2` columns you need to specify the path fore _R1_ and _R2_
  files respectively. If you have single paired reads, leave `fastq_2` column empty.
  if you have more file for the same sample, specify all the single/pair files using
  the same _sample name_: this pipeline will append all reads belonging to the
  same sample before calling _trimgalore_
- `--genome_fasta`: (required) path to genome (FASTA) file. If file is compressed,
  index calculation will be forced even if provided by CLI

There are also additional pipeline parameters that can be provided:

- `--genome_fasta_fai`: path to fasta index file (skip fasta index step)
- `--genome_bwa_index`: path to genome bwa index directory (skip bwa index step)
- `--save_bam`: (bool, def. false) save markduplicated bam files with their indexes
  in results folder
- `--save_trimmed`: (bool, def. false) save trimmed reads in results folder
- `--save_fasta_index`: (bool, def. false) save fasta index (for reuse with this pipeline)
- `--save_bwa_index`: (bool, def. false) save bwa index (for reuse with this pipeline)

### Provide parameters as a config file

In alternative (for reproducibility purpose) you can create a custom configuration
file and provide it when calling nextflow. For example if you create a file like this

```conf
params {
  input = "<samplesheet.csv>"
  genome_fasta = "<genome_fasta>"
  outdir = "<results dir>"
}
```

Then you can call `nextflow` providing such configuration file:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  -config custom.config
```

See Nextflow [Configuration](https://www.nextflow.io/docs/latest/config.html)
documentation for more information.

## Customize configuration

When running Nextflow, Nextflow looks for a file named `nextflow.config` in the
current directory and in the script base directory (if it is not the same as the
current directory). Finally it checks for the file `$HOME/.nextflow/config`.
Please modify `params` in `nextflow.config` according your needs:

- The `params` scope defines variables which applies to the whole pipeline.
- The `profiles` scope defines variables which apply to the particular profile
  invoked with the `-profile` Nextflow parameter
- The `process` scope can define parameters applied to a single process (for example
  the number of CPUs used or the required RAM)

There are also configuration files in the `conf` folder of this repository, for
example the `conf/modules.config` file keeps the configuration for each module
used within this pipeline: you can affect a module behavior without modifying the
proper module simply changing the proper process configuration section.

### Params scope

The params scope let you to define parameters to the pipeline. You can set those
params directly in config file or passing those from command line.

### Profiles scope

By default, this pipeline is supposed to work in an environment with all required softwares
already installed. You could install and manage all required software in a conda
environment, for example, and then run the pipeline inside this environment:
**this is not the recommended way to call this pipeline**.

To better manage softwares
dependencies is better to let **nextflow manage dependencies** for you, by calling
the pipeline with a [custom profile](https://www.nextflow.io/docs/edge/config.html#config-profiles).
Nextflow profiles let nextflow to manage software dependencies and custom parameters.
Three profiles are currently defined:

- **conda**: every pipeline step will manage its requirements using conda in a
  specific environment. Conda environments are created inside `work` directory
  (but you can change this behaviour using `cacheDir` option within the conda
  scope).
- **docker**: manage requirements using docker images. You will need to be part of
  the `docker` group in order to use this profile
- **singularity**: manage requirements using singularity images. You can execute
  this profile without any permissions. `singularity` software need to be installed
  and available in your `$PATH` bash environment variable

### Process scope

The process configuration scope allows you to provide the default configuration
for all the processes. You can specify here any property described in the process
directive and the executor sections and override them.
With the `withLabel` selector, you can configure of all processes annotated with
such label selector. For example:

```text
process {
  // only processes with this label have those parameters
  withLabel: process_high {
    cpus = 4
    memory = 4.GB
  }
}
```

Will affect all the processes annotated with `process_high` label. Similarly, you
can provide additional parameters to a certain step. For example, supposing to
change the `freebayes` _ploidy_ since you are dealing with a non-diploid sample
(default): you can provide a `custom.config` file in which pass extra parameters
to freebayes:

```text
process {
    withName: FREEBAYES_MULTI {
        ext.args = '--ploidy 4'
    }
}
```

Then call `nextflow` with `-config` parameter

## Calling this pipeline with local executor

The local executor is the default executor used when calling pipeline. It spawns
pipeline processes using fork and threads using all your local CPUs and memory
available. If you need to set a limit to the resources used, enable this configuration
for local executor:

```text
executor {
  // for the local executer, I will set the maximum values of CPU and MEMORY
  $local {
    cpus = 8
    memory = '16 GB'
  }
}
```

## Calling this pipeline using pbs executor

You can change the default executor by specifying the `pbs` profile. Simply add
such profile to your command line, for example:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile pbs,singularity \
  --input "<samplesheet.csv>" --genome_fasta <genome_fasta> --outdir <results dir>
```

In alternative, you can export this environment variable:

```bash
export NXF_EXECUTOR=pbs
```

## Calling this pipeline using slurm executor

Like `pbs` executor, simply add such profile to your command line, for example:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile slurm,singularity \
  --input "<samplesheet.csv>" --genome_fasta <genome_fasta> --outdir <results dir>
```

In alternative, you can export this environment variable:

```bash
export NXF_EXECUTOR=slurm
```

## Calling this pipeline using AWS batch

This pipeline could run using the [AWS batch queue system](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html).
In order to do so, you need to configure your credentials with [aws cli](https://docs.aws.amazon.com/translate/latest/dg/setup-awscli.html):
you require to configure a _IAM_ account with permission to run _batch_, _EC2_ and _S3_.
You require also a s3 bucket in which nextflow could store and retrive data (nextflow
will make a copy of the input data and will retrieve the results from here) and
a AWS batch queue with _EC2 spot instances_ as recommended compute environment.
After that, you could launch this pipeline by providing only _awsbatch_ as profile
and the _queue name_ and the AWS _region_ with the `--awsqueue` and `--awsregion`
parameters:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile awsbatch \
  -bucket-dir s3://<s3 bucket name>/<subfolder> \
  --input '<samplesheet.csv>' --genome_fasta <genome_fasta> \
  --awsqueue <aws batch queue name> --awsregion <aws region>
```

Please see the [Amazon Cloud](https://www.nextflow.io/docs/latest/awscloud.html#)
section of nextflow documentation to get other information on nextflow and AWS
usage.

## Known issues

### Ignore sample sheet check

Rarely this pipeline could fail at `INPUT_CHECK:SAMPLESHEET_CHECK`
claiming that `given sample sheet does not appear to contain a header`,
even if your sample sheet starts with `sample,fastq_1,fastq_2` header
section. This behavior is very uncommon and described in `nf-core/tools`
issue [#1539](https://github.com/nf-core/tools/issues/1539). There's a way
to overcome this issue by providing `--has_header` parameter to
`bin/check_samplesheet.py` which let to skip the header check sections. However,
you have to ensure your samplesheet have the required `sample,fastq_1,fastq_2`
header section. You can do it by customizing the `SAMPLESHEET_CHECK` step
in the process scope, for example

```text
process {
    withName: SAMPLESHEET_CHECK {
        ext.args = '--has_header'
    }
}
```
