# Resequencing Nextflow Pipeline

<!-- markdownlint-disable MD014 -->

## Overview

This pipeline will execute a resequencing analysis by calling
[freebayes](https://github.com/freebayes/freebayes)
on reads aligned to genome using [bwa](https://bio-bwa.sourceforge.net/bwa.shtml)
_mem_. Genotypes will be called from all samples, and the resulting _VCF_
file will be normalized using [bcftools](https://samtools.github.io/bcftools/bcftools.html#norm),
which is the final output of this pipeline.

## Setting up

In order to execute this pipeline, you will need `nextflow` installed and one of this
different executors: `conda`, `singularity` and `docker`. You can choose to clone
this repository if you plan to change this pipeline according your needs:

```bash
git clone https://github.com/cnr-ibba/nf-resequencing-mem
```

The other way to running this pipeline is described in
[pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html#pipeline-sharing)
nextflow manual, and lets nextflow to download and execute the pipeline, for example:

```bash
nextflow pull cnr-ibba/nf-resequencing-mem
```

You will need also to define your credentials for private
repositories. See [SCM configuration file](https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file)
for more details.

## Customize configuration

When running Nextflow, Nextflow looks for a file named `nextflow.config` in the
current directory and in the script base directory (if it is not the same as the
current directory). Finally it checks for the file `$HOME/.nextflow/config`.
Additionally, you can provide a custom config file using the `-config` option,
which has an higher priority respect to the other configuration files.
There are tree different _scopes_, which can be used to customize this
pipeline according your needs:

- The `params` scope defines variables which applies to the whole pipeline:
  those _params_ could be the _input samplesheet_ or the genome to be used
- The `profiles` scope defines variables which apply to the particular profile
  invoked with the `-profile` Nextflow parameter
- The `process` scope can define parameters applied to a single process (for example
  the number of CPUs used or the required RAM)

All this configuration scopes can be customized in each configuration file,
and configuring a value in an higher priority file will override values
defined in lower priority configuration files.
There are also configuration files in the `conf` folder of this repository, for
example the `conf/modules.config` file keeps the configuration for each module
used within this pipeline: you can affect a module behavior without modifying the
proper module simply changing the proper process configuration section.

### Params scope

The params scope let you to define pipeline parameters. You can set those
params directly in config file or pass those from command line (which have
the highest priority and override each parameter defined in other configuration
files). These are pipeline parameters which are _mandatory_:

- `--input`: (required) specify the samplesheet CSV/TSV file where `sample,fastq_1,fastq_2`
  columns are described (see `assets/samplesheet.csv` for an example). In the
  `fastq_1` and `fastq_2` columns you need to specify the path fore _R1_ and _R2_
  files respectively. If you have single paired reads, leave `fastq_2` column empty.
  if you have more file for the same sample, specify all the single/pair files using
  the same _sample name_: this pipeline will append all reads belonging to the
  same sample before calling _trimgalore_ in the `CAT_FASTQ` process.
- `--genome_fasta`: (required) path to genome (FASTA) file. If file is compressed,
  index calculation will be forced even if provided by CLI

There are also additional pipeline parameters that can be provided and can be
used to save _intermediate results_ or to skip a particular step:

- `--genome_fasta_fai`: path to fasta index file (skip fasta index step)
- `--genome_bwa_index`: path to genome bwa index directory (skip bwa index step)
- `--save_bam`: (bool, def. false) save _markduplicated_ bam files with their indexes
  in results folder
- `--save_trimmed`: (bool, def. false) save trimmed reads in results folder
- `--save_fasta_index`: (bool, def. false) save fasta index (for reusing with this pipeline)
- `--save_bwa_index`: (bool, def. false) save bwa index (for reuse with this pipeline)
- `--save_freebayes`: (bool, def. false) save freebayes output file (not normalized!)

You can have a list of available parameters by calling:

```bash
nextflow run cnr-ibba/nf-resequencing-mem --help
```

In addition, instead of passing parameters using CLI, you can create a custom configuration
file and define each params in the _params scope_. Parameters have the same name
used within the CLI, but without the `--` prefix. For example if you create a
`custom.config` file like this

```conf
params {
  input = "<samplesheet.csv>"
  genome_fasta = "<genome_fasta>"
  outdir = "<results dir>"
  save_fasta_index = true
}
```

Then you can call `nextflow` providing such configuration file:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  -config custom.config
```

See Nextflow [Configuration](https://www.nextflow.io/docs/latest/config.html)
documentation for more information.

### Profiles scope

The profile scope is used to define configuration attributes that can be activated
in different environment or context when providing the profile _name_ using the
`-profile` option. For example, pipeline dependencies are managed using nextflow
by defining three different profiles:

- **conda**: every pipeline step will manage its requirements using conda in a
  specific environment. Conda environments are created inside `work` directory
  (but you can change this behaviour using `cacheDir` option within the conda
  scope).
- **docker**: manage requirements using docker images. You will need to be part of
  the `docker` group in order to use this profile
- **singularity**: manage requirements using singularity images. You can execute
  this profile without any permissions. `singularity` software need to be installed
  and available in your `$PATH` bash environment variable

Then there are also profiles which define the _executor_ used to submit/collect a
job, that can be customized in order to affect the pipeline in a particular environment.
Profiles can be combined when calling nextflow, for example:

```bash
nextflow run -profile singularity,slurm ...
```

will enable all the configuration to be applied on slurm executor using singularity

### Process scope

The process configuration scope allows you to provide a specific configuration
for a specific process. You can specify here any property described in the process
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

will affect all the processes annotated with `process_high` label. Similarly, you
can provide additional parameters to a certain step. For example, supposing to
change the `freebayes` _ploidy_ since you are dealing with a non-diploid sample
(default): you can provide a `custom.config` file in which pass extra parameters
to freebayes:

```text
process {
    withName: FREEBAYES_CHUNK {
        ext.args = '--ploidy 4'
    }
}
```

## Nextflow specific parameters

There are parameters which are nextflow specific. They start all with a single
`-` before the parameter name and need to be provided using the nextflow CLI.
Here are the most used parameters:

- `-resume`: recover previous attempt
- `-profile`: specify one of `docker`, `singularity` and `conda` profiles. `singularity`
  is the recommended profile in a HPC environment
- `-work-dir`: save temporary files in a default temporary directory
  (def. `$projectDir/work`)

## Calling this pipeline with local executor

The local executor is the default executor used when calling pipeline. It spawns
pipeline processes using fork and threads using all your local CPUs and memory
available. If you need to set a limit to the resources used, enable this configuration
for local executor:

```text
executor {
  // for the local executor, I will set the maximum values of CPU and MEMORY
  $local {
    cpus = 8
    memory = '16 GB'
  }
}
```

Using the `$` symbol before the executor name, let you to specify different
configuration in the same _executor_ scope. See
[scope executor](https://www.nextflow.io/docs/latest/config.html?highlight=configuration#scope-executor)
documentation page for more informations.

## Calling this pipeline using pbs executor

You can change the default executor by specifying the `pbs` profile. Simply add
such profile to your command line, for example:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile pbs,singularity \
  --input <samplesheet.csv> --genome_fasta <genome.fasta> --outdir <results dir>
```

In alternative, you can export this environment variable:

```bash
export NXF_EXECUTOR=pbs
```

## Calling this pipeline using slurm executor

Like `pbs` executor, simply add such profile to your command line, for example:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile slurm,singularity \
  --input <samplesheet.csv> --genome_fasta <genome.fasta> --outdir <results dir>
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
header section, otherwise you will loose a sample.
You can do it by customizing the `SAMPLESHEET_CHECK` step
in the process scope, for example:

```text
process {
    withName: SAMPLESHEET_CHECK {
        ext.args = '--has_header'
    }
}
```

### Tuning freebayes resources usage

The freebayes step can take a lot of time when calculating SNP with a lot of data.
This process is calculated on single process by splitting the whole genome in
regions, and calling SNP on each region before collect all results in a unique file.
You can customize the region splitting, for example by using a smaller file size
in the split process like this:

```text
process {
    withName: FREEBAYES_SPLITBAM {
        ext.args = '--target-data-size 10e6'
    }
}
```

Then there are some freebayes process that can fail by reaching time or memory
limits. Nextflow can resubmit such process increasing the required resources at
each step until `maxRetries` attempts are reached: you could increase the retry
attempts like this:

```text
process {
    withName: FREEBAYES_CHUNK {
        maxRetries = 5
    }
}
```

## Acknowledgments

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
