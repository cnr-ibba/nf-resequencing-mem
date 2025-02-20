# Resequencing Nextflow Pipeline

<!-- markdownlint-disable MD014 -->

[![GitHub Actions CI Status](https://github.com/cnr-ibba/nf-resequencing-mem/actions/workflows/ci.yml/badge.svg)](https://github.com/cnr-ibba/nf-resequencing-mem/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/cnr-ibba/nf-resequencing-mem/actions/workflows/linting.yml/badge.svg)](https://github.com/cnr-ibba/nf-resequencing-mem/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Overview

This pipeline performs resequencing analysis by calling
[freebayes](https://github.com/freebayes/freebayes) on reads aligned to the
genome using [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) _mem_.
Aligned _BAM_ files are then filtered with
[markduplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
to select reads aligned unambiguously to the genome. Genotypes are called from all samples, and the resulting _VCF_ file is normalized using [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md) and [bcftools](https://samtools.github.io/bcftools/bcftools.html#norm), which is the final output of this pipeline.

There's also an additional wiki for this project with more information about
this pipeline:
[nf-resequencing-mem wiki](https://github.com/cnr-ibba/nf-resequencing-mem/wiki)

## Setting up

In order to execute this pipeline, you will need `nextflow` installed and one of this
different executors: `conda`, `singularity` and `docker`.
The proper way to running this pipeline is described in
[pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html#pipeline-sharing)
nextflow manual, and lets nextflow to download and execute the pipeline, for example:

```bash
nextflow pull cnr-ibba/nf-resequencing-mem
```

You can specify the version of the pipeline to be downloaded by providing the
`-r` option, for example:

```bash
nextflow pull cnr-ibba/nf-resequencing-mem -r v0.6.2
```

## Testing stuff

In order to test this pipeline in your infrastructure, you can call this pipeline
on the _test_ dataset using the `test` profile:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -profile test,<your profile>
```

where `<your profile>` is one of `docker`, `singularity`, `apptainer` or `conda` profiles.
Nextflow will download the required containers (if necessary) and will execute
the pipeline on the test dataset.
For more information, please see [Profiles scope](#profiles-scope) section
of this documentation.

## Customize configuration

When running Nextflow, Nextflow looks for a file named `nextflow.config` in the
current directory and in the script base directory (if it is not the same as the
current directory). Finally it checks for the file `$HOME/.nextflow/config`, which
have lower priority than the other configuration files.
There's also a _institutional configuration_ for this pipeline, which is located
at [cnr-ibba/nf-config](https://github.com/cnr-ibba/nf-configs/blob/ibba/docs/pipeline/nf-resequencing-mem/ibba.md).
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
  columns are described (see [assets/samplesheet.csv](https://raw.githubusercontent.com/cnr-ibba/nf-resequencing-mem/refs/heads/master/assets/samplesheet.csv)
  for an example). In the
  `fastq_1` and `fastq_2` columns you need to specify the path fore _R1_ and _R2_
  files respectively. If you have single paired reads, leave `fastq_2` column empty.
  if you have more file for the same sample, specify all the single/pair files using
  the same _sample name_: this pipeline will append all reads belonging to the
  same sample before calling _trimgalore_ in the `CAT_FASTQ` process.
- `--genome_fasta`: (required) path to genome (FASTA) file. If file is compressed,
  index calculation will be forced even if provided by CLI
- `--outdir`: (required) directory in which store _results_

There are also additional pipeline parameters that can be provided and can be
used to save _intermediate results_ or to skip a particular step:

- `--genome_fasta_fai`: path to fasta index file (skip fasta index step)
- `--genome_bwa_index`: path to genome bwa index directory (skip bwa index step)
- `--save_cram`: (bool, def. false) save _markduplicated_ cram files with their indexes
  in results folder
- `--save_trimmed`: (bool, def. false) save trimmed reads in results folder
- `--save_fasta_index`: (bool, def. false) save fasta index (for reusing with this pipeline)
- `--save_bwa_index`: (bool, def. false) save bwa index (for reuse with this pipeline)
- `--save_freebayes`: (bool, def. false) save freebayes output file (not normalized!)
- `--remove_fastq_duplicates`: (bool, def. false) remove FASTQ duplicates by IDs using
  [seqkit/rmdup](https://bioinf.shenwei.me/seqkit/usage/#rmdup)
- `--save_unique_fastq`: (bool, def. false) write de-duplicated FASTQ files (requires
  `--remove_fastq_duplicates` option)
- `--ploidy`: (integer, def. 2) Sets the default ploidy for the analysis
- `--gvcf`: (bool, def. false) save output in _gVCF_ format (def. false)
- `--gvcf_chunk`: (integer) When writing gVCF output emit a record for every N bases.
- `--gvcf_dont_use_chunk`: (bool, def. false) When writing the gVCF output emit a
  record for all bases, will also route an int to `--gvcf_chunk` similar to
  `--output-mode EMIT_ALL_SITES` from _GATK_
- `--skip_normalization`: (bool, def. false) skip VCF normalization steps
- `--normalization_only`: (bool, def. false) only normalize a VCF file (skip all the
  other steps, see [Normalize a vcf file](#normalize-a-vcf-file))
- `--input_vcf`: path to a VCF file to be normalized (required when `--normalization_only` is set)
- `--input_tbi`: path to a VCF index file (required when `--normalization_only` is set)
- `--snpeff_database`: annotate the VCF file with SnpEff by providing a pre-built
  database that can be found using the `java -jar snpEff.jar databases` command.
  If the database is known to SnpEff will be downloaded and managed by the pipeline
  itself. Otherwise you can provide a custom database using `--snpeff_cachedir`
  and `--snpeff_config` parameters. See our wiki [Creating a snpEff database](https://github.com/cnr-ibba/nf-resequencing-mem/wiki/Running-the-pipeline#creating-a-snpeff-database)
  for more information.
- `--snpeff_cachedir`: SnpEff cache directory. It must contain a subdirectory with
  the same name of `--snpeff_database`, with a valid SnpEff database as a content.
  Is required when annotating with SnpEff with a custom database
- `--snpeff_config`: SnpEff custom config file. Is required **only** with a custom
  database. Needs to have the same custom database defined by `--snpeff_database` option
  (see: [Building databases](https://pcingola.github.io/SnpEff/snpeff/build_db/)
  of SnpEff documentation)

You can have a list of available parameters by calling:

```bash
nextflow run cnr-ibba/nf-resequencing-mem --help
```

In addition, instead of passing parameters using CLI, you can create a custom
configuration file and define each params in the _params scope_. According
nextflow _best practices_, parameters defined in _params scope_ should be
placed in a _json_ file, for example `params.json`:

```json
{
  "input": "<samplesheet.csv>",
  "genome_fasta": "<genome_fasta>",
  "outdir": "<results dir>",
  "save_fasta_index": true
}
```

Parameters have the same name used within the CLI, but without the `--` prefix.
Nextflow can be then called like this:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  -params-file params.json
```

There's also a way to provide parameters using a _config_ file:
according to the _nextflow_ community, this
configuration file should be used to specify all the parameters
that can't be specified using the command line, for example custom arguments
to be passed to a certain module using `ext.args` option can be defined in
a `custom.config` file like this:

```conf
process {
    withName: PICARD_MARKDUPLICATES {
        ext.args = '--TMP_DIR $TMPDIR'
    }
}
```

Then you can call `nextflow` providing such configuration file:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  -config custom.config
```

You can name this file as you prefer, but if is named `nextflow.config` it will
be automatically loaded by nextflow: is better to use a different name to avoid
confusion and explictly provide the configuration file using the `-config` option.
See Nextflow [Configuration](https://www.nextflow.io/docs/latest/config.html)
documentation for more information.

### Profiles scope

The profile scope is used to define configuration attributes that can be activated
in different environment or context when providing the profile _name_ using the
`-profile` option. For example, pipeline dependencies are managed using nextflow
by defining three different profiles:

- **conda**: every pipeline step will manage its requirements using conda in a
  specific environment. Conda environments are created inside `work` directory
  (but you can change this behavior using `cacheDir` option within the conda
  scope). This in not the recommended profile in a HPC environment, should be
  used if you aren't able to use singularity or docker profiles.
- **docker**: manage requirements using docker images. You will need to be part of
  the `docker` group in order to use this profile.
- **singularity**: manage requirements using singularity images. You can execute
  this profile without any permissions. `singularity` software need to be installed
  and available in your `$PATH` bash environment variable
- **apptainer**: manage requirements using apptainer images, which is the community
  driven continuation of Singularity, maintained under the Linux Foundation.

Then there are also profiles which define the _executor_ used to submit/collect a
job, that can be customized in order to affect the pipeline in a particular environment.
Profiles can be combined when calling nextflow, for example:

```bash
nextflow run -profile singularity,slurm ...
```

will enable all the configuration to be applied on slurm executor using singularity.
Profiles are also defined in the _institutional configuration_ file, with the aim
to describe and configure the infrastructure available by a certain institution.
For this pipeline, the base repository (`--custom_config_base`)
for the institutional configuration is located at
[cnr-ibba/nf-configs](https://github.com/cnr-ibba/nf-configs) using
[ibba](https://github.com/cnr-ibba/nf-configs/tree/ibba) as default branch. You
call any profile described by the
[cnr-ibba/nf-config/README.md](https://github.com/cnr-ibba/nf-configs/blob/ibba/README.md)
file.

### Process scope

The process configuration scope allows you to provide a specific configuration
for a specific process. For example, within this scope you can declare the
maximum number of CPUs used by a certain process or the amount of memory required,
for example:

```config
process {
  // https://www.nextflow.io/docs/edge/reference/process.html#resourcelimits
  resourceLimits = [
      cpus: 32,
      memory: 64.GB
  ]
}
```

You can specify here any property described in the process
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

> NOTE: freebayes _ploidy_ can be set using the `--ploidy` parameter, see
> [Params scope](#params-scope) section of this documentation.

## Nextflow specific parameters

There are parameters which are nextflow specific. They start all with a single
`-` before the parameter name and need to be provided using the nextflow CLI.
Here are the most used parameters:

- `-resume`: recover previous attempt
- `-profile`: specify one of `docker`, `singularity`, `apptainer` and `conda` profiles.
  `singularity` and `apptainer` are the recommended profile in a HPC environment
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
documentation page for more information.

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
You require also a _S3_ bucket in which nextflow could store and retrieve data (nextflow
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

## Normalize a vcf file

With this pipeline is it possible to perform the normalization workflow on a VCF
file, without running the whole pipeline. This is useful when you have a VCF file
that needs to be normalized, for example after a _freebayes_ run. You can call
this pipeline providing the `--normalization_only` parameter and the `--input_vcf`
and `--input_tbi` parameters:

```bash
nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> \
  --normalization_only --input_vcf <input.vcf> --input_tbi <input.tbi> \
  --genome_fasta <genome.fasta> --outdir <results dir>
```

Other provided parameters will be ignored, and the pipeline will normalize the
VCF file and will store the normalized VCF file in the `outdir` directory.
You can find more information about normalization workflow in our
[Normalizing VCF files](https://github.com/cnr-ibba/nf-resequencing-mem/wiki/Running-the-pipeline#normalizing-vcf-files) and
[Running normalization only](https://github.com/cnr-ibba/nf-resequencing-mem/wiki/Running-the-pipeline#running-normalization-only)
wiki pages

## Merge different datasets of SNPs

We have a wiki page that describes how to merge different datasets of SNPs
using this pipeline: you can find it at
[Merging different datasets of SNPs](https://github.com/cnr-ibba/nf-resequencing-mem/wiki/Merging-datasets)

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
header section, otherwise this pipeline will not work as intended.
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

The freebayes step can take a lot of time when calculating SNPs with a lot of data.
This process is calculated using multiple processes by splitting the whole genome in
regions (relying on CRAM alignment sizes), and then by calling SNPs on each region
on a single process.
In the last step, all results are collected and sorted to produce the final VCF file
(see `subworkflows/local/cram_freebayes_parallel` subworkflow for more information).
You can customize the region splitting, for example by using a greater cumulative
coverage (def. is `500_000_000`) or the minimum fragment length (def. is `10_000`)
in the split process like this:

```config
process {
    withName: FREEBAYES_SPLITCRAM {
        ext.args = '--max_coverage 100_000_000 --min_length 20_000'
    }
}
```

however, having a lot of smaller regions could be problematic if the overlapping
region shared by segments become larger than the unique region.

There are some freebayes processes that can fail by reaching time or memory
limits. According to freebayes creators (see issues
[#465](https://github.com/freebayes/freebayes/issues/465)
and [#626](https://github.com/freebayes/freebayes/issues/626)) this is unpredictable
and may depends on _coverage_, _allelic complexity_ and _sample_ or _ploidy_
number. Nextflow can resubmit such process increasing the required resources at
each step until `maxRetries` attempts are reached: you could increase the retry
attempts like this:

```config
process {
    withName: FREEBAYES_CHUNK {
        maxRetries = 10
    }
}
```

If this is not sufficient and jobs continue to fail for memory and time requirements,
you could think to skip regions with a very high coverage (`-g` option) or
try to sub-sample regions with high coverage with the `--limit-coverage` option,
like this:

```config
process {
    withName: FREEBAYES_CHUNK {
        ext.args = '--limit-coverage 25'
    }
}
```

This pipeline compute a _coverage_ step before calling freebayes, so you could
try to determine which value make sense to be used for filtering.

> **N.B.** Even with downsampling this step could require a lot of time, and time
> required for each step will be unpredictable. Using a very low coverage limit
> could affect the SNP calling process. Use this option with caution.

### FASTQ has duplicated IDs

When a `FASTQ` file have duplicated IDs, the `MarkDuplicates` fails with an error
message like this:

```text
Exception in thread "main" htsjdk.samtools.SAMException: Value was put into PairInfoMap more than once.  11: RGE200002173L1C027R02800999933
  at htsjdk.samtools.CoordinateSortedPairInfoMap.ensureSequenceLoaded(CoordinateSortedPairInfoMap.java:133)
  at htsjdk.samtools.CoordinateSortedPairInfoMap.remove(CoordinateSortedPairInfoMap.java:86)
  at picard.sam.markduplicates.util.DiskBasedReadEndsForMarkDuplicatesMap.remove(DiskBasedReadEndsForMarkDuplicatesMap.java:61)
  at picard.sam.markduplicates.MarkDuplicates.buildSortedReadEndLists(MarkDuplicates.java:571)
  at picard.sam.markduplicates.MarkDuplicates.doWork(MarkDuplicates.java:270)
  at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:308)
  at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:103)
  at picard.cmdline.PicardCommandLine.main(PicardCommandLine.java:113)
```

This is a not an issue with `MarkDuplicates` (as discussed
[here](https://gatk.broadinstitute.org/hc/en-us/community/posts/4408717387803-SAMException-Value-was-put-into-PairInfoMap-more-than-once))
but an issue at _demultiplexing_ step: the only way to deal with this problem is
to make rid of duplicated IDs using [seqkit/rmdup](https://bioinf.shenwei.me/seqkit/usage/#rmdup)
by providing the `--remove_fastq_duplicates` option.

### MarkDuplicates temporary files

Markduplicates writes temporary files into `/tmp` partition by default. If your
organization have a different location where temporary files should be stored
(ex `/scratch` or any other `$TMP` position) and your jobs are running out of
spaces, you should provide a different temporary location to `MarkDuplicates` steps,
for example:

```config
process {
    withName: PICARD_MARKDUPLICATES {
        ext.args = '--TMP_DIR $TMPDIR'
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
