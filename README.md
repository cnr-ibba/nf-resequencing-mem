
Resequencing Nextflow Pipeline
==============================

Setting up
----------

In order to execute this pipeline, you will need `nextflow` installed and one of this
different executors: `conda`, `singularity` and `docker`. You can choose to clone
this repository if you plan to change this pipeline according your needs:

```
$ git clone https://github.com/cnr-ibba/nf-resequencing-mem
```

The other way to running this pipeline is described in
[pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html#pipeline-sharing)
nextflow manual. You will need also to define your credentials for private
repositories. See [SCM configuration file](https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file)
for more details. After that, call this pipeline with:

```
$ nextflow run cnr-ibba/nf-resequencing-mem -resume -profile <your profile> --reads_path "<reads_path/*_R{1,2}_*.fastq.gz>" --genome_path <genome_path> --outdir <results dir>
```

Customize configuration
-----------------------

When running Nextflow, Nextflow looks for a file named `nextflow.config` in the
current directory and in the script base directory (if it is not the same as the
current directory). Finally it checks for the file `$HOME/.nextflow/config`.
Please modify `params` in `nextflow.config` according your needs:

* The `params` scope defines variables which applies to the whole pipeline.
* The `profiles` scope defines variables which apply to the particular profile
invoked with the `-profile` Nextflow parameter
* The `process` scope can define parameters applied to a single process (for example
the number of CPUs used or the required RAM)

### params scope

The params scope let you to define parameters to the pipeline. You can set those
params directly in config file or passing those from command line.

### profiles scope

By default, this pipeline is supposed to work in an environment with all required softwares
already installed. You could install and manage all required software in a conda
environment, for example, and then run the pipeline inside this environment:
**this is not the recommended way to call this pipeline**.

To better manage softwares
dependencies is better to let **nextflow manage dependencies** for you, by calling
the pipeline with a [custom profile](https://www.nextflow.io/docs/edge/config.html#config-profiles).
Nextflow profiles let nextflow to manage software dependencies and custom parameters.
Three profiles are currently defined:

* **conda**: every pipeline step will manage its requirements using conda in a
specific environment. Conda environments are created inside `work` directory
(but you can change this behaviour using `cacheDir` option within the conda
scope).
* **docker**: manage requirements using docker images. You will need to be part of
the `docker` group in order to use this profile
* **singularity**: manage requirements using singularity images. You can execute
this profile without any permissions. `singularity` software need to be installed
and available in your `$PATH` bash environment variable

### process scope

The process configuration scope allows you to provide the default configuration
for all the processes. You can specify here any property described in the process
directive and the executor sections and override them.
With the `withLabel` selector, you can configure of all processes annotated with
such label selector. For example:

```
process {
  // only processes with this label have those parameters
  withLabel: process_high {
    cpus = 4
    memory = 4.GB
  }
}
```

Will affect all the processes annotated with `process_high` label.

Calling this pipeline with local executor
-----------------------------------------

The local executor is the default executor used when calling pipeline. It spawns
pipeline processes using fork and threads using all your local CPUs and memory
available. If you need to set a limit to the resources used, enable this configuration
for local executor:

```
executor {
  // for the local executer, I will set the maximum values of CPU and MEMORY
  $local {
    cpus = 8
    memory = '16 GB'
  }
}
```

Calling this pipeline using pbs executor
----------------------------------------

You can change the default executor by specifying the `pbs` profile. Simply add
such profile to your command line, for example:

```
$ nextflow run cnr-ibba/nf-resequencing-mem -resume -profile pbs,singularity --reads_path "<reads_path/*_R{1,2}_*.fastq.gz>" --genome_path <genome_path> --outdir <results dir>
```
