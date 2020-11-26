
Resequencing Nextflow Pipeline
==============================

When running Nextflow, Nextflow looks for a file named `nextflow.config` in the
current directory and in the script base directory (if it is not the same as the
current directory). Finally it checks for the file `$HOME/.nextflow/config`.
Please modify `params` in `nextflow.config` according your needs:

* The `params` scope defines variables which applies to the whole pipeline.
* The `profiles` scope defines variables which apply to the particular profile
invoked with the `-profile` Nextflow parameter
* The `process` scope can define parameters applied to a single process (for example
the number of CPUs used or the required RAM)

Pipeline profiles
-----------------

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

Calling nextflow using conda
----------------------------

### Creating an environment relying on environment export file

`cacheDir`: defined in [conda scope](https://www.nextflow.io/docs/latest/config.html#scope-conda),
specifies where environment should be created or found

```
profiles {
  conda {
    process.conda = "$projectDir/environment.yml"
    conda.cacheDir = "/home/paolo/Softwares/anaconda3/envs"
  }
}
```

Calling nextflow using docker
-----------------------------

```
$ nextflow run main.nf --reads_path='data/*R{1,2}_001.fastq.gz' -resume -profile docker
```

Calling nextflow using singularity
----------------------------------

```
$ nextflow run main.nf -resume -profile singularity --reads_path="/mnt/storage/140822_M01314_0102_000000000-AA7HH/*_R{1,2}_*.fastq.gz"
```
