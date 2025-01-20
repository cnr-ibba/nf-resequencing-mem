# cnr-ibba/nf-resequencing-mem: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.6.2 - dev

- Update MultiQC module
- Move functions inside workflows
- Convert freebayes specific parameters in pipeline parameters ([#80](https://github.com/cnr-ibba/nf-resequencing-mem/issues/80))
- Update `nextflow.config`
- Support for institutional configuration
- Parallelize normalization steps by chromosomes. Merge VCF files after normalization
- Normalize VCF file with `vcfallelicprimitives` ([#76](https://github.com/cnr-ibba/nf-resequencing-mem/issues/76))
- Add `freebayes_normalized` local subworkflow
- Update `nextflow` to version `24.04.0`
- Using the `resourceLimits` directive to set the max requirements for each process
- Update CI system ([#81](https://github.com/cnr-ibba/nf-resequencing-mem/issues/81))
- Use remote files with `test` profile
- Fixed issue with `picard/markduplicates` ([#77](https://github.com/cnr-ibba/nf-resequencing-mem/issues/77))
- Display current chromosome in `samtools/depth`

### `Added`

- Add `pipeline_initialization` local subworkflow
- Add institutional configuration custom repository
- Add `bcftools/concat` process
- Add `bcftools_filltags` process
- Add `freebayes_norm` process
- Add `freebayes_normalized` local subworkflow
- Add `bcftools/sort` process

### `Fixed`

- Use remote files with `test` profile
- `picard/markduplicates` was updated to the latest release in order to work
  with `*.cram` files

### `Removed`

- Remove unsupported options in `nextflow.config`
- Remove `tabix/tabix` after the first normalization step (removing overlap)
- Remove `check_max` function

## 0.6.1 - [2024-04-11]

- Calculate `samtools/depth` on each chromosomes ([#73](https://github.com/cnr-ibba/nf-resequencing-mem/issues/73))

### `Fixed`

- Passing args to `modules/local/freebayes_splitcram`
- Calculate `samtools/depth` without _0 coverage_ regions

## 0.6.0 - [2024-04-04]

- Replace `*.bam` file format with `*.cram` ([#9](https://github.com/cnr-ibba/nf-resequencing-mem/issues/9))
- Add _Read Groups_ during the alignment step ([#57](https://github.com/cnr-ibba/nf-resequencing-mem/issues/57))
- Annotate VCF file with SnpEff ([#59](https://github.com/cnr-ibba/nf-resequencing-mem/issues/59))
- Configure MultiQC analysis ([#60](https://github.com/cnr-ibba/nf-resequencing-mem/issues/60))
- Update modules ([#64](https://github.com/cnr-ibba/nf-resequencing-mem/issues/64))

### `Added`

- Add `samtools/depth` process
- Add `freebayes_splitcram` custom module to split genome in regions relying on
  _total sample coverage_
- Add `cram_markduplicates_picard` custom local subworkflow by modifying
  `bam_markduplicates_picard` to work with `*.cram` files by default
- Add `cram_stats_samtools` custom local subworkflow by modifying
  `bam_stats_samtools` to work with `*.cram` files by default
- Add `freebayes_splitcram` local module to splice alignments regions relying
  on `samtools/depth` step
- Add `snpeff/download` module
- Add `snpeff/snpeff` module
- Add `snpeff_annotate` local subworkflow

### `Fixed`

- `freebayes_parallel` subworkflow was moved to `cram_freebayes_parallel` local
  subworkflow and was modified to deal with _total sample coverage_ and to work
  with `*.cram` files
- `picard/markduplicates` now works with `.*.cram` files
- `bwa/mem` was configured to write files as `*.cram` files
- `samtools/depth` was patched to write results with headers, with 0 coverage position
  and to compress output with gzip
- `resequencing-mem` workflow was modified in order to use local subworkflow, for
  example to deal with `samtools` and `markduplicates`
- fixed a issue when providing the `--genome_bwa_index` parameter
- `snpeff_download` was patched in order to remove the `version` parameter
- `snpeff/snpeff` module was patched to support custom database annotations and
  to compress VCF output using a mulled image with `tabix`
- the configuration file for MultiQC module was updated to simplify results, to order
  them and to support all the supported modules

### `Removed`

- Remove `cnr-ibba/bamaddrg` module
- Remove `cnr-ibba/freebayes/splitbam` module

## 0.5.2 - [2023-12-21]

- Use MultiQC with all supported tools ([#53](https://github.com/cnr-ibba/nf-resequencing-mem/issues/53))
- Update modules
- Minor fixes

### `Added`

- Add `samtools/stats` module
- Add `samtools/idxstats` module
- Add `bcftools/stats` module

## 0.5.1 - [2023-11-15]

- Updated pipeline to support the latest version of `nf-core/tools` ([#46](https://github.com/cnr-ibba/nf-resequencing-mem/issues/46))

### `Added`

- Check that reads IDs are unique (using `seqkit/rmdup` - [#47](https://github.com/cnr-ibba/nf-resequencing-mem/issues/47))

## v0.5.0 - [2022-11-21]

- Improved freebayes performances ([#41](https://github.com/cnr-ibba/nf-resequencing-mem/issues/41))
- Calling freebayes in different processes ([#44](https://github.com/cnr-ibba/nf-resequencing-mem/issues/44))

### `Fixed`

- `CAT_FASTQ` result is not saved by default
- split chromosome relying on BAM size while running freebayes
- call freebayes in distinct process for each chromosome chunk

## v0.4.3 - [2022-11-03]

- Upgrade modules ([#34](https://github.com/cnr-ibba/nf-resequencing-mem/issues/34))
- Solve linting issues ([#32](https://github.com/cnr-ibba/nf-resequencing-mem/issues/27))
- Force `check_samplesheet.py` assuming header present ([#31](https://github.com/cnr-ibba/nf-resequencing-mem/issues/31))
- Normalize VCF file ([#33](https://github.com/cnr-ibba/nf-resequencing-mem/issues/33))

### `Added`

- Add `bcftools/norm` module
- Add `tabix/tabix` module

### `Fixed`

- Force samplesheet having header in `INPUT_CHECK:SAMPLESHEET_CHECK` workflow
- Fix software dump version
- Support for `--help` option
- Fix resource limits

### `Removed`

- Freebayes result is not more published by default
- Removed unused params

## v0.4.2 - [2022-07-29]

- Improved coverage step with `samtools/coverage` ([#27](https://github.com/cnr-ibba/nf-resequencing-mem/issues/27))
- Saving indexed files as outputs

### `Added`

- Add `samtools/coverage` module

### `Removed`

- Remove `BEDTOOLS_GENOMECOV` local process

## v0.4.1 - [2022-07-28]

- Calling `freebayes` on all samples ([#11](https://github.com/cnr-ibba/nf-resequencing-mem/issues/11))
- Deal with compressed genomes ([#14](https://github.com/cnr-ibba/nf-resequencing-mem/issues/14))
- Calculate sample coverage ([#17](https://github.com/cnr-ibba/nf-resequencing-mem/issues/17))
- Provide _indexes_ as parameters and save them as results

### `Added`

- Add `freebayes/multi` process
- Add `BEDTOOLS_GENOMECOV` local process

### `Removed`

- Remove `freebayes/single` module

## v0.4.0 - [2022-07-27]

- Updating modules ([#20](https://github.com/cnr-ibba/nf-resequencing-mem/issues/20), [#4](https://github.com/cnr-ibba/nf-resequencing-mem/issues/4))
- Name samples within bam ([#18](https://github.com/cnr-ibba/nf-resequencing-mem/issues/18))
- Export _BAM_ and _trimmed fastq_ ([#16](https://github.com/cnr-ibba/nf-resequencing-mem/issues/16))
- Add _module_ and _base_ config files ([#10](https://github.com/cnr-ibba/nf-resequencing-mem/issues/10))
- Get input from _samplesheet_ ([#15](https://github.com/cnr-ibba/nf-resequencing-mem/issues/15))

### `Added`

- Add `INPUT_CHECK:SAMPLESHEET_CHECK` workflow
- Add `bamaddrg` module

### `Removed`

- Remove `samtools/sort` module (now sorting is done in `bwa/mem`)

### `Fixed`

- Track sample names in `.bam` files

## v0.3.0 - [2021-11-08]

- Support _AWS_ environment ([#2](https://github.com/cnr-ibba/nf-resequencing-mem/issues/2))
- Structure pipeline with the new nextflow template ([#7](https://github.com/cnr-ibba/nf-resequencing-mem/issues/7))

## v1.0dev - [2020-11-24]

Initial release of `cnr-ibba/nf-resequencing-mem` as nexflow implementation of
_resequencing-mem_ pipeline
