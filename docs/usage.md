# nfcore/rnaseq Usage

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
module load blic-modules
module load nextflow

nextflow_rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh38
```

This will launch the pipeline with the `legion` or `myriad` configuration profile, depending on where you submit the job from.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Main Arguments

### `-profile`
This parameter is NOT necessary as the shortcut `nextflow_rnaseq` takes care of selecting the appropiate configuration profile. But just for your information, profiles are used to give 
configuration presets for different compute environments.

* `legion`
    * A generic configuration profile to be used with the UCL cluster legion
* `myriad`
    * A generic configuration profile to be used with the UCL cluster myriad

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### Library strandedness
Three command line flags / config parameters set the library strandedness for a run:

* `--forward_stranded`
* `--reverse_stranded`
* `--unstranded`

If not set, the pipeline will be run as unstranded. Specifying `--pico` makes the pipeline run in `forward_stranded` mode.

These flags affect the commands used for several steps in the pipeline - namely HISAT2, featureCounts, RSeQC (`RPKM_saturation.py`) and StringTie:

* `--forward_stranded`
  * HISAT2: `--rna-strandness F` / `--rna-strandness FR`
  * featureCounts: `-s 1`
  * RSeQC: `-d ++,--` / `-d 1++,1--,2+-,2-+`
  * StringTie: `--fr`
* `--reverse_stranded`
  * HISAT2: `--rna-strandness R` / `--rna-strandness RF`
  * featureCounts: `-s 2`
  * RSeQC: `-d +-,-+` / `-d 1+-,1-+,2++,2--`
  * StringTie: `--rf`

## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the alignment tool instead. Developed by the same group behind the popular Tophat aligner, HISAT2 has a much smaller memory footprint.

To use HISAT2, use the parameter `--aligner hisat2` or set `params.aligner = 'hisat2'` in your config file.

Additionally, the pipeline runs [kallisto](https://pachterlab.github.io/kallisto/) for transcript quantification. As it does not need alignment, this step is failry quick, but if you wnat to skip it you cna do it using the 
`--skip_kallisto` argument.

#### Minimum alignment lengths
When using STAR as the aligner then the pipeline sets the minimum alignment length to 15 base pairs. This filters out very short read alignments that can result from soft-clipping.

To customise this threshold, use the pipeline command line parameter `--min_aln_length [num]`, where `[num]` is the minimum required alignment length in base pairs. Mapping sensitivity can be improved by increasing this value, or made more flexible by decreasing it. The default value in STAR is `0`.

## Reference Genomes

The pipeline config files come bundled with paths to the following reference genome assemblies: GRCh37, GRCh38, GRCm38. The pipeline has this aprameter set up as `false` by default, you need to specify it using the `--genome` flag.

* Human
  * `--genome GRCh37`
  * `--genome GRCh38`
* Mouse
  * `--genome GRCm38`

### `--star_index`, `--hisat2_index`, `--fasta`, `--gtf`, `--bed12`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--star_index '[path to STAR index]' \
--hisat2_index '[path to HISAT2 index]' \
--fasta '[path to Fasta reference]' \
--gtf '[path to GTF file]' \
--bed12 '[path to bed12 file]'
```

Note that only one of `--star_index` / `--hisat2_index` are needed depending on which aligner you are using (see below).

The minimum requirements are a Fasta and GTF file. If these are provided and no others, then all other reference files will be automatically generated by the pipeline.

### `--downloadFasta`, `--downloadGTF`
Instead of a path to a file, a URL can be supplied to download reference Fasta and GTF files at the start of the pipeline. A required STAR index and BED12 files will then be generated from these downloaded files.

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveAlignedIntermediates`
As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM files from STAR / HISAT2 and sorting steps.

## Adapter Trimming
If specific additional trimming is required (for example, from additional tags),
you can use any of the following command line parameters. These affect the command
used to launch TrimGalore!

### `--clip_r1 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

### `--clip_r2 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

### `--three_prime_clip_r1 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

### `--three_prime_clip_r2 [int]`
Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.


## Library Prep Presets
Some command line options are available to automatically set parameters for common RNA-seq library preparation kits.

> Note that these presets override other command line arguments. So if you specify `--pico --clip_r1 0`, the `--clip_r1` bit will be ignored.

If you have a kit that you'd like a preset added for, please let us know!

### `--pico`
Sets trimming and standedness settings for the _SMARTer Stranded Total RNA-Seq Kit - Pico Input_ kit.

Equivalent to: `--forward_stranded` `--clip_r1 3` `--three_prime_clip_r2 3`


## Skipping QC steps
The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

* `--skip_kallisto` -          Skip kallisto
* `--skip_qc` -                Skip **all QC steps**, apart from MultiQC
* `--skip_fastqc` -            Skip FastQC
* `--skip_rseqc` -             Skip RSeQC
* `--skip_genebody_coverage` - Skip calculating the genebody coverage
* `--skip_preseq` -            Skip Preseq
* `--skip_dupradar` -          Skip dupRadar (and Picard MarkDups)
* `--skip_edger` -             Skip edgeR MDS plot and heatmap
* `--skip_multiqc` -           Skip MultiQC

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

###
## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.

### `--clusterOptions`
Submit arbitrary cluster scheduler options (not available for all config profiles). For instance, you could use `--clusterOptions '-p devcore'` to run on the development node (though won't work with default process time requests).

## Stand-alone scripts
The `bin` directory contains some scripts used by the pipeline which may also be run manually:

* `gtf2bed`
  * Script used to generate the BED12 reference files used by RSeQC. Takes a `.gtf` file as input
* `dupRadar.r`
  * dupRadar script used in the _dupRadar_ pipeline process.
* `edgeR_heatmap_MDS.r`
  * edgeR script used in the _Sample Correlation_ process
* `RNA-pipeline-from-BAM.sh`
  * SLURM script used to mimic pipeline QC steps, taking an aligned BAM file as input.
  * Potentially unmaintained, use at your own risk!
