# Multilocus genotyping from amplicon sequencing
The goal of this repository is to provide an easy and powerful workflow to genotype loci from amplicon sequencing in High-Throughput Sequencing (HTS), such as from Illumina MiSeq.
* Input: amplicon reads from HTS.
* Output: individual x loci genotype matrix, plus all alleles in FASTA.

It provides a modular and customizable workflow. The core of the genotyping module relies on dada2, an R package designed for determining Amplicon Sequence Variants (ASVs) from metabarcoding.

## What it does:
* genotype amplicon reads.
* it takes as input overlapping and non-overlapping reads.
* it can genotype simultaneously samples from the same or different species.
* output format in STRUCTURE, genotype table and FASTA.

## What it does **not** do:
* current version only works with **paired-end** reads.
* current version has only been tested in **diploids**.
* some processes from current version do not run in parallele.

# Set up

Clone the current repository.
In `\test` you can find a worked example.
You should have installed in added to the $PATH the following programs (also avaible in the `code` folder):
cutadapt version 2.1 or above
mmv
The scripts are dependent on the following R packages with their dependencies (also available on the `/code` folder): `dada2, adegenet, dplyr, tidyr, seqinr`
Run the script `/code/checks.R` to confirm you have all the installed dependencies.

Move your R1 and R2 reads to `data/raw`.
Edit the Excel file `data/raw/data.xls`, or alternatively, edit the following plain text files:
`data/raw/forward`, with forward primers in FASTA format
`data/raw/reverse`, with reverse primers in FASTA format
`data/raw/match`, with loci
These data has to be sorted. That is, the forward and reverse primer from line1 in the text files should amplify locus in line1 from the loci file, and so on.

# Running

Scripts should be run in order from 1 (`code/01fasf.R`)to 5 (`code/05-fasfda.R`) in the Terminal/Console. In the first time, wait until the script has finished running and check the expected results in `data/intermediate` or `output` before running the next script.
```
sh 01.script.sh
Rscript -e 01.script.R
Rscript -e 02.script.R
Rscript -e 03.script.R
Rscript -e 04.script.R
```
At any time,  you can delete intermadiate results `rm -rf data/intermediate` and output `rm -rf output`, to re-do analysis from raw data.


## Output

Output files can be found in `output`:
* 




