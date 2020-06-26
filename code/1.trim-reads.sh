###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: demultiplex amplicon data based on loci primers
#PROJECT: amplicon-genotyping
###.............................................................................
#create symlinks
#create symlinks to raw data
for fasta in /Users/miguelcamacho/Dropbox/RattusbaluensisTrusMadi/MiSeq1908/*
do
  npath=$(basename $fasta | xargs echo data/raw | sed 's| |/|')
  ln -s $fasta $npath
done

#batch renaming after patterns in data/raw/id-match
#eg: CLKDJ_1_CATCGT~GGAGCG_[12].fastq.gz	BOR0614.#1.fastq.gz
mmv < data/raw/id-match

#remove samples with no match
rm data/raw/CLKDJ_1*
#get list with samples for which R1 or R2 are missing
ls data/raw/BOR* | sed -e 's|.[12].fastq.gz||;s|data/raw/||' | sort | uniq -c \
 | grep '^.*1 ' | awk '{print $2}' > data/intermediate/unpaired-samples
#get list with samples for which R1 or R2 are missing
ls data/raw/BOR* | sed -e 's|.[12].fastq.gz||;s|data/raw/||' | sort | uniq -c \
| grep '^.*2 ' | awk '{print $2}' > data/intermediate/samples-list
#remove unpaired samples
cat data/intermediate/unpaired-samples | while read sample
do
  rm data/raw/*$sample.[12]*
done

#dual demultiplexing of reads.
#dual demultiplexing is only enabled from cutadapt v2.1
cat data/intermediate/samples-list | while read sample
do
  cutadapt \
      -e 0.15 --no-indels \
      --discard-untrimmed \
      --pair-adapters \
      -q 10 \
      -m 150 \
      -g file:data/raw/forward-primers.fa \
      -G file:data/raw/reverse-primers.fa \
      -o data/intermediate/$sample.{name1}-{name2}.1.fastq.gz \
      -p data/intermediate/$sample.{name1}-{name2}.2.fastq.gz \
      data/raw/$sample.1.fastq.gz \
      data/raw/$sample.2.fastq.gz
done
