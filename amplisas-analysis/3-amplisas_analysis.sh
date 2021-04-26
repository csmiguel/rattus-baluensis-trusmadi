cd amplisas-analysis/
#merge F and R reads
cd
docker run -v $PWD:/workdir --rm sixthresearcher/amplisat ampliMERGE.pl \
    -i 1.fastq.gz 2.fastq.gz \
    -o merged_reads > amplisas.log

# ampliCLEAN.pl seems to collapse due to the large file size. For that reason I split the fastq in to 10 parts:
seqkit split -p 10 merged_reads.fq >> amplisas.log
rm merged_reads.fq
#clean merged reads. It can also clean demultiplexed reads.
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# if -d amplicon-data.csv  is included the processing time increases x30 because it needs to demultiplex.
for partx in ls merged_reads.fq.split/*fq
do
  outname=$(echo $partx | sed 's|^.*\(part.*\).fq|\1|')
  docker run -v $PWD:/workdir --rm sixthresearcher/amplisat ampliCLEAN.pl \
      -i $partx \
      -o filtered_reads_$outname >> amplisas.log
done
rm -rf merged_reads.fq.split
# concatenate split sequences
cat filtered_reads_part*fq > filtered_reads.fq
rm *part_0*fq
#fast analysis of amplicon sequences
# Description:
#   Fast pairwise comparison among higher depth sequences in amplicon sequencing (AS) data
#   Analyzes AS data and gives as output an Excel file where are shown the similarities among the high depth sequences
#   This results can be used to optimize the parameters for a further analysis with AmpliSAS
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl a
docker run -v $PWD:/workdir --rm sixthresearcher/amplisat ampliCHECK.pl \
    -d amplicon-data.csv \
    -i filtered_reads.fq \
    -o results_amplicheck >> amplisas.log
# loci for which there where no reads:
echo "\n\nloci for which there where no reads:" >> amplisas.log
cat amplisas.log | grep 'de-multiplexed (0 reads' | cut -d"-" -f1 | uniq -c >> amplisas.log

# run amplisSAS
docker run -v $PWD:/workdir --rm sixthresearcher/amplisat ampliSAS.pl \
-i results_amplicheck/results.xlsx \
-o results_amplisas

#remove intermediate files
rm [12].fastq.gz
rm filtered_reads.fq
rm -rf
