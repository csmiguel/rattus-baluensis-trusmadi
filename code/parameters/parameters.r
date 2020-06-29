#input data
sepe = NULL #"se" or "pe"; paired-end or single-end reads
idata = "excel" #"excel" or "txt"; is input data (primers, loci, samples) in excel or txt.

#genotype calls
ab = 0.3 #threshold for allele balance (eg. 0.3 implies that reads_min_allele/reads_max_allele > 0.3)
cov = 11 #minimum coverage for an allele in an individual not be be dropped
ploidy <- 2 #ploidy of the organism

#pop gen filtering
remove_monomorphic = "yes" #remove remove monomorphic loci
max_missing_ind = 1 #max. proportion of missing calls per individual (recommended 0.3)
max_missing_loci = 1 #max. proportion of missing calls per locus (recommended 0.3)

#output
structure_str = "yes" #if "yes", genotypes a written to a STRUCUTRE formatted file is in output/strcuture
phylogen = "yes" #if "yes", a random allele from each individual is written to a fasta file in output/phylogen.fasta
ade_genind = "yes" #if "yes", genotypes are stored in adegenet genlight format in output/genlight.Rdata
aleles_fasta = "yes" #if "yes", all alleles are written to a multifasta file with the format ">sampleID-locusID-allele[a-z]"
