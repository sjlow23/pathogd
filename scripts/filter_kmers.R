#!/usr/bin/env Rscript

#library(Biostrings)
library(data.table)
library(dplyr)
#library(stringr)
library(tidyr)

#$1 : threshold for filtering
#$2 : combined counts matrix
#$3	: output directory
#$4 : run addPAM step (yes/no)

args <- commandArgs(trailingOnly=TRUE)
threshold <- as.numeric(args[1])

## Get file paths
files <- list.files(pattern="*.count")
genomecount <- length(files)

merged <- fread(args[2], header=T, sep="\t")

## Read in query kmer file
# query_kmer <- fread(args[2], header=F, sep="\t")
# names(query_kmer) <- "kmer"


# ## Function to read in all files and add genome id as header
# add_header <- function(myfile) {
# 	mydf <- fread(myfile, header=F, sep="\t")
# 	newcol <- myfile %>%
# 		str_replace(".count", "")
# 	names(mydf) <- newcol
# 	mydf
# }

# ## Add headers to each output count file
# plan(multisession, workers=8)

# if (genomecount > 50) {
# 	counts <- list()
# 	for (i in files) {
# 		newdf <- add_header(i)
# 		counts <- append(counts, newdf)
# 	}
# 	counts_split <- split(counts, ceiling(seq_along(counts)/50))
# 	counts_final <- bind_cols(future_lapply(counts_split, bind_cols))
# } else {
# 	counts_final <- future_lapply(files, add_header)
# }


## Merge kmer query list and all count vectors (row ordering should be identical in all dfs)
#merged <- bind_cols(query_kmer, counts_final)


## Filter to only kmers present above certain threshold defined
finaldf <- merged %>%
	# mutate(count_genomes_present=rowSums(dplyr::select(., -starts_with("kmer")) !=0)) %>%
	# mutate(sum_across_genomes=rowSums(dplyr::select(., -starts_with("kmer")))) %>%
	# mutate(perc_total = count_genomes_present/genomecount*100) %>%
	# mutate(avg_copy_number = sum_across_genomes/genomecount) %>%
	filter(perc_total >= threshold) %>%
	#select(kmer, count_genomes_present, avg_copy_number, perc_total) %>%
	arrange(desc(perc_total))


if (dim(finaldf)[1] == 0) {
        stop("No kmers found at specified threshold!")
}


## If specified, add PAM to each k-mer passing defined threshold
if (args[4] == "yes") {		## this is for k=22 case
		mykmers <- finaldf %>% 
		pull(kmer)
	mykmers_vec <- c(paste("TTTA", mykmers, sep=""),
		paste("TTTC", mykmers, sep=""),
		paste("TTTG", mykmers, sep=""),
		paste("TTTT", mykmers, sep=""))	
	
	## Write kmer query list to file
	fwrite(finaldf, file=paste(args[3], "query_withpam_lookup.tsv", sep="/"), col.names=T, row.names=F, quote=F, sep="\t")
	fwrite(as.data.frame(mykmers_vec), file=paste(args[3], "query_withpam.txt", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")

} else if (args[4] == "no") {	## this would be for k=kmerpam case
	mykmers <- finaldf %>% 
		#filter(avg_copy_number <=1) %>% 
		separate(kmer, into=c("pam", "guide_to_order"), sep=4, remove=T)
		# mutate(actual_pam = case_when(pam=="TAAA" ~ "TTTA", 
		# 		      pam=="GAAA" ~ "TTTC",
		# 		      pam=="CAAA" ~ "TTTG",
		# 		      pam=="AAAA" ~ "TTTT",
		# 		      TRUE ~ "NA")) %>%
		#relocate(actual_pam, .after=pam)
	#mykmers$guide_to_order <- sapply(mykmers$kmer, function(x) as.character(reverseComplement(DNAString(x))))
	#mykmers <- mykmers %>%
	#	relocate(guide_to_order, .after=kmer)
		
	mykmers_nopam <- unique(mykmers$guide_to_order)
	fwrite(as.data.frame(mykmers_nopam), file=paste(args[3], "kmers_withpam_offtarget_query.txt", sep="/"), col.names=F, row.names=F, quote=F, sep="\t")
	fwrite(as.data.frame(mykmers),file=paste(args[3], "kmers_withpam_lookup.tsv", sep="/"), col.names=T, row.names=F, quote=F, sep="\t")

}
