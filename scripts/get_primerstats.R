#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(argparse)
#library(ggplot2)

#$1: primersearch target parsed output file
#$2: primersearch nontarget parsed output file
#$3: primersearch input file
#$4: target genomes count
#$5: non-target genomes count
#$6: primer3 combined output file
#$7: primer stats output file name
#$8: guide sensitivity output file
#$9: guide specificity output file
#$10: gene annotation file (for pangenome)
#$11: workflow (kmer or pangenome)

args <- commandArgs(trailingOnly=TRUE)

targets <- fread(args[1], header=T, sep="\t")

if (file.exists(args[2])) {
	nontargets <- fread(args[2], header=T, sep="\t")
	} else {
		cat("No primer hits to non-target genomes\n")
		nontargets <- NULL
	}


primersearch_input <- fread(args[3], header=F, sep="\t")
names(primersearch_input) <- c("primerset_uniq", "fwd_primer", "rev_primer")

primer_lengths <- primersearch_input %>%
mutate(sum_primer_lengths=nchar(fwd_primer)+nchar(rev_primer)) %>%
select(primerset_uniq, sum_primer_lengths)

## Add info on sum primer mismatches to targets and nontargets df
targets <- targets %>% 
left_join(primer_lengths, by=c("primer_set"="primerset_uniq")) %>%
mutate(sum_primer_mismatches=round(sum_primer_lengths-(sum_primer_lengths*ispcr_score/1000))) %>%
select(-sum_primer_lengths)

if (!is.null(nontargets)) {
	nontargets <- nontargets %>% 
	left_join(primer_lengths, by=c("primer_set"="primerset_uniq")) %>%
	mutate(sum_primer_mismatches=round(sum_primer_lengths-(sum_primer_lengths*ispcr_score/1000))) %>%
	select(-sum_primer_lengths)
}


targetcount <- as.numeric(args[4])
nontargetcount <- as.numeric(args[5])
primer_output <- fread(args[6], header=T, sep="\t")
guide_output <- fread(args[8], header=T, sep="\t")

if (file.exists(args[9])) {
	guide_specificity <- fread(args[9], header=T, sep="\t")
	guide_specificity <- select(guide_specificity, !contains("avg"))
	} else {
		guide_specificity <- NULL
	}

if (file.exists(args[10])) {
	annot <- fread(args[10], header=T, sep="\t")
	} else {
		annot <- NULL
	}


## Handle differently for pangenome and kmer
if (args[11] == "kmer") {
	primer_output <- primer_output %>%
	mutate(primerset_new = paste("primer", primer_set, pam, guide_set, sep="|"))
	
	if (!is.null(guide_specificity)) {
		guide_output <- guide_output %>%
		mutate(guide_target_prevalence=unique_genome/targetcount*100) %>%
		select(-unique_genome) %>% 
		left_join(guide_specificity, by="guide") %>%
		mutate(guide_offtarget_prevalence_max5mm=unique_genome/nontargetcount*100) %>%
		select(-unique_genome)
		# mutate(guide_ncpam_offtarget_prevalence_max5mm=ifelse("unique_genome_ncpam" %in% colnames(.),
		# 	unique_genome_ncpam/nontargetcount*100, NA)) %>%
		# select(if (exists('unique_genome_ncpam')) -unique_genome_ncpam else everything())
		} else {
			guide_output <- guide_output %>%
			mutate(guide_target_prevalence=unique_genome/targetcount*100) %>% 
			select(-unique_genome)
			# mutate(guide_ncpam_offtarget_prevalence_max5mm=ifelse("unique_genome_ncpam" %in% colnames(.), 
			# 	unique_genome_ncpam/nontargetcount*100, NA)) %>%
			# select(if (exists('unique_genome_ncpam')) -unique_genome_ncpam else everything())
		}

} else {   #pangenome
	primer_output <- primer_output %>%
	mutate(primerset_new = paste("primer", guide_set, primer_set, sep="|"))
	primer_output$primerset_new <- gsub("--", "|", primer_output$primerset_new) 
	primer_output$primerset_new <- gsub("-gset", "|gset", primer_output$primerset_new)

	if (!is.null(guide_specificity)) {
		guide_output <- guide_output %>%
		mutate(guide_target_prevalence=unique_genome/targetcount*100) %>%
		select(-unique_genome) %>%
		left_join(guide_specificity, by="guide") %>%
		mutate(guide_offtarget_prevalence_max5mm=unique_genome/nontargetcount*100) %>%
		select(-unique_genome)
		# mutate(guide_ncpam_offtarget_prevalence_max5mm=ifelse("unique_genome_ncpam" %in% colnames(.), 
		# 	unique_genome_ncpam/nontargetcount*100, NA)) %>%
		# select(if (exists('unique_genome_ncpam')) -unique_genome_ncpam else everything())
		} else {
			guide_output <- guide_output %>%
			mutate(guide_target_prevalence=unique_genome/targetcount*100) %>%
			select(-unique_genome)
			# mutate(guide_ncpam_offtarget_prevalence_max5mm=ifelse("unique_genome_ncpam" %in% colnames(.),
			# 	unique_genome_ncpam/nontargetcount*100, NA)) %>%
			# select(if (exists('unique_genome_ncpam')) -unique_genome_ncpam else everything())
		}
	}


## Link dereplicated primers to ori
#ori_primers <- select(primer_output, primerset_uniq, fwd_primer, rev_primer)
new_lookup <- primer_output %>% 
left_join(primersearch_input, by=c("fwd_primer"="fwd_primer", "rev_primer"="rev_primer"))


## Check if nontarget parsed file exists

# targets_summary <- targets %>% 
# 			group_by(primer_set) %>% 
# 			summarize(target_exactmatch_count = n_distinct(genome[fwd_mismatches==0 & rev_mismatches==0]),
# 				target_1mismatch_count = n_distinct(genome[fwd_mismatches<=1 & rev_mismatches<=1]),
# 				target_exactmatch = target_exactmatch_count/targetcount,
# 				target_1mismatch = target_1mismatch_count/targetcount,
# 				avg_target_amplicon_length = mean(amplicon_length)) %>%
# 				select(-ends_with("count"))


targets_summary <- targets %>% 
filter(sum_primer_mismatches<=8) %>%
group_by(primer_set) %>% 
summarize(primer_target_0mm_count = n_distinct(genome[sum_primer_mismatches==0]),
	primer_target_2mm_count = n_distinct(genome[sum_primer_mismatches<=2]),
	primer_target_prevalence_0mm = primer_target_0mm_count/targetcount*100,
	primer_target_prevalence_2mm = primer_target_2mm_count/targetcount*100,
	avg_amplicon_num_target = n()/n_distinct(genome),
	avg_amplicon_length_target = mean(amplicon_length)) %>%
select(-ends_with("count"))



if (!is.null(nontargets)) {
	nontargets_summary <- nontargets %>% 
	filter(sum_primer_mismatches<=8) %>%
	group_by(primer_set) %>% 
	summarize(primer_nontarget_0mm_count = n_distinct(genome[sum_primer_mismatches==0]),
		primer_nontarget_2mm_count = n_distinct(genome[sum_primer_mismatches<=2]),
		primer_nontarget_8mm_count = n_distinct(genome[sum_primer_mismatches<=8]),
		primer_nontarget_prevalence_0mm = primer_nontarget_0mm_count/nontargetcount*100,
		primer_nontarget_prevalence_2mm = primer_nontarget_2mm_count/nontargetcount*100,
		primer_nontarget_prevalence_8mm = primer_nontarget_8mm_count/nontargetcount*100,
		avg_amplicon_num_nontarget = n()/n_distinct(genome),
		avg_amplicon_length_nontarget = mean(amplicon_length)) %>%
	select(-ends_with("count"))
	} else {
		nontargets_summary <- NULL
	}


## Merge target and nontarget output, annotation information, cas13a compatibility

if (!is.null(nontargets_summary)) {
	combined <- targets_summary %>%
	left_join(nontargets_summary, by="primer_set")

	if (args[11] == "kmer") {
		final <- new_lookup %>%
		left_join(combined, by=c("primerset_uniq"="primer_set")) %>%
		left_join(guide_output, by=c("guide_set"="guide")) %>%
		#left_join(cas13, by=c("guide_set"="guide", "primerset_uniq"="primerset")) %>%
		filter(primer_target_prevalence_0mm!=0 | primer_target_prevalence_0mm!="") %>%
		mutate(primer_cas13a_compatibility="swap_primers") %>%
		relocate(primer_cas13a_compatibility, .after=rev_primer)
		} else {
			final <- new_lookup %>%
			left_join(combined, by=c("primerset_uniq"="primer_set")) %>%
			left_join(guide_output, by=c("guide_set"="guide")) %>%
			left_join(annot, by=c("gene"="Gene")) %>%
			rename(gene_annotation=Annotation) %>%
			relocate(gene_annotation, .after=gene) %>%
			select(-probe) %>%
			mutate(primer_cas13a_compatibility=case_when(endsWith(guide_set, "f") ~ "swap primers", TRUE ~ "yes")) %>%
			filter(primer_target_prevalence_0mm!=0 | primer_target_prevalence_0mm!="") %>%
			relocate(primer_cas13a_compatibility, .after=rev_primer)	
		}

		fwrite(final, args[7], col.names=T, row.names=F, sep="\t", quote=F)

		} else {
			if (args[11] == "kmer") {
				final <- new_lookup %>%
				left_join(targets_summary, by=c("primerset_uniq"="primer_set")) %>%
				left_join(guide_output, by=c("guide_set"="guide")) %>%
				#left_join(cas13, by=c("guide_set"="guide", "primerset_uniq"="primerset")) %>%
				filter(primer_target_prevalence_0mm!=0 | primer_target_prevalence_0mm!="") %>%
				mutate(primer_cas13a_compatibility="swap_primers") %>%
				relocate(primer_cas13a_compatibility, .after=rev_primer)
				} else {
					final <- new_lookup %>%
					left_join(targets_summary, by=c("primerset_uniq"="primer_set")) %>%
					left_join(guide_output, by=c("guide_set"="guide")) %>%
					left_join(annot, by=c("gene"="Gene")) %>%
					rename(gene_annotation=Annotation) %>%
					relocate(gene_annotation, .after=gene) %>%
					select(-probe) %>% 
					mutate(primer_cas13a_compatibility=case_when(endsWith(guide_set, "f") ~ "swap primers", TRUE ~ "yes")) %>%
					filter(primer_target_prevalence_0mm!=0 | primer_target_prevalence_0mm!="") %>%
					relocate(primer_cas13a_compatibility, .after=rev_primer)
				}
				
				fwrite(final, args[7], col.names=T, row.names=F, sep="\t", quote=F)
			}


