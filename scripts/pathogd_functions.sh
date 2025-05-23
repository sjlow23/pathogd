#!/usr/bin/env bash

# Define functions
check_assemblies_taxid() {
	# This module is only for selecting or subsampling genomes for input into one of the "user-provided genomes" module
	# Make genomes_target and genomes_offtarget directories for downloading genomes
	# Run wget script within those directories to download subsampled genomes

	cd "$OUTPUT"
	logger "Checking number of genome assemblies available for '$target_species'" 
	MYORG="$domain"
	
	if [[ "$assembly_level" == "complete" ]]; then
		MYLEVEL="Complete Genome"
	elif [[ "$assembly_level" == "chromosome" ]]; then
		MYLEVEL="Chromosome"
	elif [[ "$assembly_level" == "scaffold" ]]; then
		MYLEVEL="Scaffold"
	elif [[ "$assembly_level" == "contig" ]]; then
		MYLEVEL="Contig"
	else
		MYLEVEL="all"
	fi
	
	echo "$MYLEVEL"

	if [ ! -f "genbank.txt" ] || [ ! -f "refseq.txt" ]; then  # condition to avoid redownload of summaries for testing - can be removed if behaviour not desired
			wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/$MYORG/assembly_summary.txt --output-document=genbank.txt
			wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/$MYORG/assembly_summary.txt --output-document=refseq.txt
	fi

	if [[ -n "$target_taxid" ]]; then
		# Get all children species taxids under user-provided parent taxids
		taxonkit --data-dir $TAXONKIT_DB list --ids "$target_species_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > target_sptaxids.txt
	fi

	if [[ -n "$offtarget_taxid" ]]; then
		taxonkit --data-dir $TAXONKIT_DB list --ids "$offtarget_species_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > nontarget_sptaxids.txt
	fi

	# Subsample where necessary and generate download script
	for db in genbank refseq; do

		if [ $db == "genbank" ]; then
			mydb="GenBank"
		else
			mydb="RefSeq"
		fi

		# Target genomes
		mkdir -p ${db}
		cd ${db}

		# If using taxids
		if [[ -n "$target_taxid" ]]; then
			awk 'FNR==NR { a[$1]; next } ($7 in a)' "$OUTPUT"/target_sptaxids.txt "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > target_${db}.txt

		# If using species name
		else
			# if multiple species present
			if [[ "$target_species" == *,* ]]; then
				target_taxa=$(sed 's/,/|/g' <<< "$target_species")
				egrep -w "$target_taxa" "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > target_${db}.txt
			else
				grep -w "$target_species" "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > target_${db}.txt
				
			fi
		fi
	
		target_count=$(wc -l <target_${db}.txt)
		logger ""$target_count" '$assembly_level' "$mydb" assemblies found for target taxa"

		# Subsample target genomes to maximum of 1000 per species (using species taxid column 7) if overrepresented
		cut -f7 target_${db}.txt | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/1' | 
		awk -F "\t" '{ if ($1<1000) print $2 > "keep.txt"; else print $2 > "sample.txt" }'

		if [[ -s sample.txt || -s keep.txt ]]; then
			logger "Subsampling to a maximum of 1000 target genomes per target species"

			# Function to generate wget commands
			generate_wget_commands() {
				awk -F "/" '{ print "wget "$0"/"$NF"_genomic.fna.gz -P genomes_target" }'
			}

			if [[ -s sample.txt ]]; then
				awk '{print $1}' sample.txt | xargs -I {} sh -c "awk -F '\t' -v id={} '\$7 == id {print}' target_${db}.txt | shuf -n 1000" | cut -f20 | generate_wget_commands > wget_target_${db}.sh
			fi
			if [[ -s keep.txt ]]; then
				awk 'FNR==NR {arr[$1]=1; next} $7 in arr {print $0}' keep.txt FS="\t" target_${db}.txt | cut -f20 | generate_wget_commands >> wget_target_${db}.sh
			fi

			# Cleanup
			[[ -s sample.txt ]] && rm sample.txt
			[[ -s keep.txt ]] && rm keep.txt
		fi

		# Offtarget genomes
		# If using taxids
		if [[ -n "$offtarget_taxid" ]]; then
			awk 'FNR==NR { a[$1]; next } ($7 in a)' "$OUTPUT"/nontarget_sptaxids.txt "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > offtarget_${db}_tmp.txt
			
			# rm "$OUTPUT"/nontarget_sptaxids.txt  # fix: this removes the file in the second iteration (refseq) and break the workflow
		# if using species name
		else
			if [[ "$offtarget" == *,* ]]; then
				#offtarget_taxa=$(echo "$offtarget_taxa" | sed 's/,/|/g')
				offtarget_taxa=$(sed 's/,/|/g' <<< "$offtarget_species")
				egrep -w "$offtarget_taxa" "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > offtarget_${db}_tmp.txt
				
			else
				grep -w "$offtarget_species" "$OUTPUT"/${db}.txt | awk -F "\t" '{if ("'"$MYLEVEL"'" == "all" || $12 == "'"$MYLEVEL"'") print}' > offtarget_${db}_tmp.txt
			fi
		fi

		# Remove overlapping genomes
		grep -vwF -f target_${db}.txt offtarget_${db}_tmp.txt > offtarget_${db}.txt
		sed -i 's/\"//g' offtarget_${db}.txt
		rm offtarget_${db}_tmp.txt

		offtarget_count=$(wc -l <"offtarget_${db}.txt")
		logger ""$offtarget_count" '$assembly_level' "$mydb" assemblies found for non-target taxa"

		# Subsample offtarget genomes to maximum of 100 per species (using species_taxid column 7) if overrepresented
		cut -f7 offtarget_${db}.txt | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/1' | 
		awk -F "\t" '{ if ($1<100) print $2 > "keep.txt"; else print $2 > "sample.txt" }'

		if [[ -s sample.txt || -s keep.txt ]]; then
			logger "Subsampling to a maximum of 100 genomes per non-target species"

			# Function to generate wget commands
			generate_wget_commands() {
				awk -F "/" '{ print "wget "$0"/"$NF"_genomic.fna.gz -P genomes_offtarget" }'
			}

			if [[ -s sample.txt ]]; then
				awk '{print $1}' sample.txt | xargs -I {} sh -c "awk -F '\t' -v id={} '\$7 == id {print}' offtarget_${db}.txt | shuf -n 100" | cut -f20 | generate_wget_commands > wget_offtarget_${db}.sh
			fi
			if [[ -s keep.txt ]]; then
				awk 'FNR==NR {arr[$1]=1; next} $7 in arr {print $0}' keep.txt FS="\t" offtarget_${db}.txt | cut -f20 | generate_wget_commands >> wget_offtarget_${db}.sh
			fi

			# Cleanup: remove input files only if wget script is generated
			[[ -s sample.txt ]] && rm sample.txt
			[[ -s keep.txt ]] && rm keep.txt
			
		fi
		cd "$OUTPUT"

	done
		rm genbank.txt refseq.txt
}



download_target_genomes() {
	
	cd $OUTPUT
	mkdir -p $GENOMES_TARGET $OUT_METADATA

	if [[ $db == "genbank" ]] || [[ $db == "refseq" ]]; then
		if [[ -n "$target_taxid" ]]; then
			taxonkit --data-dir "$TAXONKIT_DB" list --ids "$target_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > target_sptaxids.txt

			logger "Downloading target genomes using target species taxids"
			species_option="--species-taxids target_sptaxids.txt"
		else
			logger "Downloading target genomes using target species name"
			species_option="--genera \"$target\""
		fi

		# Common ncbi-genome-download command
		ncbi-genome-download \
			--section "$db" \
			$species_option \
			--flat-output \
			--formats fasta \
			--output-folder "$GENOMES_TARGET" \
			--metadata-table "$OUT_METADATA/metadata_target.tsv" \
			--assembly-levels "$assembly_level" \
			--parallel "$MAX_JOBS" \
			--retries 5 \
			$domain ||
			error_exit "No $db target genomes found!"	

	elif [[ $db == "both" ]]; then
		no_genomes=0

		download_genomes() {
			local section=$1
			local species_option=$2
			local metadata_file=$3

			logger "Downloading target genomes from NCBI $section using $species_option"
			ncbi-genome-download \
				--section "$section" \
				$species_option \
				--flat-output \
				--formats fasta \
				--output-folder "$GENOMES_TARGET" \
				--metadata-table "$metadata_file" \
				--assembly-levels "$assembly_level" \
				--parallel "$MAX_JOBS" \
				--retries 5 \
				$domain || ((no_genomes++))
		}

		if [[ -n "$target_taxid" ]]; then
			taxonkit --data-dir "$TAXONKIT_DB" list --ids "$target_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > target_sptaxids.txt
			
			download_genomes "refseq" "--species-taxids target_sptaxids.txt" "$OUT_METADATA/metadata_target_refseq.tsv"
			download_genomes "genbank" "--species-taxids target_sptaxids.txt" "$OUT_METADATA/metadata_target_genbank.tsv"

			if [[ $no_genomes -eq 2 ]]; then
				error_exit "No target genomes found using species taxids!"
			fi
		else
			download_genomes "refseq" "--genera \"$target\"" "$OUT_METADATA/metadata_target_refseq.tsv"
			download_genomes "genbank" "--genera \"$target\"" "$OUT_METADATA/metadata_target_genbank.tsv"

			if [[ $no_genomes -eq 2 ]]; then
				error_exit "No target genomes found using species names!"
			fi
		fi
	fi
}


download_nontarget_genomes() {

	cd $OUTPUT
	mkdir -p $GENOMES_OFFTARGET $OUT_METADATA

	if [[ $db == "genbank" ]] || [[ $db == "refseq" ]]; then
		if [[ -n "$offtarget_taxid" ]]; then
			taxonkit --data-dir "$TAXONKIT_DB" list --ids "$offtarget_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > nontarget_sptaxids.txt

			logger "Downloading non-target genomes using non-target species taxids"
			species_option="--species-taxids nontarget_sptaxids.txt"
		else
			logger "Downloading non-target genomes using non-target species name"
			species_option="--genera \"$offtarget\""
		fi

	# Common ncbi-genome-download command
		ncbi-genome-download \
			--section "$db" \
			$species_option \
			--flat-output \
			--formats fasta \
			--output-folder "$GENOMES_OFFTARGET" \
			--metadata-table "$OUT_METADATA/metadata_offtarget.tsv" \
			--assembly-levels "$assembly_level" \
			--parallel "$MAX_JOBS" \
			--retries 5 \
			$domain ||
			error_exit "No $db non-target genomes found!"	
	
	elif [[ $db == "both" ]]; then
		no_genomes=0

		download_genomes() {
			local section=$1
			local species_option=$2
			local metadata_file=$3

			logger "Downloading non-target genomes from NCBI $section using $species_option"
			ncbi-genome-download \
				--section "$section" \
				$species_option \
				--flat-output \
				--formats fasta \
				--output-folder "$GENOMES_OFFTARGET" \
				--metadata-table "$metadata_file" \
				--assembly-levels "$assembly_level" \
				--parallel "$MAX_JOBS" \
				--retries 5 \
				$domain || ((no_genomes++))
		}

		if [[ -n "$offtarget_taxid" ]]; then
			taxonkit --data-dir "$TAXONKIT_DB" list --ids "$offtarget_taxid" -r | grep "\[species\]" | awk '{ print $1 }' > nontarget_sptaxids.txt
			
			download_genomes "refseq" "--species-taxids nontarget_sptaxids.txt" "$OUT_METADATA/metadata_offtarget_refseq.tsv"
			download_genomes "genbank" "--species-taxids nontarget_sptaxids.txt" "$OUT_METADATA/metadata_offtarget_genbank.tsv"

			if [[ $no_genomes -eq 2 ]]; then
				error_exit "No non-target genomes found using species taxids!"
			fi
		else
			download_genomes "refseq" "--genera \"$target\"" "$OUT_METADATA/metadata_offtarget_refseq.tsv"
			download_genomes "genbank" "--genera \"$target\"" "$OUT_METADATA/metadata_offtarget_genbank.tsv"

			if [[ $no_genomes -eq 2 ]]; then
				error_exit "No non-target genomes found using species names!"
			fi
		fi
	fi
}


remove_overlap_genomes() {
	## Genomes assumed to be compressed
	cd $OUTPUT
	
	logger "Removing overlapping genomes"
	local target=($(ls $GENOMES_TARGET/*genomic* | awk -F "/" '{ print $NF }' | sort))
	local offtarget=($(ls $GENOMES_OFFTARGET/*genomic* | awk -F "/" '{ print $NF }' | sort))
	
	echo ${target[@]} ${offtarget[@]} | sed 's/ /\n/g' | sort | uniq -d > remove

	awk '{ print "'"$GENOMES_OFFTARGET"'/"$0"" }' remove | xargs -r rm
	rm remove

	## Report number of target and non-target genomes downloaded
	local count_target=$(find $GENOMES_TARGET/ -type f -name "*.gz" -or -name "*.fna" | wc -l)
	logger ""$count_target" target genomes found"

	local count_offtarget=$(find $GENOMES_OFFTARGET/ -type f -name "*.gz" -or -name "*.fna" | wc -l)
	logger ""$count_offtarget" non-target genomes found"
}



decompress_genomes() {
	process_genomes() {
		local genome_type=$1
		local genomes_dir=$2
		local output_file=$3
		local log_msg=$4

		cd "$genomes_dir" || error_exit "Failed to change directory to $genomes_dir"

		logger "$log_msg"
		local count=$(find . \( -type f -o -type l \) -name "*.gz" -or -name "*.fna" | wc -l)

		[[ "$count" -ge 2 ]] || error_exit "Only 1 $genome_type genome found - cannot proceed!"

		find . -type f -name "*.gz" | parallel -j "$MAX_JOBS" gunzip

		local ori=($(ls -d *.fna))
		local new=($(printf '%s\n' "${ori[@]}" | awk -F "_" '{ if ($0 ~ /^GC/) print $1"_"$2"_genomic.fna"; else print $0 }'))

		for i in "${!ori[@]}"; do
			if [[ "${ori[i]}" != "${new[i]}" ]]; then
				mv -n "${ori[i]}" "${new[i]}"
			fi
		done

		find . -maxdepth 1 \( -type f -o -type l \) -name "*.fna" -exec basename {} _genomic.fna \; > "$output_file"
		cd "$OUTPUT" || error_exit "Failed to change directory to $OUTPUT"
	}

	if [[ $1 == "target" ]]; then
		process_genomes "target" "$GENOMES_TARGET" "$OUTPUT/genomes_target.txt" "Decompressing target genomes"

	elif [[ $1 == "offtarget" ]]; then
		process_genomes "non-target" "$GENOMES_OFFTARGET" "$OUTPUT/genomes_offtarget.txt" "Decompressing non-target genomes"
	fi
}


prepare_human_genome() {
	## Kmer counting is dependent on user kmer size input
	cd $MYDIR

	if [[ ! -d "$DATABASE" || ! -s "$DATABASE"/GCF_000001405.40_GRCh38.p14_genomic.fna ]]
	then
		mkdir -p "$DATABASE"
		cd "$DATABASE"
		wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		gunzip -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		cd $OUTPUT
	fi

	## Check if kmer list present	
	cd $DATABASE
	
	local basefile="GCF_000001405.40_GRCh38.p14_genomic_"$kmer".list"
	local extfile="GCF_000001405.40_GRCh38.p14_genomic_"$kmerpam".list"
	local ref="GCF_000001405.40_GRCh38.p14_genomic.fna"
	local prefix=$(echo "${ref%.*}")

	if [[ -s $basefile && -s $extfile ]]
	then
		logger "k-mer lists for human genome for selected 'k' previously generated"
	elif [[ ! -s $basefile && ! -s $extfile ]]
	then
		logger "Counting k="$kmer" and k="$kmerpam" in human genome"
		glistmaker $ref --num_threads "$MAX_JOBS" -w "$kmer" -o $prefix
		glistmaker $ref --num_threads "$MAX_JOBS" -w "$kmerpam" -o $prefix

	elif [[ -s $basefile && ! -s $extfile ]]
	then
		logger "Counting k="$kmerpam" in human genome"
		glistmaker $ref --num_threads "$MAX_JOBS" -w "$kmerpam" -o $prefix
	elif [[ -s $extfile && ! -s $basefile ]]
	then
		logger "Counting k="$kmer" in human genome"
		glistmaker $ref --num_threads "$MAX_JOBS" -w "$kmer" -o $prefix
	fi

	cd $OUTPUT
}


subsample_genomes() {
	cd $OUTPUT

	local oricount=$(wc -l <$OUTPUT/genomes_target.txt)

	[[ -s "$OUTPUT"/genomes_target.txt ]] || { error_exit "No target genomes found"; return; }
	[[ -s "$OUTPUT"/genomes_target_subsample.txt ]] && return

	if [[ "$oricount" -gt "$subsample" ]]
	then
		logger "Subsampling "$subsample" genomes from original list of "$oricount" genomes"
		shuf -n "$subsample" $OUTPUT/genomes_target.txt > $OUTPUT/genomes_target_subsample.txt
	elif [[ "$oricount" -le "$subsample" ]]
	then
		logger "Using all "$oricount" target genomes"
		cp "$OUTPUT"/genomes_target.txt "$OUTPUT"/genomes_target_subsample.txt
	fi
}


get_pangenome() {
	cd $OUTPUT
	mkdir -p $PANGENOME/prokka

	cd $PANGENOME

	[[ -s "$OUTPUT"/genomes_target.txt ]] || error_exit "Target genomes not found!"
	[[ -s "$OUTPUT"/genomes_offtarget.txt ]] || error_exit "Non-target genomes not found!"

	if [[ "$subsample" == "no" ]]
	then
		local myfile=$OUTPUT/genomes_target.txt
	else
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	# If conda environment name provided
	set +u
	
	if [[ "$ref_annot" == "yes" ]]
	then
		logger "Annotating genes for target species with Prokka using provided genome"
		#cp $MYDIR/$REF prokka/
		#local ref=$(basename "$MYDIR"/"$REF")
		prodigal -i $MYDIR/$REF -t prokka/ref.trn
			
		#conda activate prokka-1.14.6
		
		cat $myfile | \
		parallel -j "$MAX_PARALLEL_JOBS" 'prokka --quiet --outdir prokka/"{}"_out --cpus 8 --locustag "{}" --prefix "{}" --proteins '"$MYDIR"'/'"$REF"' --prodigaltf prokka/ref.trn '"$GENOMES_TARGET"'/"{}"_genomic.fna; \
		rm prokka/"{}"_out/proteins.*'

	elif [[ "$ref_annot" == "no" ]]
	then
		logger "Annotating genes for target species with Prokka"

		cat $myfile | \
		parallel -j "$MAX_PARALLEL_JOBS" 'prokka --quiet --outdir prokka/"{}"_out --cpus 8 --locustag "{}" --prefix "{}" '"$GENOMES_TARGET"'/"{}"_genomic.fna'
	fi

	logger "Estimating pangenome for '$target_species'"

	if [[ -n "${roary_env:-}" ]]
	then
		mamba activate $roary_env
		roary -cd 90 -s -z --mafft -p "$MAX_JOBS" -f roary prokka/*_out/*.gff
		mamba deactivate
	else
		roary -cd 90 -s -z --mafft -p "$MAX_JOBS" -f roary prokka/*_out/*.gff
		mamba deactivate
	fi
	
	set +u

	## Get list of conserved genes (at least 90% prevalence)]
	if [[ "$subsample" == "no" ]]
	then
		local cutoff=$(wc -l $OUTPUT/genomes_target.txt | awk '{ print $1*0.9 }' | awk '{ print int($1+0.5)}')
	else
		local cutoff=$(wc -l $OUTPUT/genomes_target_subsample.txt | awk '{ print $1*0.9 }' | awk '{ print int($1+0.5)}')
	fi

	csvtk filter --filter "No. isolates>="$cutoff"" roary/gene_presence_absence.csv | \
	csvtk cut -f1,3 --out-delimiter $'\t' > genes_90pct_annotation.tsv

	awk -F "\t" 'NR>1' genes_90pct_annotation.tsv | cut -f1 > conserved90pct_genes.txt

	## Extract genes for primer design
	find prokka/*_out -name "*.ffn" -exec cat {} + > allnt.fna
	cdbfasta allnt.fna

	cd $OUTPUT

}


species_specific_genes() {
	## Cluster all prokka target proteins in sets
	## Cluster all prodigal nontarget proteins in sets
	## Cluster nontarget and target representative proteins
	## Eliminate target proteins for clusters where both target and nontarget proteins are present
	## Align remaining protein-coding genes from Roary output

	mkdir -p $PANGENOME/prodigal_nontarget $PANGENOME/prodigal_nontarget/combined_batches $PANGENOME/species_specific/target $PANGENOME/species_specific/nontarget
	cd $PANGENOME

	[[ -s allnt.fna ]] || error_exit "Combined genes file "allnt.fna" not found!"

	## Get genes in non-target genomes for clustering
	logger "Annotating genes for nontarget species with Prodigal"
	
	# cat $OUTPUT/genomes_offtarget.txt | \
	# parallel -j "$MAX_PARALLEL_JOBS" 'prodigal -a prodigal_nontarget/"{}".faa -q -c -m -p single -f sco -q -i '"$GENOMES_OFFTARGET"'/"{}"_genomic.fna'
	
	counter=0
	cat $OUTPUT/genomes_offtarget.txt | xargs -n 200 | \
	while read -r genome
	do
		counter=$((counter + 1))
		echo "$genome" | \
		tr " " "\n" | \
		parallel -j "$MAX_JOBS" 'prodigal -a prodigal_nontarget/"{}".faa -q -c -m -f sco -p meta -o prodigal_nontarget/tmp -q -i '"$GENOMES_OFFTARGET"'/"{}"_genomic.fna'

		cat prodigal_nontarget/*.faa > prodigal_nontarget/combined_batches/batch_"$counter".faa
		[[ -s prodigal_nontarget/combined_batches/batch_"$counter".faa ]] && rm prodigal_nontarget/*.faa

	done

	#mkdir -p species_specific/target species_specific/nontarget

	## Combine target and nontarget proteins into two files
	cat prodigal_nontarget/combined_batches/batch_*.faa > species_specific/nontarget/allnontarget.faa
	rm -rf prodigal_nontarget

	cat prokka/*_out/*.faa > species_specific/target/alltarget.faa

	# if [[ "$subsample" == "no" ]]
	# then
	# 	for i in $(cat $OUTPUT/genomes_target.txt); do cat prokka/"$i"_out/"$i".faa >> species_specific/target/alltarget.faa; done
	# else 
	# 	for i in $(cat $OUTPUT/genomes_target_subsample.txt); do cat prokka/"$i"_out/"$i".faa >> species_specific/target/alltarget.faa; done
	# fi

	#for i in $(cat $OUTPUT/genomes_offtarget.txt); do cat prodigal_nontarget/"$i".faa >> species_specific/nontarget/allnontarget.faa; done


	#find prodigal_nontarget -type f -name "*.faa" | grep -v allnontarget | xargs -L1 rm
	#rm -rf prodigal_nontarget

	## Split protein sequences into 500k chunks & cluster
	## Nontarget
	logger "Clustering genes in nontarget genomes"

	cd $PANGENOME/species_specific/nontarget

	seqkit split2 -j "$MAX_JOBS" -s 500000 allnontarget.faa -O .
	for i in allnontarget.part*faa; do cd-hit -i "$i" -M 32000 -c 0.95 -aL 0.8 -d 1000 -T 16 -o "$i".cdhit; done
	cat allnontarget.part*.cdhit >> allnontarget_round1.cdhit

	cd-hit -i allnontarget_round1.cdhit -M 32000 -c 0.95 -aL 0.8 -d 1000 -T 16 -o allnontarget_round2.cdhit
	rm allnontarget.part*.cdhit

	## Target
	logger "Clustering genes in target genomes"

	cd $PANGENOME/species_specific/target
	
	seqkit split2 -j "$MAX_JOBS" -s 500000 alltarget.faa -O .
	for i in alltarget.part*faa; do cd-hit -i "$i" -M 32000 -c 0.95 -aL 0.8 -d 1000 -T 16 -o "$i".cdhit; done
	cat alltarget.part*.cdhit >> alltarget_round1.cdhit
	cd-hit -i alltarget_round1.cdhit -M 32000 -c 0.95 -aL 0.8 -d 1000 -T 16 -o alltarget_round2.cdhit
	rm alltarget.part*.cdhit

	## Now cluster target and nontarget together (90% AAI threshold used)
	cd $PANGENOME/species_specific
	cat target/alltarget_round2.cdhit nontarget/allnontarget_round2.cdhit > combined.cdhit

	cd-hit -i combined.cdhit -M 32000 -c 0.90 -aL 0.8 -d 1000 -T 16 -o final.cdhit

	clstr2txt.pl final.cdhit.clstr > final.cdhit.clstr.parsed

	
	## Define array of target and nontarget proteins
	logger "Retaining target species-specific genes"

	local targetprot=($(grep ">" target/alltarget_round2.cdhit | awk '{ print $1 }' | sed 's/>//g'))
	local nontargetprot=($(grep ">" nontarget/allnontarget_round2.cdhit | awk '{ print $1 }' | sed 's/>//g'))

	cat final.cdhit.clstr.parsed | csvtk summary -t -f id:uniq,id:countunique -g clstr | awk -F "\t" '$3!=1 {print $1, $2 }' OFS="\t" > cluster_summary.tsv

	## Search array elements in cluster summary file
	local clusters_target=($(grep -wFf <(printf "%s\n" "${targetprot[@]}") cluster_summary.tsv | cut -f1 | sort | uniq))
	local clusters_nontarget=($(grep -wFf <(printf "%s\n" "${nontargetprot[@]}") cluster_summary.tsv | cut -f1 | sort | uniq))

	echo ${clusters_target[@]} ${clusters_nontarget[@]} | sed 's/ /\n/g' | sort | uniq -d > clusters_overlap.txt

	## Get list of target genes present in clusters overlap to eliminate
	local remove=($(grep -wF -f clusters_overlap.txt cluster_summary.tsv | cut -f2 | sed 's/; /\n/g'))
	echo ${targetprot[@]} ${remove[@]} | sed 's/ /\n/g' | sort | uniq -d > target_remove.txt 

	cd $OUTPUT
	
}


align_species_specific() {
	## Align only species-specific genes
	cd $PANGENOME

	grep -wF -f species_specific/target_remove.txt roary/gene_presence_absence.csv | awk -F "," '{ print $1 }' | sed 's/^"//g' | sed 's/"$//g' > genes_remove.txt
	grep -vwF -f genes_remove.txt conserved90pct_genes.txt > conserved90pct_genes_specific.txt

	[[ -s conserved90pct_genes_specific.txt ]] || error_exit "No target-specific genes found!"

	rm -rf species_specific

	mkdir -p gene_alignments
	
	while read gene
	do
		newgene="`echo $gene | sed 's/\///g'`"
		grep -w "$gene" roary/clustered_proteins | cut -f2- | sed 's/\t/\n/g' > gene_alignments/"$newgene"_ids.txt
		cat gene_alignments/"$newgene"_ids.txt | cdbyank allnt.fna.cidx > gene_alignments/"$newgene".fna
		rm gene_alignments/"$newgene"_ids.txt
	done < conserved90pct_genes_specific.txt

	logger "Aligning core species-specific genes"
	cat conserved90pct_genes_specific.txt | \
	parallel -j "$MAX_PARALLEL_JOBS" 'mafft --auto --thread '"$MAX_THREADS_PER_JOB"' gene_alignments/"{}".fna > gene_alignments/"{}".aln; \
	cons -sequence gene_alignments/"{}".aln -name "{}" -outseq gene_alignments/"{}".cons;
	seqkit -is replace -p "^n+|n+$" -r "" gene_alignments/"{}".cons -o gene_alignments/"{}".cons.2;
	mv gene_alignments/"{}".cons.2 gene_alignments/"{}".cons'

	# Replace conserved90pct_genes_specific.txt with filtered genes
	cd gene_alignments
	tar -cf - *.aln --remove-files | pigz -9 -p 16 > alignments.tar.gz
	tar -cf - *.fna --remove-files | pigz -9 -p 16 > genes.tar.gz

	cd $OUTPUT

}


get_guides_nopam() {
	logger "Identifying potential guide RNAs in genes"
	cd $PANGENOME

	local min_flank_len=35
	local flank_len=150
	local guide_len="$kmer"

	## Calculate maximum template length
	local total_len=$((guide_len + 2 * flank_len))

	## Loop through the gene alignments
	while read gene; do
		gene_aln=$(grep -v ">" gene_alignments/"$gene".cons | tr -d "\n" | awk '{print $0}')
		count=0

		## Remove genes smaller than minimum length
		gene_len=${#gene_aln}
		if [[ "$gene_len" -lt "$template_len" ]]; then
			continue
		else
			template_len=$((guide_len + 2 * flank_len))
			# max 40 guide RNA candidates per gene
			max_positions=$((gene_len - flank_len - guide_len))
			step_size=$(( (max_positions + 20) / 40 ))
		fi

		## Process the original strand
		for (( i=min_flank_len+1; i <= gene_len - min_flank_len - guide_len; i+=step_size )); do
			count=$((count + 1))
			guide=${gene_aln:i:guide_len}

			if [[ "$guide" =~ [nN] ]]; then
				echo "Guide contains N- skipping.."
				continue
			fi

			if [[ $i -lt $flank_len ]]; then
				template=${gene_aln:0:flank_len*2+guide_len}
			else
				template=${gene_aln:i-flank_len:flank_len*2+guide_len-1}
			fi

			echo -e "$gene\t$guide\t$template" | \
				awk -F "\t" '$2 !~ "AAAAA|CCCCC|GGGGG|TTTTT"' | \
				awk -F "\t" '$2 ~ /A/ && $2 ~ /C/ && $2 ~ /G/ && $2 ~ /T/' | \
				awk -F "\t" '$3 !~ /nnnnn/' | \
				awk -F "\t" '{ print $0, "gset""'$count'f" }' OFS="\t" | \
				awk -F "\t" '{ print $1, $2, "NA", $1"--"$4, $3 }' OFS="\t" | \
				awk '{ print $0, match($5, $2)","'$kmer'"" }' OFS="\t" >> allgenes_lookup.tsv
		done

		## Process the reverse complement strand
		rev_gene_aln=$(echo "$gene_aln" | rev | tr 'ACGT' 'TGCA')  # Get the reverse complement
		countr=0  # Reset count for reverse complement

		for (( i=min_flank_len+1; i <= gene_len - min_flank_len - guide_len; i+=step_size )); do
			countr=$((countr + 1))
			guide=${rev_gene_aln:i:guide_len}

			if [[ "$guide" =~ [nN] ]]; then
				echo "Guide contains N- skipping.."
				continue
			fi

			if [[ $i -lt $flank_len ]]; then
				template=${rev_gene_aln:0:flank_len*2+guide_len}
			else
				template=${rev_gene_aln:i-flank_len:flank_len*2+guide_len-1}
			fi

			echo -e "$gene\t$guide\t$template" | \
				awk -F "\t" '$2 !~ "AAAAA|CCCCC|GGGGG|TTTTT"' | \
				awk -F "\t" '$2 ~ /A/ && $2 ~ /C/ && $2 ~ /G/ && $2 ~ /T/' | \
				awk -F "\t" '$3 !~ /nnnnn/' | \
				awk -F "\t" '{ print $0, "gset""'$countr'r" }' OFS="\t" | \
				awk -F "\t" '{ print $1, $2, "NA", $1"--"$4, $3 }' OFS="\t" | \
				awk '{ print $0, match($5, $2)","'$kmer'"" }' OFS="\t" >> allgenes_lookup.tsv
		done
	done < conserved90pct_genes_specific.txt

	[[ -s allgenes_lookup.tsv ]] || error_exit "Failed to generate "allgenes_lookup.tsv""

	cd $OUTPUT
}


get_guides_pam_cas12a() {

	logger "Identifying PAM sites in genes"
	cd $PANGENOME
	mkdir -p templates

	## Define canonical PAM sites
	local pam1="AAAA"
	local pam2="CAAA"
	local pam3="GAAA"
	local pam4="TAAA"
	local pam5="TTTA"
	local pam6="TTTC"
	local pam7="TTTG"
	local pam8="TTTT"

	#set +e
	while read gene
	do 
		seqkit seq -w 0 gene_alignments/"$gene".cons > templates/"$gene".template

		## Search for PAM sites on reverse strand
		for i in $pam1 $pam2 $pam3 $pam4
		do
			pam_found=false
			positions=$(grep -v ">" templates/"$gene".template | grep -aob "$i" | grep -oE '[0-9]+' | wc -l)
	
			# Check if positions are empty
			if [[ "$positions" -eq 0 ]]; then
				continue
			else
				pam_found=true
				
				grep -v ">" templates/"$gene".template | 
				# count positions excluding header
				grep -aob "$i" | 
				grep -oE '[0-9]+' | 
				awk '{ print $1-"'$kmer'"+1, $1 }' | 
				## Add set number for gene-pam combinations
				awk '$1>=50 {print "'$gene'", $1"-"$2, "'$i'", "'$gene'--'$i'" }' OFS="\t" | 
				awk -F "\t" 'BEGIN {prev=0; count=1} {if (prev==$NF) count++; else {count=1;; prev=$NF } print $0, "gset"count"r" }' OFS="\t" | 
				sed 's/\t/-/4' | 
				sed 's/--AAAA/--TTTT/; s/--CAAA/--TTTG/; s/--GAAA/--TTTC/; s/--TAAA/--TTTA/1' | 
				awk -F "\t" '{ if ($3=="AAAA") $3="TTTT";
				else if ($3=="CAAA") $3="TTTG";
				else if ($3=="GAAA") $3="TTTC";
				else if ($3=="TAAA") $3="TTTA"; print $0 }' OFS="\t" >> templates/"$gene".positions
			fi
		done

		## Search for PAM sites on forward strand
		for i in $pam5 $pam6 $pam7 $pam8
		do
			pam_found=false
			# Get length of gene
			genelength=$(($(grep -v ">" templates/"$gene".template | wc -m) - 1))
			
			positions=$(grep -v ">" templates/"$gene".template | grep -aob "$i" | grep -oE '[0-9]+' | wc -l)
	
			# Check if positions are empty
			if [[ "$positions" -eq 0 ]]; then
				continue
			else
				pam_found=true
				
				grep -v ">" templates/"$gene".template | 
				grep -aob "$i" | 
				grep -oE '[0-9]+' | 
				awk '{ print $1+5, $1+"'$kmer'"+5-1 }' | 
				## discard if exceeds gene length
				awk '$2<='"$genelength"'' |
				## Add set number for gene-pam combinations
				awk '$1>=50 {print "'$gene'", $1"-"$2, "'$i'", "'$gene'--'$i'" }' OFS="\t" | \
				awk -F "\t" 'BEGIN {prev=0; count=1} {if (prev==$NF) count++; else {count=1;; prev=$NF } print $0, "gset"count"f" }' OFS="\t" | \
				sed 's/\t/-/4' >> templates/"$gene".positions
			fi
		done


		if [[ -s templates/"$gene".positions ]]
		then
			##coordinates for seqkit
			awk '{ print $1, $2, $4, "1", "+" }' OFS="\t" templates/"$gene".positions | sed 's/-/\t/1' >> templates/"$gene".coords.bed

			## extract to tab-delimited format directly
			seqkit subseq -j "$MAX_JOBS" --bed templates/"$gene".coords.bed --up-stream 150 --down-stream 150 templates/"$gene".template >> templates/"$gene"_templates.fasta
			seqkit fx2tab -j "$MAX_JOBS" templates/"$gene"_templates.fasta | awk '{ print $2, $3 }' OFS="\t" >> templates/"$gene"_templates_lookup.tsv

			for j in `cut -f2 templates/"$gene".positions`; do grep -v ">" templates/"$gene".template | cut -c"$j" >> templates/"$gene".probes; done

			paste templates/"$gene".probes templates/"$gene".positions templates/"$gene"_templates_lookup.tsv | 
				#sed 's/-/\t/1' | 
				cut -f1,2,4,5,7 | 
				#awk '{ print $0, "'$KMER'" }' OFS="\t" | 
				awk '{ print $2, $1, $3, $4, $5 }' OFS="\t" | 
				awk -F "\t" '$2 !~ "n|N"' | 
				awk -F "\t" '$2 !~ "AAAAA|CCCCC|GGGGG|TTTTT"' | 
				awk '{ print $0, match($5, $2)","'$kmer'"" }' OFS="\t" | 
				# Remove guides that do not contain all nucleotide bases
				awk -F "\t" '$2 ~ /A/ && $2 ~ /C/ && $2 ~ /G/ && $2 ~ /T/' \
				>> allgenes_lookup.tsv
			rm templates/"$gene"*
		
		else 
			logger "PAM not present in ["$gene"]..skipping.."
			#rm templates/"$gene"*
		fi

	done < conserved90pct_genes_specific.txt

	#set -e
	rm -rf templates
	cd gene_alignments
	tar -cf - *.cons --remove-files | pigz -9 -p 16 > consensus.tar.gz

	cd "$PANGENOME"
	[[ -s allgenes_lookup.tsv ]] || error_exit "Failed to generate "allgenes_lookup.tsv""

	cd $OUTPUT
}



precheck_guides_pangenome() {
	cd $OUTPUT
	mkdir -p $PANGENOME/guides_precheck

	cd $PANGENOME
	logger "Preparing Primer3 input file for guides precheck"

	while IFS=$'\t' read gene oligo pam seqname sequence target
	do
			## Reverse translate actual guide sequence
				#oligo_revseq=$(echo "$oligo" | sed 's/^/>seq/g' | seqkit seq -r -t DNA -v -p | grep -v ">")

				{
					echo "SEQUENCE_ID=$seqname"
					echo "SEQUENCE_INTERNAL_OLIGO=$oligo"
					echo "PRIMER_TASK=check_primers"
					echo "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1"
					echo "PRIMER_PICK_INTERNAL_OLIGO=1"
					echo "PRIMER_INTERNAL_MIN_SIZE=18"
					echo "PRIMER_INTERNAL_MAX_SIZE=30"
					echo "PRIMER_INTERNAL_MIN_TM=40.0"
					echo "PRIMER_EXPLAIN_FLAG=1"
					echo "="
				} >> guides_precheck/allguides.input.p3

	done < allgenes_lookup.tsv

	## Filter out guides with bad properties using Primer3
	logger "Filtering low-quality guides"
	primer3_core < guides_precheck/allguides.input.p3 > guides_precheck/allguides.output.p3

	sed 's/^=/\n/g' guides_precheck/allguides.output.p3 |\
	awk -v RS= '/\nPRIMER_INTERNAL_NUM_RETURNED=1/' |\
	grep SEQUENCE_ID |\
	awk -F "=" '{ print $2 }' >> guides_precheck/guides_allpass.txt

	grep -wF -f guides_precheck/guides_allpass.txt allgenes_lookup.tsv > guides_precheck/allgenes_lookup_pass.tsv

	## Remove guides (+ and - strand) with hits to human or non-target species genomes
	awk -F "\t" '{ if ($4 ~ /f$/) print $2 }' guides_precheck/allgenes_lookup_pass.tsv > guides_precheck/guides_pass_fwd.txt
	awk -F "\t" '{ if ($4 ~ /r$/) print ">"$4"\n"$2 }' OFS="\t" guides_precheck/allgenes_lookup_pass.tsv > guides_precheck/guides_pass_rev.fasta
	seqkit seq -j "$MAX_JOBS" -r -t DNA -v -p guides_precheck/guides_pass_rev.fasta | grep -v ">" > guides_precheck/guides_pass_rev_revcomp.txt


	## Guides that pass for Cas12-only- with PAM sites added
	if [[ $CAS != "cas13" ]]
	then
		awk '{ print "TTTA"$0; print "TTTC"$0; print "TTTG"$0; print "TTTT"$0; }' guides_precheck/guides_pass_fwd.txt > \
			guides_precheck/guides_pass_fwd_pam.txt
		awk '{ print "TTTA"$0; print "TTTC"$0; print "TTTG"$0; print "TTTT"$0; }' guides_precheck/guides_pass_rev_revcomp.txt > \
			guides_precheck/guides_pass_rev_pam.txt
	
		## Get kmers from non-target genomes to exclude if containing adjacent PAM (only for Cas12)
		for i in fwd rev; do
			glistquery $KMER_OFFTARGET/combined_offtarget_"$kmerpam"_union.list -f guides_precheck/guides_pass_"$i"_pam.txt \
			--mismatch 2 | 
			awk -F "\t" '$2!=0 {print $1}' | 
			awk '{ $1 = substr($1, 5, "'$kmerpam'") } 1' | 
			sort | uniq > remove_nontargethits_"$i";

			## Remove guides present in the human genome with adjacent PAM
			glistquery $DATABASE/GCF_000001405.40_GRCh38.p14_genomic_"$kmerpam".list -f guides_precheck/guides_pass_"$i"_pam.txt \
			--mismatch 2 | \
			awk -F "\t" '$2!=0 {print $1}' |
			awk '{ $1 = substr($1, 5, "'$kmerpam'") } 1' | 
			sort | uniq > remove_humanhits_"$i";
		done
	else
		# for non-cas12 run	
		for i in fwd rev; do
			glistquery $KMER_OFFTARGET/combined_offtarget_"$kmer"_union.list -f guides_precheck/guides_pass_"$i"*.txt \
			--mismatch 2 | \
			awk -F "\t" '$2!=0 {print $1 }' |
			sort | uniq > remove_nontargethits_"$i";

			glistquery $DATABASE/GCF_000001405.40_GRCh38.p14_genomic_"$kmer".list -f guides_precheck/guides_pass_"$i"*.txt \
			--mismatch 2 | \
			awk -F "\t" '$2!=0 {print $1 }' |
			sort | uniq > remove_humanhits_"$i";
		done

	fi
			
	# Reverse complement rev hits
	cat remove_*hits*_rev | sort | uniq | nl | 
	sed 's/^ */>seq/g' | 
	sed 's/\t/\n/g' | 
	seqkit seq -r -t DNA -v -p |
	grep -v ">" > remove_rev.txt

	cat remove_*hits*_fwd remove_rev.txt | sort | uniq > remove

	if [[ -s remove ]] 
	then
		grep -vwF -f remove guides_precheck/allgenes_lookup_pass.tsv > guides_precheck/allgenes_lookup_pass_nohuman.tsv
		mv guides_precheck/allgenes_lookup_pass_nohuman.tsv guides_precheck/allgenes_lookup_pass.tsv
	fi

	## Check if there are guides remaining after filtering
	[[ -s guides_precheck/allgenes_lookup_pass.tsv ]] || error_exit "No guides remaining after excluding non-target hits"

	cd $OUTPUT
}



create_p3_input_pangenome() {
	## Get guide id
	## Get amplicon from 'amplicons' dir; assign as template sequence
	## Get probe sequence from 'final_results' dir
	## Split probe sequence into probe + pam 
	
	mkdir -p $PRIMERS/p3_input $PRIMERS/p3_output
	cd $PRIMERS
	
	logger "Creating Primer3 input file for a maximum of 40 guides per gene (max 10 per PAM)"

	## Limit number of guides per PAM to 10 (max 40 guides per gene)
	## Split into max 500 guides per file
	csvtk mutate -tH -f 1 $PANGENOME/guides_precheck/allgenes_lookup_pass.tsv | \
	csvtk mutate -Ht -f 3 | \
	sed 's/\t/-/7' | \
	shuf - | \
		#random subsampling
		awk -F "\t" 'seen[$7]++ < 10 {print}' | \
		csvtk sort -Ht -k 7:N | \
		#csvtk uniq -Ht -f 7 -n 5 | \
		cut -f1-6 > $PANGENOME/guides_precheck/allgenes_lookup_pass_max5.tsv
		split -l 500 $PANGENOME/guides_precheck/allgenes_lookup_pass_max5.tsv allgenes_subset -d

	## Choice of PCR or RPA primer design
	if [[ $AMP_METHOD == "pcr" ]]
	then
		primerminsize=18
		primermaxsize=25
		primeroptsize=22
		primermintm=60
		primermaxtm=75
		primermingc=40
		primermaxgc=60
	else
		primerminsize=30
		primermaxsize=35
		primeroptsize=32
		primermintm=0
		primermaxtm=100
		primermingc=30
		primermaxgc=70
	fi

		for i in $(ls allgenes_subset* | awk -F "_" '{ print $2 }')
		do
			while IFS=$'\t' read gene oligo pam seqname sequence target
			do
						#sequence=$(grep -w $seqname $OUTPUT/templates_lookup.tsv | cut -f2)
						{
							echo "SEQUENCE_ID=$seqname"
							echo "SEQUENCE_TEMPLATE=$sequence"
							echo "SEQUENCE_INTERNAL_OLIGO=$oligo"
							echo "PRIMER_TASK=pick_pcr_primers_and_hyb_probe"
							echo "PRIMER_PICK_LEFT_PRIMER=1"
							echo "PRIMER_PICK_INTERNAL_OLIGO=1"
							echo "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1"
							echo "PRIMER_PICK_RIGHT_PRIMER=1"
							echo "PRIMER_INTERNAL_MIN_SIZE=18"
							echo "PRIMER_INTERNAL_MAX_SIZE=30"
							echo "PRIMER_INTERNAL_OPT_SIZE=$kmer"
							echo "PRIMER_MIN_SIZE=$primerminsize"
							echo "PRIMER_OPT_SIZE=$primeroptsize"
							echo "PRIMER_MAX_SIZE=$primermaxsize"
							echo "PRIMER_MAX_TM=$primermaxtm"
							echo "PRIMER_MIN_TM=$primermintm"
							echo "PRIMER_INTERNAL_MIN_TM=0"
							echo "PRIMER_INTERNAL_MAX_TM=100"
							echo "PRIMER_MIN_GC=$primermingc"
							echo "PRIMER_MAX_GC=$primermaxgc"
							#echo "PRIMER_NUM_RETURN=$primercount"
							echo "PRIMER_MAX_NS_ACCEPTED=0"
							echo "SEQUENCE_TARGET=$target"
							echo "PRIMER_LOWERCASE_MASKING=1"
							echo "PRIMER_PRODUCT_SIZE_RANGE=90-300"
							echo "P3_FILE_FLAG=0"
							echo "PRIMER_EXPLAIN_FLAG=1"
							echo "="
						} >> p3_input/"$i".input.p3
						#awk '{ print "primer3_core < "'$P3_INPUT'"/"$i".input.p3 > "'$P3_OUTPUT'"/"$i".output.p3" }'
				done < allgenes_"$i"

		done

		[[ "$(ls -A p3_input)" ]] || error_exit "Failed to generate Primer3 input file(s)!"
		rm allgenes_subset*
		cd $OUTPUT
}



parse_p3_output_pangenome() {

	cd $PRIMERS/p3_output
	logger "Parsing pangenome primer design results"

	for i in *output.p3; do \
	awk -F'[=,]' '/^SEQUENCE_ID/ { guide=$2; set=1 } \

	/^PRIMER_LEFT_._SEQUENCE/ { printf guide "\t" "pset"set "\t" $2 } \
	/^PRIMER_RIGHT_._SEQUENCE/ { printf "\t" $2 } \
	/^PRIMER_INTERNAL_._SEQUENCE/ { printf "\t" $2 } \

	/^PRIMER_LEFT_.\=/ { printf "\t" $3 } \
	/^PRIMER_RIGHT_.\=/ { printf "\t" $3 }  \
	/^PRIMER_INTERNAL_.\=/ { printf "\t" $3 } \

	/^PRIMER_LEFT_._GC_PERCENT/ { printf "\t" $2 } \
	/^PRIMER_RIGHT_._GC_PERCENT/ { printf "\t" $2 } \
	/^PRIMER_INTERNAL_._GC_PERCENT/ { printf "\t" $2 } \

	/^PRIMER_PAIR_._PRODUCT_SIZE/ { printf "\t" $2 "\n"; set++ }' $i >> allgenes.p3.out; done

	[[ -s allgenes.p3.out ]] || { error_exit "Failed to parse Primer3 results!"; return; }

	rm subset*output.p3
	cut -f1 allgenes.p3.out > genes.txt
	awk -F "--" '{ print $1 }' genes.txt > genes_only.txt

		# Revseq guide sequence 
		local guideseq=$(cut -f5 allgenes.p3.out | \
			nl | \
			sed 's/^ */>seq/g' | \
			sed 's/\t/\n/g' | \
			seqkit seq -r -t DNA -v -p | \
			grep -v ">")

		## Create array from prev file; extract pam
		local pam=$(awk -F "\t" 'FNR==NR{pam_lookup[$4]=$3; next} {print pam_lookup[$1]}' $PANGENOME/guides_precheck/allgenes_lookup_pass_max5.tsv genes.txt)

		## Combine
		paste allgenes.p3.out <(echo -e "$guideseq") <(echo -e "$pam") genes_only.txt --delimiters "\t" > allgenes_p3_parsed.tsv
		rm genes_only.txt

		# Keep revcomp guide sequence for reverse strand guides
		awk -F "\t" '{ if ($1 ~ /r$/) print $1,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; \
			else print $1,$5,$14,$15,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }' OFS="\t" allgenes_p3_parsed.tsv \
		> allgenes_primers_output.tsv


		## Add header
		local header=$(printf "%s\t" \
			"guide_set" \
			"guide_seq" \
			"pam" "gene" \
			"primer_set" \
			"fwd_primer" \
			"rev_primer" \
			"probe" \
			"fwd_primer_length" \
			"rev_primer_length" \
			"guide_length" \
			"fwd_primer_gc" \
			"rev_primer_gc" \
			"guide_gc" \
			"product_size" | sed 's/\t$//g')

		sed -i "1s/^/$header\n/" allgenes_primers_output.tsv

		## Summarize multiple guides per primer
		cat allgenes_primers_output.tsv | \
		csvtk summary -t -f guide_set:uniq,guide_set:countunique \
		-g primer_set,fwd_primer,rev_primer \
		> allgenes_primers_output_grouped.tsv

		## Write guides with primers designed to file (this file is not used anywhere else )
		cut -f2 allgenes_primers_output.tsv | grep -v guide | sort | uniq > $PRIMERS/allgenes_primers.tsv

		rm genes.txt #allgenes_p3_parsed.tsv revseq_input.txt 

		cd $OUTPUT

}



guide_sensitivity() {
	# Get guides that passed Primer3
	# bbmap.sh against target genomes
	# This will tell if guide is present in genome, but not necessarily inside amplicon region

	cd $PRIMERS
	mkdir -p $PRIMERS/guide_prevalence
	logger "Calculating prevalence of guides"

	if [[ "$METHOD" == "pangenome" ]]; then
		#All possible PAMs
		if [[ "$CAS" == "cas12" ]]; then
			cut -f1-2 p3_output/allgenes_primers_output.tsv | sed '1d' | sort | uniq |
			awk '{ print $0, "TTTA"; print $0, "TTTC"; print $0, "TTTG"; print $0, "TTTT"; }' |
			awk '{ print ">"$1, "\n"$3$2 }' > guides_input.fasta
		else
			cut -f1-2 p3_output/allgenes_primers_output.tsv | sed '1d' | sort | uniq |
			awk '{ print ">"$1, "\n"$3$2 }' > guides_input.fasta
		fi

	elif [[ "$METHOD" == "kmer" ]]; then
		#All possible PAMs
		if [[ "$CAS" == "cas12" ]]; then
			cut -f1,3 OFS="\t" p3_output/allguides_primers_output.tsv | sed '1d' | sort | uniq | 
			awk '{ print $0, "TTTA"; print $0, "TTTC"; print $0, "TTTG"; print $0, "TTTT"; }' | 
			awk '{ print ">"$1, "\n"$3$2 }' > guides_input.fasta
		else
			cut -f1,3 OFS="\t" p3_output/allguides_primers_output.tsv | sed '1d' | sort | uniq |
			awk '{ print ">"$1, "\n"$3$2 }' > guides_input.fasta
		fi
	fi
	
	# Map to target genomes (only canonical PAMs)
	cd guide_prevalence
	split -l 500 $OUTPUT/genomes_target.txt subset. -d -a 3

	for file in subset.*; do
		cat $file | \
		parallel -j "$MAX_PARALLEL_JOBS" \
		'bbmap.sh in='"$PRIMERS"'/guides_input.fasta \
		ref='"$GENOMES_TARGET"'/"{}"_genomic.fna \
		nodisk \
		noheader=t \
		ambig=all \
		vslow \
		perfectmode \
		maxsites=1000000 \
		threads=8 \
		outm="{}"_wg.sam'; \
		for i in *_wg.sam; do awk -F "\t" '{ print FILENAME, $0 }' OFS="\t" $i | sed 's/_wg.sam//1' | cut -f1-2 >> allguides_hits_genome.tsv; done
		rm *_wg.sam "$file"
	done
	#rm subset.*

	# cd $PRIMERS/guide_prevalence
	# for i in *_wg.sam; do awk -F "\t" '{ print FILENAME, $0 }' OFS="\t" $i | sed 's/_wg.sam//1' | cut -f1-2 >> allguides_hits_genome.tsv; done
	sed -i '1i genome\tguide' allguides_hits_genome.tsv
	#rm *_wg.sam

	# Summarize by guide 
	cat allguides_hits_genome.tsv | csvtk summary -t -f genome:count,genome:countunique -g guide | sed '1d' | awk '{ print $1, $3, $2/$3 }' OFS="\t" > $PRIMERS/allguides_target_summary.tsv
	sed -i '1i guide\tunique_genome\tguide_average_copy_number_target' $PRIMERS/allguides_target_summary.tsv
	
	# Join results
	#cd $PRIMERS
	#csvtk join --left-join --out-tabs --tabs -k allguides_target_summary_genome.tsv allguides_target_summary_amplicon.tsv > allguides_target_summary.tsv

	[[ -s "$PRIMERS"/allguides_target_summary.tsv ]] || error_exit "Failed to generate "allguides_target_summary.tsv""
	rm allguides_hits_genome.tsv
	cd $OUTPUT
}



process_cigar() {

	local input_file="$1"

	# Process SAM file line by line, excluding header lines
	while IFS=$'\t' read -r -a fields; do
		# if [[ ${fields[0]} == @* ]]; then
		# 	continue
		# fi

		# Get the relevant columns
		read_name="${fields[0]}"
		cigar="${fields[5]}"

		# Check if the mapping is on the reverse strand
		is_reversed=false
		if (( ${fields[1]} & 16 )); then
			is_reversed=true
		fi

		# Break down the CIGAR string into individual components
		components=()
		current_component=""
		for (( i=0; i<${#cigar}; i++ )); do
			char="${cigar:$i:1}"
			if [[ $char =~ [0-9] ]]; then
				current_component+=$char
			else
				components+=("$current_component$char")
				current_component=""
			fi
		done

		 # Generate the expanded CIGAR string
		 expanded_string=""
		 num_X=0
		 for component in "${components[@]}"; do
		 	count="${component%[=X]}"
		 	op="${component: -1}"
		 	expanded_string+=$(printf "$op%.0s" $(seq 1 "$count"))
		 	if [[ $op == "X" ]]; then
		 		num_X=$((num_X + count))
		 	fi
		 done

	# Reverse the expanded string if the mapping is on the reverse strand
	if $is_reversed; then
		expanded_string=$(echo "$expanded_string" | rev)
	fi

	# Output the read name and the expanded CIGAR string
	echo -e "${read_name}\t${expanded_string}\t${num_X}"
	
	done < "$input_file"

}



guide_specificity() {
	# Get guides that passed Primer3
	# compare kmers against offtarget genomes (allow 5 mismatches)
	
	cd $PRIMERS
	export -f process_cigar

	mkdir -p $PRIMERS/guide_specificity
	cd guide_specificity

	logger "Checking specificity of guides"

	# Map to offtarget genomes
	local minid=$(echo "scale=3; 1 - 5/$kmerpam" | bc)

	split -l 500 $OUTPUT/genomes_offtarget.txt subset. -d -a 3

	for file in subset.*; do
		cat $file | \
		parallel -j "$MAX_JOBS" \
		'bbmap.sh in='"$PRIMERS"'/guides_input.fasta \
		ref='"$GENOMES_OFFTARGET"'/"{}"_genomic.fna \
		nodisk \
		noheader=t \
		ambig=all \
		vslow \
		idfilter='"$minid"' \
		nmtag=t \
		xmtag=f \
		nhtag=f \
		amtag=f \
		idtag=t \
		indelfilter=0 \
		maxsites=10000000 \
		threads=4 \
		outm="{}"_wg.sam;
		process_cigar "{}"_wg.sam | awk "{ print \"{}\", \$0 }" OFS="\t" > "{}".parsed; rm "{}"_wg.sam'
		find . -name "*.parsed" -type f -print0 | xargs -0 -n100 cat >> allguides_hits_offtarget_genomes.tsv
		find . -name "*.parsed" | xargs rm
		rm $file
	done

	##################################################################################
	[[ -s "$PRIMERS/guide_specificity/allguides_hits_offtarget_genomes.tsv" ]] || { logger "No guides mapped to off-target genomes"; return; }

	## Filter likelihood of off-target activity based on PAM and seed region
	awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $4, substr($3,1,4), substr($3,5,7), substr($3,12)}' allguides_hits_offtarget_genomes.tsv > allguides_hits_offtarget_genomes_split.tsv
	
	awk 'BEGIN{FS=OFS="\t"} { mm_pam = gsub(/X/, "X", $4); \
	mm_seed = gsub(/X/, "X", $5); \
	mm_distal = gsub(/X/, "X", $6); \
	print $0, mm_pam, mm_seed, mm_distal }' allguides_hits_offtarget_genomes_split.tsv | \
	awk -F "\t" '$3<=5' | \
	awk -F "\t" '{ if ($7==0 && $8==0) print $0, "offtarget"; else print $0, "ok" }' OFS="\t" > allguides_hits_offtarget_genomes_split.tsv2
	
	mv allguides_hits_offtarget_genomes_split.tsv2 allguides_hits_offtarget_genomes_split.tsv

	sed -i '1i genome\tguide\tmismatch_count\tcigar_pam\tcigar_seed\tcigar_distal\tmm_pam\tmm_seed\tmm_distal\tremarks' allguides_hits_offtarget_genomes_split.tsv

	awk -F "\t" '$10!="ok"' allguides_hits_offtarget_genomes_split.tsv | 
	csvtk summary -t -f genome:count,genome:countunique -g guide | sed '1d' | awk '{ print $1, $3, $2/$3 }' OFS="\t" > $PRIMERS/guide_offtarget_summary.tsv
	sed -i '1i guide\tunique_genome\tguide_avg_copy_number_offtarget' $PRIMERS/guide_offtarget_summary.tsv

	[[ -s "$PRIMERS"/guide_offtarget_summary.tsv ]] || error_exit "Failed to generate "guide_offtarget_summary.tsv""
	rm allguides_hits_offtarget_genomes.tsv

	cd $OUTPUT

}



get_kmers_offtarget() {

	logger "Counting k="$kmer" in non-target genomes"
	mkdir -p $KMER_OFFTARGET/combined_batches $KMER_OFFTARGET/output
	
	cd $KMER_OFFTARGET
	
	#Process files in batches of 20
	counter=0
	cat $OUTPUT/genomes_offtarget.txt | xargs -n 20 | \
	while read -r genome
	do
		counter=$((counter + 1))
		echo "$genome" | \
		tr " " "\n" | \
		parallel -j "$MAX_JOBS" 'glistmaker '"$GENOMES_OFFTARGET"'/"{}"_genomic.fna \
		--num_threads '"$MAX_THREADS_PER_JOB"' \
		-w '"$kmer"' \
		-o "{}"' 

		glistcompare *.list --union -o combined_batches/combined_samples_"$counter"
		[[ -s combined_batches/combined_samples_"$counter"_"$kmer"_union.list ]] && rm "$KMER_OFFTARGET"/*_"$kmer".list

	done

	logger "Running union of k="$kmer" in non-target genomes"
	if [[ $(wc -l <$OUTPUT/genomes_offtarget.txt) -gt 200 ]]
	then
		glistcompare combined_batches/combined_samples_*_"$kmer"_union.list --union -o output/combined_offtarget
	else
		mv combined_batches/combined_samples_1_"$kmer"_union.list output/combined_offtarget_"$kmer"_union.list
	fi
	[[ -s $KMER_OFFTARGET/output/combined_offtarget_"$kmer"_union.list ]] && rm -rf combined_batches
	
	#if [[ $METHOD == "kmer" ]]
	#then
		mv output/combined* $KMER_OFFTARGET
		rm -rf output
	#fi

	cd $OUTPUT
	   	
}



get_kmers_offtarget_pam() {
	# For only pangenome module, also count kmers of length kmerpam
	
	logger "Counting k="$kmerpam" in non-target genomes"
	mkdir -p $KMER_OFFTARGET/combined_batches

	cd $KMER_OFFTARGET

	counter=0
	cat $OUTPUT/genomes_offtarget.txt | xargs -n 20 | \
	while read -r genome
	do
		counter=$((counter + 1))
		echo "$genome" | \
		tr " " "\n" | \
		parallel -j "$MAX_JOBS" 'glistmaker '"$GENOMES_OFFTARGET"'/"{}"_genomic.fna \
		--num_threads '"$MAX_THREADS_PER_JOB"' \
		-w '"$kmerpam"' \
		-o "{}"' 

		## Remove list so combined_offtarget_$kmer is not combined
		glistcompare *_"$kmerpam".list --union -o combined_batches/combined_samples_"$counter"
		[[ -s combined_batches/combined_samples_"$counter"_"$kmerpam"_union.list ]] && rm "$KMER_OFFTARGET"/*"$kmerpam".list

	done

	logger "Running union of k="$kmerpam" in non-target genomes"
	
	if [[ $(wc -l <$OUTPUT/genomes_offtarget.txt) -gt 200 ]]
	then
		glistcompare combined_batches/combined_samples_*_"$kmerpam"_union.list --union -o output/combined_offtarget
	else
		mv combined_batches/combined_samples_1_"$kmerpam"_union.list output/combined_offtarget_"$kmerpam"_union.list
	fi

	[[ -s $KMER_OFFTARGET/output/combined_offtarget_"$kmerpam"_union.list ]] && rm -rf combined_batches
	
	mv output/*.list $KMER_OFFTARGET
	cd $OUTPUT

}



get_kmers_target() {
	
	mkdir -p $KMER_TARGET/combined_batches -p $KMER_TARGET/output

	if [[ "$subsample" == "no" ]]
	then
		local cutoff=$(wc -l $OUTPUT/genomes_target.txt | awk '{ print $1*0.5 }' | awk '{ print int($1+0.5)}')
		local targetcount=$(wc -l $OUTPUT/genomes_target.txt | awk '{ print $1 }')
		local myfile=$OUTPUT/genomes_target.txt
	else
		local cutoff=$(wc -l $OUTPUT/genomes_target_subsample.txt | awk '{ print $1*0.5 }' | awk '{ print int($1+0.5)}')
		local targetcount=$(wc -l $OUTPUT/genomes_target_subsample.txt | awk '{ print $1 }')
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	logger "Counting k="$kmer" in target genomes"
	
	cat $myfile |
	parallel -j "$MAX_PARALLEL_JOBS" \
	'glistmaker '"$GENOMES_TARGET"'/"{}"_genomic.fna \
	--num_threads '"$MAX_THREADS_PER_JOB"' \
	-w '"$kmer"' \
	-o '"$KMER_TARGET"'/"{}"'

	cd $KMER_TARGET

	if [[ $(wc -l <$myfile) -gt 200 ]]
	then
		
		# ls *.list > samples
		# split -l 100 samples samples. -d

		logger "Running union of k="$kmer" in target genomes"
		# cat $myfile | awk '{ print $0"_"'$kmer'".list" }' | xargs -n 200 \
		# parallel -j "$MAX_PARALLEL_JOBS" 'glistcompare "{}" \
		# --union -o combined_batches/combined_samples_"{#}"'

		awk '{ print $0"_"'$kmer'".list" }' $myfile | xargs -L 20 | \
		awk '{ print NR, $0 }' OFS="\t" | \
		awk -F "\t" '{ print "glistcompare "$2" --union -o combined_batches/combined_samples_"$1"" }' > run.sh

		cat run.sh | parallel -j "$MAX_PARALLEL_JOBS"
		rm run.sh
		# ls samples.* |
		# parallel -j "$MAX_JOBS" \
		# 'glistcompare `cat "{}" | tr "\n" " "` \
		# --union -o combined_"{}"'

		glistcompare combined_batches/combined_samples_*_"$kmer"_union.list --cutoff $cutoff --union -o combined_target
		[[ -s combined_target_"$kmer"_union.list ]] && rm -rf combined_batches
		#rm combined_samples*list samples*
		
	else
		
		logger "Running union of k="$kmer" in target genomes"
		glistcompare *.list --cutoff $cutoff --union -o combined_target
		rm -rf combined_batches
		
	fi

	cd $OUTPUT

	## Remove kmers present in off-target and human genomes; reduce kmer+pam query
	## Filter 50% (raw count) of target kmers
	## Query k22; filter 80% prevalence
	## Add PAM sequence to filtered kmer list

	logger "Comparing k="$kmer" in target and non-target+human genomes"

	## First combine off-target and human genome kmers
	glistcompare $DATABASE/GCF_000001405.40_GRCh38.p14_genomic_"$kmer".list \
	$KMER_OFFTARGET/combined_offtarget_"$kmer"_union.list \
	--union -o $KMER_OFFTARGET/combined_allnontarget

	glistcompare $KMER_TARGET/combined_target_"$kmer"_union.list \
	$KMER_OFFTARGET/combined_allnontarget_"$kmer"_union.list \
	--difference -o $KMER_TARGET/combined_target_uniq

	## First pass filtering
	logger "Filtering target only k="$kmer" by first 50% aggregate cutoff and homopolymer tracts"
	glistquery $KMER_TARGET/combined_target_uniq*diff*list |
	egrep -v "NUnique|NTotal" |
	awk '$2>="$cutoff"' |
		## remove homopolymers (5 consecutive bases)
		awk '$1 !~ /AAAAA|CCCCC|GGGGG|TTTTT/ {print}' |
		cut -f1 > $KMER_TARGET/combined_target_kmers.txt

	## Query combined target kmers against all target genomes
	logger "Querying target-unique k="$kmer" against target genomes to obtain counts"

	split -l 200 $myfile subset. -d -a 3

	for file in subset.*; do
		cat $file |
		parallel -j 10 \
		'glistquery '"$KMER_TARGET"'/"{}"_'"$kmer"'.list \
		-f '"$KMER_TARGET"'/combined_target_kmers.txt | \
		cut -f2 > '"$KMER_TARGET"'/"{}".count; \
		sed -i "1i$(echo '{}')" '"$KMER_TARGET"'/"{}".count; \
		rm '"$KMER_TARGET"'/"{}"_'"$kmer"'.list'; \
		paste -d'\t' "$KMER_TARGET"/*.count >> "$KMER_TARGET"/allcounts_"$file".txt; \
		rm "$KMER_TARGET"/*.count; done

	rm subset.*
	## Rscript to combine and perform second-pass filtering by prevalence across all target genomes
	cd $KMER_TARGET 

	## Optimized merging of count files
	#header=$(ls *.count | sed 's/.count//g' | tr "\n" "\t")

	logger "Second-round filtering of user-specified "$threshold"% prevalence across target genomes"
	printf "%s\t" \
	"kmer" \
	"count_genomes_present" \
	"sum_across_genomes" \
	"avg_copy_number" \
	"perc_total" |
	awk '{ print $0 }' |
	sed 's/\t$//g' > header.txt

	## Paste count files 200 at a time
	# ls -1 *.count | split -l 200 -d - countfiles
	# for i in countfiles*; do paste $(cat $i) > merge${i##countfiles}; done
	# paste merge* > allcounts.txt

	paste allcounts_subset* > allcounts.txt
	sed  '1ikmer' combined_target_kmers.txt > kmers.txt
	paste kmers.txt allcounts.txt > allcounts_matrix.tsv
	pigz -9 -p "$MAX_JOBS" allcounts_matrix.tsv 
	rm kmers.txt allcounts_subset*

	## Sum non-zero columns here, only store aggregate info
	sed -i '1d' allcounts.txt
	awk -F "\t" '{nz=0; sum=0; for(i=1;i<=NF;i++) {nz+=($i!=0); sum+=$i} print nz, sum, sum/nz, nz/"'$targetcount'"*100 }' OFS="\t" allcounts.txt > allcounts_sum.txt
	paste combined_target_kmers.txt allcounts_sum.txt > combined_counts.txt
	#awk -F "\t" '{ print $1, $(NF-3), $(NF-2), $(NF-1), $(NF) }' OFS="\t" > combined_counts.txt

	## Add headers
	cat header.txt combined_counts.txt > combined_counts_final.txt

	## Since pre-processing done, Rscript needs to sort
	logger "Processing k="$kmer" count data in R"

	if [[ "$CAS" == "cas12" ]]
	then
		filter_kmers.R $threshold $KMER_TARGET/combined_counts_final.txt $KMER_TARGET yes cas12
	else
		filter_kmers.R $threshold $KMER_TARGET/combined_counts_final.txt $KMER_TARGET NA cas13
	fi	

	# Remove count files to save space only if previous steps successful
	if [[ (-s $KMER_TARGET/query_withpam_lookup.tsv && -s $KMER_TARGET/query_withpam.txt) || (-s $KMER_TARGET/query_nopam_lookup.tsv && -s $KMER_TARGET/query_nopam.txt) ]] 
	then
		rm combined_counts.txt header.txt allcounts.txt allcounts_sum.txt merge* countfiles* samples.*
		rm $KMER_TARGET/combined_counts_final.txt
		#find $KMER_TARGET -type f -name "*.count" | xargs -r rm
	else
		error_exit "No target-unique k="$kmer" found at specified threshold! Exiting.."
	fi

	cd $OUTPUT

}


get_kmers_target_pam() {
	#kmerpam=$(grep -v "#" $OUTPUT/$IN | grep kmer | awk -F "|" '{ print $2+4}')

	logger "Counting k="$kmerpam" in target genomes"

	mkdir -p $KMER_TARGET_PAM/combined_batches

	if [[ "$subsample" == "no" ]]
	then
		local myfile=$OUTPUT/genomes_target.txt
	else
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	## Get kmers for guide length AND guide length+4 (for PAM length)
	## Get union for target_pam, to reduce glistquery lines for later step
	cat $myfile |
	parallel -j "$MAX_PARALLEL_JOBS" \
	'glistmaker '"$GENOMES_TARGET"'/"{}"_genomic.fna \
	--num_threads '"$MAX_THREADS_PER_JOB"' \
	-w '"$kmerpam"' \
	-o '"$KMER_TARGET_PAM"'/"{}"'

	cd $KMER_TARGET_PAM

	if [[ $(wc -l <$myfile) -gt 200 ]]
	then
		
		# ls *.list > samples
		# split -l 100 samples samples. -d

		logger "Running union of k="$kmerpam" in target genomes"
		# cat $myfile | xargs -n 200
		# parallel -j "$MAX_PARALLEL_JOBS" 'glistcompare "{}"_'"$kmerpam"'.list \
		# --union -o combined_batches/combined_samples_"{#}"'

		awk '{ print $0"_"'$kmerpam'".list" }' $myfile | xargs -L 200 | \
		awk '{ print NR, $0 }' OFS="\t" | \
		awk -F "\t" '{ print "glistcompare "$2" --union -o combined_batches/combined_samples_"$1"" }' > run.sh
		
		cat run.sh | parallel -j "$MAX_PARALLEL_JOBS"
		rm run.sh

		glistcompare combined_batches/combined_samples_*_"$kmerpam"_union.list --union -o combined_target_pam
		[[ -s combined_target_"$kmerpam"_union.list ]] && rm -rf combined_batches

		#rm combined_samples*list samples*
		
	else
		
		logger "Running union of k="$kmerpam" in target genomes"
		glistcompare *.list --union -o combined_target_pam
		rm -rf combined_batches
		
	fi
	cd $OUTPUT

}



get_uniq_target_kmer_withpam() {
	## Query target+pam list against k=26
	## Filter to k=26 where kmer is present in ≥80% target genomes
	## Get substring of k=22
	## Query against non-target genomes and exclude those present in non-target genomes
	
	## Exclude queries not present in combined_target_pam union list file (using glistquery method)
	## This includes filtering for minimum count using --minfreq option, reducing computational time (use 50% as default here)
	cd $OUTPUT

	if [[ "$subsample" == "no" ]]
	then
		local targetcount=$(wc -l <$OUTPUT/genomes_target.txt)
		local myfile=$OUTPUT/genomes_target.txt
	else
		local targetcount=$(wc -l <$OUTPUT/genomes_target_subsample.txt)
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	glistquery $KMER_TARGET_PAM/combined_target_pam_"$kmerpam"_union.list \
	--minfreq $(wc -l $myfile | awk '{ print $1*0.5 }' | awk '{ print int($1+0.5) }') \
	--queryfile $KMER_TARGET/query_withpam.txt | \
	cut -f1 \
	> $KMER_TARGET/query_withpam_filtered.txt


	## Query target genome k=26 kmers using shortlisted k=22 kmers with PAM
	logger "Query k="$kmerpam" in target genomes with PAM-added target-unique set to get counts"

	split -l 200 $myfile subset. -d -a 3

	for file in subset.*; do \
		cat $file |
		parallel -j "$MAX_JOBS" \
		'glistquery '"$KMER_TARGET_PAM"'/"{}"_'"$kmerpam"'.list \
		--queryfile '"$KMER_TARGET"'/query_withpam_filtered.txt | \
		cut -f2 > '"$KMER_TARGET_PAM"'/"{}".count; \
		sed -i "1i$(echo '{}')" '"$KMER_TARGET_PAM"'/"{}".count; \
		rm '"$KMER_TARGET_PAM"'/"{}"_'"$kmerpam"'.list'; \
		paste -d'\t' "$KMER_TARGET_PAM"/*.count >> "$KMER_TARGET_PAM"/allcounts_"$file".txt
		rm "$KMER_TARGET_PAM"/*.count; done

	rm subset.*

	## Rscript to calculate prevalence of kmer+pam across target genomes
	cd $KMER_TARGET_PAM

	## Optimized merging of count files
	printf "%s\t" \
	"kmer" \
	"count_genomes_present" \
	"sum_across_genomes" \
	"avg_copy_number" \
	"perc_total" |
	awk '{ print $0 }' |
	sed 's/\t$//g' > header.txt

	## Paste count files serially 200 at a time
	# ls -1 *.count | split -l 200 -d - countfiles
	# for i in countfiles*; do paste $(cat $i) > merge${i##countfiles}; done
	# paste merge* > allcounts.txt
	
	paste allcounts_subset* > allcounts.txt
	sed  '1ikmer' $KMER_TARGET/query_withpam_filtered.txt > kmers.txt
	paste kmers.txt allcounts.txt > allcounts_matrix.tsv
	pigz -9 -p "$MAX_JOBS" allcounts_matrix.tsv 
	rm kmers.txt allcounts_subset*

	## Sum non-zero columns here, only store aggregate info
	sed -i '1d' allcounts.txt
	awk -F "\t" '{nz=0; sum=0; for(i=1;i<=NF;i++) {nz+=($i!=0); sum+=$i} print nz, sum, sum/nz, nz/"'$targetcount'"*100 }' OFS="\t" allcounts.txt > allcounts_sum.txt
	paste $KMER_TARGET/query_withpam_filtered.txt allcounts_sum.txt > combined_counts.txt
	

	## Add headers
	cat header.txt combined_counts.txt > combined_counts_final.txt

	logger "Parsing k="$kmerpam" results in R"

	filter_kmers.R "$threshold" $KMER_TARGET_PAM/combined_counts_final.txt $KMER_TARGET_PAM no $CAS
	
	cd $OUTPUT

	## Remove target+pam count files here as no longer needed if previous steps successful
	if [[ -s $KMER_TARGET_PAM/kmers_withpam_offtarget_query.txt && -s $KMER_TARGET_PAM/kmers_withpam_lookup.tsv ]]
	then
		cd $KMER_TARGET_PAM
		#find $KMER_TARGET_PAM -type f -name "*.count" | xargs -r rm
		rm merge* countfiles* combined_counts.txt header.txt allcounts.txt allcounts_sum.txt combined_counts_final.txt
		cd $OUTPUT
	else
		logger "Parsing of target kmer+pam count files failed! Exiting.."
		error_exit "Parsing of target kmer+pam count files failed! Exiting.."
	fi

	mkdir -p $KMER_FINAL

	## Search output from Rscript above against union of off-target genomes, keep those with zero hits (allowing 2 mismatches)
	if [[ $mismatch == "yes" ]]
	then
		logger "Running second specificity filtering with 2 mismatches in non-target genomes"
		glistquery $KMER_OFFTARGET/combined_offtarget_*_union.list \
		--mismatch 2 \
		-f $KMER_TARGET_PAM/kmers_withpam_offtarget_query.txt |
		awk '$2==0' |
		egrep -v "NUnique|NTotal" |
		cut -f1 |
		# must contain all bases
		awk '/A/ && /C/ && /G/ && /T/' > $KMER_FINAL/uniq_target_kmers_final.txt
	elif [[ $mismatch == "no" ]]
		then
			logger "Running second specificity filtering without mismatch in non-target genomes"
			glistquery $KMER_OFFTARGET/combined_offtarget_*_union.list \
			-f $KMER_TARGET_PAM/kmers_withpam_offtarget_query.txt |
			awk '$2==0' |
			egrep -v "NUnique|NTotal" |
			cut -f1 |
		# must contain all bases
		awk '/A/ && /C/ && /G/ && /T/' > $KMER_FINAL/uniq_target_kmers_final.txt
	fi


	## If kmers unique to target genomes present, then generate table
	## Obtain prevalence of guide+PAM across target genomes
	if [[ -s $KMER_FINAL/uniq_target_kmers_final.txt ]]
	then
		##Merge kmer and pam field 
		logger "Generating target-unique k-mer results table"
		grep -wF -f $KMER_FINAL/uniq_target_kmers_final.txt $KMER_TARGET_PAM/kmers_withpam_lookup.tsv |\
		awk '$NF>="$threshold"' |\
		nl |\
		sed 's/^ *//g' |\
		sed 's/^/g/g' > $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv

		## Add header
		echo -e "guide\tactual_pam\tguide_sequence_to_order\tcount_genomes_present\tsum_across_genomes\tavg_copy_number\tperc_total\n$(cat $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv)" \
		> $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv 
	else
		logger "No target-unique kmers found! Exiting.."	
		error_exit "No k-mers unique to target genomes found!"
	fi

	rm $KMER_FINAL/uniq_target_kmers_final.txt $KMER_TARGET_PAM/kmers_withpam_lookup.tsv $KMER_TARGET_PAM/kmers_withpam_offtarget_query.txt
	rm $KMER_TARGET/query_kmer80pct_filter2_withpam_filtered.tsv $KMER_TARGET/query_withpam_lookup.tsv $KMER_TARGET/query_withpam.txt

	cd $OUTPUT

}


get_uniq_target_kmer_nopam() {
	## Query target+pam list against k=26
	## Filter to k=26 where kmer is present in ≥80% target genomes
	## Get substring of k=22
	## Query against non-target genomes and exclude those present in non-target genomes
	
	## Exclude queries not present in combined_target_pam union list file (using glistquery method)
	## This includes filtering for minimum count using --minfreq option, reducing computational time (use 50% as default here)
	cd $OUTPUT

	mkdir -p $KMER_FINAL

	## Search output from Rscript above against union of off-target genomes, keep those with zero hits (allowing 2 mismatches)
	if [[ $mismatch == "yes" ]]
	then
		logger "Running second specificity filtering with 2 mismatches in non-target genomes"
		glistquery $KMER_OFFTARGET/combined_offtarget_*_union.list \
		--mismatch 2 \
		-f $KMER_TARGET/query_nopam.txt |
		awk '$2==0' |
		egrep -v "NUnique|NTotal" |
		cut -f1 |
		# must contain all bases
		awk '/A/ && /C/ && /G/ && /T/' > $KMER_FINAL/uniq_target_kmers_final.txt
	elif [[ $mismatch == "no" ]]
	then
		logger "Running second specificity filtering without mismatch in non-target genomes"
		glistquery $KMER_OFFTARGET/combined_offtarget_*_union.list \
		-f $KMER_TARGET/query_nopam.txt |
		awk '$2==0' |
		egrep -v "NUnique|NTotal" |
		cut -f1 |
		# must contain all bases
		awk '/A/ && /C/ && /G/ && /T/' > $KMER_FINAL/uniq_target_kmers_final.txt
	fi

	grep -wF -f $KMER_FINAL/uniq_target_kmers_final.txt $KMER_TARGET/query_nopam_lookup.tsv |
		awk '$NF>="$threshold"' |
		shuf -n 5000 |
		nl |
		sed 's/^ *//g' |
		sed 's/^/g/g' |
		awk '{ print $1, "NA", $2, $3, $4, $5, $6 }' OFS="\t" > $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv

	## Add header
	echo -e "guide\tactual_pam\tguide_sequence_to_order\tcount_genomes_present\tsum_across_genomes\tavg_copy_number\tperc_total\n$(cat $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv)" \
	> $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv 

	[[ -s $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv ]] || (logger "No target-unique kmers found! Exiting//" && error_exit "No k-mers unique to target genomes found!")

	cd $OUTPUT
}



precheck_guides() (
	cd $OUTPUT
	## Input potential guide sequence to primer3, check if probe features are okay before extracting coordinates, clustering and aligning amplicons
	mkdir -p $KMER_FINAL/guides_precheck

	# Limit number of guides to 5000
	grep -v kmer $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv | cut -f1 > $KMER_FINAL/guides_precheck.txt

	logger "Preparing Primer3 input file for guides precheck"
	
	while read guide
	do
		if [[ $guide == "guide" ]]
		then
			continue
		fi

		oligo=$(grep -w $guide $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv | cut -f3)

		{
			echo "SEQUENCE_ID=$guide"
			echo "SEQUENCE_INTERNAL_OLIGO=$oligo"
			echo "PRIMER_TASK=check_primers"
			echo "PRIMER_PICK_INTERNAL_OLIGO=1"
			echo "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1"
			echo "PRIMER_INTERNAL_MIN_SIZE=18"
			echo "PRIMER_INTERNAL_MAX_SIZE=30"
			echo "PRIMER_INTERNAL_MIN_TM=40.0"
			echo "PRIMER_EXPLAIN_FLAG=1"
			echo "="
		} >> $KMER_FINAL/guides_precheck/allguides.input.p3

	done < $KMER_FINAL/guides_precheck.txt
	
	## Filter out guides with bad properties using Primer3
	logger "Filtering low-quality guides"
	primer3_core < $KMER_FINAL/guides_precheck/allguides.input.p3 > $KMER_FINAL/guides_precheck/allguides.output.p3
	
	sed 's/^=/\n/g' $KMER_FINAL/guides_precheck/allguides.output.p3 |\
	awk -v RS= '/\nPRIMER_INTERNAL_NUM_RETURNED=1/' |\
	grep SEQUENCE_ID |\
	awk -F "=" '{ print $2 }' >> $KMER_FINAL/guides_allpass.txt

	# Exit if no guides pass
	[[ -s $KMER_FINAL/guides_allpass.txt ]] || { error_exit "No guides found!"; return; }

	## Guides that didn't pass primer3 precheck
	grep -vwF -f $KMER_FINAL/guides_allpass.txt $KMER_FINAL/guides_precheck.txt > $KMER_FINAL/guides_precheck/guides_fail_precheck.txt

	grep -wF -f $KMER_FINAL/guides_allpass.txt $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv |\
		## sort by decreasing prevalence and decreasing copy number
		sort -k7,7rn -k6,6rn |\
		## prevalence at least 50%
		awk '$7>=50' |\
		cut -f1 > $KMER_FINAL/guides_precheck/guides_allpass_50pct.txt

	# Exit if no guides above threshold
	[[ -s $KMER_FINAL/guides_precheck/guides_allpass_50pct.txt ]] || { error_exit "No guides above threshold found!"; return; }

	cp $KMER_FINAL/guides_precheck/guides_allpass_50pct.txt $KMER_FINAL/guides_shortlisted.txt

	grep -wF -f $KMER_FINAL/guides_shortlisted.txt $OUTPUT/uniq_target_kmers_final_prevalence_shortlisted.tsv

	## Guides that pass precheck but <90% prevalence
	grep -vwF -f $KMER_FINAL/guides_precheck/guides_allpass_50pct.txt $KMER_FINAL/guides_allpass.txt > $KMER_FINAL/guides_precheck/guides_pass_precheck_lt50pct.txt
	
	logger "Generating fasta for shortlisted guides" 
	grep -wF -f $KMER_FINAL/guides_shortlisted.txt $KMER_FINAL/uniq_target_kmers_final_prevalence.tsv |
	cut -f1,3 |
	sed 's/^/>/g' |
	sed 's/\t/\n/g' > $KMER_FINAL/query_uniq_target_kmers_nohomopolymer.fasta

	rm $KMER_FINAL/guides_precheck.txt

	cd $OUTPUT

)



get_guide_coordinates() {
	mkdir -p $ALIGN_OUT

	if [[ "$subsample" == "no" ]]
	then
		local myfile=$OUTPUT/genomes_target.txt
	else
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	logger "Mapping guides to target genomes"
	cat $myfile |
	parallel -j "$MAX_PARALLEL_JOBS" \
	'bbmap.sh in='"$KMER_FINAL"'/query_uniq_target_kmers_nohomopolymer.fasta \
	ref='"$GENOMES_TARGET"'/"{}"_genomic.fna \
	nodisk \
	noheader=t \
	ambig=all \
	vslow \
	perfectmode \
	maxsites=10000000 \
	threads='"$MAX_THREADS_PER_JOB"' \
	outm='"$ALIGN_OUT"'/"{}".sam'

	## Extract coordinates based on strand
	mkdir -p $COORD_OUT

	## Primary alignments: flags 0 (forward)
	## Primary alignments: flag 16 (reverse)
	## Secondary alignments: flags 256 (forward)
	## Secondary alignments: flags 272 (reverse)

	cd $ALIGN_OUT

	## Make bed file from sam
	logger "Extracting guide coordinates"
	for i in $(cat $myfile)
	do
		cut -f1-4 $i.sam |
		awk -F "\t" '{ print $3, $4, $4+26-1, $1, "1", $2 }' OFS="\t" |
		awk -F "\t" '{ if ($6==0||$6==256) print $0, "+"; else print $0, "-" }' OFS="\t" |
		cut -f1-5,7 |
		awk '{ print $1, $(NF-4), $(NF-3), $(NF-2), $(NF-1), $(NF) }' OFS="\t" > $COORD_OUT/$i.bed; done

	find . -type f -name "*.sam" | xargs rm

	cd $OUTPUT

}



split_amplicons() {

	cd $OUTPUT
	mkdir $AMPLICON

	if [[ "$subsample" == "no" ]]
	then
		local myfile=$OUTPUT/genomes_target.txt
	else
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	## Run seqkit subseq directly on bed file to get flanking regions and extract amplicons
	logger "Extracting 150bp flanking regions of guides"
	cat $myfile |
	parallel -j "$MAX_JOBS" \
	'cat '"$GENOMES_TARGET"'/"{}"_genomic.fna |
	seqkit subseq --bed '"$COORD_OUT"'/"{}".bed --up-stream 150 --down-stream 150 -j '"$MAX_THREADS_PER_JOB"' \
	>> '"$ALIGN_OUT"'/"{}"_template.fasta; \
	sed -i "s/0 g/0__g/1" '"$ALIGN_OUT"'/"{}"_template.fasta'

	find "$ALIGN_OUT" -maxdepth 1 -type f -name "*template.fasta" -exec cat {} + > "$ALIGN_OUT"/combined_templates.fasta

	rm -rf coord_out

	## Depending on assembly method of ref genomes, there may be duplicates in header id, and
	## this will cause problems downstream when extracting using subseq, as subseq only takes 
	## leading non-space characters as sequence id
	## To fix, remove space in header (before guide name), then rename using seqkit

	#sed -i 's/0 g/0__g/1' $ALIGN_OUT/combined_templates.fasta
	seqkit rename -j "$MAX_JOBS" $ALIGN_OUT/combined_templates.fasta > $ALIGN_OUT/combined_templates_renamed.fasta
	mv $ALIGN_OUT/combined_templates_renamed.fasta $ALIGN_OUT/combined_templates.fasta
	

	## Split templates into individual guide fastas
	## Get headers, sort headers by column 2, split headers by guide, then use cdbfasta

	logger "Processing templates for alignment"
	
	if [[ -s $ALIGN_OUT/combined_templates.fasta ]]
	then
		## Split fasta before indexing
		## Keep individual template files to speed up primersearch
		#find $ALIGN_OUT -type f -name "*_template.fasta" -printf '%p\n' | xargs -r rm
		seqkit split2 -j "$MAX_THREADS_PER_JOB" --by-size 1000000 --out-dir $ALIGN_OUT -o sub $ALIGN_OUT/combined_templates.fasta
		rm $ALIGN_OUT/combined_templates.fasta

		ls $ALIGN_OUT/sub*fasta | \
		parallel -j "$MAX_JOBS" 'cdbfasta "{}"; \
		grep ">" "{}" | sed "s/>//g" > "{}".headers'
				
		## Split templates header list by guide, then extract fasta
		## Must run serially to not override file
		
		for i in $ALIGN_OUT/*.headers
		do
			#awk -F "__" '{print >> $2}' $i
			sed 's/__/\t/1' $i | awk 'BEGIN{FS=OFS="\t"} {gsub(/_/, "\t", $2)} 1' | \
			sort -k2,2 > "$i".sorted
				#awk -F "\t" '{ if (NF==3) print $1"__"$2"_"$3; else print $1"__"$2 }' > "$i".sorted
		done

	else
		error_exit "Combined templates fasta file not found!"
	fi

	cd $OUTPUT
}


process_individual_guides() {
	cd $AMPLICON

	if [[ "$subsample" == "no" ]]
	then
		local myfile=$OUTPUT/genomes_target.txt
	else
		local myfile=$OUTPUT/genomes_target_subsample.txt
	fi

	## Minimum number of sequences required to keep cluster (80% of target genome count)
	if [[ "$domain" == "bacteria" ]]; then
		local cluster_cutoff=$(wc -l $myfile | awk '{ print $1*0.8 }')
		local seq_id=0.9
	else
		local cluster_cutoff=$(wc -l $myfile | awk '{ print $1*0.3 }')
		local seq_id=0.8
	fi


	# Get list of guide names
	find "$ALIGN_OUT" -type f -name "*.sorted" -exec awk '{ print $2 }' {} + | sort | uniq > list

	split -l 1000 list list. -d 

	# Process in batches of 1000
	logger "Clustering and aligning amplicons"
	for i in list.*; do 
		for guide in `cat "$i"`; do grep -hEw "$guide" "$ALIGN_OUT"/*.sorted >> "$i".txt; done
		awk -F "\t" '{ print >> $2 }' "$i".txt

		#Replace tab
		find . -maxdepth 1 -type f -name "g*" -print0 | parallel -0 -j "$MAX_JOBS" \
			'awk -F "\t" '\''{ if (NF == 3) print $1 "__" $2 "_" $3; else print $1 "__" $2 }'\'' {} > {}.2 && mv {}.2 {}'

		#Extract sequences
		parallel -a "$i" -j "$MAX_JOBS" 'for subset in '"$ALIGN_OUT"'/sub*.cidx; do cat "{}" | \
		cdbyank "$subset" >> "{}".fasta; rm "{}"; done'

		#Cluster sequences
		#logger "Clustering guide amplicons with cd-hit"
		parallel -a "$i" -j "$MAX_JOBS" 'cd-hit-est -i "{}".fasta -T '"$MAX_THREADS_PER_JOB"' \
		-c '"$seq_id"' -n 9 -g 1 -d 1000 -M 32000 -o "{}".cdhit; \
		make_multi_seq2.pl "{}".fasta "{}".cdhit.clstr alignments '"$cluster_cutoff"' "{}"; \
		rm "{}".fasta; rm "{}".cdhit; rm "{}".cdhit.clstr'

		#Align sequences
		#logger "Aligning clustered amplicons"
		find . -name "*.fa" -printf '%f\n' | sed 's/.fa//g' | \
		parallel -j "$MAX_JOBS" 'mafft --auto --thread 8 \
		alignments/{}.fa > alignments/{}.aln.fa; rm alignments/{}.fa'

		#Get consensus alignments
		find . -name "*.aln.fa" -printf '%f\n' | sed 's/.aln.fa//g' | \
		parallel -j "$MAX_JOBS" 'cons -sequence alignments/{}.aln.fa -name {} \
		-outseq alignments/{}.cons; rm alignments/{}.aln.fa'

		# Process consensus alignment
		find alignments/ -type f -name "*.cons" | \
 		parallel -j "$MAX_JOBS" 'seqkit -is replace -j 8 -p "^n+|n+$" -r "" {} -o {.}.2 && mv {.}.2 {}'

		find alignments -type f -name "g*_cl*.cons" -exec cat {} + | \
		awk '{if(NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' | \
		tr "\n>" "\t\n>" | \
		grep -v "^$" | \
		awk '{ print $1, $1, $2 }' OFS="\t" |  \
		sed 's/_/\t/2' | \
		cut -f1,2,4 > "$i".alltemplate.cons

		find alignments -type f -name "*.cons" -delete
		rm "$i" "$i".txt
	done

	rm -rf alignments
	cd $OUTPUT
}



create_p3_input() {
	## Get guide id
	# Get amplicon from 'amplicons' dir; assign as template sequence
	## Get probe sequence from 'kmer_results' dir
	## Split probe sequence into probe + pam

	#kmer=$(grep -v "#" $OUTPUT/$IN | grep kmer | awk -F "|" '{ print $2}')

	cd $OUTPUT
	mkdir -p $PRIMERS/p3_input

	logger "Generating Primer3 input files"

	## Create guide lookup table
	cd $AMPLICON

	find . -type f -name "*alltemplate.cons" -exec cat {} + |
	sort -k2b,2 > alltemplates.cons

	# Exit if no templates found
	[[ -s $AMPLICON/alltemplates.cons ]] || { error_exit "No templates found!"; return; }

	find . -type f -name "list*.cons" -delete

	cut -f2 alltemplates.cons | sort | uniq > $KMER_FINAL/guides_final.txt

	seqkit grep -j "$MAX_JOBS" -f $KMER_FINAL/guides_final.txt $KMER_FINAL/query_uniq_target_kmers_nohomopolymer.fasta | \
	tr "\n>" "\t\n>" | \
	sort -k1,1 | \
	grep -v "^$" | \
	sed 's/\t$//g' | \
	sort -k1b,1 > allguides.txt
	
	join -1 2 -2 1 -t $'\t' alltemplates.cons allguides.txt | awk '{ print $2, $4, $3 }' OFS="\t" > $AMPLICON/guides_lookup.txt

	## Add information on guide location within template sequence
	awk '{ print $0, match($3, $2)","'$kmer'"" }' OFS="\t" $AMPLICON/guides_lookup.txt > $AMPLICON/guides_lookup_final.txt

	## Remove guide directories from 'amplicons' dir
	#logger "Removing individual guide directories"
	if [[ -s $AMPLICON/guides_lookup_final.txt ]]
	then
		rm -rf alignments
		## Split guide lookup file into 500 guides per file
		split -l 500 $AMPLICON/guides_lookup_final.txt $AMPLICON/guides_subset -d
	else 
		error_exit "File guide_lookup_final.txt is missing"
	fi

	## Choice of PCR or RPA primer design
	if [[ $AMP_METHOD == "pcr" ]]
	then
		primerminsize=18
		primermaxsize=25
		primeroptsize=22
		primermintm=60
		primermaxtm=75
		primermingc=40
		primermaxgc=60
	else
		primerminsize=30
		primermaxsize=35
		primeroptsize=32
		primermintm=0
		primermaxtm=100
		primermingc=30
		primermaxgc=70
	fi

	for i in $(ls guides_subset* | awk -F "_" '{ print $2 }')
	do
		while IFS=$'\t' read guide oligo sequence target
		do
			{
				echo "SEQUENCE_ID=$guide"
				echo "SEQUENCE_TEMPLATE=$sequence"
				echo "SEQUENCE_INTERNAL_OLIGO=$oligo"
				echo "PRIMER_TASK=pick_pcr_primers_and_hyb_probe"
				echo "PRIMER_PICK_LEFT_PRIMER=1"
				echo "PRIMER_PICK_INTERNAL_OLIGO=1"
				echo "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1"
				echo "PRIMER_PICK_RIGHT_PRIMER=1"
				echo "PRIMER_INTERNAL_MIN_SIZE=18"
				echo "PRIMER_INTERNAL_MAX_SIZE=30"
				echo "PRIMER_INTERNAL_OPT_SIZE=$kmer"
				echo "PRIMER_MIN_SIZE=$primerminsize"
				echo "PRIMER_OPT_SIZE=$primeroptsize"
				echo "PRIMER_MAX_SIZE=$primermaxsize"
				echo "PRIMER_MAX_TM=$primermaxtm"
				echo "PRIMER_MIN_TM=$primermintm"
				echo "PRIMER_INTERNAL_MIN_TM=0"
				echo "PRIMER_INTERNAL_MAX_TM=100"
				echo "PRIMER_MIN_GC=$primermingc"
				echo "PRIMER_MAX_GC=$primermaxgc"
				#echo "PRIMER_NUM_RETURN=$primercount"
				echo "PRIMER_MAX_NS_ACCEPTED=0"
				echo "SEQUENCE_TARGET=$target"
				echo "PRIMER_LOWERCASE_MASKING=1"
				echo "PRIMER_PRODUCT_SIZE_RANGE=90-300"
				echo "P3_FILE_FLAG=0"
				echo "PRIMER_EXPLAIN_FLAG=1"
				echo "="
			} >> $PRIMERS/p3_input/"$i".input.p3
		done < guides_"$i"
	done

	rm alltemplates.cons allguides.txt guides_subset* guides_lookup.txt
}	



run_p3() {
	mkdir -p $PRIMERS/p3_output
	cd $PRIMERS
	ls p3_input/*.input.p3 | \
	awk -F "/" '{ print $2 }' | \
	sed 's/.input.p3//g' | \
	awk '{ print "primer3_core < p3_input/"$0".input.p3 > p3_output/"$0".output.p3" }' \
	> run_p3.sh

	if [[ -s run_p3.sh ]]
	then	
		logger "Designing primers"
		cat run_p3.sh | parallel -j "$MAX_JOBS"
	else
		error_exit "Primer3 input file not found!"
	fi
	rm run_p3.sh
	# ## Filter out guides with no primers designed
	# for i in p3_output/*.output.p3; do
	# 	sed 's/^=/\n/g' $i | awk -v RS= '/\nPRIMER_LEFT_0/ {print $0"\n" >> "'"$i"'.withprimers"; next}'; 
	# done

	cd $OUTPUT
}



parse_p3_output() {
	
	cd $PRIMERS
	logger "Parsing primer design results"

	cd p3_output

	for i in *output.p3; do \
		awk -F'[=,]' '/^SEQUENCE_ID/ { guide=$2; set=1 } \
		
		/^PRIMER_LEFT_._SEQUENCE/ { printf guide "\t" "set"set "\t" $2 } \
		/^PRIMER_RIGHT_._SEQUENCE/ { printf "\t" $2 } \
		/^PRIMER_INTERNAL_._SEQUENCE/ { printf "\t" $2 } \

		/^PRIMER_LEFT_.\=/ { printf "\t" $3 } \
		/^PRIMER_RIGHT_.\=/ { printf "\t" $3 }  \
		/^PRIMER_INTERNAL_.\=/ { printf "\t" $3 } \

		/^PRIMER_LEFT_._GC_PERCENT/ { printf "\t" $2 } \
		/^PRIMER_RIGHT_._GC_PERCENT/ { printf "\t" $2 } \
		/^PRIMER_INTERNAL_._GC_PERCENT/ { printf "\t" $2 } \

		/^PRIMER_PAIR_._PRODUCT_SIZE/ { printf "\t" $2 "\n"; set++ }' $i >> allguides.p3.out; done


	cut -f1 allguides.p3.out | awk -F "_" '{ print $1 }' > guides.txt
	
	## Associative array from prev file; extract pam
	local pam=$(awk -F "\t" 'FNR==NR{pam_lookup[$1]=$2; next} {print pam_lookup[$1]}' "$KMER_FINAL"/uniq_target_kmers_final_prevalence.tsv guides.txt)

	## Combine
	paste allguides.p3.out guides.txt <(echo -e "$pam") --delimiters "\t" > allguides_p3_parsed.tsv

	## Need to put guide seq after probe seq column, add pam

	## Add header
	local header=$(printf "%s\t" \
		"guide_set" \
		"primer_set" \
		"fwd_primer" \
		"rev_primer" \
		"guide_seq" \
		"fwd_primer_length" \
		"rev_primer_length" \
		"guide_length" \
		"fwd_primer_gc" \
		"rev_primer_gc" \
		"guide_gc" \
		"product_size" \
		"guide" "pam")
	sed -i "1s/^/$header\n/" allguides_p3_parsed.tsv

	awk -F "\t" '{ print $1,$13,$5,$14,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12 }' OFS="\t" allguides_p3_parsed.tsv \
		> allguides_primers_output.tsv

	## Summarize multiple guides per primer
	cat allguides_primers_output.tsv | \
		csvtk summary -t -f guide:uniq,guide:countunique,guide_set:uniq,guide_set:countunique -g primer_set,fwd_primer,rev_primer \
		> allguides_primers_output_grouped.tsv 

	## Write guides with primers designed to file (file not used elsewhere)
	cut -f2 allguides_primers_output.tsv | grep -v guide | sort | uniq > $KMER_FINAL/guides_shortlisted_withprimers.txt

	rm guides.txt allguides_p3_parsed.tsv
	
	cd $OUTPUT

}



run_fast_ispcr() {

	## This uses the UCSC ispcr script by Jim Kent
	## We will filter output by scores (sum mismatches in fwd+rev primers can be determined from scores)
	
	cd $PRIMERS

	## Remove target template files (for primersearch)
	find $ALIGN_OUT -type f -name "*_template.fasta" -printf '%p\n' | xargs -r rm

	local method=$1
	mkdir -p primer_check/target primer_check/nontarget 
	#mkdir -p primer_check/target_amplicons primer_check/nontarget_amplicons

	## Set headers for parsed output
	local header=$(printf "%s\t" \
		"genome" \
		"sequence" \
		"primer_set" \
		"start" \
		"end" \
		"strand" \
		"amplicon_length" \
		"ispcr_score" | \
		sed 's/\t$//g')

	## Prepare ispcr input
	if [[ $method == "pangenome" ]]
	then
		awk -F "\t" '{ print "primer|"$1"", $6, $7 }' OFS="\t" p3_output/allgenes_primers_output.tsv | \
		sed 's/--/|/1' | sed 's/-set/|set/g' | grep -v guide_set \
		> primer_check/allprimers_lookup.tsv

	elif [[ $method == "kmer" ]]
		then
			awk -F "\t" '{ print "primer|"$5"|"$4"|"$1, $6, $7 }' OFS="\t" p3_output/allguides_primers_output.tsv | \
			grep -v primer_set > primer_check/allprimers_lookup.tsv

		fi

	## Dereplicate primer sets
	cut -f2-3 primer_check/allprimers_lookup.tsv | \
	sort | uniq | nl | \
	sed 's/^ *//g' | \
	sed 's/^/primerset/1' \
	> primer_check/ispcr_input.tsv


	logger "Running in silico PCR (fast mode) on target genomes"
	cd $PRIMERS/primer_check/target

	## No need to differentiate between pangenome and kmer- run on whole genome assembly
	## Allow one mismatch for target genomes, 3 for nontarget genomes
	## 1 mismatch in target to have backups if no 100% sensitive perfect match primers available for target genomes
	## 3 mismatches in nontarget genomes to discard primers if 3 mismatches result in primer binding (make up for non-specificity of RPA)

	split -l 200 $OUTPUT/genomes_target.txt subset. -d -a 3

	for file in subset.*; do \
		cat $file | \
		parallel -j "$MAX_JOBS" \
		'isPcr -maxSize=500 \
		-out=bed \
		'"$GENOMES_TARGET"'/"{}"_genomic.fna \
		'"$PRIMERS"'/primer_check/ispcr_input.tsv \
		"{}".ispcr'; \
		for i in *.ispcr; do awk '{ print FILENAME, $0, $3-$2 }' OFS="\t" $i | sed 's/.ispcr//g' >> "$PRIMERS"/primer_check/combined_ispcr.txt; done
		rm *.ispcr
	done

	rm subset.*

	cd $PRIMERS/primer_check
	mkdir -p target_amplicons

	# Write amplicons to file (for determining cas13a gRNA strand)
	if [[ $METHOD == "kmer" ]]; then
		cat $OUTPUT/genomes_target.txt | \
		parallel -j "$MAX_JOBS" \
		'isPcr -maxSize=500 \
		-out=fa \
		'"$GENOMES_TARGET"'/"{}"_genomic.fna \
		ispcr_input.tsv \
		target_amplicons/"{}"_amp.fasta; sed -i "s/primerset/>"{}"--primerset/g" target_amplicons/"{}"_amp.fasta'

		cat target_amplicons/*_amp.fasta > amplicons_target.fasta
		rm -rf target_amplicons
		sed -i 's/^>[^ ]* \([^ ]*\).*/>\1/' amplicons_target.fasta
	fi
	
	## Combine all target output
	cd target
	#for i in *.ispcr; do awk '{ print FILENAME, $0, $3-$2 }' OFS="\t" $i | sed 's/.ispcr//g' >> combined_ispcr.txt; done
	awk -F "\t" '{ print $1, $2, $5, $3, $4, $7, $8, $6 }' OFS="\t" "$PRIMERS"/primer_check/combined_ispcr.txt > all_ontarget_parsed.tsv
	rm combined_ispcr.txt

	sed -i "1s/^/$header\n/" all_ontarget_parsed.tsv
	find . -type f -name "*.ispcr" | xargs -L1 rm

	cd $PRIMERS/primer_check

	logger "Running in silico PCR (fast mode) on non-target genomes"

	cat $OUTPUT/genomes_offtarget.txt | \
	parallel -j "$MAX_JOBS" \
	'isPcr -maxSize=3000 \
	-out=bed \
	'"$GENOMES_OFFTARGET"'/"{}"_genomic.fna \
	ispcr_input.tsv \
	nontarget/"{}".ispcr'

	# Write amplicons to file
	# cat $OUTPUT/genomes_offtarget.txt | \
	# parallel -j "$MAX_JOBS" \
	# 'isPcr -maxSize=3000 \
	# -out=fa \
	# '"$GENOMES_OFFTARGET"'/"{}"_genomic.fna \
	# ispcr_input.tsv \
	# nontarget_amplicons/"{}".amplicons'

	## Combine all nontarget output
	cd nontarget
	for i in *.ispcr; do awk '{ print FILENAME, $0, $3-$2 }' OFS="\t" $i | sed 's/.ispcr//g' >> combined_ispcr.txt; done
	awk -F "\t" '{ print $1, $2, $5, $3, $4, $7, $8, $6 }' OFS="\t" combined_ispcr.txt > all_nontarget_parsed.tsv
	rm combined_ispcr.txt
	
	sed -i "1s/^/$header\n/" all_nontarget_parsed.tsv
	find . -type f -name "*.ispcr" | xargs -L1 rm

	cd $OUTPUT

}

primer_orientation_cas13() {
	# Guide RNA should be in the same orientation as reverse primer
	# Output new column "cas13a_use"- "yes" if can be used as-is, "swap primers" if primers have to be swapped

	cd $PRIMERS/primer_check

	# Deduplicate combined target amplicons
	seqkit rmdup -P -s amplicons_target.fasta -o amplicons_target_uniq.fasta

	bbmap.sh \
	in="$PRIMERS"/guides_input.fasta \
	nodisk \
	noheader=t \
	threads="$MAX_JOBS" \
	ref=amplicons_target_uniq.fasta \
	ambig=all \
	vslow \
	perfectmode \
	maxsites=100000 \
	outm=amplicons_target_uniq.sam

	cut -f1-3 amplicons_target_uniq.sam | 
	awk '{ print $1, $2, $3 }' OFS="\t" | 
	sed 's/--/\t/1' | 
	cut -f1-2,4 | 
	sort | uniq |
	awk -F "\t" '{ if ($2==0||$2==256) print $1, $3, "swap primers"; else if ($2==16||$2==272) print $1, $3, "yes" }' OFS="\t" |
	sort | uniq > primer_cas13a_compatibility.tsv

	sed -i '1i guide\tprimerset\tprimer_cas13a_compatibility' primer_cas13a_compatibility.tsv
	#rm combined_target_amplicons.sam combined_target_amplicons_derep.fasta

	cd $OUTPUT
}


get_primerstats() {
	cd $PRIMERS
	
	local targetcount=$(wc -l <$OUTPUT/genomes_target.txt)
	local nontargetcount=$(wc -l <$OUTPUT/genomes_offtarget.txt)
	
	local targetcheck=primer_check/target
	local nontargetcheck=primer_check/nontarget

	if [[ $METHOD == "pangenome" ]]
	then
		local method="pangenome"
		local myfile=$PRIMERS/p3_output/allgenes_primers_output.tsv
		local annot=$PANGENOME/genes_90pct_annotation.tsv
		local cas13="NULL"
	elif [[ $METHOD == "kmer" ]]
	then
		local method="kmer"
		local myfile=$PRIMERS/p3_output/allguides_primers_output.tsv
		local annot="NULL"
		local cas13=$PRIMERS/primer_check/primer_cas13a_compatibility.tsv
	fi
	

	if [[ -s "$nontargetcheck"/all_nontarget_parsed.tsv ]]
	then
		
		if [[ -s guide_offtarget_summary.tsv ]]
		then
			logger "Calculating primer sensitivity and specificity (with guide offtarget hits)"
			offtarget="$nontargetcheck"/all_nontarget_parsed.tsv
			offtarget_summary=$PRIMERS/guide_offtarget_summary.tsv
		else
			logger "Calculating primer sensitivity and specificity (no guide offtarget hits)"
			offtarget="$nontargetcheck"/all_nontarget_parsed.tsv
			offtarget_summary=NULL
		fi
	
	else
		echo "No amplified sequences from non-target genomes"
		logger "Calculating primer sensitivity"
		if [[ -s guide_offtarget_summary.tsv ]]
		then
			offtarget=NULL
			offtarget_summary=$PRIMERS/guide_offtarget_summary.tsv
		else
			offtarget=NULL
			offtarget_summary=NULL
		fi
	fi

	get_primerstats.R "$targetcheck"/all_ontarget_parsed.tsv \
		$offtarget \
		$PRIMERS/primer_check/ispcr_input.tsv \
		$targetcount \
		$nontargetcount \
		$myfile \
		$OUTPUT/pathogd_primers_stats.tsv \
		$PRIMERS/allguides_target_summary.tsv \
		$offtarget_summary \
		$annot \
		$cas13 \
		$method
}


clean_files() {
	
	logger "Cleaning up intermediate files"
	local method=$1
	
	if [[ $method == "pangenome" ]]
	then
		cd $PANGENOME
		rm -rf species_specific guides_precheck
		rm -f remove* genes_remove.txt conserved*.txt 
		cd $OUTPUT
		rm -rf kmers_offtarget human_genome
		rm -rf primers/guide_prevalence
		rm -rf primers/p3_input 
		#primers/p3_output
		
	elif [[ $method == "kmer" ]]
	then
		cd $OUTPUT
		rm -rf kmers_offtarget kmers_target kmers_target_pam human_genome
		rm -rf align_out amplicons
		rm -rf primers/guide_prevalence
		rm -rf primers/p3_input 
		#primers/p3_output
	fi
	cd $OUTPUT
}




#############################################################################
