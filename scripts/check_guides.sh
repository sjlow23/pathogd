#!/usr/bin/env bash

usage() {
	echo -e "\nUsage:"
	echo -e "  $0 [-i <guides input>] [-d <genome db directory>] [-o <output file>] [-t <number of cpus>]\n"
	echo -e "	-i	Input fasta file containing guide RNA sequences\n"
	echo -e "	-d	Directory containing target and non-target genome databases\n"
	echo -e	"	-o	Output file name (default: pathogd_guide_prevalence.tsv)\n"
	echo -e	"	-t 	Number of cpus (default: 8)\n"
	echo -e "	-h  Print this help menu\n"
	exit 0
}


error_exit() {
	echo "Error: $1" 1>&2
	exit 1
}


if [[ $# -lt 4 ]]
then
	usage
	exit 0
else
	while getopts "i:d:o:t:h" opt
	do
		case $opt in
			i) IN=$OPTARG;;
			d) GENOMEDIR=$OPTARG;;
			o) OUTPUT=$OPTARG;;
			t) CPU=$OPTARG;;
			h) usage;;
			\?) echo -e "\nInvalid option: -$OPTARG"
				usage;;
			:) echo -e "\nOption -$OPTARG requires an argument!\n"
				exit 0;;
		esac
	done
	shift $(( OPTIND - 1 ))
fi


# Exit if mandatory arguments missing

if [[ -z "$IN" ]]; then
    echo -e "\nInput file \"-i\" with guide sequences must be provided!\n"
    usage
fi


if [[ -z "$GENOMEDIR" ]]; then
    echo -e "\nDirectory \"-d\" containing target and non-target genomes must be specified!\n"
    usage
fi

if [[ -z "$OUTPUT" ]]; then
	OUTPUT="pathogd_guides_prevalence.tsv"
fi


if [[ -z "$CPU" ]]; then
	CPU=8
fi




# Define working and output directories
MYDIR=$PWD
GENOMEDIR=$(basename $GENOMEDIR)
OUTDIR=$MYDIR/$GENOMEDIR/guide_stats
IN=$(readlink -f $IN)
IN_PAM=""${IN%.*}"_pam.fasta"
LOG=$OUTDIR/pathogd.log


logger() {
	echo "[$(date +"%F %T")] INFO: $1" >> $LOG
}


## Check target or non-target database directories present and not empty
cd $MYDIR/$GENOMEDIR
mkdir $OUTDIR


if [[ "$(ls -A ./genomes_target)" && "$(ls -A ./genomes_offtarget)" ]]
then
	logger "Both target and non-target genome databases present"
	find genomes_target \( -type f -o -type l \) -name "*.fna" -o -name "*.fasta" | awk -F "/" '{ print $2 }' | sed 's/_genomic.fna//g' > $OUTDIR/genomes_target.txt
	find genomes_offtarget \( -type f -o -type l \) -name "*.fna" -o -name "*.fasta" | awk -F "/" '{ print $2 }' | sed 's/_genomic.fna//g' > $OUTDIR/genomes_offtarget.txt
	DB="both"
elif [[ "$(ls -A ./genomes_target)" && ! "$(ls -A ./genomes_offtarget)" ]]
then
	logger "Only target genome database present"
	logger "Running workflow for target genomes only"
	find genomes_target \( -type f -o -type l \) -name "*.fna" -o -name "*.fasta" | awk -F "/" '{ print $2 }' | sed 's/_genomic.fna//g' > $OUTDIR/genomes_target.txt
	DB="target_only"
elif [[ ! "$(ls -A ./genomes_target)" && "$(ls -A ./genomes_offtarget)" ]]
then
	logger "Only non-target genome database present"
	logger "Running workflow for non-target genomes only"
	find genomes_offtarget \( -type f -o -type l \) -name "*.fna" -o -name "*.fasta" | awk -F "/" '{ print $2 }' | sed 's/_genomic.fna//g' > $OUTDIR/genomes_offtarget.txt
	DB="offtarget_only"
else
	logger "No database found! Please provide either target or non-target genomes"
	error_exit "No database(s) found..Exiting.."
fi


## Get number of guide sequences
guidecount=$(grep -c ">" $IN)


## Generate input sequences with canonical PAM
awk '{ print "TTTA"$0; print "TTTC"$0; print "TTTG"$0; print "TTTT"$0; }'

seqkit fx2tab $IN | 
	awk '{ print $1"--TTTA", "TTTA"$2;
	print $1"--TTTC", "TTTC"$2; 
	print $1"--TTTG", "TTTG"$2;
	print $1"--TTTT", "TTTT"$2 }' OFS="\t" |
	seqkit tab2fx > "$IN_PAM"


## Function to run guide mapping
process_cigar() {

	local input_file="$1"

	# Process SAM file line by line, excluding header lines
	while IFS=$'\t' read -r -a fields; do
		# if [[ ${fields[0]} == @* ]]; then
		# 	continue
		# fi

		# Get bounds tag
		boundstag="${fields[-1]}"

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
    	 num_S=0
    	 for component in "${components[@]}"; do
    	 	countx="${component%[=X]}"
    	 	counts="${component%[=S]}"
    	 	count="${component%[=X|S]}"
    	 	op="${component: -1}"
    	 	expanded_string+=$(printf "$op%.0s" $(seq 1 "$count"))
    	 	if [[ $op == "X" ]]; then
    	 		num_X=$((num_X + countx))
    	 	fi
    	 	if [[ $op == "S" ]]; then
    	 		num_S=$((num_S + counts))
    	 	fi
    	 done

	# Reverse the expanded string if the mapping is on the reverse strand
	if $is_reversed; then
		expanded_string=$(echo "$expanded_string" | rev)
	fi

	# Output the read name and the expanded CIGAR string
	echo -e "${read_name}\t${expanded_string}\t${num_X}\t${num_S}\t${boundstag}"
	
	done < "$input_file"

}


## Get number of mismatches
#awk '{for(i=1;i<=NF;i++) if($i~"^NM:i") {print $i}}'


map_guides() {

	cd $OUTDIR

	export -f process_cigar

	if [[ $1 == "target" ]]
	then
		local mydir="$OUTDIR"/target
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_target
		local genomelist="$OUTDIR"/genomes_target.txt
		local genomecount=$(cat "$genomelist" | wc -l)
		
		logger "Mapping "$guidecount" guides to "$genomecount" target genomes"
		local outfile="$mydir/target_guides_mapped.tsv"
		
	elif [[ $1 == "nontarget" ]]
	then
		local mydir="$OUTDIR"/offtarget
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_offtarget
		local genomelist="$OUTDIR"/genomes_offtarget.txt
		local genomecount=$(cat "$genomelist" | wc -l)
		
		logger "Mapping "$guidecount" guides to "$genomecount" non-target genomes"
		local outfile="$mydir/nontarget_guides_mapped.tsv"
	fi

	mkdir -p "$mydir"
	cd "$mydir"

	split -l 200 $genomelist subset. -d -a 3

	for file in subset.*; do \
		cat $file | \
		parallel -j "$CPU" \
		'bbmap.sh in='"$IN_PAM"' \
		ref='"$genomedir"'/"{}"_genomic.fna \
		nodisk \
		noheader=t \
		ambig=all \
		vslow \
		subfilter=1 \
		editfilter=1 \
		indelfilter=0 \
		lengthtag=t \
		nmtag=t \
		nhtag=f \
		amtag=f \
		boundstag=t \
		maxsites=1000000 \
		threads=8 \
		outm="{}"_wg.sam; \
		awk -F "\t" '\''$6 !~ /M/'\'' "{}"_wg.sam > "{}"_clean.sam; \
		process_cigar "{}"_clean.sam | awk "{ print \"{}\", \$0 }" OFS="\t" > "{}".parsed;
		rm "{}"_wg.sam "{}"_clean.sam'

		find . -name "*.parsed" -type f -print0 | xargs -0 -n100 cat >> "$outfile"
		find . -name "*.parsed" | xargs rm
		rm $file
		
	done

	cd $OUTDIR

	if [[ -s $outfile ]]; then
		logger "Parsing guide hits to "$1" genomes"
		#sed -i '1i genome\tguide\tcigar\tmismatch_count\tsoftclip_count\tbounds_tag' $outfile

		# Remove if contain soft-clip and not at edge of contig
		awk -F "\t" '{ if ($5==0 || ($5!=0 && $6=="XB:Z:O")) print $1, $2, $3, $4, $5 }' OFS="\t" $outfile > "$outfile".2

		# Filter by allowed mismatches
		# Max 2 mismatches (both must be PAM-distal)
		awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $4, substr($3,1,4), substr($3,5,7), substr($3,12)}' "$outfile".2 > tmp
		awk 'BEGIN{FS=OFS="\t"} { mm_pam = gsub(/X/, "X", $4); 
		mm_seed = gsub(/X/, "X", $5); 
		mm_distal = gsub(/X/, "X", $6); print $0, mm_pam, mm_seed, mm_distal }' tmp |
		## Filter alignments with 2 or lesser mismatches
		awk -F "\t" '$3<=2 && $7==0 && $8==0' > $outfile
		rm tmp "$outfile".2

		sed -i '1i genome\tguide\tmismatch_count\tcigar_pam\tcigar_seed\tcigar_distal\tmm_pam\tmm_seed\tmm_distal' $outfile

		# Summarize by guide (for 0, 1 and 2 mismatches)
		for i in 0 1 2; do
			csvtk filter -t -f "mismatch_count<="$i"" $outfile |
			# separate guide name from PAM before calculating prevalence
			csvtk sep -t -f 2 -s '--' -n guide,pam -N 2 --na NA -R |
			csvtk summary -t -i -f "genome:countunique,genome:count" -g guide |
			csvtk mutate2 -n guide_prevalence_mm"$i" -t -e '$genome:countunique/'"$genomecount"'*100' |
			# calculate copy number using only genomes with hits, not original input genome count
			csvtk mutate2 -n guide_avg_copy_number_mm"$i" -t -e '$genome:count/$genome:countunique' |
			csvtk round -t -f guide_prevalence_mm"$i" -n 1 |
			csvtk round -t -f guide_avg_copy_number_mm"$i" -n 1 | 
			csvtk cut -t -f -genome:countunique,-genome:count > $OUTDIR/"$1"_mm"$i"_$OUTPUT; done

		# Join if file not empty
		find . -size 0 -name "*mm*.tsv" -type f -print -delete

		csvtk join -t -f guide $OUTDIR/"$1"_mm*_$OUTPUT > $OUTDIR/"$1"_$OUTPUT
		#csvtk rename -t -f 2,4,6 -n guide_prev_mm0,guide_prev_mm1,guide_prev_mm2

		rm $OUTDIR/"$1"_mm*_$OUTPUT

	else
		logger "No guide hits to "$1" genomes"
		#rm -rf "$mydir"
	fi

	cd $MYDIR

}




if [[ $DB == "both" ]]; then
	map_guides target
	map_guides nontarget
	logger "Guide mapping complete"

elif [[ $DB == "target_only" ]]; then
	map_guides target
	logger "Guide mapping complete"

elif [[ $DB == "offtarget_only" ]]; then
	map_guides nontarget
	logger "Guide mapping complete"
fi


