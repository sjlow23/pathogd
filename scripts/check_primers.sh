#!/usr/bin/env bash

usage() {
	echo -e "\nUsage:"
	echo -e "  $0 [-i <primers input>] [-d <genome db directory>] [-a <get_amplicons>] [-o <output file>] [-t <number of cpus>]\n"
	echo -e "	-i	Input file containing three columns in the following order: primer id, forward primer sequence, reverse primer sequence\n"
	echo -e "	-d	Directory containing target and non-target genome databases\n"
	echo -e "	-a  Whether to write amplicons to file (yes/no)\n"
	echo -e	"	-o	Output file name (default: pathogd_primer_prevalence.tsv)\n"
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
	while getopts "i:d:a:o:t:h" opt
	do
		case $opt in
			i) IN=$OPTARG;;
			d) GENOMEDIR=$OPTARG;;
			a) AMPLICON=$OPTARG;;
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
    echo -e "\nInput file \"-i\" with primer sequences must be provided!\n"
    usage
fi


if [[ -z "$GENOMEDIR" ]]; then
    echo -e "\nDirectory \"-d\" containing target and non-target genomes must be specified!\n"
    usage
fi


if [[ -z "$AMPLICON" ]]; then
	AMPLICON="no"
fi


if [[ -z "$OUTPUT" ]]; then
	OUTPUT="pathogd_primer_prevalence.tsv"
fi


if [[ -z "$CPU" ]]; then
	CPU=8
fi




# Define working and output directories
MYDIR=$PWD
GENOMEDIR=$(basename $GENOMEDIR)
OUTDIR=$MYDIR/$GENOMEDIR/primer_stats
IN=$(readlink -f $IN)
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


## Get primer lengths
awk '{ print $0, length($2), length($3) }' $IN | awk '{ print $1, $4, $5 }' OFS="\t" > $OUTDIR/primer_lengths.txt
echo -e "primer_set\tfwd_primer_length\trev_primer_length\n($(cat $OUTDIR/primer_lengths.txt)" > $OUTDIR/primer_lengths.txt



## Function to run isPcr

run_ispcr() {

	if [[ $1 == "target" ]]
	then
		local mydir="$OUTDIR"/target
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_target
		local genomelist="$OUTDIR"/genomes_target.txt
		local genomecount=$(cat "$genomelist" | wc -l)
		logger "Running in silico PCR on "$genomecount" target genomes"
		local maxampliconsize=500
		local outfile="$mydir/target_combined_ispcr.tsv"
	elif [[ $1 == "non-target" ]]
	then
		local mydir="$OUTDIR"/offtarget
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_offtarget
		local genomelist="$OUTDIR"/genomes_offtarget.txt
		local genomecount=$(cat "$genomelist" | wc -l)
		logger "Running in silico PCR on "$genomecount" non-target genomes"
		local maxampliconsize=5000
		local outfile="$mydir/nontarget_combined_ispcr.tsv"
	fi

	## Set headers for parsed output
	header=$(printf "%s\t" \
		"genome" \
		"sequence" \
		"primer_set" \
		"start" \
		"end" \
		"strand" \
		"amplicon_length" \
		"ispcr_score" | \
		sed 's/\t$//g')

	mkdir -p "$mydir"
	cd "$mydir"

	split -l 200 $genomelist subset. -d -a 3

	for file in subset.*; do \
		cat $file | \
		parallel -j "$CPU" \
		'isPcr -maxSize='"$maxampliconsize"' \
		-out=bed \
		'"$genomedir"'/"{}"_genomic.fna \
		'"$IN"' \
		"{}".ispcr'; \
		for i in *.ispcr; do awk '{ print FILENAME, $0, $3-$2 }' OFS="\t" $i | 
			sed 's/.ispcr//g' >> $outfile; done
		rm *.ispcr $file
	done

	cd $OUTDIR

	if [[ -s $outfile ]]; then
		logger "Parsing primer hits to "$1" genomes"
		awk -F "\t" '{ print $1, $2, $5, $3, $4, $7, $8, $6 }' OFS="\t" $outfile > ispcr_"$1".tsv
		sed -i "1s/^/$header\n/" ispcr_"$1".tsv
		csvtk join -f primer_set ispcr_"$1".tsv primer_lengths.txt -t > tmp

		## Add mismatches
		csvtk mutate2 -n sum_primer_lengths -t -e '$fwd_primer_length+$rev_primer_length' -w 0 tmp | 
		csvtk mutate2 -n sum_primer_mismatches -t -e '$sum_primer_lengths-($sum_primer_lengths*$ispcr_score/1000)' |
		csvtk round -t -f sum_primer_mismatches -n 0 > ispcr_"$1".tsv

		## Get prevalence statistics
		## Group by primer set
		for i in 0 1 2 3 4 5 6; do
			csvtk filter -t -f "sum_primer_mismatches<="$i"" ispcr_target.tsv | 
			csvtk summary -t -i -f "genome:countunique" -g primer_set |
			csvtk mutate2 -n perc_mm"$i" -t -e '$genome:countunique/'"$genomecount"'*100' |
			csvtk round -t -f perc_mm"$i" -n 1 |
			csvtk cut -t -f -genome:countunique > $OUTDIR/"$1"_mm"$i".tsv; done
		csvtk join -t -f primer_set $OUTDIR/"$1"_mm*.tsv | 
		csvtk rename -t -f 2-8 -n perc_mm0,perc_mm1,perc_mm2,perc_mm3,perc_mm4,perc_mm5,perc_mm6 > $OUTDIR/"$1"_$OUTPUT

		rm -rf "$mydir" tmp "$1"_mm*.tsv

	else
		logger "No primer hits to "$1" genomes"
		rm -rf "$mydir"
	fi

	cd $MYDIR

}


get_amplicons() {

	cd $OUTDIR
	if [[ $1 == "target" ]]
	then
		logger "Getting amplicons for target genomes"
		local mydir="$OUTDIR"/target
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_target
		local genomelist="$OUTDIR"/genomes_target.txt
		local maxampliconsize=500
		local outfile="$OUTDIR/target_amplicons.fna"
	elif [[ $1 == "non-target" ]]
	then
		logger "Getting amplicons for non-target genomes"
		local mydir="$OUTDIR"/offtarget
		local genomedir="$MYDIR"/"$GENOMEDIR"/genomes_offtarget
		local genomelist="$OUTDIR"/genomes_offtarget.txt
		local maxampliconsize=5000
		local outfile="$OUTDIR/nontarget_amplicons.fna"
	fi
	
	mkdir -p "$mydir"
	cd "$mydir"

	split -l 200 $genomelist subset. -d -a 3

	for file in subset.*; do \
		cat $file | \
		parallel -j "$CPU" \
		'isPcr -maxSize='"$maxampliconsize"' \
		-out=fa \
		'"$genomedir"'/"{}"_genomic.fna \
		'"$IN"' \
		"{}".fasta'; \
		for i in *.fasta; do cat $i >> $outfile; done
		rm *.fasta $file
	done

	cd $OUTDIR
	rm -rf "$mydir" $genomelist

	[[ -s $outfile ]] || (rm $outfile; logger "No amplicons obtained from "$1" genomes")

}


clean_files() {
	cd $OUTDIR
	[[ -s genomes_*.txt ]] && rm genomes_*.txt
	rm primer_lengths.txt
}



if [[ $AMPLICON == "yes" ]]
then
	if [[ $DB == "both" ]]; then
		run_ispcr target
		run_ispcr non-target
		get_amplicons target
		get_amplicons non-target
	elif [[ $DB == "target_only" ]]; then
		run_ispcr target
		get_amplicons target
	elif [[ $DB == "offtarget_only" ]]; then
		run_ispcr non-target
		get_amplicons non-target
	fi
	clean_files
	logger "In silico PCR complete"

else
	if [[ $DB == "both" ]]; then
		run_ispcr target
		run_ispcr non-target
	elif [[ $DB == "target_only" ]]; then
		run_ispcr target
	elif [[ $DB == "offtarget_only" ]]; then
		run_ispcr non-target
	fi
	clean_files
	logger "In silico PCR complete"
fi