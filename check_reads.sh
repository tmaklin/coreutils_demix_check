#!/bin/bash
## Input:
##        --abundances: abundances file from mSWEEP   (required)
##        --cluster: which cluster to check           (required)
##        --reference: output from setup_reference.sh (required)
##        --fwd: forward reads from mGEMS             (required)
##        --rev: reverse reads from mGEMS             (required)
##        --threads: number of threads                (default: 1)
##        --tmpdir: tmp directory                     (default: working directory)
##        --bufsize: buffer size for `sort`           (default: 4000M)
##        --verbose: echo commands run                (default: silent)
##        --help: print this message

set -euo pipefail

arg_parser() {
    ## Default values
    abundances=""   ## Exit if not supplied
    cluster=""
    reference=""
    fwd=""
    rev=""
    threads=1     ## 1 thread
    bufsize=4000"M"  ## 4000MB
    tmpdir="."    ## Working directory
    verbose=false ## Silent

    # Parse values
    while [ $# -gt 0 ]; do
        if [[ $1 == *"--"* ]]; then
	    if [ "$1" == "--verbose" ]; then
		verbose="true"
	    elif [ "$1" == "--help" ]; then
		echo -e "check_reads.sh:\n\t--abundances: abundances file from mSWEEP   (required)\n\t--cluster: which cluster to check           (required)\n\t--reference: output from setup_reference.sh (required)\n\t--fwd: forward reads from mGEMS             (required)\n\t--rev: reverse reads from mGEMS             (required)\n\t--threads: number of threads                (default: 1)\n\t--tmpdir: tmp directory                     (default: working directory)\n\t--bufsize: buffer size for sort             (default: 4000M)\n\t--verbose: echo commands run                (default: silent)\n\t--help: print this message"
		exit 0
	    else
		param="${1/--/}"
		declare -g $param="$2"
	    fi
        fi
        shift
    done

    if [ "$abundances" == "" ]; then
	echo "--abundances was not supplied!"
	exit 1
    fi


    if [ "$cluster" == "" ]; then
	echo "--cluster was not supplied!"
	exit 1
    fi


    if [ "$reference" == "" ]; then
	echo "--reference was not supplied!"
	exit 1
    fi


    if [ "$fwd" == "" ]; then
	echo "--fwd was not supplied!"
	exit 1
    fi


    if [ "$rev" == "" ]; then
	echo "--rev was not supplied!"
	exit 1
    fi
}

arg_parser $@

if [ "$verbose" == "true" ]; then
    set -x
fi

export LC_ALL=C
export TMPDIR=$3

abun_in=$abundances
nthreads=$threads
refdir=$reference

abundance=$(grep "$cluster[[:space:]]" $abun_in | cut -f2)

seqtk_res=$(seqtk fqchk $fwd)

total_bases=$(echo $(echo $seqtk_res | grep -o "ALL[[:space:]][0-9]*" | grep -o "[0-9]*"))
total_bases=$(( total_bases * 2 )) ## paired-end data

read_len=$(echo $(echo $seqtk_res | grep -o "avg_len:[[:space:]][0-9]*" | grep -o "[0-9]*"))
read_file_lines=$(zcat --force $fwd | wc -l)
read_count=$(( read_file_lines/4 ))

cluster_avg_len=$(grep "$cluster[[:space:]]" $refdir/ref_clu_comp.tsv | cut -f3)
coverage=$(printf %.2f $(echo "$total_bases/$cluster_avg_len" | bc -l))
coverage_final=$coverage

if (( $(echo "$coverage > 100" | bc -l) )); then
    echo "subsampling"

    r1=$tmpdir/$cluster"_subsampled_1.fastq.gz"
    r2=$tmpdir/$cluster"_subsampled_2.fastq.gz"

    sub_read_count=$(( read_count/(coverage/100) ))
    seqtk sample -s 11 $fwd $sub_read_count | pigz -p $nthreads > $r1
    seqtk sample -s 11 $rev $sub_read_count | pigz -p $nthreads > $r2

    subsampled=1
    coverage_final=100
    notes="Warning: high coverage ($coverage) so $sub_read_count reads were subsampled to reduce the coverage to 100"
elif (( $(echo "$coverage < 10" | bc -l) )); then
    notes="Warning: low coverage ($coverage) so distances may not be accurate"
    subsampled=0
else
    subsampled=0
    notes=""
fi
r1=$fwd
r2=$rev

m=$(printf %.0f $(echo "$coverage_final/10 + 2" | bc -l))

clu_sketch=$tmpdir/$cluster".msh"
mash sketch -p $nthreads -s 10000 -r -m $m -I $cluster -C - -o $clu_sketch $r1 $r2 2> $tmpdir/mash_check.log

sorted_ref_info=$tmpdir/ref_info-$RANDOM".sorted.tsv"
sed '1d' $refdir/ref_info.tsv | sort -T $tmpdir -S $bufsize --parallel=$nthreads > $sorted_ref_info
dist_res=$(mash dist -p $nthreads $refdir/ref.msh $clu_sketch \
	       | sort -T $tmpdir -S $bufsize --parallel=$nthreads \
	       | join -1 1 -2 3 - $sorted_ref_info \
	       | grep "$cluster$" \
	       | datamash min 3 median 3 -t ' ')

mindis=$(echo $dist_res | cut -f1 -d' ')
meddis=$(echo $dist_res | cut -f2 -d' ')

clu_info=$(grep "$cluster[[:space:]]" $refdir"/ref_clu_thr.tsv")
clu_thr=$(echo $clu_info | cut -f3 -d ' ')
clu_dis_same_max=$(echo $clu_info | cut -f4 -d' ')
clu_dis_same_med_all=$(echo $clu_info | cut -f5 -d' ')
clu_dis_diff_med_all=$(echo $clu_info | cut -f6 -d' ')

med_vs_same=$(echo "$meddis - $clu_dis_same_med_all" | bc -l)
med_vs_diff=$(echo "$meddis - $clu_dis_diff_med_all" | bc -l)

if (( $(echo "$mindis < $clu_dis_same_max" | bc -l) )); then
    score="1"
elif (( $(echo "$mindis < $clu_thr" | bc -l) )); then
    score="2"
elif (( $(echo "$med_vs_same < $med_vs_diff" | bc -l) )); then
    score="3"
else
    score="4"
fi

echo -e $cluster'\t'$abundance'\t'$cluster_avg_len'\t'$score'\t'$read_count'\t'$total_bases'\t'$coverage'\t'$subsampled'\t'$coverage_final'\t'$notes

if (( $(echo "$subsampled" | bc -l) )); then
    rm $r1
    rm $r2
fi

rm $clu_sketch
rm $sorted_ref_info
