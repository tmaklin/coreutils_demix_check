#!/bin/bash
## Input: $1: abundances file from mSWEEP
##        $2: number of threads
##        $3: tmp directory
##        $4: buffer size for `sort`
##        $5: demix_check index folder
##        $6: forward strand
##        $7: reverse strand

##set -euxo pipefail

export LC_ALL=C
export TMPDIR=$3

abun_in=$1
nthreads=$2
tmpdir=$3
bufsize=$4
refdir=$5
fwd=$6
rev=$7

cluster=$(echo $fwd | sed 's/_1[.]fastq[.]gz//g')
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
mash sketch -p $nthreads -s 10000 -r -m $m -I $cluster -C - -o $clu_sketch $r1 $r2 2> /dev/null

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
