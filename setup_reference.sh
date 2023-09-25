#!/bin/bash
## Input:
##        --ref_info: ref_info file (format: tab separated columns ID, cluster, assembly_path).
##        --threads: number of threads      (default: 1)
##        --tmpdir: tmp directory           (default: working directory)
##        --bufsize: buffer size for `sort` (default: 4000M)
##        --verbose: echo commands run      (default: silent)
##        --help: print this message

set -euo pipefail

export LOCALE=C
export LANG=C

arg_parser() {
    ## Default values
    ref_info=""   ## Exit if not supplied
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
		echo -e "setup_reference.sh:\n\t--ref_info: ref_info file         (format: tab separated columns ID, cluster, assembly_path).\n\t--threads: number of threads      (default: 1)\n\t--tmpdir: tmp directory           (default: working directory)\n\t--bufsize: buffer size for sort   (default: 4000M)\n\t--verbose: echo commands run      (default: silent)\n\t--help: print this message"
		exit 0
	    else
		param="${1/--/}"
		declare -g $param="$2"
	    fi
        fi
        shift
    done

    if [ "$ref_info" == "" ]; then
	echo "--ref_info was not supplied!"
	exit 1
    fi
}

run_seqtk() {
    line=$1
    id=$(echo $line | cut -f1,2 -d' ');
    path=$(echo $line | cut -f3 -d' ');
    res=$(seqtk comp $path | cut -f2-13 | awk -F "[[:space:]]" '{{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}}END{{print}}');
    echo $id' '$res | tr ' ' '\t';
}
export -f run_seqtk

arg_parser $@

if [ "$verbose" == "true" ]; then
    set -x
fi

if [ "$verbose" == "true" ]; then
    set -x
fi

if command -v pigz > /dev/null; then
    compress_cmd="pigz -p $threads"
else
    compress_cmd="gzip"
fi

export LC_ALL=C
export TMPDIR=$tmpdir

tmpdir=$tmpdir

paths=ref_paths.txt
cut -f3 $ref_info | sed '1d' > $paths
mash sketch -p $threads -s 10000 -o ref -l $paths 2> $tmpdir/mash_setup.log
mash dist -p $threads ref.msh ref.msh | $compress_cmd > ref_msh_dis.tsv.gz

touch ref_comp.tsv
touch ref_clu.txt
touch ref_clu.tsv
echo -e "chr\tcluster\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > ref_comp.tsv
echo -e "id\tcluster" > ref_clu.tsv
cut -f1,2 $ref_info >> ref_clu.tsv
cut -f2 $ref_info | sed '1d' > ref_clu.txt

if command -v parallel > /dev/null; then
    parallel --will-cite -j $threads 'run_seqtk {}' < <(sed '1d' $ref_info) >> ref_comp.tsv
else
    while read line; do
	run_seqtk "$line"
    done < <(sed '1d' $ref_info)
fi

touch ref_clu_comp.tsv
echo -e "cluster\tn\tlength_ave\tlength_min\tlength_max" > ref_clu_comp.tsv
sed '1d' ref_comp.tsv | datamash --sort --group 2 count 2 mean 3 min 3 max 3 >> ref_clu_comp.tsv

## Sort the ref info for `join`
sed '1d' $ref_info | sort --parallel=$threads -S $bufsize -T $tmpdir -k3 > ref_info.sorted.tsv

echo -e "ref_id\tmet_id\tdistance\thashes\tss\tp\tref_cluster\tmet_cluster\tcategory" > ref_msh_dis_clu.tsv
gunzip -c --force ref_msh_dis.tsv.gz \
    | awk -F"[\t]" '$1!=$2 { print $0 }' \
    | sort --parallel=$threads -S $bufsize -T $tmpdir -k1 \
    | join -1 1 -2 3 - ref_info.sorted.tsv \
    | tr ' ' '\t' | sort --parallel=$threads -S $bufsize -T $tmpdir -k2 \
    | join -1 2 -2 3 - ref_info.sorted.tsv \
    | tr ' ' '\t' \
    | sed 's/\([0-9][0-9]*\)[/]\([0-9][0-9]*\)/\1\t\2/g' \
    | awk '{ print $9 "\t" $7 "\t" $3 "\t" $5 "\t" $6 "\t" $4 "\t" $10 "\t" $8}' \
    | awk -F"[\t]" '{print $0, "\t" ($7==$8?"same":"different") "_cluster"}' \
    | tr -d ' ' \
    | sed 's/_cluster/ cluster/g' >> ref_msh_dis_clu.tsv

within_dis=$tmpdir/tmp_within_dis-$RANDOM".tsv"

cat ref_msh_dis_clu.tsv | sed '1d' \
    | awk -F"[\t]" '$7==$8 { print $0 }' \
    | datamash --sort --group 7 median 3 min 3 max 3 > $within_dis

tmp_thr=$tmpdir/tmp_clu_thr-$RANDOM".tsv"

cat ref_msh_dis_clu.tsv | sed '1d' \
    | awk -F"[\t]" '$7!=$8 { print $0 }' \
    | datamash --sort --group 7 median 3 min 3 max 3\
    | join -a 1 -o 1.1,1.2,1.3,1.4,2.2,2.3,2.4 -e 0 - $within_dis | tr ' ' '\t' | sort --parallel=$threads -S $bufsize -T $tmpdir -k1 \
    | awk -v SF="0.20" '{printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$7*(1 + SF))"\n" }' > $tmp_thr

t_min=$(cut -f2 $tmp_thr | awk -v SF="0.20" '{printf($1*SF)"\n" }' | datamash --sort median 1)
t_med=$(cut -f8 $tmp_thr | grep -v "^0$" | datamash --sort median 1)
if (( $(echo "$t_med < $t_min" | bc -l) )); then
    t_comp=$t_min
else
    t_comp=$t_med
fi

sorted_clu_comp=$tmpdir/sorted_clu_comp.tsv
sed '1d' ref_clu_comp.tsv | sort --parallel=$threads -S $bufsize -T $tmpdir > $sorted_clu_comp
echo -e "cluster\tn\tthreshold\tdis_same_max\tdis_same_med_all\tdis_diff_med_all" > ref_clu_thr.tsv
cat $tmp_thr \
    | paste - <(cut -f8 $tmp_thr | awk -v SF="$t_comp" '{print ($1 < SF ? SF:$1) }') \
    | awk '{ print $1 "\t" $9 "\t" $7}' \
    | paste - <(same_med_all=$(cut -f5 $tmp_thr | grep -v "^0$" | datamash --sort median 1); yes "$same_med_all" | head -n $(wc -l $tmp_thr | cut -f1 -d' ')) <(diff_med_all=$(cut -f2 $tmp_thr | grep -v "^0$" | datamash --sort median 1); yes "$diff_med_all" | head -n $(wc -l $tmp_thr | cut -f1 -d' ')) \
    | join -1 1 -2 1 $sorted_clu_comp - \
    | tr ' ' '\t' \
    | cut -f1,2,6-9 >> ref_clu_thr.tsv

rm $sorted_clu_comp
rm $tmp_thr

$compress_cmd --force ref_msh_dis_clu.tsv
