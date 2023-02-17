#!/bin/bash
## Input: $1: ref_info file (format: tab separated columns ID, cluster, assembly_path).
##        $2: number of threads

##set -euxo pipefail

paths=ref_paths.txt
cut -f3 $1 | sed '1d' > $paths
mash sketch -p $2 -s 10000 -o ref -l $paths 2> /dev/null
mash dist -p $2 ref.msh ref.msh | pigz -p $2 > ref_msh_dis.tsv.gz

touch ref_comp.tsv
touch ref_clu.txt
touch ref_clu.tsv
echo -e "chr\tcluster\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > ref_comp.tsv
echo -e "id\tcluster" > ref_clu.tsv
while read line; do
    id=$(echo $line | cut -f1,2 -d' ')
    path=$(echo $line | cut -f3 -d' ')
    res=$(seqtk comp $path | cut -f2-13 | awk -F "[[:space:]]" '{{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}}END{{print}}')
    echo $id' '$res | tr ' ' '\t'
    echo $id | cut -f2 -d' ' >> ref_clu.txt
    echo $id | tr ' ' '\t' >> ref_clu.tsv
done < <(sed '1d' $1) >> ref_comp.tsv

touch ref_clu_comp.tsv
echo -e "cluster\tn\tlength_ave\tlength_min\tlength_max" > ref_clu_comp.tsv
sed '1d' ref_comp.tsv | datamash --sort --group 2 count 2 mean 3 min 3 max 3 >> ref_clu_comp.tsv

## Sort the ref info for `join`
sort ref_info.tsv > ref_info.sorted.tsv

echo -e "ref_id\tmet_id\tdistance\thashes\tss\tp\tref_cluster\tmet_cluster\tcategory" > ref_msh_dis_clu.tsv
zcat --force ref_msh_dis.tsv.gz \
    | awk -F"[\t]" '$1!=$2 { print $0 }' \
    | sort -k1 \
    | join -1 1 -2 3 - ref_info.sorted.tsv \
    | tr ' ' '\t' | sort -k2 \
    | join -1 2 -2 3 - ref_info.sorted.tsv \
    | tr ' ' '\t' \
    | sed 's/\([0-9][0-9]*\)[/]\([0-9][0-9]*\)/\1\t\2/g' \
    | awk '{ print $9 "\t" $7 "\t" $3 "\t" $5 "\t" $6 "\t" $4 "\t" $10 "\t" $8}' \
    | awk -F"[\t]" '{print $0, "\t" ($7==$8?"same":"different") "_cluster"}' \
    | tr -d ' ' \
    | sed 's/_cluster/ cluster/g' >> ref_msh_dis_clu.tsv


cat ref_msh_dis_clu.tsv | sed '1d' \
    | awk -F"[\t]" '$7==$8 { print $0 }' \
    | datamash --sort --group 7 median 3 min 3 max 3 > within_dis.tsv

cat ref_msh_dis_clu.tsv | sed '1d' \
    | awk -F"[\t]" '$7!=$8 { print $0 }' \
    | datamash --sort --group 7 median 3 min 3 max 3\
    | join -a 1 -o 1.1,1.2,1.3,1.4,2.2,2.3,2.4 -e 0 - within_dis.tsv | tr ' ' '\t' | sort -k1 \
    | awk -v SF="0.20" '{printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$7*(1 + SF))"\n" }' > tmp_clu_thr.tsv

t_min=$(cut -f2 tmp_clu_thr.tsv | awk -v SF="0.20" '{printf($1*SF)"\n" }' | datamash --sort median 1)
t_med=$(cut -f8 tmp_clu_thr.tsv | grep -v "^0$" | datamash --sort median 1)
if (( $(echo "$t_med < $t_min" | bc -l) )); then
    t_comp=$t_min
else
    t_comp=$t_med
fi

echo -e "cluster\tn\tthreshold\tdis_same_max\tdis_same_med_all\tdis_diff_med_all" > ref_clu_thr.tsv
cat tmp_clu_thr.tsv \
    | paste - <(cut -f8 tmp_clu_thr.tsv | awk -v SF="$t_comp" '{print ($1 < SF ? SF:$1) }') \
    | awk '{ print $1 "\t" $9 "\t" $7}' \
    | paste - <(same_med_all=$(cut -f5 tmp_clu_thr.tsv | grep -v "^0$" | datamash --sort median 1); yes "$same_med_all" | head -n $(wc -l tmp_clu_thr.tsv | cut -f1 -d' ')) <(diff_med_all=$(cut -f2 tmp_clu_thr.tsv | grep -v "^0$" | datamash --sort median 1); yes "$diff_med_all" | head -n $(wc -l tmp_clu_thr.tsv | cut -f1 -d' ')) \
    | join -1 1 -2 1 ref_clu_comp.tsv - \
    | tr ' ' '\t' \
    | cut -f1,2,6-9 >> ref_clu_thr.tsv

rm tmp_clu_thr.tsv

pigz --force -p $2 ref_msh_dis_clu.tsv
