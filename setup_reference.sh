#!/bin/bash
## Input: $1: ref_info file (format: tab separated columns ID, cluster, assembly_path).
##        $2: number of threads

##set -euxo pipefail

paths=tmp_paths_$RANDOM".txt"
cut -f3 $1 | sed '1d' > $paths
mash sketch -p $2 -s 10000 -o ref -l $paths 2> /dev/null
rm $paths
mash dist -p $2 ref.msh ref.msh | pigz -11 -p $2 > ref_msh_dis.tsv.gz

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

ref_path=tmp_1_$RANDOM".txt"
met_path=tmp_2_$RANDOM".txt"
distance=tmp_3_$RANDOM".txt"
prob=tmp_4_$RANDOM".txt"
hashes=tmp_5_$RANDOM".txt"
ss=tmp_6_$RANDOM".txt"

zcat --force ref_msh_dis.tsv.gz | cut -f1 > $ref_path
zcat --force ref_msh_dis.tsv.gz | cut -f2 > $met_path
zcat --force ref_msh_dis.tsv.gz | cut -f3 > $distance
zcat --force ref_msh_dis.tsv.gz | cut -f4 > $prob
zcat --force ref_msh_dis.tsv.gz | cut -f5 | cut -f1 -d'/' > $hashes
zcat --force ref_msh_dis.tsv.gz | cut -f5 | cut -f2 -d'/' > $ss

ref_id=tmp_7_$RANDOM".txt"
met_id=tmp_8_$RANDOM".txt"

ref_cluster=tmp_9_$RANDOM".txt"

while read line; do
    info=$(grep $line $1 | head -n1)
    id=$(echo $info | cut -f1 -d' ')
    cluster=$(echo $info | cut -f2 -d' ')
    echo $id >> $ref_id
    echo $cluster >> $ref_cluster
done < $ref_path

met_cluster=tmp_10_$RANDOM".txt"

while read line; do
    info=$(grep $line $1 | head -n1)
    id=$(echo $info | cut -f1 -d' ')
    cluster=$(echo $info | cut -f2 -d' ')
    echo $id >> $met_id
    echo $cluster >> $met_cluster
done < $met_path

category=tmp_11_$RANDOM".txt"

while IFS= read -r ref_cl && IFS= read -r met_cl <&3; do
    if [[ "$ref_cl" == "$met_cl" ]]; then
	echo "same cluster"
    else
	echo "different cluster"
    fi
done < $ref_cluster 3< $met_cluster > $category

echo -e "ref_id\tmet_id\tdistance\thashes\tss\tp\tref_cluster\tmet_cluster\tcategory" > ref_msh_dis_clu.tsv
paste $ref_id $met_id $distance $hashes $ss $prob $ref_cluster $met_cluster $category >> ref_msh_dis_clu.tsv

Rscript get_thresholds.R

pigz --force -11 -p $2 ref_msh_dis_clu.tsv

rm $ref_path
rm $met_path
rm $distance
rm $prob
rm $hashes
rm $ss
rm $ref_id
rm $met_id
rm $ref_cluster
rm $met_cluster
rm $category
