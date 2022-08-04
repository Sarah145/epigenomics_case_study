#! /bin/bash

celltype_dir=./cluster_beds
cell_files=($(ls ${celltype_dir}))

for i in ${cell_files[@]}
        do
                cd $celltype_dir
                sample=$(echo "${i%%.*}")
                echo "Working on ${sample}"
                out_dir=cluster/${sample}
                findMotifsGenome.pl $i hg38 $out_dir -size given -S 10
        done
