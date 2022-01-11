#!/bin/bash

IN=/share/afsc/snowcrab-OA-2022/raw-data/concat

FILES=$(ls ${IN}/*.fastq | awk -F "/" '{print $NF}')

for file in ${FILES}
do
    echo ${file}
    gzip -c ${IN}/${file} > ${IN}/${file}.gz
done
