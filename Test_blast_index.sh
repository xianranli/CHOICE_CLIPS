#!/bin/bash


for file_bsz in B97	CML103	CML228	CML247	CML277	CML322	CML333	CML52	CML69	HP301	Il14H	Ki11	Ki3	Ky21	M162W	M37W	Mo18W	Ms71	NC350	NC358	Oh43	Oh7B	P39	Tx303	Tzi8


do

echo "processing $file_bsz ..."

cd $file_bsz


for m in chr2	chr3	chr4	chr5	chr6	chr7	chr8	chr9	chr10

do

rm "$file_bsz"_"$m".fa.gz.fai "$file_bsz"_"$m".fa.gz.gzi "$file_bsz"_"$m".fa.n*

bgzip -d "$file_bsz"_"$m".fa.gz

makeblastdb -in "$file_bsz"_"$m".fa -dbtype nucl

bgzip "$file_bsz"_"$m".fa

samtools faidx "$file_bsz"_"$m".fa.gz

done


cd ..

done





