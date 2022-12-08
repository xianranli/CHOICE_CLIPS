#!/bin/bash


for file_bsz in jagger	julius	lancer	landmark	mace	mattis	norin61	stanley	IWGSC	spelta	

do

echo "processing $file_bsz ..."

cd $file_bsz


for m in chr1A	chr1B	chr1D	chr2A	chr2B	chr2D	chr3A	chr3B	chr3D	chr4A	chr4B	chr4D	chr5A	chr5B	chr5D	chr6A	chr6B	chr6D	chr7A	chr7B	chr7D

do

rm "$file_bsz"_"$m".fa.gz.fai "$file_bsz"_"$m".fa.gz.gzi "$file_bsz"_"$m".fa.n*

bgzip -d "$file_bsz"_"$m".fa.gz

makeblastdb -in "$file_bsz"_"$m".fa -dbtype nucl

bgzip "$file_bsz"_"$m".fa

samtools faidx "$file_bsz"_"$m".fa.gz

done


cd ..

done





