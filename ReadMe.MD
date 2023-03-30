#### Introduction
_CHOICE_ and _CLIPS_ are two unsupervised learning mechanisms used in [BRIDGEcereal](https://bridgecereal.scinet.usda.gov/), which is developed to survey and graph large indel based haplotype graph from publicly accessable pan-genomes for a gene of interest .

A large indel, if presents, is indicated by being flanked by two High-scoring Segment Pairs (HSPs) for a pair of segments identified by [Blastn](https://www.ncbi.nlm.nih.gov/books/NBK279690/). A Pan-genome with mulitple assemblies enables the possibility to identify large indels from a gene of interest through the desired All-vs-All comparison.

However, two main challenges are:
1. Identification of the segment containing the ortholog from each assembly.
2. Presentation of identified indels from a large number of assembiles.

_CHOICE_ and _CLIPS_ are developed to tackle these two challenges, respectively.

#### Algorithm illustration
_CHOICE_, Clustering HSPs for Ortholog Identification via Coordinates and Equivalence, identifies and extracts the segment harboring the ortholog from each assembly.
![CHOICES](./figs/CHOICE.png)

CLIPS, Clustering via Large-Indel Permuted Slopes, groups these segments sharing the same set of indels into haplotypes. 
![CLIPS](./figs/CLIPS.png)


#### Demonstration
(R libaray *dendextend* and *data.table* are required to go through the demonstration)

TraesCS4A02G058900 in chromosome 4A of **Chinese Spring** is used as the demonstration. 

CDS of TraesCS4A02G058900 is used to blast agains chr4A of each wheat assembly (11 in total) via Blastn with -outfmt 6. The blastn output is reorgzined as [TraesCS4A02G058900_syn.txt](./demo/TraesCS4A02G058900_syn.txt), which is the input for _CHOICE_. _CHOICE_ which will determine the region of each assembly contain the ortholog. For example, the coordinates for variety **spelta** are from 52,231,267 to 52,245,777, while are from 52,859,043 to 52,873,560 for variety **lancer**. The coordinates will be then used to exact the corresponding segment out from each assembile. Two output files (CHOICE_summary.txt and Filtered_By_CHOICE.txt) will be also saved.

These segments are then blasted against each other to identify HSP via Blastn (-outfmt 6). The output file [TraesCS4A02G058900_seg-Self_out_m6.txt](./demo/TraesCS4A02G058900_seg-Self_out_m6.txt) is then fed into _CLIPS_. The clusting tree suggests two main haplotypes for this genes, one is represented by **Chinese Spring**, while all other varieties share another haplotype. One out file "slope_matrix_output_from_CLIPS.txt" will be also saved. 

#### Citation
Zhang B, Huang H, Tibbs-Cortes LE, Vanous A, Zhang Z, Sanguinet K, Garland-Campbell KA, Yu J, Li X. [Streamline unsupervised machine learning to survey and graph indel-based haplotypes from pan-genomes](https://doi.org/10.1101/2023.02.11.527743). bioRxiv. 2023. doi: 10.1101/2023.02.11.527743