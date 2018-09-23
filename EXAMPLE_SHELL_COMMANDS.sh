#!/usr/bin/bash


#CHANGE THE FOLLOWING FILEPATHS TO MATCH YOUR OWN DIRECTORY STRUCTURE
#A MINIMUM NUMBER OF OPTIONS ARE BEING SPECIFIED FOR GREATER SIMPLICITY

analysisDir=/oasis/tscc/scratch/srohrbac/research/pipelinetest2/
fastqDir=/oasis/tscc/scratch/srohrbac/research/pipelinetest2/Fastq/
species=mm10
mapidx=/home/srohrbac/research/reference_files/mappingRef/mm10_index/mm10.index
pipelineDir=/home/srohrbac/research/scripts/GITtmp/somaticCNVpipeline
samtools=/home/k4zhang/softwares/samtools-0.1.19/samtools


cd $pipelineDir

python somaticCNVcalling preprocess $fastqDir -5 30 -l 36
python somaticCNVcalling map -m $samtools $fastqDir $mapidx
python somaticCNVcalling count $analysisDir $species
python somaticCNVcalling segment $analysisDir $species
python somaticCNVcalling interpret $analysisDir $species


