#!/usr/bin/bash


#CHANGE THE FOLLOWING FILEPATHS TO MATCH YOUR OWN DIRECTORY STRUCTURE
#A MINIMUM NUMBER OF OPTIONS ARE BEING SPECIFIED FOR GREATER SIMPLICITY

analysisDir=/oasis/tscc/scratch/srohrbac/research/pipelinetest2/
fastqDir=/oasis/tscc/scratch/srohrbac/research/pipelinetest2/Fastq/
species=mm10
mapidx=/home/srohrbac/research/reference_files/mappingRef/mm10_index/mm10.index
pipelineDir=/home/srohrbac/research/scripts/GITtmp/somaticCNVpipeline



cd $pipelineDir

python somaticCNVcalling preprocess $fastqDir -5 30 -l 36
python somaticCNVcalling map $fastqDir $mapidx
python somaticCNVcalling count $analysisDir $species
python somaticCNVcalling $analysisDir $species
python somaticCNVcalling $analysisDir $species


