#!/usr/bin/python
import os
import shutil

from map import mapfile
import common










def runAll(args):
	
	print('\n\n\nYou have requested to map fastq files')
	print('\tWARNING:')
	print('\t\tPLEASE MAKE SURE YOU ARE USING')
	print('\t\t\tBowtie v1 and Samtools v0.1.19')
	print('\t\t\tBowtie v1 mapping indexes for either mm10 or hg38')
	print('\n')

	
	
	#set up environment#
	args.FastqDirectory = common.fixDirName(args.FastqDirectory)
	
	samDir = os.path.dirname(args.FastqDirectory[:-1]) + '/Sam/'
	if args.output:
		samDir = common.fixDirName(args.output)
	
	statsDir = os.path.dirname(args.FastqDirectory[:-1]) + '/PipelineStats/'
	if args.statdir:
		statsDir = common.fixDirName(args.statdir)

	tempDir = args.FastqDirectory + 'Temp/'

	for i in [samDir, tempDir, statsDir]:
		common.makeDir(i)		

	fastqFiles = common.getSampleList(args.FastqDirectory, args.samples, 'fastq')
	
	

	#run multiprocessing of all mapping commands#
	argList = [(x, args.MapIndex, args.trim, statsDir, tempDir, samDir, args.bowtie, args.samtools) for x in fastqFiles]
	common.daemon(mapfile.runOne, argList, 'map fastq files', cpuPerProcess=8)


	
	#remove all temporary files#
	shutil.rmtree(tempDir[:-1])
  
  
	
  	print('\nMapping complete\n\n\n')

	
	
	
	
	
	
	
	
	
