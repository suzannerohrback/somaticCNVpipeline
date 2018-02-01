#!/usr/bin/python
import os
import shutil

from map import mapfile
import common










def runAll(args):
	
	print('\n\n\nYou have requested to map fastq files')
	print('\tWARNING:')
	print('\t\tIF USING ANY REFERENCES OTHER THAN THOSE I PROVIDE I CANNOT GUARANTEE RESULT ACCURACY')
	print('\n')

	
	
	#set up environment#
	args.FastqDirectory = common.fixDirName(args.FastqDirectory)
	
	if not args.output:
		samDir = os.path.dirname(args.FastqDirectory) + '/' + Sam + '/'
	else:
		samDir = common.fixDirName(args.output)
	
	tempDir = args.FastqDirectory + 'Temp/'
	statsDir = args.FastqDirectory + 'MapStats/'
	for i in [samDir, tempDir, statsDir]:
		common.makeDir(i)		

	fastqFiles = common.getSampleList(args.FastqDirectory, args.samples, 'fastq')
	
	

	#run multiprocessing of all mapping commands#
	argList = [(x, args.species, args.trim, statsDir, tempDir, samDir) for x in fastqFiles]
		
	common.daemon(mapfile.runOne, argList, 'map fastq files', cpuPerProcess=8)


	
	#remove all temporary files#
	shutil.rmtree(tempDir[:-1])
  
  
	
  	print('\nMapping complete\n\n\n')

	
	
	
	
	
	
	
	
	
