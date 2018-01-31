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

		

	if not args.samples:
		fastqFiles = [ x for x in os.listdir(args.FastqDirectory) if 'fastq' in x.split('.')[-2:] ]
	else:
		fastqFiles = common.importSampleList(args.samples)		
	fastqFiles = [args.FastqDirectory + x for x in fastqFiles]	
	
	
	
  

	#run multiprocessing of all mapping commands#
	argList = [(x, args.species, args.trim, statsDir, tempDir, samDir) for x in fastqFiles]
		
	common.daemon(mapfile.runOne, argList, 'map fastq files', cpuPerProcess=8)

	

	
	
	#remove all temporary files#
	shutil.rmtree(tempDir[:-1])
  
  
	
  	print('Mapping complete\n\n\n')

	
	
	
	
	
	
	
	
	
