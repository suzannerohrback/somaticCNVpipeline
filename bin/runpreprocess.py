#!/usr/bin/python
import os

from preprocess import trimfile
from common import daemon










def runAll(args):
	
	print('\n\n\nYou have requested preprocess (trim) fastq files')
	print('\tWARNING:')
	print('\t\tIF USING ANY LENGTH OTHER THAN 36 BP, REFERENCE FILES ARE NOT SUPPORTED FOR DOWNSTREAM PROCESSING')
	print('\n')
	

	#ensure fastq directory is in the expected format#
	if args.FastqDirectory[-1] != '/':
		args.FastqDirectory += '/'
	
	
	
	#get list of fastq files to process (depending on args.samples)
	if os.path.exists(args.samples):
		fastqFiles = []
		with open(args.samples, 'r') as IN:
			for x in IN:
				fastqFiles.append(x.rstrip())

	else:
		fastqFiles = [ x for x in os.listdir(args.FastqDirectory) if 'fastq' in x.split('.')[-2:] ]
	
	fastqFiles = [args.FastqDirectory + x for x in fastqFiles]
	
	
	
	#use the daemon to and preprocessing code to trim all fastq files with parallel processing
	if args.remove and not os.path.exists(args.FastqDirectory + '/FullLength/'):
		os.mkdir(args.FastqDirectory + '/FullLength/')
		
	if args.remove:
		argList = [(x, args.trim5, args.length, remove=True,) for x in fastqFiles]
	else:
		argList = [(x, args.trim5, args.length,) for x in fastqFiles]
		
	daemon(trimfile.preprocessOne, argList, 'trim sequencing reads to desired length', cpuPerProcess=1)
		
	
	
	print('Pre-processing complete\n\n\n')

