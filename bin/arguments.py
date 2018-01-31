#!/usr/bin/python
import argparse










###basic parser for parent help statement###
def parentArgs():

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''\
	Suzanne's pipeline to identify somatic CNVs from single-cell whole-genome sequencing data
	=========================================================================================

	You must specify a function to perform:
	  *preprocess (trim fastq reads to the appropriate length)
	  *[More functions coming soon]
		''')

	parser.print_help()

	raise SystemExit









###interpret arguments needed to perform preprocessing of fastq files###
def preprocessArgs():

	parser = argparse.ArgumentParser(description='Trim fastq reads to the appropriate length')

	#required arguments#
	parser.add_argument('FastqDirectory', 
		help = 'The path to the folder that contains fastq files to be processed')
	
	#optional arguments#
	parser.add_arugment('-5', '--trim5', metavar='X', type=int, default=0,
		help = "Number of 5' bases to trim from fastq reads")
	parser.add_arugment('-l', '--length', metavar='X', type=int, default=36,
		help = 'The desired read length')
	parser.add_arugment('-r', '--remove', action='store_true'
		help = 'Set this flag if you want to delete the full length fastq files')
	parser.add_argument('-s', '--samples', metavar='/path/to/sample_list.txt', default=False,
		help='Path to a text file containing a list of fastq files to be processed')










###interpret arguments needed to perform mapping of fastq files###
def mapArgs():
	
	parser = argparse.ArgumentParser(description='Map fastq files to the appropriate reference genome')

	#required arguments#
	parser.add_argument('FastqDirectory', 
		help = 'The path to the folder that contains fastq files to be processed')
	parser.add_argument('species', choices=['hg38', 'mm10'], 
		help = 'The genome build of the species being assessed')
	
	#optional arguments#
	parser.add_arugment('-t', '--trim', metavar='X', nargs=2, type=int, default=[0, 0],
		help = "Number of 5' and 3' bases to trim from fastq reads during mapping")
	parser.add_arugment('-o', '--output', metavar='/path/to/output_directory/', default=False,
		help = 'A filepath to the desired directory where you would like sam files saved, if not in the same parent directory as the fastq files')
	parser.add_argument('-s', '--samples', metavar='/path/to/sample_list.txt', default=False,
		help='Path to a text file containing a list of fastq files to be processed')
	
	

	





###interpret arguments needed to perform counting of unique.sam files###
def countArgs():
	
	parser = argparse.ArgumentParser(description='Count the reads per genomic bin from unique sam files')

	#required arguments#
	parser.add_argument('SamDirectory', 
		help = 'The path to the folder that contains unique.sam files to be processed')
	parser.add_argument('species', choices=['hg38', 'mm10'], 
		help = 'The genome build of the species being assessed')
	
	#optional arguments#
	parser.add_arugment('-o', '--output', metavar='/path/to/output_directory/', default=False,
		help = 'A filepath to the desired directory where you would like bincount.txt files saved, if not in the same parent directory as the sam files')
	parser.add_argument('-s', '--samples', metavar='/path/to/sample_list.txt', default=False,
		help='Path to a text file containing a list of unique.sam files to be processed')

	
	
	
	
	
	
	
	
	
def fullParser(input):

	functionDict =	{ 
			'-h': parentArgs,
			'--help': parentArgs,
			'preprocess': preprocessArgs, 
			'map': mapArgs,
			'count': countArgs,
			}

	
	
	if input == []:
		parentArgs()
		
	if input not in functionDict.keys():
		return input[0], False

	
	
	parser = functionDict[input[0]]()
	args = parser.parse_args(input[1:])

	return input[0], args









