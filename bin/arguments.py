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
	parser.add_arugment('-r', '--remove', action='store_false'
		help = 'Set this flag if you want to delete the full length fastq files')










def fullParser(input):

	functionDict =	{ 
			'-h': parentArgs,
			'--help': parentArgs,
			'preprocess': preprocessArgs, 
			}

	if input == []:
		parentArgs()

	parser = functionDict[input[0]]()
	args = parser.parse_args(input[1:])

	return input[0], args
