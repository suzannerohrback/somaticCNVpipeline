#!/usr/bin/python















def runAll(args):

	print('\n\n\nYou have requested to normalize and segment bincounts files')
	print('\tWARNING:')
	print('\t\tIF USING ANY REFERENCES OTHER THAN THOSE I PROVIDE I CANNOT GUARANTEE RESULT ACCURACY')
	print('\n')
 
	printText = '\nTHIS FUNCTION IS NOT YET COMPLETE, SORRY\n'
	print(printText)
	raise SystemExit
	
	
	#Set up environment#
	#make folders, get list of sample names, etc
	
	
	
	
	
	#Run normalization#
	#running on a subset of samples? gc only or method too? 
	#if method, make the method reference
	#use multiprocessing daemon to run for each sample in parallel
	
	
	
	print('\nNormalization complete\n\n\n')
	
	if args.normalizeonly:
		return 0
	


	
	#Run segmentation#
	#write matlab script
	#run matlab script
	#run for all samples in parallel with multiprocessing daemon
	
	print('\nSegmentation complete\n\n\n')

	
	
	
	
#	parser.add_argument('CountDirectory', 
#		help = 'The path to the folder that contains bincounts.txt files to be processed')
#	parser.add_argument('species', choices=['hg38', 'mm10'], 
#		help = 'The genome build of the species being assessed')
#	parser.add_arugment('-o', '--output', metavar='/path/to/output_directory/', default=False,
#		help = 'A filepath to the desired directory where you would like lowess.bincounts.txt and segments.txt files saved, if not in the same parent directory as the bincounts files')
#	parser.add_argument('-i', '--infofile', metavar='/path/to/sample.info.txt', default=False,
#		help='Path to a .txt file containing information about the samples to be processed (unique name, amplification method, number of cells)\n\tIf not all are identical')
#	parser.add_argument('-c', '--columns', metavar='X X X', default=[0, 1, 2], type=int, nargs=3,
#		help='The zero-indexed locations of the columns to import from the infofile in the order: name, method, cell number (if not the first 3 columns)')
#	parser.add_arugment('-g', '--gconly', action='store_true'
#		help = 'Set this flag if you only want GC-correction to be performed during normalization')
#	parser.add_argument('-s', '--samples', metavar='/path/to/sample_list.txt', default=False,
#		help='Path to a file containing a list of unique.sam files to be processed\n\tsample names only, no path or file extension needed')

	
