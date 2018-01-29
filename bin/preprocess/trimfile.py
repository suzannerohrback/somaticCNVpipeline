#!/usr/bin/python
import os
import gzip










#open input and output files#
def openFiles(file):
	
	if file[-5:] == 'fastq':
		IN = open(file, 'r')
		outfilename = file[:-5] + 'k' + str(length) + '.fastq'
		
	elif file[-8:] == 'fastq.gz':
		IN = gzip.open(file, 'rb')
		outfilename = file[:-8] + 'k' + str(length) + '.fastq'
		
	else:
		errorText = '\nERROR: unrecognized file type to be processing, please ensure only fastq files are being submitted to this pipeline'
		errorText += '\n\t' + file + '\n'
		print(errorText)
		raise SystemExit
	
	OUT = open(outfilename, 'w')
	
	return IN, OUT

	

	
	
	
	
	
	
	
#do the trimming (and print out changes)#
def trimOne(IN, OUT, trim, length):
	lineCount = 0
	includeCount = 0
	includeTest = True
	for x in IN:
		lineCount += 1
		
		#name line, save information#
		if (lineCount + 3) % 4 == 0:
			lines = x
			
		#sequence line, test trimming ability#
		elif (lineCount + 2) % 4 == 0:
			aline = x.rstrip()
			if len(aline) >= trim + length:
				includeTest = True
				newLine = aline[trim:(trim + length)]
				lines += newLine + '\n'
			else:
				includeTest = False
				
		#name line 2, save information#				
		elif (lineCount + 1) % 4 == 0:
			lines += x
			
		#quality line, write data to file if read is long enough#
		elif (lineCount) % 4 == 0:
			if includeTest:
				newLine = x[trim:(trim + length)]
				lines += newLine + '\n'
				OUT.write(lines)
				includeCount += 1
			includeTest = True
			
	IN.close()
	OUT.close()
	
	printText = '\t\tMaintained ' + str(includeCount) + ' of ' + str(lineCount) + ' fastq reads for ' + file.split('/')[-1] 
	print(printText)
	
	
	
	
	
	
	
	
	
	
def preprocessOne(file, trim, length, remove=False):
	IN, OUT = openFiles(file)
	trimOne(IN, OUT, trim, length)
	
	#remove or move full length fastq#
	if remove:
		os.remove(file)
	else:
		newLoc = os.path.dirname(file) + '/FullLength/' + os.path.basename(file)
		os.rename(file, newLoc)
	
	
	
	
	
