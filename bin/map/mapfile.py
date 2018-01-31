#!/usr/bin/python
import sys
import os
import inspect
import subprocess as sub
import shlex

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import common
import config as cfg
	
	
	
	
	
	
	
	
	
	
def runCommand(cmd, outfile=False):
	cmd = shlex.split(cmd)
	
	if outfile:
		if os.path.exists(outfile):
			stdout = open(outfile, 'a')
			stdout.flush()
		else:
			stdout = open(outfile, 'w')
		p = sub.Popen(cmd, stdout=stdout, stderr=sub.STDOUT)
		
	else:
		p = sub.Popen(cmd)
		
	p.wait()
	
	return 0










def runOne(fastqFile, species, trim, statsDir, tempDir, samDir):
  
	#get environment prepared#
	mapVars = cfg.Map()

 	fastqFile = common.zipping(fastqFile)
	
	sampleName = os.path.basename(fastqFile)[:-6]
  
	statFile = statsDir + sampleName + '.map.stats.txt'

  
  
	#run bowtie#
	cmd =	mapVars.bowtieOptions +	\
			[
			'-5', str(trim[0]), 
			'-3', str(trim[1]), 
			mapVars.indexDict[species], 
			fastqFile, 
			tempDir + sampleName + '.sam'
			]
	cmd = ' '.join(cmd)
	runCommand(cmd, outfile=statFile)
	
	
	
	#sam to bam#
	cmd =	[
			mapVars.samtools, 'view', '-bS', 
		  	'-o', tempDir + sampleName + '.bam',
		 	tempDir + sampleName + '.sam'
			]
	cmd = ' '.join(cmd)
	runCommand(cmd)
	

	
	#sort bam#
	cmd =	[
			mapVars.samtools, 'sort', 
			'-T', tempDir + sampleName + '.sorted',
			'-o', tempDir + sampleName + '.sorted.bam',
		 	tempDir + sampleName + '.bam'
			]
	cmd = ' '.join(cmd)
	runCommand(cmd)
	
	
	
	#remove duplicates#
	cmd =	[
			mapVars.samtools, 'rmdup', '-s', 
			tempDir + sampleName + '.sorted.bam',
			tempDir + sampleName + '.unique.bam'
			]
	cmd = ' '.join(cmd)
	runCommand(cmd, outfile=statFile)

	
	
	#bam to sam#
	cmd =	[
			mapVars.samtools, 'view', '-h', 
		  	'-o', samDir + sampleName + '.unique.sam',
		 	tempDir + sampleName + '.unique.bam'
			]
	cmd = ' '.join(cmd)
	runCommand(cmd)


	
	printText = '\t\tFinished mapping and removing PCR duplicates for ' + os.path.basename(fastqFile)
	print(printText)




	
