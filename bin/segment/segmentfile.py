#!/usr/bin/python
import sys
import os
import inspect
import subprocess as sub
import shlex

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
import config as cfg










def writeMatlabScript(sample, species, tempDir, lowessDir, segmentDir):
	segVars = cfg.Segment()
	
#	scriptFile = tempDir + sample + '.m'
#	OUT = open(scriptFile, 'w')
	matlabName = ''.join(sample.split('_'))
	matlabName = ''.join(matlabName.split('-'))
	matlabName = ''.join(matlabName.split('.'))
	scriptFile = tempDir + matlabName + '.m'
	OUT = open(scriptFile, 'w')

	
	
	OUT.write('%Sample specific variable definitions\n')
	
	OUT.write(str("refFile = '" + segVars.binDict[species] + "';\n"))
	OUT.write(str("binFile = '" + lowessDir + sample + ".lowess.txt';\n"))
	OUT.write(str("saveFile = '" + segmentDir + sample + ".segments.txt';\n"))
	OUT.write(str("chromNum = " + str(segVars.chromNumDict[species]) + ";\n"))
	OUT.write(str("alpha = " + str(segVars.CBSalpha) + ";\n"))
#	OUT.write(str("refFile = '" + refFile + "';\n"))
#	OUT.write(str("binFile = '" + lowessDir + sample + lowessExt + "';\n"))
#	OUT.write(str("saveFile = '" + segmentDir + sample + ".segments.txt';\n"))
#	OUT.write(str("chromNum = " + str(chromNumDict[species]) + ";\n"))
#	OUT.write(str("alpha = " + str(CBSalpha) + ";\n"))

	OUT.write('\n\n\n\n\n%Generic processing code\n')
	
	
	
	IN = open(segVars.matlabBase, 'r')
	
	for x in IN:
		OUT.write(x)
	
	OUT.close()
	IN.close()
	
	return matlabName










def segmentOne(sample, species, tempDir, lowessDir, segmentDir):
	#write matlab script
	scriptName = writeMatlabScript(sample, species, tempDir, lowessDir, segmentDir)
	
	os.chdir(tempDir)
	print os.getcwd()
	
	#run matlab script
	stdoutFile = tempDir + sample + '.stderr.txt'
	stdout = open(stdoutFile, 'w')
	
	cmd = 'matlab -nodisplay -r ' + scriptName
	cmd = shlex.split(cmd)
	
	p = sub.Popen(cmd, stdout = stdout, stderr = sub.STDOUT)
	p.wait()
	
	stdout.close()
	
	#delete intermediate files
#	os.remove(tempDir + scriptName + '.m')
#	os.remove(stdoutFile)
	#THIS SHOULD BE DONE IN runsegment.py WITH SHUTIL, LIKE FOR runmap.py
	
	
	
	
