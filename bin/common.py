#!/usr/bin/python
import os
import multiprocessing as mp
import subprocessing as sub
import shlex









def fixDirName(dirpath):
	if dirpath[-1] != '/':
		dirpath += '/'
	return dirpath





def makeDir(dirpath):
	if not os.path.exists(dirpath):
		os.mkdir(dirpath)
	return 0










def importSampleList(infile):
	if os.path.exists(infile):
		files = []
		with open(infile, 'r') as IN:
			for x in IN:
				files.append(x.rstrip())
	else:
		errorText = '\nERROR: the specified sample name file does not exist, please fix\n\t' + infile + '\n'
		print(errorText)
		raise SystemExit
	
	if len(files) == 0:
		errorText = '\nERROR: The sample name file does not contain any sample names, please fix\n'
		print(errorText)
		raise SystemExit
		
	return files
		
		
		







###daemon to run multiprocessing and parallelize tasks###
def daemon(target, argList, name, cpuPerProcess=1):
	print( str( '\t' + str(len(argList)) + ' processes to run to ' + name ) )
	numCPU = mp.cpu_count()
	numWorkers = min( [int(numCPU / cpuPerProcess), len(argList)] )


	pool = mp.Pool(numWorkers)
	processes = [pool.apply_async(target, args=x) for x in argList]
	pool.close()

	for i,j in enumerate(processes):
		j.wait()

		if not j.successful():
			pool.terminate()

			print '\n\n\nprocessing failed, getting traceback now...'
	
			p = mp.Process(target=target, args=argList[i])
			p.start()
			p.join()
	#	else:
	#		print( str( '\t\t' + str(i+1) + ' of ' + str(len(argList)) + ' processes complete' ) )

	print( str( '\tAll processing to ' + name + ' complete\n' ) )

	
	
	
	
	
	
	
	
	
def zipping(filepath, gunzip=True):
	if filepath.split('.')[-1] != 'gz' and gunzip:
		return filepath
	elif filepath.split('.')[-1] == 'gz' and not gunzip:
		return filepath
	
	if gunzip:
		cmd = 'gunzip ' + filepath
		fixname = filepath[:-3]
	else:
		cmd = 'gzip ' + filepath
		fixname = filepath + '.gz'
	
	cmd = shlex.split(cmd)
	p = sub.Popen(cmd)
	p.wait()
	
	return fixname
	
	
	
	
	
	
	
	
	
	
