#!/usr/bin/python
import multiprocessing as mp










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
		else:
			print( str( '\t\t' + str(i+1) + ' of ' + str(len(argList)) + ' processes complete' ) )

	print( str( '\tAll processing to ' + name + ' complete\n' ) )
