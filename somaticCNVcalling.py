#!/usr/bin/python
import sys

import arguments





if __name__ == '__main__':
	objective, args=arguments.fullParser(sys.argv[1:])
	
	if objective == 'preprocess':
		import preprocess
		runpreprocess.runAll(args)
