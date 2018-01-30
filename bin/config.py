#!/usr/bin/python
import os
import inspect


"""
This code can be modified to suit your specific environment
  Such as file locations
  Execultable locations
  
PLEASE DO NOT SEND ME PULL REQUESTS INCLUDING ANY MODIFICATIONS YOU MAKE TO THIS FILE, I WILL REJECT THEM
  If I have missed adding in a specific reference, let me know by message
"""


class Map:
	
	def __init__(self):
		self.currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
		
		self.indexDict =	{
							'hg38': self.currentdir + '/reference/hg38_index/hg38.index',
							'mm10': self.currentdir + '/referencemm10_index/mm10.index',
							}
		
		self.samtools = '/home/k4zhang/softwares/samtools-0.1.19/samtools'
		
###		self.bowtie = 



