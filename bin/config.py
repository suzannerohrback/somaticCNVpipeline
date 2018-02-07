#!/usr/bin/python
import os
import inspect
import sys


"""
This code can be modified to suit your specific environment
  Such as reference file locations
  Software locations
  
PLEASE DO NOT SEND ME PULL REQUESTS INCLUDING ANY MODIFICATIONS YOU MAKE TO THIS FILE, I WILL ALMOST CERTAINLY REJECT THEM
  If I have missed adding in a specific reference, let me know by message
"""










class Map:
	
	def __init__(self):
		self.currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
		
		self.indexDict =	{
							'hg38': self.currentdir + '/reference/hg38_index/hg38.index',
							'mm10': self.currentdir + '/reference/mm10_index/mm10.index',
							}
		
		#THIS MUST BE VERSION 0.1.19 OR THE PIPELINE MAY NOT WORK#
		self.samtools = '/home/k4zhang/softwares/samtools-0.1.19/samtools'
		
		self.bowtieOptions =	[
								'bowtie',
								 '-p 8',
								 '-S', 
								 '-t',
								 '-n 2',
								 '-e 70',
								 '-m 1',
								 '--phred33-quals',
								 '--best',
								 '--strata',
								 '--chunkmbs 200',
								]




		
		
		
		
		
		
class Count:
	
	def __init__(self):
		self.currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
		
		self.binDict =	{
						'hg38': self.currentdir + '/reference/hg38.varbin.gc.content.25k.bowtie.k36.txt',
						'mm10': self.currentdir + '/reference/mm10.varbin.gc.content.25k.bowtie.k36.txt',
						}
		
		self.chromDict =	{
							'hg38': self.currentdir + '/reference/hg38.chrom.sizes.txt',
							'mm10': self.currentdir + '/reference/mm10.chrom.sizes.txt',
							}






		
		
		
class Segment:
	
	def __init__(self):
		self.currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
		
		self.binDict =	{
						'hg38': self.currentdir + '/reference/hg38.varbin.fullRef.25k.bowtie.k36.txt',
						'mm10': self.currentdir + '/reference/mm10.varbin.fullRef.25k.bowtie.k36.txt',
						}
			
		self.statsmodelsLocation = '/home/srohrbac/research/packages/statsmodels-0.6.1'
		sys.path.append(self.statsmodelsLocation)
		
		self.chromNumDict = {
							'hg38': '25',
							'mm10': '21',
							}
		
		self.matlabBase = self.currentdir + '/reference/matlab.base.txt'
		self.CBSalpha = 0.01
		
		
		
		
		
		
