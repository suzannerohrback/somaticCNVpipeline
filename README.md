# somaticCNVpipeline
Pipeline for somatic CNV detection from single-cell whole-genome sequencing data of non-cancerous mouse or human tissue.

Please see the files in the README folder for more detailed information on how to use this code


========================================================================================================================
INSTALLATION
========================================================================================================================

1) Software Dependancies -- install these separate programs if you do not currently have them
	a) Bowtie 1 
		http://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie
	b) Samtools v0.1.19
		https://sourceforge.net/projects/samtools/files/samtools/0.1.19/
	c) Statsmodels (python package)
		For the best results, enter the following at the command line
			pip install --user statsmodels
	d) Matlab
		And this should be accessible at the command line
		Only required for the CBS step of the segment function
		
2) Create mapping reference(s) for mm10 and/or hg38, if you do not currently have them
	http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer
	I recommend follwing the strategy outlined in the following publication for the most accurate results 
		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069701/	

3) Obtain this repository
	At the command line, enter the following
		cd /desired/directory/for/code (ex. cd ~/GitRepos/)
		git clone https://github.com/suzannerohrback/somaticCNVpipeline

4) If updates have been made to the code
	At the command line, enter the following
		cd /desired/directory/for/code/somaticCNVpipeline
		git fetch
		git merge origin/master



========================================================================================================================
Command Line
========================================================================================================================

The following commands assumes you have cd'd the the directory the somaticCNVcalling script is located
If you have added the /path/to/somaticCNVpipeline folder to your PATH variable 
	It is not necessary to cd to that folder
	You do not need to specify python when calling the function

Get a help message:
	python somaticCNVcalling -h
	python somaticCNVcalling --help
	python somaticCNVcalling [function] -h
	python somaticCNVcalling [function] --help
	
General usage (see function README files for further information):
	python somaticCNVcalling preprocess [options]
	python somaticCNVcalling map [options]
	python somaticCNVcalling count [options]
	python somaticCNVcalling segment [options]
	python somaticCNVcalling interpret [options]

See EXAMPLE_SHELL_COMMANDS.sh for actual examples


========================================================================================================================
CITATION
========================================================================================================================
Feel free to use this code, but please provide credit if utilized for a publication
Our paper has not yet been published, so currently cite as a URL:
  Rohrback, Suzanne. Somatic CNV Analysis Pipeline. https://github.com/suzannerohrback/somaticCNVpipeline/. Accession Date.
  ###PAPER WAS RECENTLY ACCEPTED TO PNAS, FULL CITATION COMING SOON###
  


========================================================================================================================
DISCLAIMER
========================================================================================================================
Not all of the code in this repository will work on all systems "out of the box".
  That is because one purpose of this repository is to support continuity within the Jerold Chun lab.
  Portions of the code that are most likely to require modification or independent creation include:
    Mapping
    Bin location reference
    Segmentation
  There are multiple strategies for each of these actions, and it is beyond the scope of this repo to be universally compatible.





  

