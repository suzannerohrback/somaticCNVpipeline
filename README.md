# somaticCNVpipeline
Pipeline for somatic CNV detection from single-cell whole-genome sequencing data of non-cancerous mouse or human tissue.



========================================================================================================================
CITATION
========================================================================================================================
Feel free to use this code, but please provide credit if utilized for a publication
Our paper has not yet been published, so currently cite as a URL:
  Rohrback, Suzanne. Somatic CNV Analysis Pipeline. https://github.com/suzannerohrback/somaticCNVpipeline/. <Accession Date>.
  


========================================================================================================================
DISCLAIMER
========================================================================================================================
Not all of the code in this repository will work on all systems "out of the box"
  That is because one purpose of this repository is to support continuity within the Jerold Chun lab
  Portions of the code that are most likely to require modification or independent creation include
    Mapping
    Bin location reference
    Segmentation
  There are multiple strategies for each of these actions, and it is beyond the scope of this repo to be universally compatible
This document will aim to indicate which portions of the code are less universal 



========================================================================================================================
INSTALLATION
========================================================================================================================
Currently nothing fancy is supported (sorry, I'm still learning)
Individual users are responsible for ensuring all dependencies are met

1, Download the entire repo/folder and put in the desired location
2, Add this filepath to your PATH variable
  e.g. in ~/.bashrc
    PATH=$PATH:/path/to/repository/foldername
    separate from other existing paths in that variable by a : (no spaces)
3, Add the filepath to the bin folder to your PYTHONPATH variable
  e.g. in ~/.bashrc
    PYTHONPATH=/path/to/repository/foldername/bin
4, Make sure the parent code file is executable
   e.g. chmod +x /path/to/repository/foldername/parentCode
   
   
   
========================================================================================================================
PRE-PROCESSING
========================================================================================================================
Code for this section is not yet prepared for release
   


========================================================================================================================
MAPPING
========================================================================================================================
Code for this section is not yet prepared for release
UNLESS YOU ARE ON THE TSCC CLUSTER THIS WILL REQUIRE MODIFICATIONS TO THE SAMTOOLS FILEPATH (v0.1.19 required)




========================================================================================================================
COUNTING
========================================================================================================================
Code for this section is not yet prepared for release
  However, it will be nearly identical to that created by Baslan, T. et al. Nat Protoc. 7, 1024. 2012.
    Please cite this publication as well if using this section of the code
    I did make one modification to this script to correct an error that prevented counting reads in the last bin of the genome
THIS SHOULD BE UNIVERSALLY FUNCTIONAL (assuming you are either using the same bin boundaries, or have replaced that reference file)



========================================================================================================================
NORMALIZING
========================================================================================================================
Code for this section is not yet prepared for release
THIS SHOULD BE UNIVERSALLY FUNCTIONAL (assuming access to the statsmodels package)
  However, I have had to add in corrections regularly for unusual datasets, so please report issues



========================================================================================================================
SEGMENTING
========================================================================================================================
Code for this section is not yet prepared for release
THIS CODE REQUIRES COMMAND LINE ACCESS TO MATLAB



========================================================================================================================
QC TESTING
========================================================================================================================
Code for this section is not yet prepared for release
THIS SHOULD BE UNIVERSALLY FUNCTIONAL



========================================================================================================================
CNV CALLING
========================================================================================================================
Code for this section is not yet prepared for release
THIS SHOULD BE UNIVERSALLY FUNCTIONAL



========================================================================================================================
TEXT FILES
========================================================================================================================
Code for this section is not yet prepared for release
THIS SHOULD BE UNIVERSALLY FUNCTIONAL
Two types of files are included
  1, Reference files used at different processing stages
  2, Test data files from each stage of processing to enable independent troubleshooting



========================================================================================================================
CNV ANALYSIS
========================================================================================================================
Code for this section is not yet prepared for release
THIS IS NOT INTEGRATED INTO THE PIPELINE
  Just a "data dump" from Suzanne's PhD work
  It may require substantial modifications to become universally functional
  Included primarily for examples of the forms initial analysis can take
  
  
  
========================================================================================================================
CUTOFF DEVELOPMENT
========================================================================================================================
Code for this section is not yet prepared for release
THIS IS NOT INTEGRATED INTO THE PIPELINE
Includes code used to perform the analysis that allowed me to create cutoffs for both
  1, Quality control requirements (contraction algorithm)
  2, CNV call post-hoc filtering (bootstrapped one-class support vector machines)
The reference files I used for performing these analyses are included as well 
  But you desire to do so you will need to write your own parent script to submit them for processing
  
  

========================================================================================================================
SIMULATIONS
========================================================================================================================
Code for this section is not yet prepared for release
THIS IS NOT INTEGRATED INTO THE PIPELINE
Includes code used to simulate data for both
  1, Identifying the copy number state dependence of noise
  2, Creating false positive CNV calls
  3, Creating true positive CNV calls
The data files I used are NOT included here, as there are hundreds of them
  A selection of the BioSamples in PRJNA415480
  

