oncoIMPACT

-----------
What is it?
-----------
Cancer-omics Data Integration pipeline



------------------
The latest version
------------------
Details of the latest version can be found at our SourceForge site under https://sourceforge.net/projects/oncoimpact/.



------------
Installation
------------
Download and uncompress the latest version
	tar -zxvf oncoIMPACT.tgz



-----
Usage
-----
Before running the pipeline, please ensure that you have the following files ready.

__Input file format__
-CNV data: a tab-separated file with samples as columns and genes as rows. Each entry of this matrix must hold either a value of -1/0/1. A value of -1 would indicate a DELETION event while a value of 1 will indicate an AMPLIFICATION event, and a value of 0 if CNV data is not available for that particular sample.

-SNP data: a tab-separated file with samples as columns and genes as rows. Each entry of this matrix must hold either the binary values of 0/1. A value of 1 would indicate the presence of a snp in the gene of that particular sample while a value of 0 would indicate otherwise.

-EXP data: a tab-separated file with samples as columns and genes as rows. Each entry of this matrix would represent the log fold-change value.


__Config file__
You will need to create a config file for your specific project. The config file needs to contain the following parameters (you may refer to the sample config file provided with the scripts under <Scripts Directory/sampleConfig.cfg>).

	-outDir: Full path to destination folder. Two sub-folders will be created in this directory: INCOMPLETE_SAMPLES (for incomplete samples) and COMPLETE_SAMPLES (for complete samples)

    -scriptDir: Full path to folder where oncoIMPACT is installed

    -numThreads: Number of threads to use

    -cnv: Full path to cnv data matrix

    -exp: Full path to expression data matrix

    -snp: Full path to snp data matrix
    
    -testMode: Boolean flag to toggle test mode (valid options: 0 / 1)

When you are ready to run the oncoIMPACT pipeline, simply enter the following command
oncoIMPACT.pl <path to config file> <fraction of samples used during parameters estimation>
	


-----------
Sample data
-----------
We have provided two sample datasets for to test the oncoIMPACT pipeline which you may download from the website https://sourceforge.net/projects/oncoimpact/files;
-Glioblastoma.tgz: Glioblastoma TCGA dataset
-Ovarian.tgz: Ovarian cancer TCGA dataset

To run OncoIMPACT on those dataset: oncoIMPACT.pl <path to config file> 0.2


----------
Change Log
----------
v0.9.2:
- NEW: Option in configuration file to run oncoIMPACT in test mode which performs the simulation with fewer iterations and fixed seed. In this mode, oncoIMPACT should complete in less than 2 hours using a single thread.
- NEW: Sanity checks to ensure validity of parameters provided by user
- FIX: Improved disk space utilization
- FIX: Improved compatibility with Mac OS.


v0.9.1:
- oncoIMPACT will now avoid reproducing the input files if COMPLETE_SAMPLES folder exists
- fix for bugs introduced in last version



---------
Licensing
---------
The MIT License (MIT)
Copyright (c) 2014 Genome Institute of Singapore



--------
Contacts
--------
If you want community-driven support, please visit the forum at https://sourceforge.net/p/oncoimpact/discussion/.

If you have a bug to report, you may raise a ticket at https://sourceforge.net/p/oncoimpact/tickets/.

If you have other questions or feedback, you may direct them to Burton Chia (chiakhb@gis.a-star.edu.sg).
