-----------
What is it?
-----------
OncoIMPACT is a first-in-class algorithmic framework that nominates patient-specific driver genes by integratively modeling genomic mutations (point, structural and copy-number) and the resulting perturbations in transcriptional programs via defined molecular networks.

OncoIMPACT is configured to run in two modes: (1) a database mode that allows it to determine parameter settings from the data sets provided and (2) a discovery mode where information in the provided database is used to predict driver genes for each sample in an additional data set (which can be the same as the one used to create the database).

------------
Installation
------------
Download and uncompress the latest version
	tar -zxvf oncoIMPACT.tgz
    chmod +x oncoIMPACT_v0.9.3/*.pl

-----
Usage
-----
Before running the pipeline, please ensure that you have the following files ready.


__Input file format__

- CNV_data: a tab-separated file with samples as columns and genes as rows. Each entry of this matrix must hold either a value of -1/0/1. A value of -1 would indicate a DELETION event while a value of 1 will indicate an AMPLIFICATION event, and a value of 0 if the is not affected by copy number gene.

- SNP_data (point mutations and short indels): a tab-separated file with samples as columns and genes as rows. Each entry of this matrix must hold either the binary values of 0/1. A value of 1 would indicate the presence of a snp in the gene of that particular sample while a value of 0 would indicate otherwise.

- EXP_data: a tab-separated file with samples as columns and genes as rows. Each entry of this matrix would represent the log2 fold-change of the gene expression of the tumor compared to a normal control.


__Config file__
You will need to create a config file for your specific project. The config file needs to contain the following parameters (you may refer to the sample config file provided with the scripts under <_Scripts Directory/sampleConfig.cfg_\>).

- outDir: Full path to destination folder

- scriptDir: Full path to folder where oncoIMPACT is installed

- numThreads: Number of threads to use

- cnv: Full path to cnv data matrix

- exp: Full path to expression data matrix

- dataType: Flag for expression data type. Vvalid options: ARRAY (default), RNA_SEQ

- snp: Full path to snp data matrix

- dataBase: Full path to the pre-computed database (implies discovery mode only)

- databaseExport: Full path where the database will be exported (implies database + discovery mode)

- testMode: Boolean flag to toggle test mode (valid options: 0 / 1)


When you are ready to run the OncoIMPACT pipeline, simply enter the following command
	oncoIMPACT.pl <path to config file\>

-----------
Sample data
-----------
We have provided 1 sample datasets for to test the OncoIMPACT pipeline;

- Glioblastoma.tgz: Glioblastoma (GBM) TCGA dataset

- Glioblastoma_database.tgz: pre-computed database of the TCGA GBM dataset

- single_patient.tgz: data for a single single GBM patient (to test the discovery mode)


----------
Change Log
----------

__v0.9.3__:

- NEW: Enable the processing of RNA-seq data

- NEW: Enable the construction of databases

- NEW: Enable the discovery mode using a pre-computed database


__v0.9.2__:

- NEW: Option in configuration file to run OncoIMPACT in test mode which performs the simulation with fewer iterations and fixed seed. In this mode, OncoIMPACT should complete in less than 2 hours using a single thread.

- NEW: Sanity checks to ensure validity of parameters provided by user

- FIX: Improved disk space utilization

- FIX: Improved compatibility with Mac OS.


__v0.9.1__:

- oncoIMPACT will now avoid reproducing the input files if COMPLETE_SAMPLES folder exists

- fix for bugs introduced in last version



---------
Licensing
---------
The MIT License (MIT)
Copyright (c) 2014 Genome Institute of Singapore



-----------------
Citing OncoIMPACT
-----------------
Bertrand _et. al._. __Patient-specific driver gene prediction and risk assessment through integrated network analysis of cancer omics profiles.__ [_Nucleic Acid Research_ 2015, 43 (7): e44](http://nar.oxfordjournals.org/content/43/7/e44.long)



--------
Contacts
--------
If you want community-driven support, please visit the forum at https://sourceforge.net/p/oncoimpact/discussion/.

If you have a bug to report, you may raise a ticket at https://sourceforge.net/p/oncoimpact/tickets/.

If you have other questions or feedback, you may direct them to Burton Chia (<chiakhb@gis.a-star.edu.sg>).
