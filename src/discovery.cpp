/*
 * discovery.cpp
 *
 *  Created on: May 18, 2015
 *      Author: nok
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

#include "header/utilities.h"
#include "header/input.h"
#include "header/explained_genes.h"
#include "header/phenotype_genes.h"
#include "header/driver_genes.h"
#include "header/merge_and_trim.h"
#include "header/results.h"
#include "header/impact_scores.h"
#include "header/data_structures.h"


#ifdef _WIN32
	#include <direct.h>
	#include <windows.h>
	#include <direct.h>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
#endif

using namespace std;

int discovery(string outDir, string networkFilename, string expFilename, string snpFilename, string cnvFilename,
		string benchmarkGeneListFilename, string dbPath, int numThreads, string cancerType, bool noFoldchangeCutoff){

//	//create sub-directory to store output files
	#if defined(_WIN32)	//_WIN32 - Defined for applications for Win32 and Win64
//		_mkdir((outDir + "/sensitive").c_str());
//		_mkdir((outDir + "/stringent").c_str());
//		_mkdir((outDir + "/sensitive/samples").c_str());
//		_mkdir((outDir + "/stringent/samples").c_str());
		_mkdir((outDir + "/samples").c_str());
	#else
//		mkdir((outDir + "/sensitive").c_str(), 0777); // notice that 777 is different than 0777
//		mkdir((outDir + "/stringent").c_str(), 0777); // notice that 777 is different than 0777
//		mkdir((outDir + "/sensitive/samples").c_str(), 0777); // notice that 777 is different than 0777
//		mkdir((outDir + "/stringent/samples").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/samples").c_str(), 0777);
	#endif

	/*
	 * Read all input files of input samples
	 */

	//Read Network File
	cout << "reading network file and create mapping <id, geneSymbol> ..."	<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	readNetwork(networkFilename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId, '\t');
	int totalGenes = network.size();
	int totalGenesUpDown = totalGenes * 2;


	//Read gene expression
	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix ..." << endl;
	//sample id mapping
	vector<string> sampleIdToName;
	vector<int> genesEx; // gene ids ; size = # DE genes (consider only genes that are in the network)
	vector<string> geneSymbols;
	GeneExpression geneExpression;	//all info of gene expression
	geneExpression.genes = &genesEx;
	geneExpression.matrix = &originalGeneExpressionMatrix;

	readGeneExpression(expFilename.c_str(), &geneExpression, '\t',
			&geneSymbolToId, &sampleIdToName);

	int totalInputSamples = originalGeneExpressionMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesEx = originalGeneExpressionMatrix.size();

	cout << "\ttotal genes in gene expression matrix is " << numGenesEx << endl;
	cout << "\ttotal samples in gene expression matrix is " << totalInputSamples
			<< endl;

	//Read SNP (row = gene, col = samples)
	cout << "reading point mutation matrix ..." << endl;
	vector<int> genesPointMut;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalPointMutationsMatrix;
	PointMutations pointMutations;	//all info of point mutation
	pointMutations.genes = &genesPointMut;
	pointMutations.matrix = &originalPointMutationsMatrix;
	readPointMutations(snpFilename.c_str(), &pointMutations, '\t', &geneSymbolToId);

	totalInputSamples = originalPointMutationsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesPointMut = originalPointMutationsMatrix.size();

	cout << "\ttotal genes in point mutation matrix is " << numGenesPointMut << endl;
	cout << "\ttotal samples in point mutation matrix is " << totalInputSamples
			<< endl;

	//Read CNV
	cout << "reading CNV matrix ..." << endl;
	vector<int> genesCNV;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalCNVsMatrix;
	CopyNumberVariation CNVs;	//all info of CNV
	CNVs.genes = &genesCNV;
	CNVs.matrix = &originalCNVsMatrix;
	readCopyNumberVariation(cnvFilename.c_str(), &CNVs, '\t', &geneSymbolToId);

	totalInputSamples = originalCNVsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesCNV = originalCNVsMatrix.size();

	cout << "\ttotal genes in CNV matrix is " << numGenesCNV << endl;
	cout << "\ttotal samples in CNV matrix is " << totalInputSamples
			<< endl;

	/*
	 * Combining point mutation and CNV matrix
	 */

	//create a list of mutated genes, which contain point mutation or CNV or both
	vector<int> genesMut; // gene ids ; size = # mutated genes
	vector<bool>* isMutated = new vector<bool>(totalGenes, false);
	//1. check point mutation
	for (int i = 0; i < numGenesPointMut; ++i) {
		isMutated->at(genesPointMut[i]) = true;
	}
	//2. check CNV
	for (int i = 0; i < numGenesCNV; ++i) {
		isMutated->at(genesCNV[i]) = true;
	}
	//create a list of all mutated genes
	for (int i = 0; i < totalGenes; ++i) {
		if(isMutated->at(i)){
			genesMut.push_back(i);
		}
	}
	delete isMutated;

	//combine point mutation matrix and CNV matrix (row = gene, column = sample)
	cout << "combining point mutation matrix and CNV matrix ..." << endl;

	int numGenesMut = genesMut.size();
	TIntegerMatrix originalMutationMatrix;
	//for each mutated gene i
	for (int i = 0; i < numGenesMut; ++i) {
		int currentGeneId = genesMut[i];

		//find an index of a current gene in point mutation matrix
		int pmRowIndex = findIndex(&genesPointMut, currentGeneId);

		//find an index of a current gene in CNV matrix
		int cnvRowIndex = findIndex(&genesCNV, currentGeneId);

		vector<int> hasAMutatedGene(totalInputSamples);
		//for each sample j
		for (int j = 0; j < totalInputSamples; ++j) {
			if(pmRowIndex >= 0 and originalPointMutationsMatrix[pmRowIndex][j] != 0){
				hasAMutatedGene[j] = 1;
			}
			if(cnvRowIndex >= 0 and originalCNVsMatrix[cnvRowIndex][j] != 0){
				hasAMutatedGene[j] = 1;
			}
		}
		originalMutationMatrix.push_back(hasAMutatedGene);
	}

	Mutations mutations;	//all info of mutations (point mutation and CNV)
	mutations.genes = &genesMut;
	mutations.matrix = &originalMutationMatrix;

	totalInputSamples = originalMutationMatrix[0].size();

	cout << "\ttotal genes in mutation matrix is " << numGenesMut << endl;
	cout << "\ttotal samples in mutation matrix is " << totalInputSamples
			<< endl;

	/*
	 * Rename input samples to avoid using the same name as samples in the database
	 */

	for (int si = 0; si < totalInputSamples; ++si) {
		sampleIdToName[si] = sampleIdToName[si] + "_INPUT";
		//cout << sampleIdToName[si] << endl;
	}

	cout << "read files from database ..." << endl;

	/*
	 * Read parameters from file
	 */

	int L = 0;
	int D = 0;
	double F = 0.0;

	ifstream inParameterFile;
	string parameterFilename = dbPath + "/" + cancerType + "/JS.dat";	//F D L

	//read only the first line of the file
	inParameterFile.open(parameterFilename.c_str(), ifstream::in);
	if (inParameterFile.is_open()) {
		if (inParameterFile.good()) {
			string line;
			if (!getline(inParameterFile, line)){
				cerr << "Error reading parameter file (JS.dat)\n";
				return 0;
			}

			//getline(inParameterFile, line);	//skip the header
			//cout << line << endl;
			istringstream lineStream(line);

			int ci = 0;
			while (lineStream) {	//for each column (parameter)
				string param;
				if (!getline(lineStream, param, '\t'))
					break;
				if(ci == 0){	//F
					F = atof(param.c_str());
				}else if(ci == 1){	//D
					D = atoi(param.c_str());
				}else if(ci == 2){	//L
					L = atoi(param.c_str());
				}
				ci++;
			}

		}
		inParameterFile.close();
	} else {
		cerr << "Error opening parameter file (JS.dat)\n";
	}

	cout << "\tParameters (L,D,F) are set to " << L << ", " << D << ", " << F << endl;

	//TODO if the data type does not match, do not use cut off value for fold-change (F)
	if(noFoldchangeCutoff){
		F = 0.0;
		cout << "\tThe parameters (L,D,F) are set to " << L << ", " << D << ", " << F << " because the data types do not match" << endl;
	}


	/*
	 * Read phenotype genes from file
	 */

	string phenotypeFileName = dbPath + "/" + cancerType + "/PHENO.dat";
	vector<bool> isPhenotypeGenesUpDown(totalGenesUpDown, false);
	vector<int> phenotypeGeneIdsUpDown;	// phenotype gene ids
	//call function in input.h
	readPhenotypeGenesFromFile(phenotypeFileName.c_str(), &phenotypeGeneIdsUpDown, &geneSymbolToId);
	int numPhenotypeGene = phenotypeGeneIdsUpDown.size();
	cout << "\tThere are " << numPhenotypeGene << " phenotype genes" << endl;
	for (int i = 0; i < numPhenotypeGene; ++i) {
		isPhenotypeGenesUpDown[phenotypeGeneIdsUpDown[i]] = true;
	}

//	/*
//	 * Read driver list with the statistics
//	 */
//
//	vector<DriverGeneFromFile> driverGenesFromFileSensitive(totalGenes);
//	string driverFilenameSensitive = dbPath + "/sensitive/driver_list.txt";
//	readDriverGenesFromFile(driverFilenameSensitive.c_str(), &driverGenesFromFileSensitive, &geneSymbolToId);
//	cout << "Read drivers (sensitive) from file\n";
//
//	vector<DriverGeneFromFile> driverGenesFromFileStringent(totalGenes);
//	string driverFilenameStringent = dbPath + "/stringent/driver_list.txt";
//	readDriverGenesFromFile(driverFilenameStringent.c_str(), &driverGenesFromFileStringent, &geneSymbolToId);
//	cout << "Read drivers (stringent) from file\n";


	/*
	 * Read mutated gene list with their statistics
	 */

	vector<MutatedGeneFromFile> mutatedGenesFromFile(totalGenes);
	string mutatedGeneFilenameStringent = dbPath + "/" + cancerType + "/ALTERATION.dat";
	readMutatedGenesFromFile(mutatedGeneFilenameStringent.c_str(), &mutatedGenesFromFile, &geneSymbolToId);
	cout << "\tRead mutated genes from file (stringent) \n";

	/*
	 * Read cancer benchmark gene list
	 */

	vector<int> cancerBenchmarkGenes;
	vector<string> cancerBenchmarkGeneNames;
	readBenchmarkGeneList(benchmarkGeneListFilename, &cancerBenchmarkGenes, &geneSymbolToId, &cancerBenchmarkGeneNames);
	vector<bool> isCancerBenchmarkGenes(totalGenes);
	int numBenchmarkGenes = cancerBenchmarkGenes.size();
	for (int i = 0; i < numBenchmarkGenes; ++i) {
		isCancerBenchmarkGenes[cancerBenchmarkGenes[i]] = true;
	}
	cout << "\tRead " << numBenchmarkGenes << " cancer genes from Cancer Gene Census\n";

	/*
	 * Read drug-gene association from GDSC
	 */

	//read drug list
	map<int, string> drugIdToName;
	map<string, int> drugNameToId;
	string drugsListFilename = dbPath + "/GDSC_drug_list.txt";
	readDrugsListFromFile(drugsListFilename.c_str(), &drugIdToName, &drugNameToId);
	int numDrugs = drugIdToName.size();
	cout << "\tRead " << numDrugs << " drugs from GDSC\n";

	//read gene-drug assoc
	vector< vector<int> > geneDrugsAssocList(totalGenes);
	string drugsGeneAssocFilename = dbPath + "/GDSC_drug_gene_assoc.txt";
	readGenesDrugsAssocFromFile(drugsGeneAssocFilename.c_str(), &geneDrugsAssocList, &geneSymbolToId);
	int numGenesAssociatedWithDrugs = 0;
	for(int gi = 0; gi < totalGenes; gi++){
		int drugCount = geneDrugsAssocList[gi].size();
		if(drugCount > 0){
			numGenesAssociatedWithDrugs++;
//			cout << geneIdToSymbol[gi] << "\t" << drugCount << endl;
		}
	}
	cout << "\tRead " << numGenesAssociatedWithDrugs << " genes associated with drugs from GDSC\n";

	/*
	 * Construct modules for input samples
	 */

	cout << "constructing modules for input samples ...\n";

	//A vector for collecting mutated genes and the corresponding explained genes for all samples (sample, mutated genes, explained genes)
	vector< vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;
	vector< vector<int> > mutatedGeneIdsListReal;	//to save a list of mutated gene ids for each sample

	//add modules for input samples
	for (int si = 0; si < totalInputSamples; ++si) {

		//get gene expression of a current sample
		vector<double> sampleGeneExpression(totalGenes);// to save expression of of all genes in the network
		getGeneExpressionFromSampleId(&originalGeneExpressionMatrix, &genesEx,
				&sampleGeneExpression, si);

		//find mutated genes of a current sample
		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, si, &genesMut);
		int numMutatedGenes = mutatedGeneIds.size();

		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes(totalGenes); //contains explained genes of each driver

		//to tell whether each gene is an explained gene in this sample i
		//initialize the mutatedAndExplainedGenes for all genes
		for (int gi = 0; gi < totalGenes; ++gi) {
			mutatedAndExplainedGenes[gi].isExplainedGenesUpDown = new vector<bool>(totalGenes * 2, false);
		}

		for (int mi = 0; mi < numMutatedGenes; ++mi) {	// for each mutated genes
			int mutatedGeneId = mutatedGeneIds[mi];
			MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes[mutatedGeneId];

			//BFS for explained genes of the current mutated gene
			BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(&network, mutatedGeneId, L, D, F,
					meg->isExplainedGenesUpDown, &sampleGeneExpression, si, &geneIdToSymbol, &geneSymbolToId);
		}

		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
		mutatedGeneIdsListReal.push_back(mutatedGeneIds);

//		cout << sampleIdToName[si] << " has " << numMutatedGenes << endl;
//		cout << "added modules of " << sampleIdToName[si] << " to a list" << endl;

	}

	cout << "\ttotal input sample = " << sampleIdToName.size() << endl;

	/*
	 * Read MODULE file from database
	 */

	cout << "reading modules from database ..." << endl;

	string moduleFilename = dbPath + "/" + cancerType + "/MODULE.dat";

	//get id and name of samples (id of dbSamples start from totalSamples, the number of input samples);
	map<string, int> sampleNameToId;
	readModulesFromFile(&moduleFilename, &sampleIdToName, &sampleNameToId, &geneIdToSymbol, &geneSymbolToId,
			&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal);
	cout << "\ttotal sample = " << sampleIdToName.size() << endl;

	//update size of samples
	int totalSamples = sampleIdToName.size();

	/*
	 * Print the modules of pooled samples
	 */

//	//OUTPUT: print all modules in all samples (as original)
//	vector<string>* outputStr = new vector<string>;
//	for (int i = 0; i < totalSamples; ++i) {
//		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
//		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
//
//		//for each mutated genes
//		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
//			int currentMutatedGeneId = mutatedGeneIds[j];
//			string str = sampleIdToName[i] + "\t" + geneIdToSymbol[currentMutatedGeneId] + "\t";
//			vector<bool>* isExplainedGenesUpDownForAMutatedGene = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
//			int numMember = 0;
//			for (int k = 0; k < totalGenesUpDown; ++k) {
//				if(isExplainedGenesUpDownForAMutatedGene->at(k)){
//					if(k < totalGenes){
//						str += geneIdToSymbol[k] + "_UP" + ";";
//						numMember++;
//					}else{
//						str += geneIdToSymbol[k-totalGenes] + "_DOWN" + ";";
//						numMember++;
//					}
//				}
//			}
//			if(numMember > 0){
//				outputStr->push_back(str);
//			}
//		}
//	}
//	string outModulefilename = outDir + "/MODULE.dat";
//	writeStrVector(outModulefilename.c_str(), outputStr);
//	delete outputStr;

	/*
	 * STRINGENT ONLY
	 */

	cout << "STRINGENT mode ...\n";

	int mode = 1;

	{
		//collect all mutated genes in all samples => use mutatedGeneIdsListReal (samples, mutated genes, explained genes)
		list<int> mutatedGeneIdsList;	//mutated gene ids
		vector<bool> isMutatedGenes(totalGenes, false);
		getAllMutatedGenes(&mutatedGeneIdsListReal, &isMutatedGenes, &mutatedGeneIdsList);
		cout << "\ttotal number of mutated genes is " << mutatedGeneIdsList.size() << endl;

		vector< list<Module> > modulesListOfAllSamples(totalSamples);	//to save all modules in all samples

		//create bipartite graph ( mutated gene --- phenotype gene ). This is done at sample level, so have to remember sample id.
		//BipartiteEdge: each mutated gene contains a pair of (phenotype gene id, sample id)
		vector<BipartiteEdge>* bipartiteGraph = new vector<BipartiteEdge>(totalGenes);
		createBipartiteGraph(&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
				&isPhenotypeGenesUpDown, bipartiteGraph, &geneIdToSymbol);

		vector<DriverGene> driverGenes;	//driver gene id with sample ids

		cout << "\tperforming greedy minimum set cover algorithm ...\n";
		cout << "\tfinding driver genes ...\n";
		findDriverGenes(bipartiteGraph, &mutatedGeneIdsList, &driverGenes);
		delete bipartiteGraph;

		int numDriverGenes = driverGenes.size();
		cout << "\ttotal driver genes = " << numDriverGenes << endl;

		cout << "\tgenerating modules for all samples ...\n";

		//merge modules and trim explained genes for each sample
		//use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genesUpDown)
		findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
				&mutatedGeneIdsListReal, &isPhenotypeGenesUpDown, &driverGenes, &phenotypeGeneIdsUpDown, mode);

		cout << "\ttrimming explained genes for all samples ...\n";
		trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D, &geneIdToSymbol);

		//write only final module of input sample
		cout << "\twriting final module to FINAL_MODULE.dat ...\n";
		string outFinalModuleFilenameStringent = outDir + "/FINAL_MODULE.dat";
		saveModulesOfInputSamples(&modulesListOfAllSamples, &mutatedAndExplainedGenesListReal, outFinalModuleFilenameStringent,
				&geneIdToSymbol, &sampleIdToName, totalInputSamples);

		cout << "\tcalculating IMPACT scores for all input samples ...\n";
		vector< vector<Driver> > driversOfAllSamples(totalInputSamples);
		string outDriverOfAllSamplesDirName = outDir + "/samples";

		//calculated impact score for the input samples
		calculateImpactScoresForAllInputSamples(totalInputSamples, &modulesListOfAllSamples, &driversOfAllSamples,
				&originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol,
				&sampleIdToName);

		cout << "\taggregating IMPACT scores across all samples ...\n";
		vector<double> driverAggregatedScores(totalGenes, 0);
		vector<int> driversFrequency(totalGenes, 0);
		aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol);

		cout << "\tprinting impact scores for all samples ...\n";
		printSampleDriverListForInputSamples(totalInputSamples, &driversOfAllSamples,
				outDriverOfAllSamplesDirName, &geneIdToSymbol, &sampleIdToName,
				&mutatedGenesFromFile, &cancerBenchmarkGeneNames,
				&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&drugIdToName, &geneDrugsAssocList);


		cout << "\tprinting aggregated impact scores ...\n";
		string outDriverListfilename = outDir + "/driver_list.dat";
		printAggregatedDriverListForInputSamples(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &mutatedGenesFromFile, &cancerBenchmarkGeneNames, &drugIdToName, &geneDrugsAssocList);

	}


	/*
	 * Clean-up
	 */

	//delete the vector<int>* mutatedAndExplainedGenesListReal
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		for (int gi = 0; gi < totalGenes; ++gi) {
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[gi].isExplainedGenesUpDown;
			delete isExplainedGenesUpDown;
		}
	}

	return 0;
}


