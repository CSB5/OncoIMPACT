//============================================================================
// Name        : tester.cpp
// Author      : Nok C Suphavilai
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <set>
#include <omp.h>
#include "header/utilities.h"
#include "header/input.h"
#include "header/sampling.h"
#include "header/explained_genes.h"
#include "header/phenotype_genes.h"
#include "header/driver_genes.h"
#include "header/merge_and_trim.h"
#include "header/results.h"
#include "header/parameters.h"
#include "header/impact_scores.h"

//#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main() {

	/*
	 * Start timer
	 */

	clock_t begin_time = clock();
	clock_t step_time = clock();

	/*
	 * TODO Read configuration from a file
	 */

	//0 = sensitive, 1 = stringent
	int mode = -1;
	int numThreads = 10;

	map<string, string> conf;

	/*
	 * Read network file into adjacency list
	 */

	cout << "reading network file and create mapping <id, geneSymbol> ..."
			<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	string networkFilename = "network_FIsInGene_041709.txt";
	readNetwork(networkFilename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId,
			'\t');
	int totalGenes = network.size();
	int totalGenesUpDown = totalGenes * 2;

	/*
	 * Read gene expression matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of samples is the same as in mutation matrix (also the same order)
	 */

	//sample id mapping
	vector<string> sampleIdToName;

	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix ..." << endl;
	string expFilename = "GBM/EXPR.txt";
	vector<int> genesEx; // gene ids ; size = # DE genes (consider only genes that are in the network)
	vector<string> geneSymbols;
	GeneExpression geneExpression;	//all info of gene expression
	geneExpression.genes = &genesEx;
	geneExpression.matrix = &originalGeneExpressionMatrix;
	readGeneExpression(expFilename.c_str(), &geneExpression, '\t',
			&geneSymbolToId, &sampleIdToName);

	int totalSamples = originalGeneExpressionMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesEx = originalGeneExpressionMatrix.size();

	cout << "\ttotal genes in gene expression matrix is " << numGenesEx << endl;
	cout << "\ttotal samples in gene expression matrix is " << totalSamples
			<< endl;

	/*
	 * Read point mutation matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of samples is the same as in gene expression matrix (also the same order)
	 */

	cout << "reading point mutation matrix ..." << endl;
	string snpFilename = "GBM/SNP.txt";
	vector<int> genesPointMut;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalPointMutationsMatrix;
	PointMutations pointMutations;	//all info of point mutation
	pointMutations.genes = &genesPointMut;
	pointMutations.matrix = &originalPointMutationsMatrix;
	readPointMutations(snpFilename.c_str(), &pointMutations, '\t', &geneSymbolToId);

	totalSamples = originalPointMutationsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesPointMut = originalPointMutationsMatrix.size();

	cout << "\ttotal genes in point mutation matrix is " << numGenesPointMut << endl;
	cout << "\ttotal samples in point mutation matrix is " << totalSamples
			<< endl;

	/*
	 * Read CNV matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of sample is the same as in gene expression matrix (including order)
	 */

	cout << "reading CNV matrix ..." << endl;
	string cnvFilename = "GBM/CNV.txt";
	vector<int> genesCNV;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalCNVsMatrix;
	CopyNumberVariation CNVs;	//all info of CNV
	CNVs.genes = &genesCNV;
	CNVs.matrix = &originalCNVsMatrix;
	readCopyNumberVariation(cnvFilename.c_str(), &CNVs, '\t', &geneSymbolToId);

	totalSamples = originalCNVsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesCNV = originalCNVsMatrix.size();

	cout << "\ttotal genes in CNV matrix is " << numGenesCNV << endl;
	cout << "\ttotal samples in CNV matrix is " << totalSamples
			<< endl;

	/*
	 * Combining point mutation and CNV matrix
	 */

	//TODO (next version) refactoring the code e.g.separate the combining code to a function in input.cpp
	//create a list of mutated genes, which contain point mutation or CNV or both
	vector<int> genesMut; // gene ids ; size = # mutated genes
	vector<bool>* isMutated = new vector<bool>(totalGenes);
	//initialize isMutated, pointMutationCount, CNVCount
	for (int i = 0; i < totalGenes; ++i) {
		isMutated->at(i) = false;
	}
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

		vector<int> hasAMutatedGene(totalSamples);
		//for each sample j
		for (int j = 0; j < totalSamples; ++j) {
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

	totalSamples = originalMutationMatrix[0].size();

	cout << "\ttotal genes in mutation matrix is " << numGenesMut << endl;
	cout << "\ttotal samples in mutation matrix is " << totalSamples
			<< endl;

	/*
	 * Calculate JS divergence for each set of parameters
	 */

//	//Note: gene expression and mutation have the same set of samples
//	//If number of samples is < 50, use all samples to tune the parameters
//	int numSamples = 50;
//	if(totalSamples < 50){
//		numSamples = totalSamples;
//	}
//
//	cout << "\t\tDONE reading and preparing data (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";
//
//	step_time = clock();
//	cout << "tuning parameters by using " << numSamples << " randomly choosing samples ..."
//			<< endl;
//
//	cout << "computing JS divergence for all parameters (L,D,F) ... " << endl;
//
//	//TODO initialize all the parameters to be tested
////	int LsVal[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
//	int LsVal[] = {16, 18, 20};	//fewer parameters for testing
//	vector<int> Ls(LsVal, LsVal + sizeof LsVal / sizeof LsVal[0]);
//	int numLs = Ls.size();
////	int DsVal[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
//	int DsVal[] = {60, 65, 70};	//fewer parameters for testing
//	vector<int> Ds(DsVal, DsVal + sizeof DsVal / sizeof DsVal[0]);
//	int numDs = Ds.size();	//fewer parameters for testing
////	double FsVal[] = {1, 1.5, 2, 2.5, 3};
//	double FsVal[] = {2, 2.5};	//fewer parameters for testing
//	vector<double> Fs(FsVal, FsVal + sizeof FsVal / sizeof FsVal[0]);
//	int numFs = Fs.size();
//
//	//TODO debug the js divegence calculation (the results are quite different from the original)
//	int numCombinations = numLs * numDs * numFs;
//	//initialize the vector to save the divergence of each set of parameters
//	vector<JSDivergence> jsDivergences(numCombinations);
//
//	int numPermutationsForJSDivergence = 100;
//	//TODO parameter findind: the result is quite different from the original
//	findParameters(&jsDivergences, &Ls, &Ds, &Fs, totalGenes, &geneExpression, &mutations, &network, numSamples, numPermutationsForJSDivergence, &geneSymbolToId, numThreads);
//
//	//write the JS divergence result to a file
//	string outJSFilename = "output/parameters.dat";
//	//save the JS divergence results
//	saveJSDivergences(&jsDivergences, outJSFilename);
//
//	cout << "\t\tDONE finding parameters (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";
//
//	//choose the best parameters
//	JSDivergence maxJs;
//	findMaximumJsDivergence(&jsDivergences, &maxJs);
//
//	cout << "the maximum divergence is " << maxJs.divergence << " when L, D, F = " << maxJs.L << ", " << maxJs.D << ", " << maxJs.F << endl;
//
//	//set the L D F to maxJs
//	int L = maxJs.L;
//	int D = maxJs.D;
//	double F = maxJs.F;

	//[DEBUG]
	int L = 20;
	int D = 65;
	double F = 2.5;

	/*
	 * Find phenotype genes
	 */

	step_time = clock();
	cout << "finding phenotype genes using (L,D,F) = (" << L << "," << D << "," << F << ") ..." << endl;

	//A vector for collecting mutated genes and the corresponding explained genes for all samples (sample, mutated genes, explained genes)
	vector< vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;
	//Note: size of MutatedAndExplianedGenes = totalGenes, so index = mutated gene id, value = a list of explained gene ids
	vector< vector<int> > mutatedGeneIdsListReal;	//to save a list of mutated gene ids for each sample

	cout << "\tgetting explained genes frequency of the real samples ...\n";

	for (int i = 0; i < totalSamples; ++i) {

		//get gene expression of a current sample
		vector<double> sampleGeneExpression(totalGenes);// to save expression of of all genes in the network
		getGeneExpressionFromSampleId(&originalGeneExpressionMatrix, &genesEx,
				&sampleGeneExpression, i);

		//find mutated genes of a current sample
		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
				&genesMut);
		int numMutatedGenes = mutatedGeneIds.size();

		//find explained genes of a current sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes(totalGenes); //contains explained genes of each driver
		for (int j = 0; j < numMutatedGenes; ++j) {	// for each mutated genes
			int mutatedGeneId = mutatedGeneIds[j];
			MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes[mutatedGeneId];
			meg->isExplainedGenesUpDown = new vector<bool>(totalGenes * 2);	//to tell whether each gene is an explained gene in this sample i
			//BFS for explained genes of the current mutated gene
			BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(&network, mutatedGeneId, L, D, F,
					meg->isExplainedGenesUpDown, &sampleGeneExpression, i, &geneSymbolToId);
		}

		//add explained genes to the list of all samples
		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
		mutatedGeneIdsListReal.push_back(mutatedGeneIds);

	}


	// get the frequency of genes in the real dataset
	vector<bool> isExplainedGenesUpDown(totalGenesUpDown);
	vector<int> explainedGenesFrequencyRealUpDown(totalGenes * 2);	// this will be used to compare with the frequency from random dataset
	addFrequncyForRealDataset(&explainedGenesFrequencyRealUpDown,
			&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
			&isExplainedGenesUpDown);

	//OUTPUT: print all modules in all samples (mutatedAndExplainedGenesListReal) (for Cytoscape)
	vector<string>* outputStr = new vector<string>;
	outputStr->push_back("sample_id\tgene_symbol\ttype");
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];

		//for each mutated genes
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			outputStr->push_back(intToStr(i) + "\t" + geneIdToSymbol[currentMutatedGeneId] + "\t" + "MUTATED");
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(isExplainedGenesUpDown->at(k)){
					if(k < totalGenes){
						outputStr->push_back(sampleIdToName[i] + "\t" + geneIdToSymbol[k] + "_UP\t" + "EXPLAINED");
					}else{
						outputStr->push_back(sampleIdToName[i] + "\t" + geneIdToSymbol[k-totalGenes] + "_DOWN\t" + "EXPLAINED");
					}
				}
			}
		}
	}
	string outModuleCysfilename = "output/original_modules_cys.tsv";
	writeStrVector(outModuleCysfilename.c_str(), outputStr);
	delete outputStr;
	//OUTPUT: print all modules in all samples (as original)
	outputStr = new vector<string>;
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];

		//for each mutated genes
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			string str = sampleIdToName[i] + "\t" + geneIdToSymbol[currentMutatedGeneId] + "\t";
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
			int numMember = 0;
			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(isExplainedGenesUpDown->at(k)){
					if(k < totalGenes){
						str += geneIdToSymbol[k] + "_UP" + ";";
						numMember++;
					}else{
						str += geneIdToSymbol[k-totalGenes] + "_DOWN" + ";";
						numMember++;
					}
				}
			}
			if(numMember > 0){
				outputStr->push_back(str);
			}
		}
	}
	string outModulefilename = "output/MODULE.dat";
	writeStrVector(outModulefilename.c_str(), outputStr);
	delete outputStr;

	//TODO debug the finding phenotype part (the results are quite different from the original)
	//create a vector for counting the number of times (for each gene) the random samples have greater frequency than the real samples
	vector<int> geneFrequencyGreaterThanRealFrequencyCounter(totalGenesUpDown);
	for (int i = 0; i < totalGenesUpDown; ++i) {
		geneFrequencyGreaterThanRealFrequencyCounter[i] = 0;
	}

	//the following have to be done 500-1000 times to generate the null distribution
	int numPermutationsForPhenotypeGene = 500;

	cout << "\tcreating null distribution (using " << numPermutationsForPhenotypeGene << " permutations) ... ";

	int progress = 1;
	int interval = numPermutationsForPhenotypeGene / 100;

	for (int r = 0; r < numPermutationsForPhenotypeGene; ++r) {
		// a list for explained genes of each sample
		vector< vector<bool> > explainedGenesFrequencyUpDownRandom;
		explainedGenesFrequencyUpDownRandom.resize( totalSamples , vector<bool>( totalGenesUpDown , false) );

		//print progression
		if (r % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

		//permute the gene labels of each sample independently
		omp_set_num_threads(numThreads);
		#pragma omp parallel for
		for (int i = 0; i < totalSamples; ++i) {

			//get gene expression of a current sample
			vector<double> sampleGeneExpression(totalGenes); // to save expression of of all genes in the network
			getGeneExpressionFromSampleId(&originalGeneExpressionMatrix,
					&genesEx, &sampleGeneExpression, i);

			//permute labels
			vector<int> permutedGeneLabelsMut;
			permuteGeneLabels(&genesMut, &permutedGeneLabelsMut);
			//find mutated genes of a current sample
			vector<int> mutatedGeneIds; // to store gene id of mutated genes
			getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
					&permutedGeneLabelsMut);

			//find explained genes of a current sample
			getExplainedGenesIdOnlyUpDownIncludingMutatedGene(&explainedGenesFrequencyUpDownRandom[i], &network,
					&sampleGeneExpression, &mutatedGeneIds, L, D, F, &geneSymbolToId);

		}

		//for each gene, count the frequency that exceed the real frequency
		countGeneFrequencyGreaterThanRealFrequency(&geneFrequencyGreaterThanRealFrequencyCounter,
				&explainedGenesFrequencyUpDownRandom, &explainedGenesFrequencyRealUpDown);
	}

	cout << endl; //for print progression

	//collect phenotype genes
	vector<bool> isPhenotypeGenesUpDown(totalGenesUpDown);
	for (int i = 0; i < totalGenesUpDown; ++i) {
		isPhenotypeGenesUpDown[i] = false;
	}
	vector<int> phenotypeGeneIdsUpDown;	// phenotype gene ids
	vector<double> pValues(totalGenesUpDown);

	findPhenotypeGenesUsingCounter(&isPhenotypeGenesUpDown, &phenotypeGeneIdsUpDown, &pValues, &explainedGenesFrequencyRealUpDown,
			&geneFrequencyGreaterThanRealFrequencyCounter, &isExplainedGenesUpDown, numPermutationsForPhenotypeGene, totalSamples, &geneIdToSymbol);

	cout << "\tDONE finding phenotype genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";
	step_time = clock();	//update the clock
	cout << "\tthere are " << phenotypeGeneIdsUpDown.size() << " phenotype genes" << endl;


	//[DEBUG] use phenotype genes list from previous version
//	readGenesListUpDown("original_phenotype_genes.txt", &phenotypeGeneIds, &geneSymbolToId);
//	int numPhenotypeGene = phenotypeGeneIds.size();
//	cout << "use phonotype genes list from previous version (" << numPhenotypeGene << " genes)" << endl;
//	for (int i = 0; i < numPhenotypeGene; ++i) {
//		isPhenotypeGenesUpDown[phenotypeGeneIds[i]] = true;
//	}

	//OUTPUT: print gene frequency and phenotype genes of the real dataset
	printExplinedGenesFrequencyAndPhonotype(&explainedGenesFrequencyRealUpDown, &pValues, &isPhenotypeGenesUpDown, &geneIdToSymbol, &network, &originalGeneExpressionMatrix, &genesEx, F);

	/*
	 * Find driver genes
	 */

	//[DEBUG] use driver genes list from previous version
//	vector<int> driverGeneIds;	// to save diver gene ids
//	readGenesList("original_driver_genes.txt", &driverGeneIds, &geneSymbolToId);
//	int numDriverGenesFromFile = driverGeneIds.size();
//	cout << "\ttotal driver genes = " << numDriverGenesFromFile << endl;
//	for (int i = 0; i < numDriverGenesFromFile; ++i) {
//		DriverGene driverGene;
//		driverGene.geneId = driverGeneIds[i];
//		driverGenes.push_back(driverGene);
//	}

	/*
	 * SENSITIVE MODE
	 */

	cout << "SENSITIVE mode ...\n";


	mode = 0;

	{
		//collect all mutated genes in all samples => use mutatedGeneIdsListReal (samples, mutated genes, explained genes)
		list<int> mutatedGeneIdsList;	//mutated gene ids
		vector<bool> isMutatedGenes(totalGenes);
		getAllMutatedGenes(&mutatedGeneIdsListReal, &isMutatedGenes, &mutatedGeneIdsList);
		cout << "\ttotal number of mutated genes is " << mutatedGeneIdsList.size() << endl;

		vector< list<Module> > modulesListOfAllSamples(totalSamples);	//to save all modules in all samples

		//create bipartite graph ( mutated gene --- phenotype gene ). This is done at sample level, so have to remember sample id.
		//BipartiteEdge: each mutated gene contains a pair of (phenotype gene id, sample id)
		vector<BipartiteEdge>* bipartiteGraph = new vector<BipartiteEdge>(totalGenes);
		createBipartiteGraph(&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
				&isPhenotypeGenesUpDown, bipartiteGraph, &geneIdToSymbol);

		cout << "\tperforming greedy minimum set cover algorithm ...\n";
		vector<DriverGene> driverGenes;	//driver gene id with sample ids
		cout << "\tfinding driver genes ...\n";
		findDriverGenes(bipartiteGraph, &mutatedGeneIdsList, &driverGenes);
		delete bipartiteGraph;
		cout << "\tDONE finding driver genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC)
				<< " sec)\n";
		step_time = clock();	//update the clock

		int numDriverGenes = driverGenes.size();
		cout << "\ttotal driver genes = " << numDriverGenes << endl;

		cout << "\tgenerating modules for all samples ...\n";

		//merge modules and trim explained genes for each sample
		//use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genesUpDown)
		findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
				&mutatedGeneIdsListReal, &isPhenotypeGenesUpDown, &driverGenes, &phenotypeGeneIdsUpDown, mode);

	//	string outMergedModulesCysFilename = "output/merged_modules_cys.tsv";
	//	saveModulesCytoscape(&modulesListOfAllSamples, outMergedModulesCysFilename, &geneIdToSymbol);

		cout << "\ttrimming explained genes for all samples ...\n";
		trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D, &geneIdToSymbol);

	//	string outFinalModulesCysFilename = "output/merged_modules_cys.tsv";
	//	saveModulesCytoscape(&modulesListOfAllSamples, outFinalModulesCysFilename, &geneIdToSymbol);

		cout << "\twriting final module to FINAL_MODULE.dat ...\n";
		string outFinalModuleFilenameSensitive = "output/sensitive/FINAL_MODULE.dat";
		saveModules(&modulesListOfAllSamples, &mutatedAndExplainedGenesListReal, outFinalModuleFilenameSensitive, &geneIdToSymbol, &sampleIdToName);

		cout << "\tcalculating IMPACT scores for all samples ...\n";

		{
		vector< vector<Driver> > driversOfAllSamples(totalSamples);
		calculateImpactScoresForAllSamples(&modulesListOfAllSamples, &driversOfAllSamples, &originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol, mode);

		vector<double> driverAggregatedScores(totalGenes);
		vector<int> driversFrequency(totalGenes);
		vector<int> pointMutationDriversFrequency(totalGenes);
		vector<int> deletionDriversFrequency(totalGenes);
		vector<int> amplificationDriversFrequency(totalGenes);
		vector<int> mutationFrequency(totalGenes);
		vector<int> pointMutationFrequency(totalGenes);
		vector<int> deletionFrequency(totalGenes);
		vector<int> amplificationFrequency(totalGenes);
		//initialization
		for (int i = 0; i < totalGenes; ++i) {
			driverAggregatedScores[i] = 0;
			driversFrequency[i] = 0;
			pointMutationDriversFrequency[i] = 0;
			deletionDriversFrequency[i] = 0;
			amplificationDriversFrequency[i] = 0;
			mutationFrequency[i] = 0;
			pointMutationFrequency[i] = 0;
			deletionFrequency[i] = 0;
			amplificationFrequency[i] = 0;
		}

		cout << "\taggregating IMPACT scores across all samples ...\n";
		aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol, totalGenes);

//		cout << "getting driver frequency ...\n";
		getDetailDriversFreqeuncy(&driversOfAllSamples,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&originalPointMutationsMatrix, &originalCNVsMatrix,
				&genesPointMut, &genesCNV);

//		cout << "getting mutation frequency ...\n";
		getMutationFrequency(&originalMutationMatrix, &mutationFrequency, &genesMut);
		getDetailMutationFrequency(&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

		string outSampleResultDirName = "output/sensitive/samples/";
//		cout << "printing impact scores for all samples ...\n";
		printSampleDriverList(&driversOfAllSamples, outSampleResultDirName, &geneIdToSymbol, &sampleIdToName,
				&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency);

		cout << "\tprinting aggregated impact scores ...\n";
		string outDriverListfilename = "output/sensitive/driver_list.txt";
		printAggregatedDriverList(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

		}

	}

	/*
	 * STRINGENT MODE
	 */

	cout << "STRINGENT mode ...\n";
	mode = 1;

	{
		//collect all mutated genes in all samples => use mutatedGeneIdsListReal (samples, mutated genes, explained genes)
		list<int> mutatedGeneIdsList;	//mutated gene ids
		vector<bool> isMutatedGenes(totalGenes);
		getAllMutatedGenes(&mutatedGeneIdsListReal, &isMutatedGenes, &mutatedGeneIdsList);
		cout << "\ttotal number of mutated genes is " << mutatedGeneIdsList.size() << endl;

		vector< list<Module> > modulesListOfAllSamples(totalSamples);	//to save all modules in all samples

		//create bipartite graph ( mutated gene --- phenotype gene ). This is done at sample level, so have to remember sample id.
		//BipartiteEdge: each mutated gene contains a pair of (phenotype gene id, sample id)
		vector<BipartiteEdge>* bipartiteGraph = new vector<BipartiteEdge>(totalGenes);
		createBipartiteGraph(&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
				&isPhenotypeGenesUpDown, bipartiteGraph, &geneIdToSymbol);

		cout << "\tperforming greedy minimum set cover algorithm ...\n";
		vector<DriverGene> driverGenes;	//driver gene id with sample ids
		cout << "\tfinding driver genes ...\n";
		findDriverGenes(bipartiteGraph, &mutatedGeneIdsList, &driverGenes);
		cout << "\tcalculated driver genes from stringent mode\n";
		delete bipartiteGraph;
		cout << "\tDONE finding driver genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC)
				<< " sec)\n";
		step_time = clock();	//update the clock

		int numDriverGenes = driverGenes.size();
		cout << "\ttotal driver genes = " << numDriverGenes << endl;

		cout << "\tgenerating modules for all samples ...\n";

		//merge modules and trim explained genes for each sample
		//use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genesUpDown)
		findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
				&mutatedGeneIdsListReal, &isPhenotypeGenesUpDown, &driverGenes, &phenotypeGeneIdsUpDown, mode);

	//	string outMergedModulesCysFilename = "output/merged_modules_cys.tsv";
	//	saveModulesCytoscape(&modulesListOfAllSamples, outMergedModulesCysFilename, &geneIdToSymbol);

		cout << "\ttrimming explained genes for all samples ...\n";
		trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D, &geneIdToSymbol);

	//	string outFinalModulesCysFilename = "output/merged_modules_cys.tsv";
	//	saveModulesCytoscape(&modulesListOfAllSamples, outFinalModulesCysFilename, &geneIdToSymbol);

		cout << "\twriting final module to FINAL_MODULE.dat ...\n";
		string outFinalModuleFilenameSensitive = "output/stringent/FINAL_MODULE.dat";
		saveModules(&modulesListOfAllSamples, &mutatedAndExplainedGenesListReal, outFinalModuleFilenameSensitive, &geneIdToSymbol, &sampleIdToName);

		cout << "\tcalculating IMPACT scores for all samples ...\n";

		vector< vector<Driver> > driversOfAllSamples(totalSamples);
		calculateImpactScoresForAllSamples(&modulesListOfAllSamples, &driversOfAllSamples, &originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol, mode);

		vector<double> driverAggregatedScores(totalGenes);
		vector<int> driversFrequency(totalGenes);
		vector<int> pointMutationDriversFrequency(totalGenes);
		vector<int> deletionDriversFrequency(totalGenes);
		vector<int> amplificationDriversFrequency(totalGenes);
		vector<int> mutationFrequency(totalGenes);
		vector<int> pointMutationFrequency(totalGenes);
		vector<int> deletionFrequency(totalGenes);
		vector<int> amplificationFrequency(totalGenes);
		//initialization
		for (int i = 0; i < totalGenes; ++i) {
			driverAggregatedScores[i] = 0;
			driversFrequency[i] = 0;
			pointMutationDriversFrequency[i] = 0;
			deletionDriversFrequency[i] = 0;
			amplificationDriversFrequency[i] = 0;
			mutationFrequency[i] = 0;
			pointMutationFrequency[i] = 0;
			deletionFrequency[i] = 0;
			amplificationFrequency[i] = 0;
		}

		cout << "\taggregating IMPACT scores across all samples ...\n";
		aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol, totalGenes);

//		cout << "getting driver frequency ...\n";
		getDetailDriversFreqeuncy(&driversOfAllSamples,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&originalPointMutationsMatrix, &originalCNVsMatrix,
				&genesPointMut, &genesCNV);

//		cout << "getting mutation frequency ...\n";
		getMutationFrequency(&originalMutationMatrix, &mutationFrequency, &genesMut);
		getDetailMutationFrequency(&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

		string outSampleResultDirName = "output/stringent/samples/";
//		cout << "printing impact scores for all samples ...\n";
		printSampleDriverList(&driversOfAllSamples, outSampleResultDirName, &geneIdToSymbol, &sampleIdToName,
				&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency);

		cout << "\tprinting aggregated impact scores ...\n";
		string outDriverListfilename = "output/stringent/driver_list.txt";
		printAggregatedDriverList(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

	}

	//delete the vector<int>* explainedGenesFreqency
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
			delete isExplainedGenesUpDown;
		}
	}

	/*
	 * Stop timer
	 */

	cout << "DONE (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

	return 0;
}

