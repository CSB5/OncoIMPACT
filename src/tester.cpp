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

	/*
	 * TODO Read configuration from a file
	 */

	map<string, string> conf;

	/*
	 * Read network file into adjacency list
	 */

	cout << "reading network file and create mapping <id, geneSymbol> ..."
			<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	string filename = "network_FIsInGene_041709.txt";
	readNetwork(filename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId,
			'\t');
	int totalGenes = network.size();

	/*
	 * Read gene expression matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of samples is the same as in mutation matrix (also the same order)
	 */

	//sample id mapping
	vector<string> sampleIdToName;

	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix ..." << endl;
	filename = "GBM/EXPR.txt";
	vector<int> genesEx; // gene ids ; size = # DE genes (consider only genes that are in the network)
	vector<string> geneSymbols;
	GeneExpression geneExpression;	//all info of gene expression
	geneExpression.genes = &genesEx;
	geneExpression.matrix = &originalGeneExpressionMatrix;
	readGeneExpression(filename.c_str(), &geneExpression, '\t',
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
	filename = "GBM/SNP.txt";
	vector<int> genesPointMut;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalPointMutationsMatrix;
	PointMutations pointMutations;	//all info of point mutation
	pointMutations.genes = &genesPointMut;
	pointMutations.matrix = &originalPointMutationsMatrix;
	readPointMutations(filename.c_str(), &pointMutations, '\t', &geneSymbolToId);

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
	filename = "GBM/CNV.txt";
	vector<int> genesCNV;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalCNVsMatrix;
	CopyNumberVariation CNVs;	//all info of CNV
	CNVs.genes = &genesCNV;
	CNVs.matrix = &originalCNVsMatrix;
	readCopyNumberVariation(filename.c_str(), &CNVs, '\t', &geneSymbolToId);

	totalSamples = originalCNVsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesCNV = originalCNVsMatrix.size();

	cout << "\ttotal genes in CNV matrix is " << numGenesCNV << endl;
	cout << "\ttotal samples in CNV matrix is " << totalSamples
			<< endl;

	//TODO separate the combining code to a function in input.cpp
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
	 * Randomly choose sub-sample for tuning parameters
	 */

	//Note: gene expression and mutation have the same set of samples
	//If number of samples is < 50, use all samples to tune the parameters
	int numSamples = 50;
	if(totalSamples < 50){
		numSamples = totalSamples;
	}
	cout << "tuning parameters by using " << numSamples << " randomly choosing samples ..."
			<< endl;

	cout << "computing JS divergence for all parameters (L,D,F) ... " << endl;

	//initialize all the parameters to be tested
//	int LsVal[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
	int LsVal[] = {14, 16, 18, 20};	//fewer parameters for testing
	vector<int> Ls(LsVal, LsVal + sizeof LsVal / sizeof LsVal[0]);
	int numLs = Ls.size();
//	int DsVal[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
	int DsVal[] = {55, 60, 65, 70};	//fewer parameters for testing
	vector<int> Ds(DsVal, DsVal + sizeof DsVal / sizeof DsVal[0]);
	int numDs = Ds.size();	//fewer parameters for testing
//	double FsVal[] = {1, 1.5, 2, 2.5, 3};
	double FsVal[] = {2, 2.5};	//fewer parameters for testing
	vector<double> Fs(FsVal, FsVal + sizeof FsVal / sizeof FsVal[0]);
	int numFs = Fs.size();


	int numCombinations = numLs * numDs * numFs;
	//initialize the vector to save the divergence of each set of parameters
	vector<JSDivergence> jsDivergences(numCombinations);

	findParameters(&jsDivergences, &Ls, &Ds, &Fs, totalGenes, &geneExpression, &mutations, &network, numSamples);

	//TODO write the JS divergence result to a file
	filename = "output/js_divergences";
	saveJSDivergences(&jsDivergences, filename);

	cout << "DONE tunning parameters (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";
	begin_time = clock();	//update the clock

	//choose the best parameters
	JSDivergence maxJs;
	findMaximumJsDivergence(&jsDivergences, &maxJs);

	cout << "the maximum divergence is " << maxJs.divergence << " when L, D, F = " << maxJs.L << ", " << maxJs.D << ", " << maxJs.F << endl;

	//TODO set the L D F to maxJs
	int L = 20;
	int D = 65;
	double F = 2.5;

	/*
	 * Find phenotype genes
	 */

	cout << "finding phenotype genes using (L,D,F) = (" << L << "," << D << "," << F << ") ..." << endl;

	// A vector for collecting mutated genes and the corresponding explained genes for all samples
	// sample, mutated genes, explained genes
	vector< vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;
	//Note: size of MutatedAndExplianedGenes = totalGenes, so index = mutated gene id, value = a list of explained gene ids
	vector< vector<int> > mutatedGeneIdsListReal;	//to save a list of mutated gene ids for each sample
	vector<bool> isExplainedGenes(totalGenes);

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
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes(totalGenes);
		for (int j = 0; j < numMutatedGenes; ++j) {	// for each mutated genes
			int mutatedGeneId = mutatedGeneIds[j];
			MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes[mutatedGeneId];
			meg->explainedGenesFreqency = new vector<int>(totalGenes);
			//BFS for explained genes of the current mutated gene
			BFSforExplainedGenesIdOnly(&network, mutatedGeneId, L, D, F,
					meg->explainedGenesFreqency, &sampleGeneExpression);
		}

		//add explained genes to the list of all samples
		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
		mutatedGeneIdsListReal.push_back(mutatedGeneIds);

	}

	//OUTPUT: print all modules in all samples (mutatedAndExplainedGenesListReal)
	vector<string>* outputStr = new vector<string>;
	outputStr->push_back("sample_id\tgene_symbol\ttype");
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			outputStr->push_back(intToStr(i) + "\t" + geneIdToSymbol[currentMutatedGeneId] + "\t" + "MUTATED");
			vector<int>* explainedGenesFreqency = mutatedAndExplainedGenes[currentMutatedGeneId].explainedGenesFreqency;
			for (int k = 0; k < totalGenes; ++k) {
				if(explainedGenesFreqency->at(k) > 0){
					outputStr->push_back(intToStr(i) + "\t" + geneIdToSymbol[k] + "\t" + "EXPLAINED");
				}
			}
		}
	}
	filename = "output/raw_modules.tsv";
	writeStrVector(filename.c_str(), outputStr);
	delete outputStr;

	// get the frequency of genes in the real dataset
	vector<int> genesFrequencyReal(totalGenes);
	addFrequncyForRealDataset(&genesFrequencyReal,
			&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
			&isExplainedGenes);


	//OUTPUT: print gene frequency of the real dataset
	vector<string>* outputRealGenesFrequency = new vector<string>;
	outputRealGenesFrequency->push_back("GENE\tFREQUENCY");
	for (int i = 0; i < totalGenes; ++i) {
		string str = geneIdToSymbol[i] + "\t" + intToStr(genesFrequencyReal[i]);
		outputRealGenesFrequency->push_back(str);
	}
	filename = "output/genes_frequency.tsv";
	writeStrVector(filename.c_str(), outputRealGenesFrequency);
	delete outputRealGenesFrequency;

	//create a vector for counting the number of times (for each gene)
	//the random samples have greater frequency than the real dataset
	vector<int> geneFrequencyGreaterThanRealFrequencyCounter(totalGenes);

	//the following have to be done 500-1000 times to generate the null distribution
	int round = 500;

	cout << "\tcreating null distribution (using " << round << " permutations) ... ";

	int progress = 1;
	int interval = round / 100;

	for (int r = 0; r < round; ++r) {
		// a list for explained genes of each sample
		vector<vector<int> > explainedGenesFrequencyForPhenotypeGenes;

		//print progression
		if (r % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

		// permute the gene labels of each sample independently
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
			vector<int> explainedGenesFrequency(totalGenes);
			getExplainedGenesIdOnly(&explainedGenesFrequency, &network,
					&sampleGeneExpression, &mutatedGeneIds, L, D, F);

			//add explained genes to the list of all samples
			explainedGenesFrequencyForPhenotypeGenes.push_back(explainedGenesFrequency);
		}

		//for each gene, count the frequency that exceed the real frequency
		countGeneFrequencyGreaterThanRealFrequency(&geneFrequencyGreaterThanRealFrequencyCounter,
				&explainedGenesFrequencyForPhenotypeGenes, &genesFrequencyReal);
	}

	cout << endl; //for print progression

	cout << "DONE finding phenotype genes (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";
	begin_time = clock();	//update the clock

	//collect phenotype genes
	vector<bool> isPhenotypeGenes(totalGenes);	// 1 for yes 0 for no
	vector<int> phenotypeGeneIds;	// phenotype gene ids

	findPhenotypeGenesUsingCounter(&isPhenotypeGenes, &phenotypeGeneIds, &genesFrequencyReal,
			&geneFrequencyGreaterThanRealFrequencyCounter, &isExplainedGenes, round, totalSamples, &geneIdToSymbol);

	cout << "\tthere are " << phenotypeGeneIds.size() << " phenotype genes";

	/*
	 * Find driver genes
	 */

	cout << "finding driver genes ...\n";

	// map phenotype genes to sample level (use isPhenotypeGenes to check)

	// collect all mutated genes in all samples
	// (use mutatedGeneIdsListReal (samples, mutated genes, explained genes)
	list<int> mutatedGeneIdsList;	//mutated gene ids
	vector<bool> isMutatedGenes(totalGenes);
	getAllMutatedGenes(&mutatedGeneIdsListReal, &isMutatedGenes, &mutatedGeneIdsList);
	cout << "\ttotal number of mutated genes is " << mutatedGeneIdsList.size() << endl;

	cout << "\tcreating bipartite graph ...\n";
	//create bipartite graph ( mutated gene --- phenotype gene ). This is done at sample level, so have to remember sample id.
	//BipartiteEdge: each mutated gene contains a pair of (phenotype gene id, sample id)
	vector<BipartiteEdge>* bipartiteGraph = new vector<BipartiteEdge>(totalGenes);
	createBipartiteGraph(&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal, &isPhenotypeGenes, bipartiteGraph);

	//greedy minimum set covering
	cout << "\tperforming greedy minimum set cover algorithm ...\n";
	vector<int> driverGeneIds;	// to save diver gene ids
	findDriverGenes(bipartiteGraph, &mutatedGeneIdsList, &driverGeneIds);

	int numDriverGenes = driverGeneIds.size();
	cout << "\ttotal number of driver genes = " << numDriverGenes << endl;
	filename = "output/driver_genes.tsv";
	saveGeneSymbols(filename.c_str(), &driverGeneIds, &geneIdToSymbol);

	delete bipartiteGraph;

	cout << "DONE finding driver genes (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";
	begin_time = clock();	//update the clock

	cout << "merging modules for all samples ...\n";

	//merge modules and trim explained genes for each sample
	vector<bool> isDriverGenes(totalGenes);
	for (int i = 0; i < numDriverGenes; ++i) {
		isDriverGenes[driverGeneIds[i]] = true;
	}

	vector< list<Module> > modulesListOfAllSamples(totalSamples);	//to save all modules in all samples
	findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
			&mutatedGeneIdsListReal, &isPhenotypeGenes, &isDriverGenes, &phenotypeGeneIds);

	filename = "output/merged_modules.tsv";
	saveModules(&modulesListOfAllSamples, filename, &geneIdToSymbol);

	cout << "trimming explained genes for all samples ...\n";
	//TODO FIX BUGS of trimming modules
	trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D, &geneIdToSymbol);

	//TODO print sample modules
	filename = "output/trimmed_modules.tsv";
	saveModules(&modulesListOfAllSamples, filename, &geneIdToSymbol);

	cout << "calculating IMPACT scores for all samples ...\n";

	vector< vector<Driver> > driversOfAllSamples(totalSamples);
	calculateImpactScoresForAllSamples(&modulesListOfAllSamples, &driversOfAllSamples, &originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol);

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

	cout << "aggregating IMPACT scores across all samples ...\n";
	aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol, totalGenes);

	cout << "getting driver frequency ...\n";
	getDetailDriversFreqeuncy(&driversOfAllSamples,
			&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
			&originalPointMutationsMatrix, &originalCNVsMatrix,
			&genesPointMut, &genesCNV);

	cout << "getting mutation frequency ...\n";
	getMutationFrequency(&originalMutationMatrix, &mutationFrequency, &genesMut);
	getDetailMutationFrequency(&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
			&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

	filename = "output/samples/";
	cout << "printing impact scores for all samples ...\n";
	printSampleDriverList(&driversOfAllSamples, filename, &geneIdToSymbol, &sampleIdToName,
			&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
			&driverAggregatedScores, &driversFrequency, &mutationFrequency);

	cout << "printing aggregated impact scores ...\n";
	filename = "output/drivers_aggregation.tsv";
	printAggregatedDriverList(&driverGeneIds, filename, &geneIdToSymbol, &sampleIdToName,
			&driverAggregatedScores, &driversFrequency, &mutationFrequency,
			&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
			&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

	//delete the vector<int>* explainedGenesFreqency
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			vector<int>* explainedGenesFreqency = mutatedAndExplainedGenes[currentMutatedGeneId].explainedGenesFreqency;
			delete explainedGenesFreqency;
		}
	}

	/*
	 * Stop timer
	 */

	cout << "DONE (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

	return 0;
}
