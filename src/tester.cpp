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
#include "utilities.h"
#include "input.h"
#include "sampling.h"
#include "explained_genes.h"
#include "phenotype_genes.h"
#include "driver_genes.h"

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main() {

	/*
	 * Start timer
	 */

	const clock_t begin_time = clock();

	/*
	 * Read configuration from a file
	 */

	map<string, string> conf;

	/*
	 * Read network file into adjacency list
	 */

	cout << "reading network file and create mapping <id, geneSymbol>..."
			<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	string filename = "network_FIsInGene_041709.txt";
	readNetwork(filename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId,
			'\t');
	int totalGenes = network.size();
	//printAdjacencyList(&network);
	//int toSearch = geneSymbolToId.find("EGFR")->second;
	//BFS(&network, toSearch);
	//DFS(&network, toSearch);

	/*
	 * Read gene expression matrix from file
	 * (row = gene, col = samples)
	 */

	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix..." << endl;
	filename = "GBM/EXPR.txt";
	vector<int> genesEx;
	vector<string> geneSymbols;
	GeneExpression geneExpression;
	geneExpression.genes = &genesEx;
	geneExpression.matrix = &originalGeneExpressionMatrix;
	readGeneExpression(filename.c_str(), &geneExpression, '\t',
			&geneSymbolToId);

	//traceMatrix(&originalGeneExpressionMatrix);

	int totalSamples = originalGeneExpressionMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesEx = originalGeneExpressionMatrix.size();

	cout << "\ttotal genes in gene expression matrix is " << numGenesEx << endl;
	cout << "\ttotal samples in gene expression matrix is " << totalSamples
			<< endl;

	/*
	 * Read mutation matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of sample is the same as in gene expression matrix (including order)
	 */

	cout << "reading mutation matrix..." << endl;
	filename = "GBM/SNP.txt";
	vector<int> genesMut;
	TIntegerMatrix originalMutationMatrix;
	Mutations mutations;
	mutations.genes = &genesMut;
	mutations.matrix = &originalMutationMatrix;
	readMutations(filename.c_str(), &mutations, '\t', &geneSymbolToId);

	totalSamples = originalMutationMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesMut = originalMutationMatrix.size();

	cout << "\ttotal genes in gene mutation matrix is " << numGenesMut << endl;
	cout << "\ttotal samples in gene mutation matrix is " << totalSamples
			<< endl;

	/*
	 * Sample statistics
	 */

	// TEST count differentially expressed genes on each sample
	int deGenes;
	int sId = 0;
	double minFoldChange = 2;
	deGenes = countDifferentiallyExpressedGeneForSampleId(
			&originalGeneExpressionMatrix, sId, minFoldChange);
	cout << "Sample #0 has " << deGenes << " differentially expressed genes" << endl;

	/*
	 * Randomly choose sub-sample for tuning parameters
	 */

	//Note gene expression and mutation have the same set of samples
	int numSamples = 50;
	cout << "randomly choosing " << numSamples << " samples..." << endl;

	//list of samples id to be chose
	vector<int> rrank(totalSamples);
	createPermutation(&rrank);

	//gene expression matrix
	TDoubleMatrix subGeneExpressionMatrix;
	GeneExpression subGeneExpression;
	subGeneExpression.genes = &genesEx;
	subGeneExpression.matrix = &subGeneExpressionMatrix;
	randomlyChooseSamplesDouble(&originalGeneExpressionMatrix,
			&subGeneExpressionMatrix, &rrank, numSamples);

	//mutation matrix
	TIntegerMatrix subMutationMatrix;
	Mutations subMutations;
	subMutations.genes = &genesMut;
	subMutations.matrix = &subMutationMatrix;
	randomlyChooseSamplesInteger(&originalMutationMatrix, &subMutationMatrix,
			&rrank, numSamples);

	/*
	 * Create gene label permutation for both gene expression and mutation matrix
	 * Output: permutedGeneLabelsEx and permutedGeneLabelsMut
	 */

	//1. gene expression
	vector<int> permutedGeneLabelsEx;
	permuteGeneLabels(&genesEx, &permutedGeneLabelsEx);

	//2. mutation
	vector<int> permutedGeneLabelsMut;
	permuteGeneLabels(&genesMut, &permutedGeneLabelsMut);

	/*
	 * Find explained genes
	 */

	cout << "finding explained genes..." << endl;

	/*
	 * Find explained genes: set parameters
	 */

	int minL = 2; // to be used
	int maxL = 5; // to be used
	int L = 4;
	int D = 150;
	double F = 2;

	cout << "\tcurrent parameters (L, D, F) is " << L << ", " << D << ", " << F << endl;

	/*
	 * find explained genes for real sub-sample (without gene label permutation)
	 */

	vector<vector<ExplainedGene> > explainedGenesListReal;
	int sampleId = 0; //the first sample
	for (; sampleId < numSamples; sampleId++) {
		//cout << "Sample #" << sampleId << endl;
		vector<double> sampleGeneExpression(totalGenes);	//expression of all genes in the network
		getGeneExpressionFromSampleId(subGeneExpression.matrix, &genesEx,
				&sampleGeneExpression, sampleId, &geneIdToSymbol);

		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&subMutations, &mutatedGeneIds, sampleId,
				&genesMut);

		//get all the explained genes for all mutations in the current sample
		vector<ExplainedGene> explainedGenes(totalGenes);
		getExplainedGenes(&explainedGenes, &network, &sampleGeneExpression,
				&mutatedGeneIds, 4, D, F);
		explainedGenesListReal.push_back(explainedGenes);
	}

	/*
	 * find explained genes for random sub-sample (with gene label permutation)
	 */

	vector<vector<ExplainedGene> > explainedGenesListRandom;
	sampleId = 0; //the first sample
	for (; sampleId < numSamples; sampleId++) {
		//cout << "Sample #" << sampleId << endl;
		vector<double> sampleGeneExpression(totalGenes);	// of all genes
		getGeneExpressionFromSampleId(subGeneExpression.matrix,
				&permutedGeneLabelsEx, &sampleGeneExpression, sampleId,
				&geneIdToSymbol);

		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&subMutations, &mutatedGeneIds, sampleId,
				&permutedGeneLabelsMut);
		//cout << "List of " << mutatedGeneIds.size() << " mutated genes\n";
		//printGeneSymbols(&mutatedGeneIds, &geneIdToSymbol);

		vector<ExplainedGene> explainedGenes(totalGenes);
		getExplainedGenes(&explainedGenes, &network, &sampleGeneExpression,
				&mutatedGeneIds, 4, D, F);
		explainedGenesListRandom.push_back(explainedGenes);
	}

	/*
	 * Find phenotype genes
	 */

	cout << "finding phenotype genes..." << endl;


	// for collecting mutated genes and the corresponding explained genes for all samples
	// sample, mutated genes, explained genes
	vector< vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;

	for (int i = 0; i < totalSamples; ++i) {
//		cout << "\tSample #" << i << endl;

		//get gene expression of a current sample
		vector<double> sampleGeneExpression(totalGenes);// to save expression of of all genes in the network
		getGeneExpressionFromSampleId(&originalGeneExpressionMatrix,
				&genesEx, &sampleGeneExpression, i, &geneIdToSymbol);
//		cout << "\tget gene expression" << endl;

		//find mutated genes of a current sample
		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
				&genesMut);
//		cout << "\tget mutated genes" << endl;

		//find explained genes of a current sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes;
		getMutatedAndExplainedGenes(&mutatedAndExplainedGenes, &network, &sampleGeneExpression,
				&mutatedGeneIds, 4, D, F);
//		cout << "\tget mutated genes and their corresponding explained genes" << endl;

		//add explained genes to the list of all samples
		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
	}

	//cout << "Finished explained genes for = " << mutatedAndExplainedGenesListReal.size() << endl;

	// get the frequency of genes in the real dataset
	vector<int> genesFrequencyReal(totalGenes);
	addFrequncyForRealDataset(&genesFrequencyReal, &mutatedAndExplainedGenesListReal);


	// the following have to be done 500 times to generate the null distribution
	vector<vector<int> > nullDistribution(totalGenes);
	int round = 100;

	cout << "creating null distribution for each gene...";
	int progress = 1;
	int interval = round / 100;

	for (int r = 0; r < round; ++r) {
		// a list for explained genes of each sample
		vector< vector< int > > explainedGenesListForPhenotypeGenes;

		//print progression
		if (r % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

		// permute the gene labels of each sample independently
		for (int i = 0; i < totalSamples; ++i) {

//			cout << "For round " << r << " sample #" << i << endl;

			//permute labels
			vector<int> permutedGeneLabelsMut;
			permuteGeneLabels(&genesMut, &permutedGeneLabelsMut);
			//cout << "\tpermuted gene label" << endl;

			//get gene expression of a current sample
			vector<double> sampleGeneExpression(totalGenes);// to save expression of of all genes in the network
			getGeneExpressionFromSampleId(&originalGeneExpressionMatrix,
					&genesEx, &sampleGeneExpression, i, &geneIdToSymbol);
			//cout << "\tcopied gene expression value" << endl;

			//find mutated genes of a current sample
			vector<int> mutatedGeneIds; // to store gene id of mutated genes
			getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
					&permutedGeneLabelsMut);
			//cout << "\tfound " << mutatedGeneIds.size() << " mutated genes" << endl;

			//find explained genes of a current sample
			vector<int> explainedGeneIds(totalGenes);
			getExplainedGenesOnlyId(&explainedGeneIds, &network, &sampleGeneExpression,
					&mutatedGeneIds, 4, D, F);
			//cout << "\tfound " << explainedGenes.size() << " explained genes" << endl;

			//add explained genes to the list of all samples
			explainedGenesListForPhenotypeGenes.push_back(explainedGeneIds);
		}

		// add the frequency of a current round to generate null distribution (to be added 500 times)
		addFrequencyForNullDistribution(&nullDistribution,
				&explainedGenesListForPhenotypeGenes);

	}

	cout << endl;

	// collect phenotype genes
	//consider only genes available in the gene expression matrix
	vector<bool>* isInGeneExpressionMatrix = new vector<bool>(totalGenes);
	int numGeneEx = genesEx.size();
	for (int i = 0; i < numGeneEx; ++i) {
		isInGeneExpressionMatrix->at(genesEx[i]) = true;
	}
	vector<bool> isPhenotypeGenes(totalGenes);	// 1 for yes 0 for no
	findPhenotypeGenes(&isPhenotypeGenes, &genesFrequencyReal, &nullDistribution, isInGeneExpressionMatrix);
	delete isInGeneExpressionMatrix;
	vector<int> pg;
	for (int i = 0; i < totalGenes; ++i) {
		if(isPhenotypeGenes[i] == 1){
			pg.push_back(i);
		}
	}
	cout << "There are " << pg.size() << " phenotype genes" << endl;
//	printGeneSymbols(&pg, &geneIdToSymbol);


	/*
	 * Find driver genes
	 */

	// map phenotype genes to sample level
	// just use isPhenotypeGenes to check

	// collect all mutated genes in all samples
	// just use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genes)
	vector<bool>* isAMutatedGenesInSamples = new vector<bool>(totalGenes);
	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {
		//get all mutated genes and their corresponding explained genes
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		int numMutatedGenes = mutatedAndExplainedGenes.size();

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			isAMutatedGenesInSamples->at(mutatedAndExplainedGenes[j].mutatedGeneId) = true;
		}
	}


	//create bipartite graph ( mutated gene --- phenotype gene ). This is done on sample level, so have to remember sample id.

	vector<BipartiteEdge>* bipartiteEdges = new vector<BipartiteEdge>[totalGenes];

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {
		//get all mutated genes and their corresponding explained genes
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		int numMutatedGenes = mutatedAndExplainedGenes.size();

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			vector<int> explainedGeneIds = mutatedAndExplainedGenes[j].explainedGeneIds;
			int mutatedGeneId = mutatedAndExplainedGenes[j].mutatedGeneId;
			int numExplainedGenes = explainedGeneIds.size();

			//for each explained gene k
			for (int k = 0; k < numExplainedGenes; ++k) {
				//if the explained gene k is a phenotype gene, then add to the bipartite graph
				if(isPhenotypeGenes[explainedGeneIds[k]]){
					BipartitePhenotypeNode node;
					node.phenotypeGeneId = explainedGeneIds[k];
					node.sampleId = i;
					bipartiteEdges->at(mutatedGeneId).phenotypeGeneIdsAndSampleIds.push_back(node);
				}
			}
		}
	}

	//greedy minimum set covering



	delete isAMutatedGenesInSamples;
	delete bipartiteEdges;



	/*
	 * Stop timer
	 */

	cout << (float(clock() - begin_time) / CLOCKS_PER_SEC) << " sec\n";

	return 0;
}
