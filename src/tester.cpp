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
	//cout << "Sample #0 has " << deGenes << " differentially expressed genes" << endl;

	/*
	 * Randomly choose sub-sample for tuning parameters
	 */

	//Note gene expression and mutation have the same set of samples
	int numSamples = 50;
	cout << "tuning parameters by using " << numSamples << " samples ..."
			<< endl;

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
	 * Find explained genes: parameters setting
	 */

	int minL = 2; // to be used
	int maxL = 5; // to be used
	int L = 4;
	int D = 150;
	double F = 2;

	/*
	 * Find explained genes
	 */

	cout << "\current parameters (L, D, F) is " << L << ", " << D << ", " << F
			<< endl;
	//cout << "finding explained genes ..." << endl;

	/*
	 * find explained genes for real sub-sample (without gene label permutation)
	 */

	vector<vector<ExplainedGene> > explainedGenesListReal;
	int sampleId = 0; //the first sample
	for (; sampleId < numSamples; sampleId++) {
		//cout << "Sample #" << sampleId << endl;
		vector<double> sampleGeneExpression(totalGenes);//expression of all genes in the network
		getGeneExpressionFromSampleId(subGeneExpression.matrix, &genesEx,
				&sampleGeneExpression, sampleId, &geneIdToSymbol);

		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&subMutations, &mutatedGeneIds, sampleId,
				&genesMut);

		//get all the explained genes for all mutations in the current sample
		vector<ExplainedGene> explainedGenes(totalGenes);
		getExplainedGenes(&explainedGenes, &network, &sampleGeneExpression,
				&mutatedGeneIds, L, D, F);
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
				&mutatedGeneIds, L, D, F);
		explainedGenesListRandom.push_back(explainedGenes);
	}

	/*
	 * Find phenotype genes
	 */

	cout << "finding phenotype genes ..." << endl;

	// for collecting mutated genes and the corresponding explained genes for all samples
	// sample, mutated genes, explained genes
	vector<vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;
	vector<vector<int> > mutatedGeneIdsListReal;
	vector<bool> isExplainedGenes(totalGenes);

	cout << "\tgetting explained genes frequency of the real samples ...\n";

	for (int i = 0; i < totalSamples; ++i) {

		//get gene expression of a current sample
		vector<double> sampleGeneExpression(totalGenes);// to save expression of of all genes in the network
		getGeneExpressionFromSampleId(&originalGeneExpressionMatrix, &genesEx,
				&sampleGeneExpression, i, &geneIdToSymbol);

		//find mutated genes of a current sample
		vector<int> mutatedGeneIds; // to store gene id of mutated genes
		getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
				&genesMut);

		//find explained genes of a current sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes(totalGenes);
		getMutatedAndExplainedGenes(&mutatedAndExplainedGenes, &network,
				&sampleGeneExpression, &mutatedGeneIds, 4, D, F);

		//add explained genes to the list of all samples
		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
		mutatedGeneIdsListReal.push_back(mutatedGeneIds);

	}

	// get the frequency of genes in the real dataset
	vector<int> genesFrequencyReal(totalGenes);
	addFrequncyForRealDataset(&genesFrequencyReal,
			&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
			&isExplainedGenes);

	// the following have to be done 500 times to generate the null distribution
	vector<vector<int> > nullDistribution(totalGenes);
	int round = 100;

	cout << "\tcreating null distribution (using " << round
			<< " permutations) ... ";

	int progress = 1;
	int interval = round / 100;

	for (int r = 0; r < round; ++r) {
		// a list for explained genes of each sample
		vector<vector<int> > explainedGenesListForPhenotypeGenes;

		//print progression
//		if (r % interval == 0) {
//			const string progStatus = intToStr(progress) + "%";
//			cout << progStatus << flush;
//			progress++;
//			cout << string(progStatus.length(), '\b');
//		}

		// permute the gene labels of each sample independently
		for (int i = 0; i < totalSamples; ++i) {

			//cout << "For round " << r << " sample #" << i << endl;

			//permute labels
			vector<int> permutedGeneLabelsMut;
			permuteGeneLabels(&genesMut, &permutedGeneLabelsMut);
			//cout << "\tpermuted gene label" << endl;

			//get gene expression of a current sample
			vector<double> sampleGeneExpression(totalGenes); // to save expression of of all genes in the network
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
			getExplainedGenesOnlyId(&explainedGeneIds, &network,
					&sampleGeneExpression, &mutatedGeneIds, 4, D, F);
			//cout << "\tfound " << explainedGenes.size() << " explained genes" << endl;

			//add explained genes to the list of all samples
			explainedGenesListForPhenotypeGenes.push_back(explainedGeneIds);
		}

		// add the frequency of a current round to generate null distribution (to be added 500 times)
		addFrequencyForNullDistribution(&nullDistribution,
				&explainedGenesListForPhenotypeGenes);

	}

	cout << endl;

	//collect phenotype genes
	//consider only explained genes

	vector<bool> isPhenotypeGenes(totalGenes);	// 1 for yes 0 for no
	findPhenotypeGenes(&isPhenotypeGenes, &genesFrequencyReal,
			&nullDistribution, &isExplainedGenes);

	vector<int> phenotypeGeneIds;
	for (int i = 0; i < totalGenes; ++i) {
		if (isPhenotypeGenes[i] == 1) {
			phenotypeGeneIds.push_back(i);
		}
	}
	string phenotypeGeneFileName = "phenotype_genes.txt";
	cout << "\twriting " << phenotypeGeneIds.size() << " phenotype genes to "
			<< phenotypeGeneFileName << " ..." << endl;
	saveGeneSymbols(phenotypeGeneFileName.c_str(), &phenotypeGeneIds,
			&geneIdToSymbol);

	/*
	 * Find driver genes
	 */

	cout << "finding driver genes ...\n";

	// map phenotype genes to sample level
	// just use isPhenotypeGenes to check

	// collect all mutated genes in all samples
	// just use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genes)
	vector<bool> isMutatedGenes(totalGenes);
	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {

		//get all mutated genes and their corresponding explained genes
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
		int numMutatedGenes = mutatedGeneIds.size();

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			isMutatedGenes[mutatedGeneIds[j]] = true;
		}
	}

	cout << "\tcreating bipartite graph ...\n";
	//create bipartite graph ( mutated gene --- phenotype gene ). This is done on sample level, so have to remember sample id.

	vector<BipartiteEdge>* bipartiteEdges = new vector<BipartiteEdge>(totalGenes);

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {
		cout << "sample #" << i << endl;

		//get all mutated genes and their corresponding explained genes
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
				mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];
		int numMutatedGenes = mutatedGeneIds.size();

		cout << "# mutated genes = " << numMutatedGenes << endl;

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			vector<int> explainedGenesFrequency =
					mutatedAndExplainedGenes[mutatedGeneIds[j]].explainedGenesFreqency;

			for (int k = 0; k < totalGenes; ++k) {
				if(explainedGenesFrequency[k] > 0){
					cout << k;
					if(isPhenotypeGenes[k]){
						cout << " is phenotype gene\n";
						BipartitePhenotypeNode node;
						node.phenotypeGeneId = k;
						node.sampleId = i;
						bipartiteEdges->at(k).phenotypeGeneIdsAndSampleIds.push_back(node);
					}else{
						cout << " is not phenotype gene\n";
					}
				}
			}
		}

		if(i == 4)
			break;
	}

	//greedy minimum set covering
	cout << "\tperforming greedy minimum set cover algorithm ...\n";

	int numPhenotypeGenes = 0;
	for (int i = 0; i < totalGenes; ++i) {
		if(isMutatedGenes[i]){
			numPhenotypeGenes += bipartiteEdges->at(i).phenotypeGeneIdsAndSampleIds.size();
		}
	}

	cout << "\t\ttotal number of edges is " << numPhenotypeGenes << endl;
	//until all the phenotype genes are covered
//	while(){
//
//	}

//	delete bipartiteEdges;

	/*
	 * Stop timer
	 */

	cout << "DONE (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

	return 0;
}
