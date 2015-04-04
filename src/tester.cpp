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

	//TODO read CNV matrix and combine with point mutation matrix

	/*
	 * Sample statistics
	 */

	// TEST count differentially expressed genes on each sample
//	int deGenes;
//	int sId = 0;
//	double minFoldChange = 2;
//	deGenes = countDifferentiallyExpressedGeneForSampleId(
//			&originalGeneExpressionMatrix, sId, minFoldChange);
	//cout << "Sample #0 has " << deGenes << " differentially expressed genes" << endl;

	/*
	 * Randomly choose sub-sample for tuning parameters
	 */

	//Note gene expression and mutation have the same set of samples
	//If number of samples is < 50, use all samples to tune the parameters
	int numSamples = 50;
	if(totalSamples < 50){
		numSamples = totalSamples;
	}
	cout << "tuning parameters by using " << numSamples << " samples ..."
			<< endl;

	//TODO check if repermutation and resampling are need

	//Create gene label permutation for both gene expression and mutation matrix
	//1. gene expression
	vector<int> permutedGeneLabelsEx;
	permuteGeneLabels(&genesEx, &permutedGeneLabelsEx);
	//2. mutation
	vector<int> permutedGeneLabelsMut;
	permuteGeneLabels(&genesMut, &permutedGeneLabelsMut);

	//list of samples id to be chose
	vector<int> rrank(totalSamples);
	createPermutation(&rrank);

	//gene expression submatrix
	TDoubleMatrix subGeneExpressionMatrix;
	GeneExpression subGeneExpression;
	subGeneExpression.genes = &genesEx;
	subGeneExpression.matrix = &subGeneExpressionMatrix;
	randomlyChooseSamplesDouble(&originalGeneExpressionMatrix,
			&subGeneExpressionMatrix, &rrank, numSamples);

	//mutation submatrix
	TIntegerMatrix subMutationMatrix;
	Mutations subMutations;
	subMutations.genes = &genesMut;
	subMutations.matrix = &subMutationMatrix;
	randomlyChooseSamplesInteger(&originalMutationMatrix, &subMutationMatrix,
			&rrank, numSamples);

	//parameters setting

	cout << "computing JS divergence for all parameters (L,D,F) ... ";

//	int minL = 2; // to be used
//	int maxL = 5; // to be used

//	int Ls[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
	int Ls[] = {14, 16, 18, 20};
	int numLs = sizeof(Ls)/sizeof(*Ls);
//	int Ds[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
	int Ds[] = {55, 60, 65, 70};
	int numDs = sizeof(Ds)/sizeof(*Ds);
//	double Fs[] = {1, 1.5, 2, 2.5, 3};
	double Fs[] = {2, 2.5};
	int numFs = sizeof(Fs)/sizeof(*Fs);

	//TODO separate the following loop to parameters.cpp

	int numCombinations = numLs * numDs * numFs;
	vector<JSDivergence>* jsDivergences = new vector<JSDivergence>(numCombinations);

	int count = 0;	//count number of combinations
	for (int li = 0; li < numLs; ++li) {
		for (int di = 0; di < numDs; ++di) {
			for (int fi = 0; fi < numFs; ++fi) {

				int L = Ls[li];
				int D = Ds[di];
				double F = Fs[fi];

				cout << "\tcurrent parameters (L, D, F) is " << L << ", " << D << ", " << F << endl;
				//save values
				jsDivergences->at(count).L = L;
				jsDivergences->at(count).D = D;
				jsDivergences->at(count).F = F;

				//calculated JS divergence

				//TODO 100 iterations to compute JS divergence
				vector<int> realDistributionAll;
				vector<int> randomDistributionAll;

				int round = 100;
				for (int i = 0; i < round; ++i) {

					//find explained genes for real sub-sample (without gene label permutation)

					vector<int> realDistribution(totalGenes);
					int sampleId = 0; //the first sample
					for (; sampleId < numSamples; sampleId++) {
						//cout << "Sample #" << sampleId << endl;
						vector<double> sampleGeneExpression(totalGenes); //expression of all genes in the network
						getGeneExpressionFromSampleId(subGeneExpression.matrix,
								&genesEx, &sampleGeneExpression, sampleId,
								&geneIdToSymbol);

						vector<int> mutatedGeneIds; // to store gene id of mutated genes
						getMutatedGeneIdsFromSampleId(&subMutations,
								&mutatedGeneIds, sampleId, &genesMut);

						vector<int> explainedGenesFrequency(totalGenes);
						getExplainedGenesOnlyId(&explainedGenesFrequency,
								&network, &sampleGeneExpression,
								&mutatedGeneIds, L, D, F);

						//update real distribution
						for (int j = 0; j < totalGenes; ++j) {
							if (explainedGenesFrequency[j] > 0) {
								realDistribution[j]++;
							}
						}

					}
					realDistributionAll.insert(realDistributionAll.end(),
							realDistribution.begin(), realDistribution.end());

					//find explained genes for random sub-sample (with gene label permutation)

					vector<int> randomDistribution(totalGenes);
					sampleId = 0; //the first sample
					for (; sampleId < numSamples; sampleId++) {
						//cout << "Sample #" << sampleId << endl;
						vector<double> sampleGeneExpression(totalGenes);// of all genes
						getGeneExpressionFromSampleId(subGeneExpression.matrix,
								&permutedGeneLabelsEx, &sampleGeneExpression,
								sampleId, &geneIdToSymbol);

						vector<int> mutatedGeneIds; // to store gene id of mutated genes
						getMutatedGeneIdsFromSampleId(&subMutations,
								&mutatedGeneIds, sampleId,
								&permutedGeneLabelsMut);

						vector<int> explainedGenesFrequency(totalGenes);
						getExplainedGenesOnlyId(&explainedGenesFrequency,
								&network, &sampleGeneExpression,
								&mutatedGeneIds, L, D, F);

						//update random distribution
						for (int j = 0; j < totalGenes; ++j) {
							if (explainedGenesFrequency[j] > 0) {
								randomDistribution[j]++;
							}
						}

					}
					randomDistributionAll.insert(randomDistributionAll.end(),
							randomDistribution.begin(),
							randomDistribution.end());

				}

				jsDivergences->at(count).divergence = 0;
				count++;

			}
		}
	}

	cout << "DONE tunning parameters (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

//	jsDivergences* maxJs = new JSDivergence;
//	findMaximumJsDivergence(jsDivergences, maxJs);

	//TODO choose the best parameters
	int L = 16;
	int D = 65;
	double F = 2.5;

	delete jsDivergences;

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
				&sampleGeneExpression, &mutatedGeneIds, L, D, F);

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
	vector< vector<int> > nullDistribution(totalGenes);
	int round = 500;

	cout << "\tcreating null distribution (using " << round << " permutations) ... ";

	int progress = 1;
	int interval = round / 100;

	for (int r = 0; r < round; ++r) {
		// a list for explained genes of each sample
		vector<vector<int> > explainedGenesListForPhenotypeGenes;

		//print progression
		if (r % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

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
					&sampleGeneExpression, &mutatedGeneIds, L, D, F);
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
	vector<int> phenotypeGeneIds;
	findPhenotypeGenes(&isPhenotypeGenes, &phenotypeGeneIds, &genesFrequencyReal,
			&nullDistribution, &isExplainedGenes);
	string phenotypeGeneFileName = "phenotype_genes.txt";
	saveGeneSymbols(phenotypeGeneFileName.c_str(), &phenotypeGeneIds, &geneIdToSymbol);

	/*
	 * Find driver genes
	 */

	cout << "finding driver genes ...\n";

	// map phenotype genes to sample level (use isPhenotypeGenes to check)

	// collect all mutated genes in all samples
	// (use mutatedGeneIdsListReal (samples, mutated genes, explained genes)
	list<int> mutatedGeneIdsList;
	vector<bool> isMutatedGenes(totalGenes);
	getAllMutatedGenes(&mutatedGeneIdsListReal, &isMutatedGenes, &mutatedGeneIdsList);
	cout << "\ttotal number of mutated genes is " << mutatedGeneIdsList.size() << endl;

	cout << "\tcreating bipartite graph ...\n";
	//create bipartite graph ( mutated gene --- phenotype gene ). This is done on sample level, so have to remember sample id.
	//Edge: each mutated gene contains (phenotype gene id, sample id)
	vector<BipartiteEdge>* bipartiteEdges = new vector<BipartiteEdge>(totalGenes);
	createBipartiteGraph(&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal, &isPhenotypeGenes, bipartiteEdges);

	int numPhenotypeGenes = 0;
	for (int i = 0; i < totalGenes; ++i) {
		if(isMutatedGenes[i]){
			numPhenotypeGenes += bipartiteEdges->at(i).phenotypeGeneIdsAndSampleIds.size();
		}
	}

	//cout << "\t\ttotal number of edges is " << numPhenotypeGenes << endl;

	//greedy minimum set covering
	cout << "\tperforming greedy minimum set cover algorithm ...\n";
	vector<int> driverGeneIds;
	findDriverGenes(bipartiteEdges, &mutatedGeneIdsList, &driverGeneIds);

	int numDriverGenes = driverGeneIds.size();
	cout << "\ttotal number of driver genes = " << numDriverGenes << endl;
	filename = "driver_genes.txt";
	saveGeneSymbols(filename.c_str(), &driverGeneIds, &geneIdToSymbol);

	delete bipartiteEdges;

	cout << "merging modules for all samples ...\n";

	//merge modules and trim explained genes for each sample
	vector<bool> isDriverGenes(totalGenes);
	for (int i = 0; i < numDriverGenes; ++i) {
		isDriverGenes[driverGeneIds[i]] = true;
	}

	vector< list<Module> > modulesListOfAllSamples(totalSamples);
	findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
			&mutatedGeneIdsListReal, &isPhenotypeGenes, &isDriverGenes, &phenotypeGeneIds);

	cout << "trimming explained genes for all samples ...\n";

	//TODO trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D);

	filename = "modules.txt";
	saveModules(&modulesListOfAllSamples, filename, &geneIdToSymbol);


	/*
	 * Stop timer
	 */

	cout << "DONE (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

	return 0;
}
