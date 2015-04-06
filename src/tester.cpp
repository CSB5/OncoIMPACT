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

//#include <windows.h>
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

	int totalSamples = originalGeneExpressionMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesEx = originalGeneExpressionMatrix.size();

	cout << "\ttotal genes in gene expression matrix is " << numGenesEx << endl;
	cout << "\ttotal samples in gene expression matrix is " << totalSamples
			<< endl;

	/*
	 * Read point mutation matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of sample is the same as in gene expression matrix (including order)
	 */

	cout << "reading point mutation matrix..." << endl;
	filename = "GBM/SNP.txt";
	vector<int> genesPointMut;
	TIntegerMatrix originalPointMutationsMatrix;
	PointMutations pointMutations;
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

	cout << "reading CNV matrix..." << endl;
	filename = "GBM/CNV.txt";
	vector<int> genesCNV;
	TIntegerMatrix originalCNVsMatrix;
	CopyNumberVariation CNVs;
	CNVs.genes = &genesCNV;
	CNVs.matrix = &originalCNVsMatrix;
	readCopyNumberVariation(filename.c_str(), &CNVs, '\t', &geneSymbolToId);

	totalSamples = originalCNVsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesCNV = originalCNVsMatrix.size();

	cout << "\ttotal genes in CNV matrix is " << numGenesCNV << endl;
	cout << "\ttotal samples in CNV matrix is " << totalSamples
			<< endl;


	//combine lists of point mutaton and CNV genes
	vector<int> genesMut;
	vector<bool>* isMutated = new vector<bool>(totalGenes);
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

		//find a current gene in point mutation matrix
		int pmRowIndex = findIndex(&genesPointMut, currentGeneId);

		//find a current gene in CNV matrix
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

	Mutations mutations;
	mutations.genes = &genesMut;
	mutations.matrix = &originalMutationMatrix;

	totalSamples = originalMutationMatrix[0].size(); // = originalMutationMatrix.size()

	cout << "\ttotal genes in mutatation matrix is " << numGenesMut << endl;
	cout << "\ttotal samples in mutatation matrix is " << totalSamples
			<< endl;

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

	//TODO do we need to resampling the samples for every round?
	//list of samples id to be used for tuning the parameters
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

	cout << "computing JS divergence for all parameters (L,D,F) ... " << endl;

//	int minL = 2; // to be used
//	int maxL = 5; // to be used

//	int Ls[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
	int LsVal[] = {14, 16, 18, 20};
	vector<int> Ls(LsVal, LsVal + sizeof LsVal / sizeof LsVal[0]);
	int numLs = Ls.size();
//	int Ds[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
	int DsVal[] = {55, 60, 65, 70};
	vector<int> Ds(DsVal, DsVal + sizeof DsVal / sizeof DsVal[0]);
	int numDs = Ds.size();
//	double Fs[] = {1, 1.5, 2, 2.5, 3};
	double FsVal[] = {2, 2.5};
	vector<double> Fs(FsVal, FsVal + sizeof FsVal / sizeof FsVal[0]);
	int numFs = Fs.size();


	int numCombinations = numLs * numDs * numFs;
	vector<JSDivergence>* jsDivergences = new vector<JSDivergence>(numCombinations);

	//findParameters(jsDivergences, &Ls, &Ds, &Fs, totalGenes, &subGeneExpression, &subMutations, &network);

	cout << "DONE tunning parameters (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";


	//choose the best parameters
	//JSDivergence maxJs;
	//findMaximumJsDivergence(jsDivergences, &maxJs);

	//cout << "the maximum divergence is " << maxJs.divergence << " when L, D, F = " << maxJs.L << ", " << maxJs.D << ", " << maxJs.F << endl;

	//TODO set the L D F to maxJs
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
				&sampleGeneExpression, i);

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
					&genesEx, &sampleGeneExpression, i);
			//cout << "\tcopied gene expression value" << endl;

			//find mutated genes of a current sample
			vector<int> mutatedGeneIds; // to store gene id of mutated genes
			getMutatedGeneIdsFromSampleId(&mutations, &mutatedGeneIds, i,
					&permutedGeneLabelsMut);
			//cout << "\tfound " << mutatedGeneIds.size() << " mutated genes" << endl;

			//find explained genes of a current sample
			vector<int> explainedGeneIds(totalGenes);
			getExplainedGenesIdOnly(&explainedGeneIds, &network,
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
