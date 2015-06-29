//============================================================================
// Name        : tester.cpp
// Author      : Nok C Suphavilai
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <set>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
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

int database(string outDir, string networkFilename, string expFilename, string snpFilename, string cnvFilename,
		string benchmarkGeneListFilename, int numThreads, int mode) {

	/*
	 * System configuration and timer
	 */

	clock_t begin_time = clock();
	clock_t step_time = clock();

	/*
	 * TODO Read configuration from a file
	 */

//	string outDir = "output";
//	string networkFilename = networkFilenameIn.substr(0, networkFilenameIn.length() - 1);
//	string expFilename = "GBM/EXPR.txt";
//	string snpFilename = "GBM/SNP.txt";
//	string cnvFilename = "GBM/CNV.txt";
//	string benchmarkGeneListFilename = "cancer_gene_census.csv";
//	//0 = sensitive, 1 = stringent
//	int mode = 0;
//	int numThreads = 10;
	omp_set_num_threads(numThreads);

	string dbPath = "";
	string dbPathExport = "";


	//create sub-directory to store output files
	#if defined(_WIN32)	//_WIN32 - Defined for applications for Win32 and Win64
		_mkdir((outDir + "/sensitive").c_str());
		_mkdir((outDir + "/stringent").c_str());
		_mkdir((outDir + "/sensitive/samples").c_str());
		_mkdir((outDir + "/stringent/samples").c_str());
		_mkdir((outDir + "/test_param").c_str());
	#else
		mkdir((outDir + "/sensitive").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/stringent").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/sensitive/samples").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/stringent/samples").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/test_param").c_str(), 0777); // notice that 777 is different than 0777
	#endif

	//parameters of analysis
	int numPermutationsForJSDivergence = 100;
	int numRandomDatasetForPhenotype = 500;
	//srand ( unsigned ( std::time(0) ) );	TODO uncomment this if want the different result for each run
	//srand (8);	//TODO uncomment this if want the fixed seed

	//create log file for each run
	string logFilename;
	logFilename = outDir + "/run.log";
	ofstream outLogStream;

	outLogStream.open(logFilename.c_str());
	writeToLogFile(&outLogStream, "Start Running oncoIMPACT 1.0 (using " + intToStr(numThreads) + " threads)");

	/*
	 * Read network file into adjacency list
	 */

	cout << "reading network file and create mapping <id, geneSymbol> ..."	<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	if(readNetwork(networkFilename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId, '\t') != 0){
		cout << "Error reading network file" << endl;
		return 1;
	}
	int totalGenes = network.size();
	int totalGenesUpDown = totalGenes * 2;
	writeToLogFile(&outLogStream, "Read network from " + networkFilename);

	/*
	 * Read gene expression matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of samples is the same as in mutation matrix (also the same order)
	 */

	//sample id mapping
	vector<string> sampleIdToName;

	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix ..." << endl;
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
	writeToLogFile(&outLogStream, "Read gene expression matrix from " + expFilename);

	/*
	 * Read point mutation matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of samples is the same as in gene expression matrix (also the same order)
	 */

	cout << "reading point mutation matrix ..." << endl;
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
	writeToLogFile(&outLogStream, "Read SNP from " + snpFilename);

	/*
	 * Read CNV matrix from file
	 * (row = gene, col = samples)
	 * Note: the set of sample is the same as in gene expression matrix (including order)
	 */

	cout << "reading CNV matrix ..." << endl;
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
	writeToLogFile(&outLogStream, "Read CNV from " + cnvFilename);

	/*
	 * Combining point mutation and CNV matrix
	 */

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
	writeToLogFile(&outLogStream, "SNP and CNV matrix are combined");

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

	/*
	 * Calculate JS divergence for each set of PARAMETERS (L,D,F)
	 */

	//Note: gene expression and mutation have the same set of samples
	//If number of samples is < 50, use all samples to tune the parameters
	int numSamples = 50;
	if(totalSamples < 50){
		numSamples = totalSamples;
	}

	cout << "\tDONE reading and preparing data (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";

	step_time = clock();
	cout << "tuning parameters by using " << numSamples << " randomly choosing samples ..."
			<< endl;

	cout << "computing JS divergence for all parameters (L,D,F) ... " << endl;

//	int LsVal[] = {16, 18, 20};	//fewer parameters for testing
//	int DsVal[] = {60, 65, 70};	//fewer parameters for testing
//	double FsVal[] = {2, 2.5};	//fewer parameters for testing

	//TODO initialize all the parameters to be tested
	int LsVal[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
	vector<int> Ls(LsVal, LsVal + sizeof LsVal / sizeof LsVal[0]);
	int numLs = Ls.size();
	int DsVal[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
	vector<int> Ds(DsVal, DsVal + sizeof DsVal / sizeof DsVal[0]);
	int numDs = Ds.size();	//fewer parameters for testing
	double FsVal[] = {1, 1.5, 2, 2.5, 3};
	vector<double> Fs(FsVal, FsVal + sizeof FsVal / sizeof FsVal[0]);
	if(mode == 1){	//if RNA_SEQ, fix fold-change to 1
		Fs.resize(1, 1.0);
	}
	int numFs = Fs.size();

	int numCombinations = numLs * numDs * numFs;
	//initialize the vector to save the divergence of each set of parameters
	vector<JSDivergence> jsDivergences(numCombinations);


//	//TEST permutation function
//	vector<string> testOut;
//	string testOutFilename = "permutedGenesEx.txt";
//	vector<int> permutedGeneLabelsEx;
//	permuteGeneLabels(&genesEx, &permutedGeneLabelsEx);
//	for (unsigned x = 0; x < permutedGeneLabelsEx.size(); ++x) {
//		testOut.push_back(geneIdToSymbol[genesEx[x]] + "\t" + geneIdToSymbol[permutedGeneLabelsEx[x]]);
//	}
//	writeStrVector(testOutFilename.c_str(), &testOut);
//	vector<string> testOut2;
//	string testOutFilename2 = "permutedGenesMut.txt";
//	vector<int> permutedGeneLabelsMut;
//	permutedGeneLabelsUsingAllGeneInNetwork(&genesMut, &permutedGeneLabelsMut, totalGenes);
//	for (unsigned x = 0; x < permutedGeneLabelsMut.size(); ++x) {
//		testOut2.push_back(geneIdToSymbol[genesMut[x]] + "\t" + geneIdToSymbol[permutedGeneLabelsMut[x]]);
//	}
//	writeStrVector(testOutFilename2.c_str(), &testOut2);
//	return 0;

	//TODO parameter finding: for js divegence calculation the result is quite different from the original
	findParameters(&jsDivergences, &Ls, &Ds, &Fs, totalGenes, &geneExpression, &mutations, &network, numSamples,
			numPermutationsForJSDivergence, &geneIdToSymbol, &geneSymbolToId, numThreads);

	//write the JS divergence result to a file
	string outJSFilename = outDir + "/parameters.dat";
	//save the JS divergence results
	//TODO sort the divergence first
	saveJSDivergences(&jsDivergences, outJSFilename);
	writeToLogFile(&outLogStream, "Write parameters results in " + outJSFilename);

	cout << "\tDONE finding parameters (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";

	//choose the best parameters
	JSDivergence maxJs;
	findMaximumJsDivergence(&jsDivergences, &maxJs);

	cout << "the maximum divergence is " << maxJs.divergence << " when L, D, F = " << maxJs.L << ", " << maxJs.D << ", " << maxJs.F << endl;
	writeToLogFile(&outLogStream, "The maximum JS divergence is " + doubleToStr(maxJs.divergence, 5) +
			" corresponding to (L,D,F) = (" + intToStr(maxJs.L) + "," + intToStr(maxJs.D) + "," + doubleToStr(maxJs.F, 1) + ")");

	//set the L D F to maxJs
	int L = maxJs.L;
	int D = maxJs.D;
	double F = maxJs.F;

	//[DEBUG]
//	int L = 12; //old 20
//	int D = 65;
//	double F = 2.5;

	/*
	 * Find PHENOTYPE GENES
	 */

	step_time = clock();
	cout << "finding phenotype genes using (L,D,F) = (" << L << "," << D << "," << F << ") ..." << endl;

	//A vector for collecting mutated genes and the corresponding explained genes for all samples (sample, mutated genes, explained genes)
	vector< vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal;
	//Note: size of MutatedAndExplianedGenes = totalGenes, so index = mutated gene id, value = a list of explained gene ids
	vector< vector<int> > mutatedGeneIdsListReal;	//to save a list of mutated gene ids for each sample

	cout << "\tgetting explained genes frequency of the REAL samples ...\n";
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

		//to tell whether each gene is an explained gene in this sample i
		//initialize the mutatedAndExplainedGenes for all genes
		for (int gi = 0; gi < totalGenes; ++gi) {
			mutatedAndExplainedGenes[gi].isExplainedGenesUpDown = new vector<bool>(totalGenes * 2, false);
		}

		for (int j = 0; j < numMutatedGenes; ++j) {	// for each mutated genes
			int mutatedGeneId = mutatedGeneIds[j];
			MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes[mutatedGeneId];
//			meg->isExplainedGenesUpDown = new vector<bool>(totalGenes * 2);	//already init above
			//BFS for explained genes of the current mutated gene
			BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(&network, mutatedGeneId, L, D, F,
					meg->isExplainedGenesUpDown, &sampleGeneExpression, i, &geneIdToSymbol, &geneSymbolToId);
		}

		//add explained genes to the list of all samples
		mutatedAndExplainedGenesListReal.push_back(mutatedAndExplainedGenes);
		mutatedGeneIdsListReal.push_back(mutatedGeneIds);

	}

	// get the frequency of genes in the real dataset
	vector<bool> isExplainedGenesUpDown(totalGenesUpDown, false);
	vector<int> explainedGenesFrequencyRealUpDown(totalGenes * 2, 0);	// this will be used to compare with the frequency from random dataset
	addFrequncyForRealDataset(&explainedGenesFrequencyRealUpDown,
			&mutatedAndExplainedGenesListReal, &mutatedGeneIdsListReal,
			&isExplainedGenesUpDown);

	/*
	 * Print original modules file
	 */

	//OUTPUT: print all modules in all samples (as original)
	vector<string>* outputStr = new vector<string>;
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal[i];

		//for each mutated genes
		for (unsigned int j = 0; j < mutatedGeneIds.size(); ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			string str = sampleIdToName[i] + "\t" + geneIdToSymbol[currentMutatedGeneId] + "\t";
			vector<bool>* isExplainedGenesUpDownForAMutatedGene = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
			int numMember = 0;
			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(isExplainedGenesUpDownForAMutatedGene->at(k)){
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
	string outModulefilename = outDir + "/MODULE.dat";
	writeStrVector(outModulefilename.c_str(), outputStr);
	delete outputStr;
	writeToLogFile(&outLogStream, "Save calculated gene modules to " + outModulefilename);

	cout << "\tDONE calculated real distribution (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";
	step_time = clock();

	writeToLogFile(&outLogStream, "Start finding phenotype genes");

	/*
	 * Precalculation for random datasets
	 */

	cout << "Start finding phenotype genes\n";
	cout << "\tcalculating look up tables for finding phenotype genes ... ";

	//SPEED UP: for each gene in the network, if it is mutated, find which genes are explained for it in a given sample
	//cannot directly construct the bool matrix because there is not enough memory
	vector< vector< vector<int> > > explainedGeneIdsForAMutatedGeneInASample(totalGenes, vector< vector<int> >(totalSamples));

	int progress = 1;
	int interval = totalSamples / 100;

	//for each sample, create a list of explained genes
	for (int sampleId = 0; sampleId < totalSamples; ++sampleId) {
		//print progression
		if (sampleId % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

		//get gene expression of a current sample
		vector<double> sampleGeneExpression(totalGenes); // to save expression of of all genes in the network
		getGeneExpressionFromSampleId(&originalGeneExpressionMatrix,
				&genesEx, &sampleGeneExpression, sampleId);

		//for each (possible mutated) gene, fint the corresponding explained genes for the current sample
		#pragma omp parallel for
		for (int geneId = 0; geneId < totalGenes; ++geneId) {
			vector<bool> isExplainedUpDown(totalGenesUpDown, false);
			BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(&network, geneId, L, D, F,
					&isExplainedUpDown, &sampleGeneExpression, sampleId, &geneIdToSymbol, &geneSymbolToId);
			//add explained genes
			for (int i = 0; i < totalGenesUpDown; ++i) {
				if(isExplainedUpDown[i]){
					explainedGeneIdsForAMutatedGeneInASample[geneId][sampleId].push_back(i);
				}
			}
		}
	}

	cout << endl;	//for progression

	cout << "\tDONE precalculated tables for SPEED UP (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";


	//create a vector for counting the number of times (for each gene) the random samples have greater frequency than the real samples
	vector<int> geneFrequencyGreaterThanRealFrequencyCounter(totalGenesUpDown, 0);
	//create a vector for counting the number of times the ramdom samples have greater frequency than the real samples

	step_time = clock();

	cout << "\tcreating null distribution (using " << numRandomDatasetForPhenotype << " random data sets) ... ";

	progress = 1;
	interval = numRandomDatasetForPhenotype / 100;

	#pragma omp parallel for
	for (int r = 0; r < numRandomDatasetForPhenotype; ++r) {

		//print progression
		if (r % interval == 0) {
			const string progStatus = intToStr(progress) + "%";
			cout << progStatus << flush;
			progress++;
			cout << string(progStatus.length(), '\b');
		}

		//construct matrix: row = samples col = permuted mutated gene ids
		vector< vector<int> > permutedGeneIdsForAll(totalSamples, vector<int>(numGenesMut, 0));
		vector<int> allGeneIds(totalGenes);	//to permute all genes in the network
		for (int si = 0; si < totalSamples; ++si) {

			for (int i = 0; i < totalGenes; i++) {
				allGeneIds[i] = i;
			}
			random_shuffle(allGeneIds.begin(), allGeneIds.end());

			for (int mi = 0; mi < numGenesMut; ++mi) {
				permutedGeneIdsForAll[si][mi] = allGeneIds[mi];
			}
		}


		vector< vector<bool> > isExplainedGenesUpDownRandom(totalSamples, vector<bool>(totalGenesUpDown, false));

		//permute the gene labels of each sample independently (this have been done in precomputed part)
		for (int si = 0; si < totalSamples; ++si) {


			//for each mutated gene
			for (int mi = 0; mi < numGenesMut; ++mi) {
				//if a current gene is mutated in this sample
				if (mutations.matrix->at(mi)[si] != 0) {
					//use the precomputed permuted gene label
					int currentMutatedGeneId = permutedGeneIdsForAll[si][mi];
					//get the explained gene of the current mutated gene of the current sample
					vector<int> explainedGensIds = explainedGeneIdsForAMutatedGeneInASample[currentMutatedGeneId][si];

					int numExplainedGenes = explainedGensIds.size();
					for (int ei = 0; ei < numExplainedGenes; ++ei) {
						isExplainedGenesUpDownRandom[si][explainedGensIds[ei]] = true;
					}
				}
			}

		}	//end for each sample


		//count the frequency that exceed the real frequency (This is of the old version where all explained genes are not consider independently
		countGeneFrequencyGreaterThanRealFrequency(&geneFrequencyGreaterThanRealFrequencyCounter,
				&isExplainedGenesUpDownRandom, &explainedGenesFrequencyRealUpDown);

	}	//end loop of 500 rounds

	cout << endl; //for print progression

	vector<bool> isPhenotypeGenesUpDown(totalGenesUpDown, false);
	vector<int> phenotypeGeneIdsUpDown;	// phenotype gene ids
	vector<double> pValues(totalGenesUpDown, 0.0);

	findPhenotypeGenesUsingCounter(&isPhenotypeGenesUpDown, &phenotypeGeneIdsUpDown, &pValues, &explainedGenesFrequencyRealUpDown,
			&geneFrequencyGreaterThanRealFrequencyCounter, &isExplainedGenesUpDown, numRandomDatasetForPhenotype, totalSamples, &geneIdToSymbol);

	cout << "\tDONE finding phenotype genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC) << " sec)\n";
	step_time = clock();	//update the clock
	cout << "\tthere are " << phenotypeGeneIdsUpDown.size() << " phenotype genes" << endl;

//	//[DEBUG] use phenotype genes list from previous version
//	readGenesListUpDown("original_phenotype_genes.txt", &phenotypeGeneIdsUpDown, &geneSymbolToId);
//	int numPhenotypeGene = phenotypeGeneIdsUpDown.size();
//	cout << "use phonotype genes list from previous version (" << numPhenotypeGene << " genes)" << endl;
//	for (int i = 0; i < numPhenotypeGene; ++i) {
//		isPhenotypeGenesUpDown[phenotypeGeneIdsUpDown[i]] = true;
//	}
//	//End [DEBUG]

	//OUTPUT: print gene frequency and phenotype genes of the real dataset
	string outExplainedAndPhenotypeFilename = outDir + "/exp_gene_freq.dat";
	printExplinedGenesFrequencyAndPhonotype(&explainedGenesFrequencyRealUpDown, &pValues, &isPhenotypeGenesUpDown, &geneIdToSymbol, &network,
			&originalGeneExpressionMatrix, &genesEx, F, outExplainedAndPhenotypeFilename);
	writeToLogFile(&outLogStream, "Save explained and phenotype gene statistics to " + outExplainedAndPhenotypeFilename);

	/*
	 * Find driver genes
	 */

	/*
	 * SENSITIVE MODE
	 */

	cout << "SENSITIVE mode ...\n";
	writeToLogFile(&outLogStream, "Start sensitive mode");

	mode = 0;

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
		cout << "\tDONE finding driver genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC)
				<< " sec)\n";
		step_time = clock();	//update the clock

//		//[DEBUG] use driver genes list from previous version
//		vector<int> driverGeneIds;	// to save diver gene ids
//		readGenesList("original_driver_genes.txt", &driverGeneIds, &geneSymbolToId);
//		int numDriverGenesFromFile = driverGeneIds.size();
//		for (int i = 0; i < numDriverGenesFromFile; ++i) {
//			DriverGene driverGene;
//			driverGene.geneId = driverGeneIds[i];
//			driverGenes.push_back(driverGene);
//		}
//		//End [DEBUG]

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
		string outFinalModuleFilenameSensitive = outDir + "/sensitive/FINAL_MODULE.dat";
		saveModules(&modulesListOfAllSamples, &mutatedAndExplainedGenesListReal, outFinalModuleFilenameSensitive, &geneIdToSymbol, &sampleIdToName);
		writeToLogFile(&outLogStream, "Save final module of sensitive mode to " + outFinalModuleFilenameSensitive);

		cout << "\tcalculating IMPACT scores for all samples ...\n";
		writeToLogFile(&outLogStream, "calculate IMPACT scores for all samples");

		vector< vector<Driver> > driversOfAllSamples(totalSamples);
		string outDriverOfAllSamplesFilename;
		if(mode == 0){
			outDriverOfAllSamplesFilename = outDir + "/sensitive/driver_all_samples.dat";
		}else{
			outDriverOfAllSamplesFilename = outDir + "/stringent/driver_all_samples.dat";
		}
		calculateImpactScoresForAllSamples(&modulesListOfAllSamples, &driversOfAllSamples,
				&originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol, outDriverOfAllSamplesFilename, &sampleIdToName);

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
		aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol);

//		cout << "getting driver frequency ...\n";
		getDetailDriversFreqeuncy(&driversOfAllSamples,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&originalPointMutationsMatrix, &originalCNVsMatrix,
				&genesPointMut, &genesCNV);

//		cout << "getting mutation frequency ...\n";
		getMutationFrequency(&originalMutationMatrix, &mutationFrequency, &genesMut);
		getDetailMutationFrequency(&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

		string outSampleResultDirName = outDir + "/sensitive/samples/";
//		cout << "printing impact scores for all samples ...\n";
		printSampleDriverList(&driversOfAllSamples, outSampleResultDirName, &geneIdToSymbol, &sampleIdToName,
				&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency, &isCancerBenchmarkGenes);

		cout << "\tprinting aggregated impact scores ...\n";
		string outDriverListfilename = outDir + "/sensitive/driver_list.txt";
		printAggregatedDriverList(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency, &isCancerBenchmarkGenes);

		writeToLogFile(&outLogStream, "Save final driver gene list to " + outDriverListfilename);

	}

	/*
	 * STRINGENT MODE
	 */

	cout << "STRINGENT mode ...\n";
	writeToLogFile(&outLogStream, "Start stringent mode");
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

		vector<DriverGene> driverGenes;	//driver gene id with sample ids

		cout << "\tperforming greedy minimum set cover algorithm ...\n";
		cout << "\tfinding driver genes ...\n";
		findDriverGenes(bipartiteGraph, &mutatedGeneIdsList, &driverGenes);
		delete bipartiteGraph;
		cout << "\tDONE finding driver genes (" << (float(clock() - step_time) / CLOCKS_PER_SEC)
				<< " sec)\n";
		step_time = clock();	//update the clock

//		//[DEBUG] use driver genes list from previous version
//		vector<int> driverGeneIds;	// to save diver gene ids
//		readGenesList("original_driver_genes.txt", &driverGeneIds, &geneSymbolToId);
//		int numDriverGenesFromFile = driverGeneIds.size();
//		for (int i = 0; i < numDriverGenesFromFile; ++i) {
//			DriverGene driverGene;
//			driverGene.geneId = driverGeneIds[i];
//			driverGenes.push_back(driverGene);
//		}
//		End [DEBUG]

		int numDriverGenes = driverGenes.size();
		cout << "\ttotal driver genes = " << numDriverGenes << endl;

		cout << "\tgenerating modules for all samples ...\n";

		//merge modules and trim explained genes for each sample
		//use mutatedAndExplainedGenesListReal (samples, mutated genes, explained genesUpDown)
		findModulesInAllSamples(&mutatedAndExplainedGenesListReal, &modulesListOfAllSamples,
				&mutatedGeneIdsListReal, &isPhenotypeGenesUpDown, &driverGenes, &phenotypeGeneIdsUpDown, mode);

		cout << "\ttrimming explained genes for all samples ...\n";
		trimSomeExplainedGenes(&modulesListOfAllSamples, &network, L, D, &geneIdToSymbol);

		cout << "\twriting final module to FINAL_MODULE.dat ...\n";
		string outFinalModuleFilenameSensitive = outDir + "/stringent/FINAL_MODULE.dat";
		saveModules(&modulesListOfAllSamples, &mutatedAndExplainedGenesListReal, outFinalModuleFilenameSensitive, &geneIdToSymbol, &sampleIdToName);
		writeToLogFile(&outLogStream, "Save final module of stringent mode to " + outFinalModuleFilenameSensitive);

		cout << "\tcalculating IMPACT scores for all samples ...\n";
		writeToLogFile(&outLogStream, "calculate IMPACT scores for all samples");

		vector< vector<Driver> > driversOfAllSamples(totalSamples);
		string outDriverOfAllSamplesFilename;
		if(mode == 0){
			outDriverOfAllSamplesFilename = outDir + "/stringent/driver_all_samples.dat";
		}else{
			outDriverOfAllSamplesFilename = outDir + "/stringent/driver_all_samples.dat";
		}
		calculateImpactScoresForAllSamples(&modulesListOfAllSamples, &driversOfAllSamples,
				&originalGeneExpressionMatrix, &genesEx, totalGenes, F, &geneIdToSymbol, outDriverOfAllSamplesFilename, &sampleIdToName);

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
		aggregateDriversAcrossSamples(&driversOfAllSamples, &driverAggregatedScores, &driversFrequency, &geneIdToSymbol);

//		cout << "getting driver frequency ...\n";
		getDetailDriversFreqeuncy(&driversOfAllSamples,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&originalPointMutationsMatrix, &originalCNVsMatrix,
				&genesPointMut, &genesCNV);

//		cout << "getting mutation frequency ...\n";
		getMutationFrequency(&originalMutationMatrix, &mutationFrequency, &genesMut);
		getDetailMutationFrequency(&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency);

		string outSampleResultDirName = outDir + "/stringent/samples/";
//		cout << "printing impact scores for all samples ...\n";
		printSampleDriverList(&driversOfAllSamples, outSampleResultDirName, &geneIdToSymbol, &sampleIdToName,
				&originalPointMutationsMatrix, &originalCNVsMatrix, &genesPointMut, &genesCNV,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency, &isCancerBenchmarkGenes);

		cout << "\tprinting aggregated impact scores ...\n";
		string outDriverListfilename = outDir + "/stringent/driver_list.txt";
		printAggregatedDriverList(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency, &isCancerBenchmarkGenes);
		//print in the outDir
		outDriverListfilename = outDir + "/driver_list.txt";
		printAggregatedDriverList(&driverGenes, outDriverListfilename, &geneIdToSymbol, &sampleIdToName,
				&driverAggregatedScores, &driversFrequency, &mutationFrequency,
				&pointMutationDriversFrequency, &deletionDriversFrequency, &amplificationDriversFrequency,
				&pointMutationFrequency, &deletionFrequency, &amplificationFrequency, &isCancerBenchmarkGenes);

		writeToLogFile(&outLogStream, "Save final driver gene list to " + outDriverListfilename);

	} //end STRINGENT mode

	//delete the vector<int>* mutatedAndExplainedGenesListReal
	for (int i = 0; i < totalSamples; ++i) {
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal[i];
		for (int gi = 0; gi < totalGenes; ++gi) {
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[gi].isExplainedGenesUpDown;
			delete isExplainedGenesUpDown;
		}
	}

	/*
	 * Stop timer
	 */

	cout << "DONE (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)
			<< " sec)\n";

	//close the log file
	outLogStream.close();

	return 0;
}

