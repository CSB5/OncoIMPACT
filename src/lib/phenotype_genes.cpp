/*
 * phenotype_genes.cpp
 *
 *  Created on: 28 Mar, 2015
 *      Author: Nok
 */

#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include "../header/phenotype_genes.h"
#include "../header/utilities.h"

//input: explainedGenesListForPhenotype = a list of explained genes of all samples
//output: nullDistribution = a list contain the frequency of each round (of 500 round)
/*void addFrequencyForNullDistribution(vector<vector<int> >* nullDistribution,
		vector<vector<int> >* explainedGenesListForPhenotypeGenes) {
	int totalGenes = nullDistribution->size(); // number of genes in the network
	int numSamples = explainedGenesListForPhenotypeGenes->size();
	vector<int> frequency(totalGenes);

	//for each sample i in a current round
	for (int i = 0; i < numSamples; ++i) {
		vector<int> explainedGeneIds = explainedGenesListForPhenotypeGenes->at(
				i);

		//for each gene in sample i
		for (int j = 0; j < totalGenes; ++j) {
			if (explainedGeneIds[j] > 0) {
				frequency[j]++;	//count the number of samples
			}
		}
	}

	//add a frequency of a current round to the null
	for (int i = 0; i < totalGenes; ++i) {
		nullDistribution->at(i).push_back(frequency[i]);
	}
}*/

void countGeneFrequencyGreaterThanRealFrequency(vector<int>* geneFrequencyGreaterThanRealFrequencyCounter,
		vector< vector<int> >* explainedGenesFrequencyForPhenotypeGenes, vector<int>* genesFrequencyReal){
	int totalGenes = geneFrequencyGreaterThanRealFrequencyCounter->size(); // number of genes in the network
	int numSamples = explainedGenesFrequencyForPhenotypeGenes->size();

	vector<int> frequencies(totalGenes);

	//for each sample i
	for (int i = 0; i < numSamples; ++i) {
		vector<int> explainedGenesFrequency = explainedGenesFrequencyForPhenotypeGenes->at(i);

		for (int j = 0; j < totalGenes; ++j) {
			if(explainedGenesFrequency[j] > 0){
				frequencies[j]++;	//count the number of samples
			}
		}
	}


	for (int i = 0; i < totalGenes; ++i) {
		//if the frequency is greater than the real frequency
		if(frequencies[i] > genesFrequencyReal->at(i)){
			geneFrequencyGreaterThanRealFrequencyCounter->at(i)++;
		}
	}

}

void addFrequncyForRealDataset(vector<int>* genesFrequency,	vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isExplainedGenes) {
	int totalGenes = genesFrequency->size();
	int numSamples = mutatedAndExplainedGenesListReal->size();

	//for each sample i in a current round
	for (int i = 0; i < numSamples; ++i) {

		// mutated gene and its corresponding genes for each sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
				mutatedAndExplainedGenesListReal->at(i);
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal->at(i);

		// combine list of explained genes
		vector<bool> explainedGeneIdsAll(totalGenes); //for all mutated genes of this sample
		combineListOfExplainedGenes(&mutatedAndExplainedGenes, &mutatedGeneIds,
				&explainedGeneIdsAll, totalGenes);

		//copy the value to genesFrequency (# samples in which each gene is explained)
		//for each unique explained gene in sample i
		for (int j = 0; j < totalGenes; ++j) {
			if (explainedGeneIdsAll[j]) {
				genesFrequency->at(j)++;
				isExplainedGenes->at(j) = true;
			}
		}
	}

	vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
					mutatedAndExplainedGenesListReal->at(0);
	int numMut = mutatedAndExplainedGenes.size();
	for (int l = 0; l < numMut; ++l) {
		vector<int>* explainedGenesFreqency = mutatedAndExplainedGenes[l].explainedGenesFreqency;
	}

}

void combineListOfExplainedGenes(
		vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, vector<int>* mutatedGeneIds,
		vector<bool>* explainedGenesAll, int totalGenes) {

	int numGenesMut = mutatedGeneIds->size();

	//for each mutated gene
	for (int i = 0; i < numGenesMut; ++i) {
		vector<int>* explainedGenesFreqency =
				mutatedAndExplainedGenes->at(mutatedGeneIds->at(i)).explainedGenesFreqency;
		//for each corresponding gene
		for (int j = 0; j < totalGenes; ++j) {
			if (explainedGenesFreqency->at(j) > 0) {
				explainedGenesAll->at(j) = true;
			}
		}
	}

}

void findPhenotypeGenes(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds,
		vector<int>* genesFrequency, vector<vector<int> >* nullDistribution,
		vector<bool>* isExplainedGenes) {
	int totalGenes = genesFrequency->size();
	int round = nullDistribution->at(0).size();

	//TO DO use statistical test to find the phenotype genes (created the findPhenotypeGenesUsingCounter)
	//TEST simple use percentile
	int p95 = ceil(round / 100 * 95);

	//for each gene in the network
	for (int i = 0; i < totalGenes; ++i) {
		vector<int> distribution = nullDistribution->at(i);
		sort(distribution.begin(), distribution.end());
		int cutoff = distribution[p95];
		if (isExplainedGenes->at(i))
//			cout << "cutoff for gene " << i << " is " << cutoff << endl;
			if (isExplainedGenes->at(i)
					and genesFrequency->at(i) > cutoff) {
				isPhenotypeGenes->at(i) = true; // mark that this gene is a phenotype genes
				phenotypeGeneIds->push_back(i);
				//cout << i << endl;
			}
	}
}

void findPhenotypeGenesUsingCounter(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds,
		vector<int>* genesFrequencyReal, vector<int>* geneFrequencyGreaterThanRealFrequencyCounter,
		vector<bool>* isExplainedGenes, int round, int totalSamples, vector<string>* geneIdToSymbol){
	int totalGenes = genesFrequencyReal->size();
	vector<double> pValues(totalGenes);

	//for each gene in the network
	for (int i = 0; i < totalGenes; ++i) {
		if (isExplainedGenes->at(i)){	//is the explained genes in the real samples
			pValues[i] = 1.0 * geneFrequencyGreaterThanRealFrequencyCounter->at(i) / round;
			if(pValues[i] == 0){
				isPhenotypeGenes->at(i) = true;
				phenotypeGeneIds->push_back(i);
			}
		}
	}

	//OUTPUT: print phenotype genes and their p-values
	vector<string>* output = new vector<string>;
	output->push_back("GENE\tP_VALUE\tIS_PHENOTYPE");
	for (int i = 0; i < totalGenes; ++i) {
		if (isExplainedGenes->at(i)){
			string str = geneIdToSymbol->at(i) + "\t" + doubleToStr(pValues[i]);
			if(isPhenotypeGenes->at(i)){
				str = str + "\t" + "1";
			}else{
				str = str + "\t" + "0";
			}
			output->push_back(str);
		}
	}
	string filename = "output/phenotype_genes.tsv";
	writeStrVector(filename.c_str(), output);
	delete output;
}

void printPhenotypeGenes(vector<bool>* isPhenotypeGenes, string phenotypeGeneFileName, vector<string>* geneIdToSymbol){
	int totalGenes = isPhenotypeGenes->size();
	vector<int> phenotypeGeneIds;
	for (int i = 0; i < totalGenes; ++i) {
		if (isPhenotypeGenes->at(i) == 1) {
			phenotypeGeneIds.push_back(i);
		}
	}

	cout << "\twriting " << phenotypeGeneIds.size() << " phenotype genes to "
			<< phenotypeGeneFileName << " ..." << endl;
	saveGeneSymbols(phenotypeGeneFileName.c_str(), &phenotypeGeneIds,
			geneIdToSymbol);
}

