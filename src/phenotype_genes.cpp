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
#include "phenotype_genes.h"

//input: explainedGenesListForPhenotype = a list of explained genes of all samples
//output: nullDistribution = a list contain the frequency of each round (of 500 round)
void addFrequencyForNullDistribution(vector< vector<int> >* nullDistribution, vector< vector<int> >* explainedGenesListForPhenotypeGenes){
	int totalGenes = nullDistribution->size();	// number of genes in the network
	int numSamples = explainedGenesListForPhenotypeGenes->size();
	vector<int> frequency(totalGenes);

	//for each sample i in a current round
	for (int i = 0; i < numSamples; ++i) {
		vector<int> explainedGeneIds = explainedGenesListForPhenotypeGenes->at(i);

		//for each gene in sample i
		for (int j = 0; j < totalGenes; ++j) {
			if(explainedGeneIds[j] > 0){
				frequency[j]++;	//count the number of samples
			}
		}
	}

	//add a frequency of a current round to the null
	for (int i = 0; i < totalGenes; ++i) {
		nullDistribution->at(i).push_back(frequency[i]);
	}
}

void addFrequncyForRealDataset(vector<int>* genesFrequency,
		vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal) {
	int totalGenes = genesFrequency->size();

	int numSamples = mutatedAndExplainedGenesListReal->size();

	//for each sample i in a current round
	for (int i = 0; i < numSamples; ++i) {

		// mutated gene and its corresponding genes for each sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
				mutatedAndExplainedGenesListReal->at(i);

		// combine list of explained genes
		vector<int> explainedGeneIdsAll(totalGenes); //for all mutated genes of this sample
		combineListOfExplainedGenes(&mutatedAndExplainedGenes,
				&explainedGeneIdsAll, totalGenes);

		//copy the value to genesFrequency (# samples in which each gene is explained)
		//for each unique explained gene in sample i
		for (int j = 0; j < totalGenes; ++j) {
			if (explainedGeneIdsAll[j] > 0) {
				genesFrequency->at(j)++;}
			}
		}
	}

void combineListOfExplainedGenes(
		vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes,
		vector<int>* explainedGenesAll, int totalGenes) {

	int numGenesMut = mutatedAndExplainedGenes->size();

	//for each mutated gene
	for (int i = 0; i < numGenesMut; ++i) {
		vector<int> explainedGeneIds =
				mutatedAndExplainedGenes->at(i).explainedGeneIds;
		//for each corresponding gene
		for (int j = 0; j < totalGenes; ++j) {
			if (explainedGeneIds[j] > 0) {
				explainedGenesAll->at(j) = 1;
			}
		}
	}

}

void findPhenotypeGenes(vector<bool>* isPhenotypeGenes,
		vector<int>* genesFrequency, vector<vector<int> >* nullDistribution, vector<bool>* isInGeneExpressionMatrix) {
	int totalGenes = genesFrequency->size();
	int round = nullDistribution->at(0).size();

	//TEST simple use percentile
	int p95 = ceil(round / 100 * 95);

	//for each gene in the network
	for (int i = 0; i < totalGenes; ++i) {
		vector<int> distribution = nullDistribution->at(i);
		sort(distribution.begin(), distribution.end());
		int cutoff = distribution[p95];
		if(isInGeneExpressionMatrix->at(i))
//			cout << "cutoff for gene " << i << " is " << cutoff << endl;
		if (isInGeneExpressionMatrix->at(i) and genesFrequency->at(i) > cutoff) {
			isPhenotypeGenes->at(i) = true; // mark that this gene is a phenotype genes
			//cout << i << endl;
		}
	}
}

