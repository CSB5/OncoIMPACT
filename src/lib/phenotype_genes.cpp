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
		vector< vector<bool> >* explainedGenesFrequencyUpDownRandom, vector<int>* genesFrequencyReal){
	int totalGenesUpDown = geneFrequencyGreaterThanRealFrequencyCounter->size(); // number of genes in the network
	int numSamples = explainedGenesFrequencyUpDownRandom->size();

	vector<int> frequencies(totalGenesUpDown);	//for counting the number of samples that have the explained genes

	//for each sample i
	for (int i = 0; i < numSamples; ++i) {
		vector<bool> explainedGenesFrequencyUpDownOfASample = explainedGenesFrequencyUpDownRandom->at(i);
		for (int j = 0; j < totalGenesUpDown; ++j) {
			if(explainedGenesFrequencyUpDownOfASample[j]){
				frequencies[j]++;	//count the number of samples
			}
		}
	}

	for (int i = 0; i < totalGenesUpDown; ++i) {
		//if the frequency is greater than the real frequency
		if(frequencies[i] > genesFrequencyReal->at(i)){	//TODO check > or >=
			geneFrequencyGreaterThanRealFrequencyCounter->at(i)++;
		}
	}

}

void addFrequncyForRealDataset(vector<int>* genesFrequencyRealUpDown, vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isExplainedGenes) {
	int totalGenesUpDown = genesFrequencyRealUpDown->size();
	int totalGenes = totalGenesUpDown / 2;
	for (int i = 0; i < totalGenesUpDown; ++i) {
		genesFrequencyRealUpDown->at(i) = 0;
	}

	int numSamples = mutatedAndExplainedGenesListReal->size();

	//for each sample i in a current round
	for (int i = 0; i < numSamples; ++i) {

		// mutated gene and its corresponding genes for each sample
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
				mutatedAndExplainedGenesListReal->at(i);
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal->at(i);

		//combine list of explained genes
		vector<bool> isExplainedGenesForAll(totalGenesUpDown); //for all mutated genes of this sample
		for (int i = 0; i < totalGenesUpDown; ++i) {
			isExplainedGenesForAll[i] = false;
		}
		//set the isExplainedGenesForAll to be true if the gene is explained in this sample i
		combineListOfExplainedGenes(&mutatedAndExplainedGenes, &mutatedGeneIds,
				&isExplainedGenesForAll, totalGenesUpDown);

		//increase frequency in genesFrequency (contains # samples in which each gene is explained)
		//for each unique explained gene in sample i
		for (int j = 0; j < totalGenesUpDown; ++j) {
			if (isExplainedGenesForAll[j]) {
				genesFrequencyRealUpDown->at(j)++;
				//CHANGED
				if(j < totalGenes){
					isExplainedGenes->at(j) = true;
				}else{
					isExplainedGenes->at(j-totalGenes) = true;
				}
			}
		}
	}

}

void combineListOfExplainedGenes(
		vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, vector<int>* mutatedGeneIds,
		vector<bool>* isExplainedGenesForAll, int totalGenesUpDown) {

	int numGenesMut = mutatedGeneIds->size();

	//for each mutated gene
	for (int i = 0; i < numGenesMut; ++i) {
		vector<bool>* explainedGenesFreqency =
				mutatedAndExplainedGenes->at(mutatedGeneIds->at(i)).isExplainedGenesUpDown;
		//for each corresponding gene
		for (int j = 0; j < totalGenesUpDown; ++j) {
			if (explainedGenesFreqency->at(j)) {
				isExplainedGenesForAll->at(j) = true;
			}
		}
	}

}

//void findPhenotypeGenes(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds,
//		vector<int>* genesFrequency, vector<vector<int> >* nullDistribution,
//		vector<bool>* isExplainedGenes) {
//	int totalGenes = genesFrequency->size();
//	int round = nullDistribution->at(0).size();
//
//	//TO DO use statistical test to find the phenotype genes (created the findPhenotypeGenesUsingCounter)
//	//TEST simple use percentile
//	int p95 = ceil(round / 100 * 95);
//
//	//for each gene in the network
//	for (int i = 0; i < totalGenes; ++i) {
//		vector<int> distribution = nullDistribution->at(i);
//		sort(distribution.begin(), distribution.end());
//		int cutoff = distribution[p95];
//		if (isExplainedGenes->at(i))
////			cout << "cutoff for gene " << i << " is " << cutoff << endl;
//			if (isExplainedGenes->at(i)
//					and genesFrequency->at(i) > cutoff) {
//				isPhenotypeGenes->at(i) = true; // mark that this gene is a phenotype genes
//				phenotypeGeneIds->push_back(i);
//				//cout << i << endl;
//			}
//	}
//}

void findPhenotypeGenesUsingCounter(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds, vector<double>* pValues,
		vector<int>* genesFrequencyRealUpDown, vector<int>* geneFrequencyGreaterThanRealFrequencyCounter,
		vector<bool>* isExplainedGenes, int rounds, int totalSamples, vector<string>* geneIdToSymbol){

	int totalGenesUpDown = genesFrequencyRealUpDown->size();
	int totalGenes = isExplainedGenes->size();

	for (int i = 0; i < totalGenes; ++i) {
		isPhenotypeGenes->at(i) = false;
		pValues->at(i) = 0;
		pValues->at(i+totalGenes) = 0;
	}

	//for each gene in the network
	for (int i = 0; i < totalGenesUpDown; ++i) {
		//print only the explained genes
		if(genesFrequencyRealUpDown->at(i) > 0){
			if (i < totalGenes) {
				//up regulated
				pValues->at(i) = 1.0 * geneFrequencyGreaterThanRealFrequencyCounter->at(i) / rounds;
				cout << geneIdToSymbol->at(i) << "_UP has p-value = " << pValues->at(i) << endl;
				if (pValues->at(i) == 0) {
					isPhenotypeGenes->at(i) = true;
					//do not need to check since the up regulated genes always be considered first
					phenotypeGeneIds->push_back(i);
				}
			} else {
				//down regulated
				pValues->at(i) = 1.0 * geneFrequencyGreaterThanRealFrequencyCounter->at(i) / rounds;
				cout << geneIdToSymbol->at(i-totalGenes) << "_DOWN has p-value = " << pValues->at(i) << endl;
				if (pValues->at(i) == 0) {
					isPhenotypeGenes->at(i - totalGenes) = true;
					//check if the gene i is already added (in the up reguated part) to the vector of phenotypeGeneIds
					if (!isPhenotypeGenes->at(i-totalGenes)) {
						phenotypeGeneIds->push_back(i - totalGenes);
					}
				}
			}
		}
	}
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

