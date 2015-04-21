/*
 * phenotype_genes.h
 *
 *  Created on: 28 Mar, 2015
 *      Author: Nok
 */

#ifndef PHENOTYPE_GENES_H_
#define PHENOTYPE_GENES_H_

#include "utilities.h"
#include "explained_genes.h"

//void addFrequencyForNullDistribution(vector< vector<int> >* nullDistribution, vector< vector<int> >* explainedGenesListForPhenotypeGenes);
void countGeneFrequencyGreaterThanRealFrequency(vector<int>* geneFrequencyGreaterThanRealFrequencyCounter,
		vector< vector<bool> >* explainedGenesFrequencyUpDownRandom, vector<int>* genesFrequencyReal);

void addFrequncyForRealDataset(vector<int>* genesFrequency,	vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isExplainedGenes);

void combineListOfExplainedGenes(
		vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, vector<int>* mutatedGeneIds,
		vector<bool>* explainedGenesAll, int totalGenes);

void findPhenotypeGenes(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds,
		vector<int>* genesFrequency, vector<vector<int> >* nullDistribution,
		vector<bool>* isExplainedGenes);

void findPhenotypeGenesUsingCounter(vector<bool>* isPhenotypeGenes, vector<int>* phenotypeGeneIds, vector<double>* pValues,
		vector<int>* genesFrequencyRealUpDown, vector<int>* geneFrequencyGreaterThanRealFrequencyCounter,
		vector<bool>* isExplainedGenes, int round, int totalSamples, vector<string>* geneIdToSymbol);

void printPhenotypeGenes(vector<bool>* isPhenotypeGenes, string phenotypeGeneFileName, vector<string>* geneIdToSymbol);

#endif /* PHENOTYPE_GENES_H_ */
