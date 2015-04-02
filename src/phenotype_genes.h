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

void addFrequencyForNullDistribution(vector< vector<int> >* nullDistribution, vector< vector<int> >* explainedGenesListForPhenotype);

void addFrequncyForRealDataset(vector<int>* genesFrequency,	vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isExplainedGenes);

void combineListOfExplainedGenes(
		vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, vector<int>* mutatedGeneIds,
		vector<bool>* explainedGenesAll, int totalGenes);

void findPhenotypeGenes(vector<bool>* phenotypeGeneIds, vector<int>* genesFrequency,
		vector< vector<int> >* nullDistribution, vector<bool>* isExplainedGenes);

#endif /* PHENOTYPE_GENES_H_ */
