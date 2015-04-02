/*
 * explained_genes.h
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#ifndef EXPLAINED_GENES_H_
#define EXPLAINED_GENES_H_

#include "utilities.h"
#include "input.h"

/*
 * Structure for explained genes
 */

struct ExplainedGene{
	//int id; this is equal to the index of the vector
	double expression;
	int degree;
	int frequency; //save frequency of this gene for each sample
};

struct MutatedAndExplianedGenes{
	int mutatedGeneId;
	vector<int> explainedGeneIds;
};

/*
 * Functions for finding explained genes
 */

// return sampleExpression
void getGeneExpressionFromSampleId(TDoubleMatrix* geneExpressionMatrix, vector<int>* genesEx,
		vector<double>* sampleExpression, int sampleId, vector<string>* geneIdToSymbol);

// return mutatedGeneIds
void getMutatedGeneIdsFromSampleId(Mutations* mutations,
		vector<int>* mutatedGeneIds, int sampleId, vector<int>* genesMut);

// return explainedGenes
void getExplainedGenes(vector<ExplainedGene>* explainedGenes, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, int F);

void getExplainedGenesOnlyId(vector<int>* explainedGeneIds, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, int F);

void getMutatedAndExplainedGenes(vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, TIntAdjList* network,
		vector<double>* sampleGeneExpression, vector<int>* mutatedGeneIds, int L, int D, double F);

// BFS search for finding explained genes of geneId
void BFSforExplainedGenes(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<ExplainedGene>* explainedGenes, vector<double>* sampleGeneExpression)  ;

void BFSforExplainedGenesIdOnly(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<int>* explainedGenes, vector<double>* sampleGeneExpression);

// return genesCount: # samples containing explained genes
void getGeneFrequencyOfSamples(vector<int>* genesCount, int totalGenes,
		vector< vector<ExplainedGene> >* explainedGenesList);

#endif /* EXPLAINED_GENES_H_ */
