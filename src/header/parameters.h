/*
 * parameters.h
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "utilities.h"
#include "../header/input.h"
#include "../header/utilities.h"
#include "../header/explained_genes.h"

struct JSDivergence{
	int L;
	int D;
	double F;
	double divergence;
};

void findParameters(vector<JSDivergence>* jsDivergences, vector<int>* Ls,
		vector<int>* Ds, vector<double>* Fs, int totalGenes,
		GeneExpression* geneExpression, Mutations* mutations,
		TIntAdjList* network, int numSamples, int numDatasets,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId, int numThreads);

double calculateJSDivergence(const vector<vector<int> >* realDistributionAll,
		const vector<vector<int> >* randomDistributionAll, int numSamples, //int sumOfNumDeregulatedGenes,
		vector< vector<bool> >* isDeregulatedGensAll, int L, int D, double F);
double getMedianNumberOfDeregulatedGenes(TDoubleMatrix* geneExpressionMatrix, double F);
void findMaximumJsDivergence(vector<JSDivergence>* jsDivergences, JSDivergence* maxJs);

int countNumberOfDeregulatedGenes(TDoubleMatrix* geneExpressionMatrix, double F, vector<bool>* isDeregulatedGens, vector<int>* geneEx);

#endif /* PARAMETERS_H_ */
