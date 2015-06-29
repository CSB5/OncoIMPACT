/*
 * impact_scores.h
 *
 *  Created on: Apr 14, 2015
 *      Author: nok
 */

#ifndef IMPACT_SCORES_H_
#define IMPACT_SCORES_H_

#include "../header/merge_and_trim.h"
#include "../header/input.h"

struct Driver{
	int geneId;
	int sampleId;
	int moduleId;
	double impactScore;
	double aggregatedImpactScore;
	bool isDeregulated;
};

void calculateImpactScoresForAllSamples(vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<Driver> >* drivers, TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* GenesEx,
		int totalGenes, double F, vector<string>* geneIdToSymbol, string filname, vector<string>* sampleIdToName);
void calculateImpactScoresForAllInputSamples(int totalInputSamples, vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<Driver> >* driversOfAllSamples, TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* GenesEx,
		int totalGenes, double F, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName);

void aggregateDriversAcrossSamples(vector< vector<Driver> >* driversOfAllSamples, vector<double>* driverAggregatedScores,
		vector<int>* driversFrequency, vector<string>* geneIdToSymbol);

void getDetailDriversFreqeuncy(vector< vector<Driver> >* driversOfAllSamples,
		vector<int>* pointMutationDriversFrequency, vector<int>* deletionDriversFrequency, vector<int>* amplificationDriversFrequency,
		TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix,
		vector<int>* genesPointMut, vector<int>* genesCNV);

#endif /* IMPACT_SCORES_H_ */
