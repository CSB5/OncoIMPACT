/*
 * results.h
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include "utilities.h"
#include "merge_and_trim.h"
#include "impact_scores.h"

struct SampleDriver{
	string gene;
	string type;
	double impactScore;
	double aggregatedImpactScore;
	double driverFrequency;
	double mutationFrequency;
	string cancerCensus;
	string panCancer;
};

struct AggregatedDriver{
	string gene;
	double driverFrequency;
	double driverPointMutationFrequency;
	double driverDeletionFrequency;
	double driverAmplificationFrequency;
	string cancerCensus;
	string panCancer;
	double aggregatedImpactScore;
	double mutationFrequency;
	double pointMutationFrequency;
	double deletionFrequency;
	double amplificationFrequency;
};

void saveModules(vector< list<Module> > * modulesListOfAllSamples, string filename, vector<string>* geneIdToSymbol);

void printSampleDriverList(vector< vector<Driver> >* driversOfAllSamples, string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut, vector<int>* genesCNV,
		vector<double>* driverAggregatedScores, vector<int>* driversFrequency, vector<int>* mutationFrequency);
bool sortByImpactScore(const SampleDriver& first, const SampleDriver& second);

string getDriverType(int driverGeneId, int sampleId, TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix,
		vector<int>* genesPointMut, vector<int>* genesCNV);

void getMutationFrequency(TIntegerMatrix* originalMutationMatrix, vector<int>* mutationFrequency, vector<int>* genesMut);
void getDetailMutationFrequency(TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut, vector<int>* genesCNV,
		vector<int>* pointMutationFrequency, vector<int>* deletionFrequency, vector<int>* amplificationFrequency);

void printAggregatedDriverList(vector<int>* driverGeneIds, string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		vector<double>* driverAggregatedScores, vector<int>* driversFrequency, vector<int>* mutationFrequency,
		vector<int>* pointMutationDriversFrequency, vector<int>* deletionDriversFrequency, vector<int>* amplificationDriversFrequency,
		vector<int>* pointMutationFrequency, vector<int>* deletionFrequency, vector<int>* amplificationFrequency);
bool sortByAggregatedImpactScore(const AggregatedDriver& first, const AggregatedDriver& second);

#endif /* RESULTS_H_ */
