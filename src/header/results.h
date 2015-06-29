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
#include "parameters.h"

struct SampleDriver{
	int geneId;
	string gene;
	string type;
	double impactScore;
	double databaseImpactScore;
	double databaseDriverFrequency;
	double databaseMutationFrequency;
	bool isCancerCensus;
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

struct AggregatedDriverForInputSample{
	int geneId;
	string gene;
	double aggregatedImpactScore;
	double dbImpactScore;
	double dbDriverFrequency;
	double dbDriverPointMutationFrequency;
	double dbDriverDeletionFrequency;
	double dbDriverAmplificationFrequency;
	double dbMutationFrequency;
	double dbPointMutationFrequency;
	double dbDeletionFrequency;
	double dbAmplificationFrequency;
	string cancerCensus;
};

struct ExplainedGeneDetail{
	string gene;
	int degree;
	int numSampleDeregulated;
	int numSampleExplained;
	bool isPhenotype;
	double pValue;
};

void saveModules(vector<list<Module> > * modulesListOfAllSamples, vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName);
void saveModulesOfInputSamples(vector<list<Module> > * modulesListOfAllSamples, vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName, int numInputSamples);
void saveModulesCytoscape(vector<list<Module> > * modulesListOfAllSamples,
		string filename, vector<string>* geneIdToSymbol);

void printSampleDriverList(vector< vector<Driver> >* driversOfAllSamples, string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut, vector<int>* genesCNV,
		vector<double>* driverAggregatedScores, vector<int>* driversFrequency, vector<int>* mutationFrequency, vector<bool>* isCancerBenchmarkGenes);
void printSampleDriverListForInputSamples(int totalInputSamples, vector<vector<Driver> >* driversOfAllSamples,
		string pathname, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		vector<MutatedGeneFromFile>* mutatedGeneFromFile, vector<string>* cancerBenchmarkGeneNames,
		TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut, vector<int>* genesCNV,
		map<int, string>* drugIdToName, vector< vector<int> >* geneDrugsAssocList);

bool sortByImpactScore(const SampleDriver& first, const SampleDriver& second);

string getDriverType(int driverGeneId, int sampleId, TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix,
		vector<int>* genesPointMut, vector<int>* genesCNV);

void getMutationFrequency(TIntegerMatrix* originalMutationMatrix, vector<int>* mutationFrequency, vector<int>* genesMut);
void getDetailMutationFrequency(TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut, vector<int>* genesCNV,
		vector<int>* pointMutationFrequency, vector<int>* deletionFrequency, vector<int>* amplificationFrequency);

void printAggregatedDriverList(vector<DriverGene>* driverGenes, string filename, vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		vector<double>* driverAggregatedScores, vector<int>* driversFrequency, vector<int>* mutationFrequency,
		vector<int>* pointMutationDriversFrequency, vector<int>* deletionDriversFrequency, vector<int>* amplificationDriversFrequency,
		vector<int>* pointMutationFrequency, vector<int>* deletionFrequency, vector<int>* amplificationFrequency, vector<bool>* isCancerBenchmarkGenes);
void printAggregatedDriverListForInputSamples(vector<DriverGene>* driverGenes, string filename,
		vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		vector<double>* driverAggregatedScores, vector<MutatedGeneFromFile>* mutatedGeneFromFile,
		vector<string>* cancerBenchmarkGeneNames, map<int, string>* drugIdToName, vector< vector<int> >* geneDrugsAssocList);
bool sortByAggregatedImpactScore(const AggregatedDriverForInputSample& first,
		const AggregatedDriverForInputSample& second);

void saveJSDivergences(vector<JSDivergence>* jsDivergences, string filename);

void printExplinedGenesFrequencyAndPhonotype(vector<int>* explainedGenesFrequencyRealUpDown, vector<double>* pValues, vector<bool>* isPhenotypeGenes,
		vector<string>* geneIdToSymbol, TIntAdjList* network, TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* genesEx, double F, string filename);
int getNumSamplesOfDeregulatedGene(TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* genesEx, double F, int geneIdUpDown, int totalGenesUpDown);

#endif /* RESULTS_H_ */
