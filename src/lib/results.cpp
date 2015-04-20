/*
 * results.cpp
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#include <iostream>
#include <string>
#include "../header/results.h"
#include "../header/utilities.h"

void saveModules(vector<list<Module> > * modulesListOfAllSamples,
		string filename, vector<string>* geneIdToSymbol) {
	int totalSamples = modulesListOfAllSamples->size();
	vector<string> outputStr;
	outputStr.push_back("GENE\tGENE_TYPE\tSAMPLE_ID\tMODULE_ID");

	//get a list of nodes (gene symbol, type)
	vector<int> geneIds;
	for (int i = 0; i < totalSamples; ++i) {

		list<Module> modulesList = modulesListOfAllSamples->at(i);

		//for each module
		int j = 0;
		for (list<Module>::iterator it = modulesList.begin();
				it != modulesList.end(); it++) {

			//get driver genes
			for (list<int>::iterator g = it->driverGeneIds.begin();
					g != it->driverGeneIds.end(); g++) {
				string str = "" + geneIdToSymbol->at(*g) + "\tDRIVER\t"
						+ intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			//get phenotype genes
			for (list<int>::iterator g = it->phenotypeGeneIds.begin();
					g != it->phenotypeGeneIds.end(); g++) {
				string str = "" + geneIdToSymbol->at(*g) + "\tPHENOTYPE\t"
						+ intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			//get explained genes
			for (list<int>::iterator g = it->explainedGeneIds.begin();
					g != it->explainedGeneIds.end(); g++) {
				string str = "" + geneIdToSymbol->at(*g) + "\tEXPLAINED\t"
						+ intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			j++;
		}
	}

	//get a list of edges (gene symbol, gene symbol)
	//do not need now because the original network file can be used

	writeStrVector(filename.c_str(), &outputStr);
}

void printSampleDriverList(vector<vector<Driver> >* driversOfAllSamples,
		string pathname, vector<string>* geneIdToSymbol,
		vector<string>* sampleIdToName,
		TIntegerMatrix* originaloriginalPointMutationsMatrix,
		TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut,
		vector<int>* genesCNV, vector<double>* driverAggregatedScores,
		vector<int>* driversFrequency, vector<int>* mutationFrequency) {
	//TODO create directory function for linux and windows

	int totalSamples = driversOfAllSamples->size();

	for (int i = 0; i < totalSamples; ++i) {
		vector<Driver> drivers = driversOfAllSamples->at(i);
		int numDrivers = drivers.size();

		list<SampleDriver> sampleDriversList;

		//for each driver
		for (int j = 0; j < numDrivers; ++j) {
			SampleDriver driver;
			driver.gene = geneIdToSymbol->at(drivers[j].geneId);
			driver.type = getDriverType(drivers[j].geneId, i,
					originaloriginalPointMutationsMatrix, originalCNVsMatrix,
					genesPointMut, genesCNV);
			driver.impactScore = drivers[j].impactScore;
			driver.aggregatedImpactScore = driverAggregatedScores->at(
					drivers[j].geneId);
			driver.driverFrequency = 1.0
					* driversFrequency->at(drivers[j].geneId) / totalSamples;
			driver.mutationFrequency = 1.0
					* mutationFrequency->at(drivers[j].geneId) / totalSamples;
			driver.cancerCensus = "NA";
			driver.panCancer = "NA";

			sampleDriversList.push_back(driver);
		}

		string filename = pathname + sampleIdToName->at(i) + ".tsv";
		vector<string> outputStr;
		outputStr.push_back(
				"GENE\tTYPE\tSAMPLE_IMPACT\tDATA_SET_IMPACT\tDRIVER_FREQUENCY\tMUTATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER");

		sampleDriversList.sort(sortByImpactScore);
		for (list<SampleDriver>::iterator it = sampleDriversList.begin();
				it != sampleDriversList.end(); it++) {
			string str = it->gene + "\t" + it->type + "\t"
					+ doubleToStr(it->impactScore) + "\t"
					+ doubleToStr(it->aggregatedImpactScore) + "\t"
					+ doubleToStr(it->driverFrequency) + "\t"
					+ doubleToStr(it->mutationFrequency) + "\t"
					+ it->cancerCensus + "\t" + it->panCancer;
			outputStr.push_back(str);
		}

		writeStrVector(filename.c_str(), &outputStr);

	}
}

//Usage: sort(sampleDriversList.begin(), sampleDriversList.end(), sortByImpactScore);
bool sortByImpactScore(const SampleDriver& first, const SampleDriver& second) {
	if (first.impactScore > second.impactScore) {
		return true;
	} else {
		return false;
	}
}

void printAggregatedDriverList(vector<int>* driverGeneIds, string filename,
		vector<string>* geneIdToSymbol, vector<string>* sampleIdToName,
		vector<double>* driverAggregatedScores, vector<int>* driversFrequency,
		vector<int>* mutationFrequency,
		vector<int>* pointMutationDriversFrequency,
		vector<int>* deletionDriversFrequency,
		vector<int>* amplificationDriversFrequency,
		vector<int>* pointMutationFrequency, vector<int>* deletionFrequency,
		vector<int>* amplificationFrequency) {

	vector<string> outputStr;
	outputStr.push_back(
			"GENE\tDRIVER_FREQUENCY\tDRIVER_SNV_FREQUENCY\tDRIVER_DELTION_FREQUENCY\tDRIVER_AMPLIFICATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER\t"
					"IMPACT\tMUTATION_FREQUENCY\tSNV_FREQUENCY\tDELTION_FREQUENCY\tAMPLIFICATION_FREQUENCY");
	//TODO for aggregated results: outputStr.push_back("GENE\tDRIVER_FREQUENCY\tDRIVER_SNV_FREQUENCY\tDRIVER_DELTION_FREQUENCYDRIVER_AMPLIFICATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER\tIMPACT\tMUTATION_FREQUENCY\tSNV_FREQUENCY\tDELTION_FREQUENCY\tAMPLIFICATION_FREQUENCY");

	int totalSamples = sampleIdToName->size();
	int totalDrivers = driverGeneIds->size();

	list<AggregatedDriver> aggregatedDriversList;

	for (int i = 0; i < totalDrivers; ++i) {

		AggregatedDriver driver;
		int currentDriverGeneId = driverGeneIds->at(i);

		driver.gene = geneIdToSymbol->at(currentDriverGeneId);
		driver.driverFrequency = 1.0 * driversFrequency->at(currentDriverGeneId)
				/ totalSamples;
		driver.driverPointMutationFrequency = 1.0
				* pointMutationDriversFrequency->at(currentDriverGeneId)
				/ totalSamples;
		driver.driverDeletionFrequency = 1.0
				* deletionDriversFrequency->at(currentDriverGeneId)
				/ totalSamples;
		driver.driverAmplificationFrequency = 1.0
				* amplificationDriversFrequency->at(currentDriverGeneId)
				/ totalSamples;
		driver.cancerCensus = "NA";
		driver.panCancer = "NA";
		driver.aggregatedImpactScore = driverAggregatedScores->at(
				currentDriverGeneId);
		driver.mutationFrequency = 1.0
				* mutationFrequency->at(currentDriverGeneId) / totalSamples;

		//TODO create function to find these values
		driver.pointMutationFrequency = 1.0
				* pointMutationFrequency->at(currentDriverGeneId)
				/ totalSamples;
		driver.deletionFrequency = 1.0
				* deletionFrequency->at(currentDriverGeneId) / totalSamples;
		driver.amplificationFrequency = 1.0
				* amplificationFrequency->at(currentDriverGeneId)
				/ totalSamples;

		aggregatedDriversList.push_back(driver);

	}

	aggregatedDriversList.sort(sortByAggregatedImpactScore);

	for (list<AggregatedDriver>::iterator it = aggregatedDriversList.begin();
			it != aggregatedDriversList.end(); it++) {
		string str = it->gene + "\t" + doubleToStr(it->driverFrequency) + "\t"
				+ doubleToStr(it->driverPointMutationFrequency) + "\t"
				+ doubleToStr(it->driverDeletionFrequency) + "\t"
				+ doubleToStr(it->driverAmplificationFrequency) + "\t"
				+ it->cancerCensus + "\t" + it->panCancer + "\t"
				+ doubleToStr(it->aggregatedImpactScore) + "\t"
				+ doubleToStr(it->mutationFrequency) + "\t"
				+ doubleToStr(it->pointMutationFrequency) + "\t"
				+ doubleToStr(it->deletionFrequency) + "\t"
				+ doubleToStr(it->amplificationFrequency);
		outputStr.push_back(str);
	}

	writeStrVector(filename.c_str(), &outputStr);

}

bool sortByAggregatedImpactScore(const AggregatedDriver& first,
		const AggregatedDriver& second) {
	if (first.aggregatedImpactScore > second.aggregatedImpactScore) {
		return true;
	} else {
		return false;
	}
}

string getDriverType(int driverGeneId, int sampleId,
		TIntegerMatrix* originalPointMutationsMatrix,
		TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut,
		vector<int>* genesCNV) {
	bool hasPointMut;
	bool hasCNV;

	if (findIndex(genesPointMut, driverGeneId) != -1) {
		hasPointMut = true;
	} else {
		hasPointMut = false;
	}

	if (findIndex(genesCNV, driverGeneId) != -1) {
		hasCNV = true;
	} else {
		hasCNV = false;
	}

	int pointMutation = 0;
	if (hasPointMut) {
		pointMutation = originalPointMutationsMatrix->at(
				findIndex(genesPointMut, driverGeneId))[sampleId];
	}

	int cnv = 0;
	if (hasCNV) {
		cnv =
				originalCNVsMatrix->at(findIndex(genesCNV, driverGeneId))[sampleId];
	}

	string typeStr;

	if (cnv != 0) {
		if (cnv > 0) {
			typeStr = "AMPL";
			if (pointMutation > 0) {
				typeStr += "_MUT";
			}
		} else {
			typeStr = "DEL";
			if (pointMutation > 0) {
				typeStr += "_MUT";
			}
		}
	} else {
		if (pointMutation > 0) {
			typeStr += "MUT";
		}
	}

	return typeStr;
}

void getMutationFrequency(TIntegerMatrix* originalMutationMatrix,
		vector<int>* mutationFrequency, vector<int>* genesMut) {
	int totalSamples = originalMutationMatrix->at(0).size();
	int totalMutatedGenes = genesMut->size();

	for (int i = 0; i < totalMutatedGenes; ++i) {

		int currentMutatedGeneId = genesMut->at(i);
		int currentRowIdInMutationMatrix = findIndex(genesMut,
				currentMutatedGeneId);

		for (int j = 0; j < totalSamples; ++j) {
			if (originalMutationMatrix->at(currentRowIdInMutationMatrix)[j]
					> 0) {
				mutationFrequency->at(currentMutatedGeneId)++;}
			}
		}
	}

void getDetailMutationFrequency(TIntegerMatrix* originalPointMutationsMatrix,
		TIntegerMatrix* originalCNVsMatrix, vector<int>* genesPointMut,
		vector<int>* genesCNV, vector<int>* pointMutationFrequency,
		vector<int>* deletionFrequency, vector<int>* amplificationFrequency) {
	int totalSamples = originalCNVsMatrix->at(0).size();

	int totalPointMutationGenes = genesPointMut->size();
	for (int i = 0; i < totalPointMutationGenes; ++i) {
		int currentGeneId = genesPointMut->at(i);
		int currentRowIdInMutationMatrix = findIndex(genesPointMut,
				currentGeneId);

		for (int j = 0; j < totalSamples; ++j) {
			if (originalPointMutationsMatrix->at(currentRowIdInMutationMatrix)[j] > 0) {
				pointMutationFrequency->at(currentGeneId)++;}
			}
		}

	int totalCNVGenes = genesCNV->size();
	for (int i = 0; i < totalCNVGenes; ++i) {
		int currentGeneId = genesCNV->at(i);
		int currentRowIdInMutationMatrix = findIndex(genesCNV, currentGeneId);

		for (int j = 0; j < totalSamples; ++j) {
			if (originalCNVsMatrix->at(currentRowIdInMutationMatrix)[j] != 0) {
				if (originalCNVsMatrix->at(currentRowIdInMutationMatrix)[j] > 0) {
					amplificationFrequency->at(currentGeneId)++;
				}
				else {
					deletionFrequency->at(currentGeneId)++;
				}
			}
		}
	}
}

void saveJSDivergences(vector<JSDivergence>* jsDivergences, string filename){
	vector<string> outputStr;
	outputStr.push_back("L\tD\tF\tJS Divergence");

	int total = jsDivergences->size();
	int D, L;
	double F, divergence;
	for (int i = 0; i < total; ++i) {
		D = jsDivergences->at(i).D;
		L = jsDivergences->at(i).L;
		F = jsDivergences->at(i).F;
		divergence = jsDivergences->at(i).divergence;

		string str = intToStr(L) + "\t" + intToStr(D) + "\t" + doubleToStr(F) + "\t" + doubleToStr(divergence);
		outputStr.push_back(str);
	}

	writeStrVector(filename.c_str(), &outputStr);
}
