/*
 * impact_scores.cpp
 *
 *  Created on: Apr 14, 2015
 *      Author: nok
 */

#include "../header/impact_scores.h"
#include "../header/utilities.h"
#include "../header/input.h"
#include <cmath>
#include <iostream>

void calculateImpactScoresForAllSamples(vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<Driver> >* driversOfAllSamples, TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* GenesEx,
		int totalGenes, double F, vector<string>* geneIdToSymbol, string filename, vector<string>* sampleIdToName){

	int totalSamples = modulesListOfAllSamples->size();

	//map the gene id to the row id in the gene expression matrix
	vector<int> rowId(totalGenes);
	for (int i = 0; i < totalGenes; ++i) {
		rowId[i] = findIndex(GenesEx, i);
	}

	//OUTPUT: print drivers and impact scores for all samples
	vector<string>* outputDrivers = new vector<string>;
	outputDrivers->push_back("SAMPLE_ID\tDRIVER\tIMPACT_SCORE\tIS_DEREGULATED\tMODULE_SIZE\tNUM_DRIVERS");

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {

		list<Module> modules = modulesListOfAllSamples->at(i);

		//for each module
		for (list<Module>::iterator it = modules.begin(); it != modules.end(); it++) {
			Module module = *it;
			double score = 0;

			int moduleSize = 0;
			int numDrivers = module.driverGeneIds.size();

			//sum up the fold change

			for(list<int>::iterator git = module.explainedGeneIdsUpDown.begin(); git != module.explainedGeneIdsUpDown.end(); git++){
				int currentExplainedGeneId = *git;
				if(currentExplainedGeneId < totalGenes){	//up
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentExplainedGeneId])[i]);
				}else{										//down
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentExplainedGeneId - totalGenes])[i]);
				}
				moduleSize++;
			}

			for(list<int>::iterator git = module.phenotypeGeneIdsUpDown.begin(); git != module.phenotypeGeneIdsUpDown.end(); git++){
				int currentPhenotypeGeneId = *git;
				if(currentPhenotypeGeneId < totalGenes){	//up
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentPhenotypeGeneId])[i]);
				}else{										//down
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentPhenotypeGeneId - totalGenes])[i]);
				}
				moduleSize++;
			}

			//save the drivers and their scores
			for(list<int>::iterator git = module.driverGeneIds.begin(); git != module.driverGeneIds.end(); git++){
				Driver driver;
				driver.geneId = *git;
				driver.sampleId = i;
				driver.impactScore = score;

				//OUTPUT: print drivers and impact scores for all samples (cont.)
				string str = sampleIdToName->at(i) + "\t" + geneIdToSymbol->at(*git) + "\t" + doubleToStr(score, 3) + "\t";

				//check if the driver gene is also a deregulated gene
				if(rowId[*git] != -1 and fabs(originalGeneExpressionMatrix->at(rowId[*git])[i]) >= F){
					driver.isDeregulated = true;
					str = str + "1\t" + intToStr(moduleSize) + "\t" + intToStr(numDrivers);
				}else{
					driver.isDeregulated = false;
					str = str + "0\t" + intToStr(moduleSize) + "\t" + intToStr(numDrivers);
				}

				driversOfAllSamples->at(i).push_back(driver);

				outputDrivers->push_back(str);
			}

		}	//end for each module

	}

	//OUTPUT: print drivers and impact scores for all samples (cont.)
	writeStrVector(filename.c_str(), outputDrivers);
	delete outputDrivers;
}

void calculateImpactScoresForAllInputSamples(int totalInputSamples, vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<Driver> >* driversOfAllSamples, TDoubleMatrix* originalGeneExpressionMatrix, vector<int>* GenesEx,
		int totalGenes, double F, vector<string>* geneIdToSymbol, string filename, vector<string>* sampleIdToName){

	//map the gene id to the row id in the gene expression matrix
	vector<int> rowId(totalGenes);
	for (int i = 0; i < totalGenes; ++i) {
		rowId[i] = findIndex(GenesEx, i);
	}

	//OUTPUT: print drivers and impact scores for all samples
	vector<string>* outputDrivers = new vector<string>;
	outputDrivers->push_back("SAMPLE_ID\tDRIVER\tIMPACT_SCORE\tIS_DEREGULATED\tMODULE_SIZE\tNUM_DRIVERS_IN_MODULE");

	//for each sample i
	for (int i = 0; i < totalInputSamples; ++i) {

		list<Module> modules = modulesListOfAllSamples->at(i);

		//for each module
		for (list<Module>::iterator it = modules.begin(); it != modules.end(); it++) {
			Module module = *it;
			double score = 0;

			int moduleSize = 0;
			int numDrivers = module.driverGeneIds.size();

			//sum up the fold change

			for(list<int>::iterator git = module.explainedGeneIdsUpDown.begin(); git != module.explainedGeneIdsUpDown.end(); git++){
				int currentExplainedGeneId = *git;
				if(currentExplainedGeneId < totalGenes){	//up
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentExplainedGeneId])[i]);
				}else{										//down
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentExplainedGeneId - totalGenes])[i]);
				}
				moduleSize++;
			}

			for(list<int>::iterator git = module.phenotypeGeneIdsUpDown.begin(); git != module.phenotypeGeneIdsUpDown.end(); git++){
				int currentPhenotypeGeneId = *git;
				if(currentPhenotypeGeneId < totalGenes){	//up
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentPhenotypeGeneId])[i]);
				}else{										//down
					score += fabs(originalGeneExpressionMatrix->at(rowId[currentPhenotypeGeneId - totalGenes])[i]);
				}
				moduleSize++;
			}

			//save the drivers and their scores
			for(list<int>::iterator git = module.driverGeneIds.begin(); git != module.driverGeneIds.end(); git++){
				Driver driver;
				driver.geneId = *git;
				driver.sampleId = i;
				driver.impactScore = score;

				//OUTPUT: print drivers and impact scores for all samples (cont.)
				string str = sampleIdToName->at(i) + "\t" + geneIdToSymbol->at(*git) + "\t" + doubleToStr(score, 3) + "\t";

				//check if the driver gene is also a deregulated gene
				if(rowId[*git] != -1 and fabs(originalGeneExpressionMatrix->at(rowId[*git])[i]) >= F){
					driver.isDeregulated = true;
					str = str + "1\t" + intToStr(moduleSize) + "\t" + intToStr(numDrivers);
				}else{
					driver.isDeregulated = false;
					str = str + "0\t" + intToStr(moduleSize) + "\t" + intToStr(numDrivers);
				}

				driversOfAllSamples->at(i).push_back(driver);

				outputDrivers->push_back(str);
			}

		}	//end for each module

	}

	//OUTPUT: print drivers and impact scores for all samples (cont.)
	writeStrVector(filename.c_str(), outputDrivers);
	delete outputDrivers;
}

void aggregateDriversAcrossSamples(vector< vector<Driver> >* driversOfAllSamples, vector<double>* driverAggregatedScores,
		vector<int>* driversFrequency, vector<string>* geneIdToSymbol, int totalGenes){
	vector<double> sumImpact(totalGenes);
	vector<bool> isDriver(totalGenes);

	int totalSamples = driversOfAllSamples->size();

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {
		vector<Driver> drivers = driversOfAllSamples->at(i);

		//for each driver j of the sample i
		for (unsigned int j = 0; j < drivers.size(); ++j) {
			isDriver[drivers[j].geneId] = true;
			driversFrequency->at(drivers[j].geneId)++;
			sumImpact[drivers[j].geneId] += drivers[j].impactScore;
		}
	}

	for (int i = 0; i < totalGenes; ++i) {
		driverAggregatedScores->at(i) = 0;	//initialization

		if(isDriver[i]){
			double aggregatedScore = sumImpact[i]/totalSamples;
			driverAggregatedScores->at(i) = aggregatedScore;
		}
	}
}

void getDetailDriversFreqeuncy(vector< vector<Driver> >* driversOfAllSamples,
		vector<int>* pointMutationDriversFrequency, vector<int>* deletionDriversFrequency, vector<int>* amplificationDriversFrequency,
		TIntegerMatrix* originalPointMutationsMatrix, TIntegerMatrix* originalCNVsMatrix,
		vector<int>* genesPointMut, vector<int>* genesCNV){

	int totalSamples = driversOfAllSamples->size();

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {
		vector<Driver> drivers = driversOfAllSamples->at(i);

		//for each driver j of the sample i
		for (unsigned int j = 0; j < drivers.size(); ++j) {
			int driverGeneId = drivers[j].geneId;

			//check point mutation, deletion, amplification
			bool hasPointMut;
			bool hasCNV;

			if(findIndex(genesPointMut, driverGeneId) != -1){
				hasPointMut = true;
			}else{
				hasPointMut = false;
			}

			if(findIndex(genesCNV, driverGeneId) != -1){
				hasCNV = true;
			}else{
				hasCNV = false;
			}

			int pointMutation = 0;
			if(hasPointMut){
				pointMutation = originalPointMutationsMatrix->at(findIndex(genesPointMut, driverGeneId))[i];
			}

			int cnv = 0;
			if(hasCNV){
				cnv = originalCNVsMatrix->at(findIndex(genesCNV, driverGeneId))[i];
			}

			//if driver j is a point mutation driver in sample i
			if(pointMutation > 0){
				pointMutationDriversFrequency->at(driverGeneId)++;
			}

			//if driver j is a deletion driver in sample i
			if(cnv > 0){
				amplificationDriversFrequency->at(driverGeneId)++;
			}

			//if driver j is a amplification driver in sample i
			if(cnv < 0){
				deletionDriversFrequency->at(driverGeneId)++;
			}

		}
	}
}



