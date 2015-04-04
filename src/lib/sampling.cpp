/*
 * sampling.cpp
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#include "../header/sampling.h"

void randomlyChooseSamplesDouble(TDoubleMatrix* original, TDoubleMatrix* subSamples, vector<int>* rrank, int numSamples){
	int totalGenes = original->size();

	for (int i = 0; i < totalGenes; ++i) {
		vector<double> row;
		for (int j = 0; j < numSamples; ++j) {
			row.push_back(original->at(i)[rrank->at(j)]);
		}
		subSamples->push_back(row);
	}
}


void randomlyChooseSamplesInteger(TIntegerMatrix* original, TIntegerMatrix* subSamples, vector<int>* rrank, int numSamples){
	int totalGenes = original->size();

	for (int i = 0; i < totalGenes; ++i) {
		vector<int> row;
		for (int j = 0; j < numSamples; ++j) {
			row.push_back(original->at(i)[rrank->at(j)]);
		}
		subSamples->push_back(row);
	}
}

