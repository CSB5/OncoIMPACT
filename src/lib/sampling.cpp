/*
 * sampling.cpp
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#include "../header/sampling.h"
#include <iostream>

void randomlyChooseSamplesDouble(TDoubleMatrix* original, TDoubleMatrix* subSamples, vector<int>* rrank, int numSamples){
	int totalGenes = original->size();

	for (int i = 0; i < totalGenes; ++i) {
		vector<double> rowDouble;
		for (int j = 0; j < numSamples; ++j) {
			rowDouble.push_back(original->at(i)[rrank->at(j)]);
		}
		subSamples->push_back(rowDouble);
	}
}


void randomlyChooseSamplesInt(const TIntegerMatrix* original, TIntegerMatrix* subSamples, const vector<int>* rrank, int numSamples){
	int totalGenes = original->size();

	for (int i = 0; i < totalGenes; ++i) {
		vector<int> rowInt;
		for (int j = 0; j < numSamples; ++j) {
			rowInt.push_back(original->at(i)[rrank->at(j)]);
		}
		subSamples->push_back(rowInt);
	}

}

