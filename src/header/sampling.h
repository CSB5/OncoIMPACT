/*
 * sampling.h
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#ifndef SAMPLING_H_
#define SAMPLING_H_

#include "utilities.h"

/*
 * Functions for sampling data
 */

void randomlyChooseSamplesDouble(TDoubleMatrix* original, TDoubleMatrix* subSamples, vector<int>* rrank, int numSamples);
void randomlyChooseSamplesInt(const TIntegerMatrix* original, TIntegerMatrix* subSamples, const vector<int>* rrank, int numSamples);

#endif /* SAMPLING_H_ */
