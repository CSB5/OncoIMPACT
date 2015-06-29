/*
 * discovery.h
 *
 *  Created on: May 18, 2015
 *      Author: nok
 */

#ifndef ANNOTATOR_H_
#define ANNOTATOR_H_

#include <string>
#include "header/utilities.h"
#include "header/annotator_input_output.h"
#include "header/annotator_calculator.h"
#include "header/input.h"

int annotator(string mSigDbPath, string moduleFileName, string outputPrefix, double cutoff, int top, string geneListFileName) ;


#endif /* DISCOVERY_H_ */
