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

void saveModules(vector< list<Module> > * modulesListOfAllSamples, string filename, vector<string>* geneIdToSymbol);
string geneIdsToNodesList(const vector<string>& geneSymbols);



#endif /* RESULTS_H_ */
