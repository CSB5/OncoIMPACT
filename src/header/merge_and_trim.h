/*
 * merge_and_trim.h
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#ifndef MERGE_AND_TRIM_H_
#define MERGE_AND_TRIM_H_

#include <list>
#include <algorithm>
#include "utilities.h"
#include "explained_genes.h"
#include "driver_genes.h"

struct Module{
	int moduleId; //on sample level
	list<int> driverGeneIds;
	list<int> phenotypeGeneIdsUpDown;
	list<int> explainedGeneIdsUpDown;
};

//merge module moduleId into the currentModuleId and deleted moduleId
void mergeModules(int currentModuleId, int moduleId, list<Module>* modulesList);


void findModulesInAllSamples(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenesUpDown, vector<DriverGene>* driverGenes, vector<int>* phenotypeGeneIds, int mode);
void findModulesInAllInputSamples(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenesUpDown, vector<DriverGene>* driverGenes, vector<int>* phenotypeGeneIds, int mode);


void trimSomeExplainedGenes(vector< list<Module> >* modulesListOfAllSamples, TIntAdjList* network, int L, int D, vector<string>* geneIdToSymbol);

void trimModule(Module* module, TIntAdjList* network, int L, int D, vector<string>* geneIdToSymbol);
void findShortestPath(int geneId, vector<int>* shortestPathToDeiverGenes, vector<int>* shortestPathToPhenotypeGenes,
		vector<int>* driverGeneIdsOfTheShortestPath, vector<int>* phenotypGeneIdsOfTheShortestPath,
		vector<bool>* isDriverGeneInThisModule, vector<bool>* isPhenotypeGeneInThisModule, TIntAdjList* network,
		int D, vector<bool>* isGeneInThisModule, vector<string>* geneIdToSymbol);

#endif /* MERGE_AND_TRIM_H_ */
