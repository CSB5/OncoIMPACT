/*
 * merge_and_trim.cpp
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#include <iostream>
#include <algorithm>
#include <queue>
#include "../header/merge_and_trim.h"

//merge module moduleId into the currentModuleId and deleted moduleId
void mergeModules(int currentModuleId, int moduleId, list<Module>* modulesList){

	//find currentModuleId
	list<Module>::iterator currentRef;
	for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
		if(it->moduleId == currentModuleId){
			currentRef = it;
			break;
		}
	}

	//find moduleId to be merged
	list<Module>::iterator toBeMergedRef;
	for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
		if(it->moduleId == moduleId){
			toBeMergedRef = it;
			break;
		}
	}

	//merge divers genes
	for(list<int>::iterator it = toBeMergedRef->driverGeneIds.begin(); it != toBeMergedRef->driverGeneIds.end(); it++){
		currentRef->driverGeneIds.push_back(*it);
	}

	//merge phenotype genes
	for(list<int>::iterator it = toBeMergedRef->phenotypeGeneIds.begin(); it != toBeMergedRef->phenotypeGeneIds.end(); it++){
		currentRef->phenotypeGeneIds.push_back(*it);
	}
	currentRef->phenotypeGeneIds.unique();

	//merge explained genes
	for(list<int>::iterator it = toBeMergedRef->explainedGeneIds.begin(); it != toBeMergedRef->explainedGeneIds.end(); it++){
		currentRef->explainedGeneIds.push_back(*it);
	}
	currentRef->explainedGeneIds.unique();

	//delete toBeMergedRef after merging
	modulesList->erase(toBeMergedRef);

}

void findModulesInAllSamples(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< list<Module> >* modulesListOfAllSamples,
		vector<vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenes, vector<bool>* isDriverGenes, vector<int>* phenotypeGeneIds){
	int totalSamples = mutatedGeneIdsListReal->size();
	int totalGenes = isPhenotypeGenes->size();

	for (int i = 0; i < totalSamples; ++i) {

		//get a list of modules of sample i
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal->at(i);

		//get a list of driver gene ids of sample i
		list<int> driverGeneIdsForASample;
		copy(mutatedGeneIdsListReal->at(i).begin(), mutatedGeneIdsListReal->at(i).end(), std::back_inserter(driverGeneIdsForASample));

		//delete modules in which the mutated gene in not a driver
		list<int>::iterator it = driverGeneIdsForASample.begin();
		while(it != driverGeneIdsForASample.end()){
			int currentMutatedGeneId = *it;
			if(isDriverGenes->at(currentMutatedGeneId)){
				it++;	//go to next element
			}else{
				it = driverGeneIdsForASample.erase(it);	//remove from the list (the pointer automatically go to the next element)
			}
		}

		//create module (prepare for merging and trimming)
		int numModules = driverGeneIdsForASample.size();

		//cout << "sample #" << i << " has " << numModules << " modules\n";

		vector<Module> modules(numModules);
		int j = 0;
		for(list<int>::iterator it = driverGeneIdsForASample.begin(); it != driverGeneIdsForASample.end(); it++, j++){

			int currentMutatedGeneId = *it;
			vector<int>* explainedGenesFreqency = mutatedAndExplainedGenes[currentMutatedGeneId].explainedGenesFreqency;

			modules[j].moduleId = j;
			modules[j].driverGeneIds.push_back(currentMutatedGeneId);
			for (int k = 0; k < totalGenes; ++k) {
				if(explainedGenesFreqency->at(k) > 0){
					if(isPhenotypeGenes->at(k)){
						modules[j].phenotypeGeneIds.push_back(k);
					}else{
						modules[j].explainedGeneIds.push_back(k);
					}
				}
			}
		}

		list<Module>* modulesList = &modulesListOfAllSamples->at(i);
		copy(modules.begin(), modules.end(), std::back_inserter(*modulesList));


		//TODO merge modules that share phenotype genes
		//for each phenotype gene, find which modules contain it
		int numPhenotypeGenes = phenotypeGeneIds->size();
		for (int j = 0; j < numPhenotypeGenes; ++j) {

			int currentPhenotypeGeneId = phenotypeGeneIds->at(j);
			vector<int> moduleIdsToMerge;
			bool isPhenotypeGeneInThisSample = false;

			//for each module, find if it contains phenotype gene currentPhenotypeGeneId
			for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
				list<int> phenotypeGeneIdsOfAModule = it->phenotypeGeneIds;
				bool found = (find(phenotypeGeneIdsOfAModule.begin(), phenotypeGeneIdsOfAModule.end(), currentPhenotypeGeneId) != phenotypeGeneIdsOfAModule.end());
				if(found){ // phenotypeGeneIds is a phenotype genes in this module of sample i
					moduleIdsToMerge.push_back(it->moduleId);
					isPhenotypeGeneInThisSample = true;
				}
			}

			if (isPhenotypeGeneInThisSample) {	//then merge modules
				int numModulesToMerge = moduleIdsToMerge.size();
				if (numModulesToMerge > 1) {
					int currentModuleId = moduleIdsToMerge[0];
					for (int k = 1; k < numModulesToMerge; ++k) {
						int moduleId = moduleIdsToMerge[k];
						//merge module moduleId into the currentModuleId and deleted moduleId
						mergeModules(currentModuleId, moduleId, modulesList);
					}
				}	//else do not need to merge
			}
		}

		//cout << "sample #" << i << " has " << modulesListOfAllSamples->at(i).size() << " modules after merging\n";

	}
}

void trimSomeExplainedGenes(vector< list<Module> >* modulesListOfAllSamples, TIntAdjList* network, int L, int D){
	int numSamples = modulesListOfAllSamples->size();
	//for each sample i
	for (int i = 0; i < numSamples; ++i) {
		list<Module>* modulesList = &modulesListOfAllSamples->at(i);
		//for each module
		for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
			Module* currentModule = &(*it);
				trimModule(currentModule, network, L, D);
		}

		if(i == 1)
			break; //TEST for fist two samples
	}


}

void trimModule(Module* module, TIntAdjList* network, int L, int D){
	int totalGenes = network->size();

	//find degree of all explained genes
	list<int>* explainedGeneIdsList = &module->explainedGeneIds;
	vector<int> explainedGensIds;
	copy( explainedGeneIdsList->begin(), explainedGeneIdsList->end(), explainedGensIds.begin() );
	int numExplainedGenes = explainedGensIds.size();
	vector<bool> isExplainedGeneInThisModule(totalGenes);
	for (int i = 0; i < numExplainedGenes; ++i) {
		isExplainedGeneInThisModule[explainedGensIds[i]] = true;
	}

	vector<bool> isPhenotypeGeneInThisModule(totalGenes);
	list<int>* phenotypeGeneIds = &module->phenotypeGeneIds;
	for(list<int>::iterator it = phenotypeGeneIds->begin(); it != phenotypeGeneIds->end(); it++){
		isPhenotypeGeneInThisModule[*it] = true;
	}

	vector<bool> isDriverGeneInThisModule(totalGenes);
	list<int>* driverGeneIds = &module->driverGeneIds;
	for(list<int>::iterator it = driverGeneIds->begin(); it != driverGeneIds->end(); it++){
		isDriverGeneInThisModule[*it] = true;
	}

	vector<bool> isGeneInThisModule(totalGenes);
	for (int i = 0; i < totalGenes; ++i) {
		if(isExplainedGeneInThisModule[i] or
				isPhenotypeGeneInThisModule[i] or
				isDriverGeneInThisModule[i])
			isGeneInThisModule[i] = true;
	}

	//for each explained gene i find the shortest path to any driver gene and any phenotype genes
	vector<int> shortestPathToDeiverGenes;
	vector<int> shortestPathToPhenotypeGenes;
	for (int i = 0; i < numExplainedGenes; ++i) {
		int currentExplainedGeneId = explainedGensIds[i];
		findShortestPath(currentExplainedGeneId, &shortestPathToDeiverGenes, &shortestPathToPhenotypeGenes,
				&isDriverGeneInThisModule, &isPhenotypeGeneInThisModule, network, D, &isGeneInThisModule);
	}

	//for each explained gene, check if it belong to at least one path (with length < L) between a mutated gene and a phenotype gene
	vector<int> geneIdToBeRemoved;
	for (int i = 0; i < numExplainedGenes; ++i) {
		int currentExplainedGeneId = explainedGensIds[i];
		int disToDriver = shortestPathToDeiverGenes[currentExplainedGeneId];
		int disToPhenotype = shortestPathToPhenotypeGenes[currentExplainedGeneId];
		if(disToDriver < 0 or disToPhenotype < 0){
			geneIdToBeRemoved.push_back(currentExplainedGeneId);
		}else if(disToDriver + disToPhenotype > L){
			geneIdToBeRemoved.push_back(currentExplainedGeneId);
		}
	}

	int numDelete = geneIdToBeRemoved.size();
	for(int i = 0; i < numDelete; i++){
		list<int>::iterator dit = find(explainedGeneIdsList->begin(), explainedGeneIdsList->end(), geneIdToBeRemoved[i]);
		explainedGeneIdsList->erase(dit);
	}
}

void findShortestPath(int geneId, vector<int>* shortestPathToDeiverGenes, vector<int>* shortestPathToPhenotypeGenes,
		vector<bool>* isDriverGeneInThisModule, vector<bool>* isPhenotypeGeneInThisModule, TIntAdjList* network,
		int D, vector<bool>* isGeneInThisModule){
	// Mark all the vertices as not visited
	int totalGenes = network->size();
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	bool foundDriver = false;
	bool foundPhenotype = false;
	int distantToDriver = -1;
	int distantToPhenotype = -1;

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes); // store level of each nodes
	// Mark the current node as visited
	q.push(geneId);
	visited[geneId] = true;
	int currentGeneId;
	int currentLevel = 0;

	//consider all nodes
	while (!q.empty()) {
		//read the root node
		currentGeneId = q.front();
		q.pop();
		currentLevel = levels[currentGeneId];

		if(foundDriver and foundPhenotype){
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			int currentNeighborGeneId = (*network)[currentGeneId][j];
			if (isGeneInThisModule->at(currentNeighborGeneId) and !visited[currentNeighborGeneId]) {
				// check conditions before push to the queue
				// |fold change| > F

				if(!foundDriver and isDriverGeneInThisModule->at(currentNeighborGeneId)){
					distantToDriver = levels[currentGeneId] + 1;
					foundDriver = true;
				}

				if(!foundPhenotype and isPhenotypeGeneInThisModule->at(currentNeighborGeneId)){
					distantToPhenotype = levels[currentGeneId] + 1;
					foundPhenotype = true;
				}

				int degree = getNodeDegree(network, currentNeighborGeneId);
				if (degree <= D and !isDriverGeneInThisModule->at(currentNeighborGeneId) and !isPhenotypeGeneInThisModule->at(currentNeighborGeneId)) {
					// push explained genes to queue
					q.push(currentNeighborGeneId);
					// save level
					levels[currentGeneId] = currentLevel + 1;
				}

				visited[currentNeighborGeneId] = true;
			}
		}
	}


	delete[] visited;

	shortestPathToDeiverGenes->push_back(distantToDriver);
	shortestPathToPhenotypeGenes->push_back(distantToPhenotype);

}

