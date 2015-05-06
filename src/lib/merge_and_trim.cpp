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
		//TODO check this
//		bool found = find(toBeMergedRef->driverGeneIds.begin(), toBeMergedRef->driverGeneIds.end(), *it) != toBeMergedRef->driverGeneIds.end();
//		if(!found){
			currentRef->phenotypeGeneIds.push_back(*it);
//		}
	}
	currentRef->phenotypeGeneIds.sort();
	currentRef->phenotypeGeneIds.unique();

	//merge explained genes
	for(list<int>::iterator it = toBeMergedRef->explainedGeneIds.begin(); it != toBeMergedRef->explainedGeneIds.end(); it++){
		//TODO check this
//		bool found = find(toBeMergedRef->driverGeneIds.begin(), toBeMergedRef->driverGeneIds.end(), *it) != toBeMergedRef->driverGeneIds.end();
//		if(!found){
			currentRef->explainedGeneIds.push_back(*it);
//		}
	}
	currentRef->explainedGeneIds.sort();
	currentRef->explainedGeneIds.unique();

	//delete toBeMergedRef after merging
	modulesList->erase(toBeMergedRef);

}

//void findModulesInAllSamples(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< list<Module> >* modulesListOfAllSamples,
//		vector<vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenes, vector<bool>* isDriverGenes, vector<int>* phenotypeGeneIds){
void findModulesInAllSamples(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< list<Module> >* modulesListOfAllSamples,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenes, vector<bool>* isDriverGenes, vector<int>* phenotypeGeneIds){
	int totalSamples = mutatedAndExplainedGenesListReal->size();
	int totalGenes = isPhenotypeGenes->size();
	int totalGenesUpDown = totalGenes * 2;

	for (int i = 0; i < totalSamples; ++i) {

		//get a list of modules of sample i (this vector has size of totalGene)
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes = mutatedAndExplainedGenesListReal->at(i);

		//get a list of driver gene ids of sample i
		list<int> driverGeneIdsForASample;	//this initially contain all the mutated genes of the current sample
		vector<int> mutatedGeneIdsList = mutatedGeneIdsListReal->at(i);


		int mode = 0;	//0 = sensitive, 1 = stringent
		if(mode == 0){
			//TODO SENSITIVE MODE
			for (unsigned j = 0; j < mutatedGeneIdsList.size(); ++j) {
				int currentMutatedGeneId = mutatedGeneIdsList[j];
				if (isDriverGenes->at(currentMutatedGeneId)) { //current mutated gene is a driver in at least one sample

					//check if current mutated gene is connected to at least one phenotype in this sample
					bool isConnectedWithPhenotype = false;
					vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
					for (int k = 0; k < totalGenesUpDown; ++k) {
						if(isExplainedGenesUpDown->at(k)){
							if(k < totalGenes){
								if(isPhenotypeGenes->at(k)){
									isConnectedWithPhenotype = true;
								}
							}else{
								if(isPhenotypeGenes->at(k-totalGenes)){
									isConnectedWithPhenotype = true;
								}
							}
						}
					}
					//then the current mutated gene is a possible driver
					if(isConnectedWithPhenotype){
						driverGeneIdsForASample.push_back(currentMutatedGeneId);
					}
				}
			}
		}else{
			//TODO STRINGENT MODE
		}

		//create module (prepare for merging and trimming)
		list<Module> modules;
		int j = 0;
		for(list<int>::iterator it = driverGeneIdsForASample.begin(); it != driverGeneIdsForASample.end(); it++, j++){

			int currentDriverGeneId = *it;

			//get the list of explained genes of the current driver gene in sample i
			vector<bool>* isExplainedGenesUpDown = mutatedAndExplainedGenes[currentDriverGeneId].isExplainedGenesUpDown;

			Module currentModule;

			currentModule.moduleId = j;
			currentModule.driverGeneIds.push_back(currentDriverGeneId);

			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(k < totalGenes){	//up
					if(isExplainedGenesUpDown->at(k)){	//k is an explained gene
						if(isPhenotypeGenes->at(k)){
							currentModule.phenotypeGeneIds.push_back(k);
						}else{
							currentModule.explainedGeneIds.push_back(k);
						}
					}
				}else{				//down
					if(isExplainedGenesUpDown->at(k)){	//k is an explained gene
						if(isPhenotypeGenes->at(k - totalGenes)){
							currentModule.phenotypeGeneIds.push_back(k - totalGenes);
						}else{
							currentModule.explainedGeneIds.push_back(k - totalGenes);
						}
					}
				}
			}

			//add the module to module list
			modules.push_back(currentModule);
		}

		//cout << "sample #" << i << " has " << modules.size() << " driver genes\n";

		//delete the drivers that do not connect to any phenotype genes in this sample
		//on sample level, if the driver genes do not connect to any phenotype genes on the sample, they are not drivers
		for(list<Module>::iterator it = modules.begin(); it != modules.end(); ){
			Module currentModule  = *it;
			if(currentModule.phenotypeGeneIds.size() == 0){		//there is no phenotype gene in this module
				it = modules.erase(it);
			}else{												//there are at least one phenotype gene in this module
				it++;
			}
		}

		//cout << "sample #" << i << " has " << modules.size() << " modules before merging\n";

		//add modules data to the modulesListOfAllSamples
		list<Module>* modulesList = &modulesListOfAllSamples->at(i);
		copy(modules.begin(), modules.end(), std::back_inserter(*modulesList));

		//merge modules that share phenotype genes
		//for each phenotype gene, find which modules contain it
		int numPhenotypeGenes = phenotypeGeneIds->size();
		for (int j = 0; j < numPhenotypeGenes; ++j) {

			int currentPhenotypeGeneId = phenotypeGeneIds->at(j);
			vector<int> moduleIdsToMerge;
			bool isPhenotypeGeneInThisSample = false;

			//for each module, find if it contains phenotype gene currentPhenotypeGeneId
			for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
				list<int> phenotypeGeneIdsOfAModule = it->phenotypeGeneIds;
				bool found = (find(phenotypeGeneIdsOfAModule.begin(), phenotypeGeneIdsOfAModule.end(), currentPhenotypeGeneId)
						!= phenotypeGeneIdsOfAModule.end());
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
				}	//no sharing phenotype gene, do not need to merge
			}
		}

		//cout << "sample #" << i << " has " << modulesListOfAllSamples->at(i).size() << " modules after merging\n";

	}
}

void trimSomeExplainedGenes(vector< list<Module> >* modulesListOfAllSamples, TIntAdjList* network, int L, int D, vector<string>* geneIdToSymbol){
	int numSamples = modulesListOfAllSamples->size();

	//for each sample i
	for (int i = 0; i < numSamples; ++i) {
		list<Module>* modulesList = &modulesListOfAllSamples->at(i);
		//for each module
		for(list<Module>::iterator it = modulesList->begin(); it != modulesList->end(); it++){
			Module* currentModule = &(*it);
			trimModule(currentModule, network, L, D, geneIdToSymbol);
		}
	}
}

void trimModule(Module* module, TIntAdjList* network, int L, int D, vector<string>* geneIdToSymbol){
	int totalGenes = network->size();

	//get a list of explained gene ids
	list<int>* explainedGeneIdsList = &module->explainedGeneIds;

	//create a bool vector for telling which genes are explained in the current module of sample i
	vector<bool> isExplainedGeneInThisModule(totalGenes);
	for (list<int>::iterator it = explainedGeneIdsList->begin(); it != explainedGeneIdsList->end(); it++) {
		isExplainedGeneInThisModule[*it] = true;
	}

	//get a list of phenotype gene ids and create a bool vector
	vector<bool> isPhenotypeGeneInThisModule(totalGenes);
	list<int>* phenotypeGeneIds = &module->phenotypeGeneIds;
	for(list<int>::iterator it = phenotypeGeneIds->begin(); it != phenotypeGeneIds->end(); it++){
		isPhenotypeGeneInThisModule[*it] = true;
	}

	//get a list of driver gene ids and create a bool vector
	vector<bool> isDriverGeneInThisModule(totalGenes);
	list<int>* driverGeneIds = &module->driverGeneIds;
	for(list<int>::iterator it = driverGeneIds->begin(); it != driverGeneIds->end(); it++){
		isDriverGeneInThisModule[*it] = true;
	}

	//get a list of gene ids belong to this module
	int numGenes = 0;
	vector<bool> isGeneInThisModule(totalGenes);
	for (int i = 0; i < totalGenes; ++i) {
		if(isExplainedGeneInThisModule[i] or
				isPhenotypeGeneInThisModule[i] or
				isDriverGeneInThisModule[i]){
			isGeneInThisModule[i] = true;
			numGenes++;
		}
	}
	//cout << "\tthere are " << numGenes << " genes in this module" << endl;

	//for each explained gene, find the shortest path to any driver gene and any phenotype genes
	vector<int> shortestPathToDriverGenes;		//has the same size as explainedGeneIdsList
	vector<int> shortestPathToPhenotypeGenes;	//has the same size as explainedGeneIdsList
	for (list<int>::iterator it = explainedGeneIdsList->begin(); it != explainedGeneIdsList->end(); it++) {
		int currentExplainedGeneId = *it;
		//cout << "finding the shortest paths for explained gene " << geneIdToSymbol->at(currentExplainedGeneId) << endl;
		findShortestPath(currentExplainedGeneId, &shortestPathToDriverGenes, &shortestPathToPhenotypeGenes,
				&isDriverGeneInThisModule, &isPhenotypeGeneInThisModule, network, D, &isGeneInThisModule, geneIdToSymbol);
	}

	//for each explained gene, check if it belong to at least one path (with length =< L) between a mutated gene and a phenotype gene
	vector<int> geneIdToBeRemoved;
	int i = 0;
	for (list<int>::iterator it = explainedGeneIdsList->begin(); it != explainedGeneIdsList->end(); it++, i++) {
		int currentExplainedGeneId = *it;
		int disToDriver = shortestPathToDriverGenes[i];
		int disToPhenotype = shortestPathToPhenotypeGenes[i];
		if(disToDriver < 0 or disToPhenotype < 0){
			geneIdToBeRemoved.push_back(currentExplainedGeneId);
		}else if(disToDriver + disToPhenotype > L){
			geneIdToBeRemoved.push_back(currentExplainedGeneId);
		}
	}

	//TODO 2 BUGS: 1) not trim some gene 2) disconnected module
	int numDelete = geneIdToBeRemoved.size();
	for(int i = 0; i < numDelete; i++){
		list<int>::iterator dit = find(explainedGeneIdsList->begin(), explainedGeneIdsList->end(), geneIdToBeRemoved[i]);
		explainedGeneIdsList->erase(dit);
		//cout << "deleted " << geneIdToSymbol->at(geneIdToBeRemoved[i]) << endl;
	}
}


void findShortestPath(int geneId, vector<int>* shortestPathToDeiverGenes, vector<int>* shortestPathToPhenotypeGenes,
		vector<bool>* isDriverGeneInThisModule, vector<bool>* isPhenotypeGeneInThisModule, TIntAdjList* network,
		int D, vector<bool>* isGeneInThisModule, vector<string>* geneIdToSymbol){
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

		//stop BFS after both driver and phenotype are found
		if(foundDriver and foundPhenotype){
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {

			int currentNeighborGeneId = adj[j];

			//if the current neighbor gene is in this module and is not visited
			if (isGeneInThisModule->at(currentNeighborGeneId) and !visited[currentNeighborGeneId]) {
				//check conditions before push to the queue

				//if a driver is not found yet and the current neighbor gene is a driver
				if(!foundDriver and isDriverGeneInThisModule->at(currentNeighborGeneId)){
					distantToDriver = levels[currentGeneId] + 1;
					foundDriver = true;
				}

				//if a phenotype is not found yet and the current neighbor gene is a phenotype
				if(!foundPhenotype and isPhenotypeGeneInThisModule->at(currentNeighborGeneId)){
					distantToPhenotype = levels[currentGeneId] + 1;
					foundPhenotype = true;
				}

				int degree = getNodeDegree(network, currentNeighborGeneId);
				//do not allow the search to go through the genes that is a driver, is a phenotype, has degree > D
				if (degree <= D and !isDriverGeneInThisModule->at(currentNeighborGeneId) and !isPhenotypeGeneInThisModule->at(currentNeighborGeneId)) {
					//push explained genes to queue
					q.push(currentNeighborGeneId);
					//save level
					levels[currentNeighborGeneId] = currentLevel + 1;
				}

				visited[currentNeighborGeneId] = true;
			}
		}
	}


	delete[] visited;

//	cout << "\t\tdistance from " << geneIdToSymbol->at(geneId) << " to a driver is " << distantToDriver << endl;
//	cout << "\t\tdistance from " << geneIdToSymbol->at(geneId) << " to a phenotype is " << distantToPhenotype << endl;

	shortestPathToDeiverGenes->push_back(distantToDriver);
	shortestPathToPhenotypeGenes->push_back(distantToPhenotype);

}

