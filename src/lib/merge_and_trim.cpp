/*
 * merge_and_trim.cpp
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#include <iostream>
#include <algorithm>
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
		//use vector<vector<MutatedAndExplianedGenes> > mutatedAndExplainedGenesListReal
		//use vector<vector<int> > mutatedGeneIdsListReal
		//use vector<bool> isPhenotypeGenes and isDriverGenes
		//use vector<int> phenotypeGeneIds

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
			vector<int> explainedGenesFreqency = mutatedAndExplainedGenes[currentMutatedGeneId].explainedGenesFreqency;

			modules[j].moduleId = j;
			modules[j].driverGeneIds.push_back(currentMutatedGeneId);
			for (int k = 0; k < totalGenes; ++k) {
				if(explainedGenesFreqency[k] > 0){
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

}

