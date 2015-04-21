/*
 * driver_genes.cpp
 *
 *  Created on: 28 Mar, 2015
 *      Author: Nok
 */

#include <algorithm>
#include <set>
#include <iostream>
#include "../header/driver_genes.h"

void getAllMutatedGenes(vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isMutatedGenes, list<int>* mutatedGeneIdsList){
	int totalSamples = mutatedGeneIdsListReal->size();
	set<int> uniqueMutatedGeneIdsOfAllSamples;
	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {

		//get all mutated genes and their corresponding explained genes
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal->at(i);
		int numMutatedGenes = mutatedGeneIds.size();

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			isMutatedGenes->at(mutatedGeneIds[j]) = true;
			uniqueMutatedGeneIdsOfAllSamples.insert(mutatedGeneIds[j]);
		}
	}
	mutatedGeneIdsList->assign(uniqueMutatedGeneIdsOfAllSamples.begin(), uniqueMutatedGeneIdsOfAllSamples.end());
}

void createBipartiteGraph(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenes,
		vector<BipartiteEdge>* bipartiteEdges){
	int totalSamples = mutatedAndExplainedGenesListReal->size();
	int totalGenes = isPhenotypeGenes->size();

	//for each sample i
	for (int i = 0; i < totalSamples; ++i) {

		//get all mutated genes and their corresponding explained genes
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes =
				mutatedAndExplainedGenesListReal->at(i);
		vector<int> mutatedGeneIds = mutatedGeneIdsListReal->at(i);
		int numMutatedGenes = mutatedGeneIds.size();

		//for each mutated gene j
		for (int j = 0; j < numMutatedGenes; ++j) {
			int currentMutatedGeneId = mutatedGeneIds[j];
			vector<bool>* isExplainedGenesUpDown =
					mutatedAndExplainedGenes[currentMutatedGeneId].isExplainedGenesUpDown;
			int totalGenesUpDown = isExplainedGenesUpDown->size();

			//convert to non UpDown
			vector<bool> isExplianedGenes(totalGenes);
			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(isExplainedGenesUpDown->at(k)){
					if(k < totalGenes){
						isExplianedGenes[k] = true;
					}else{
						isExplianedGenes[k-totalGenes] = true;
					}
				}

			}


			for (int k = 0; k < totalGenes; ++k) {
				if(isExplianedGenes[k] > 0){ // gene k is explained in sample i
					if(isPhenotypeGenes->at(k)){		// gene k is a phenotype gene
						BipartitePhenotypeNode node;
						node.phenotypeGeneId = k;
						node.sampleId = i;
						bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.push_back(node);
					}
				}
			}
		}
	}
}

void findDriverGenes(vector<BipartiteEdge>* bipartiteEdges, list<int>* mutatedGeneIdsList, vector<int>* driverGeneIds){

	bool coveredAll = false;
	//until all the phenotype genes are covered
	while(!coveredAll){

		//find a mutated genes that cover the maximum number of genes (count = 2 if the same genes is phenotype in 2 samples)
		int max = 0;
		int maxGeneId = -1;
		list<int>::iterator maxRef;
		list<BipartitePhenotypeNode> coveredPhenotypeGens;

		list<int>::iterator it = mutatedGeneIdsList->begin();
		while (it != mutatedGeneIdsList->end()){
			int currentMutatedGeneId = *it;
			int numEdges = bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.size();
			//TODO if numEdges are equal, choose the one that cover more patients
			if(numEdges > max){
				max = numEdges;
				maxGeneId = currentMutatedGeneId;
				maxRef = it;
			}
			if(numEdges > 0){
				it++;
			}else{	//deleted the mutated genes that do not connect to any phenotype gene
				it = mutatedGeneIdsList->erase(it); //erase the current element and return the pointer to next element
			}
		}

//		cout << "gene " << maxGeneId << " have the max # edges = " << max << endl;
		coveredPhenotypeGens = bipartiteEdges->at(maxGeneId).phenotypeGeneIdsAndSampleIds;	// a list of covered genes to be deleted

		driverGeneIds->push_back(maxGeneId);

		//deleted mutated gene from list
		mutatedGeneIdsList->erase(maxRef);

		//deleted all edges connecting to the covered phenotype genes
		//for each mutated genes
		for (list<int>::iterator mutIt = mutatedGeneIdsList->begin(); mutIt != mutatedGeneIdsList->end(); mutIt++){
			int currentMutatedGeneId = *mutIt;
			list<BipartitePhenotypeNode>* currentNodes = &bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds;

			//for each covered genes
			for (list<BipartitePhenotypeNode>::iterator coverIt = coveredPhenotypeGens.begin(); coverIt != coveredPhenotypeGens.end(); coverIt++){
				int currentPhenotypeGeneId = coverIt->phenotypeGeneId;
				int currentSampleId = coverIt->sampleId;

				list<BipartitePhenotypeNode>::iterator nodeIt = currentNodes->begin();
				while(nodeIt != currentNodes->end()){
					if(nodeIt->sampleId == currentSampleId and nodeIt->phenotypeGeneId == currentPhenotypeGeneId){	//already covered, then delete
						nodeIt = currentNodes->erase(nodeIt); //erase the current element and return the pointer to next element
					}else{
						nodeIt++;
					}
				}
			}

		}

		//check if all phenotypes gene in all samples have been covered (set coveredAll = true)
		//also delete mutated genes from list if their phenotype are alredy covered
		int numNotCovered = 0;
		it = mutatedGeneIdsList->begin();
		while (it != mutatedGeneIdsList->end()){
			int currentMutatedGeneId = *it;
			int numEdges = bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.size();
			if(numEdges > 0){
				numNotCovered += numEdges;
				it++;
			}else{	//deleted the mutated genes that do not connect to any phenotype gene
				it = mutatedGeneIdsList->erase(it); //erase the current element and return the pointer to next element
			}
		}
		if(numNotCovered == 0){
			coveredAll = true; //done
		}

	}
}

