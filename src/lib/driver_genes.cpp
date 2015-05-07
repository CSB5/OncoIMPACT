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
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenesUpDown,
		vector<BipartiteEdge>* bipartiteEdges, vector<string>* geneIdToSymbol){
	int totalSamples = mutatedAndExplainedGenesListReal->size();
	int totalGenesUpDown = isPhenotypeGenesUpDown->size();
	int totalGenes = totalGenesUpDown / 2;


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

			//for each explained gene k
			for (int k = 0; k < totalGenesUpDown; ++k) {
				if(isExplainedGenesUpDown->at(k)){ 			// gene k is explained in sample i

					if(isPhenotypeGenesUpDown->at(k)){
						BipartitePhenotypeNode node;
						node.phenotypeGeneIdUpDown = k;		//this is fine because each sample has either up or down
						node.sampleId = i;
						bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.push_back(node);
					}
//					if(k < totalGenes){	//upregulated
//						if(isPhenotypeGenes->at(k)){		// gene k is a phenotype gene (also consider up and down)
//							BipartitePhenotypeNode node;
//							node.phenotypeGeneIdUpDown = k;		//this is find because each sample has either up or down (not both)
//							node.sampleId = i;
//							bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.push_back(node);
//						}
//					}else{				//downregulated
//						if(isPhenotypeGenes->at(k-totalGenes)){		// gene k-totalGenes is a phenotype gene
//							BipartitePhenotypeNode node;
//							node.phenotypeGeneIdUpDown = k;
//							node.sampleId = i;
//							bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds.push_back(node);
//						}
//					}
				}
			}
		}
	}
}

void findDriverGenes(vector<BipartiteEdge>* bipartiteEdges, list<int>* mutatedGeneIdsList, vector<DriverGene>* driverGenes){

	bool coveredAll = false;
	//until all the phenotype genes are covered
	while(!coveredAll){

		//find a mutated genes that cover the maximum number of genes (e.g. count = 2 if the same genes is phenotype in 2 samples)
		int maxCovered = 0;
		vector<int> maxGeneIds;

		vector<list<int>::iterator> maxRefs;	//pointer to the mutated gene that cover maximum number of phenotype genes
		vector<int> coveredSampleCount;
		vector< vector<int> > sampleIdsOfMaxs;

		list<int>::iterator it = mutatedGeneIdsList->begin();

		//loop for finding the mutated gene that cover maximum number of phenotype genes
		while (it != mutatedGeneIdsList->end()){
			int currentMutatedGeneId = *it;
			list<BipartitePhenotypeNode> edgesCovered = bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds;
			int numEdgesCovered = edgesCovered.size();	//# (phenotype gene id, sample id)

			//count number of samples
			set<int> uniqueSampleIds;
			for (list<BipartitePhenotypeNode>::iterator eit = edgesCovered.begin(); eit != edgesCovered.end(); ++eit) {
				uniqueSampleIds.insert(eit->sampleId);
			}
			int numSampleCovered = uniqueSampleIds.size();
//			cout << "# samples in set = " << numSampleCovered << endl;
//			cout << "# samples in vector = " << samplesIds.size() << endl;

			if(numEdgesCovered >= maxCovered){
				maxCovered = numEdgesCovered;
				maxGeneIds.push_back(currentMutatedGeneId);
				maxRefs.push_back(it);

				coveredSampleCount.push_back(numSampleCovered);
				vector<int> samplesIds(uniqueSampleIds.begin(), uniqueSampleIds.end());
				sampleIdsOfMaxs.push_back(samplesIds);
			}

			if(numEdgesCovered > 0){
				it++;
			}else{	//deleted the mutated genes that do not connect to any phenotype gene
				it = mutatedGeneIdsList->erase(it); //erase the current element and return the pointer to next element
			}
		}


		int maxGeneId = -1;
		list<int>::iterator maxRef;
		vector<int> sampleIdsOfMax;

		//if there are more than one mutated genes that have the same number of covers (CHECK with the original code first)
		if(maxGeneIds.size() > 1){
			//choose the one that cover more samples first
			int maxSample = 0;
			int numGenes = maxGeneIds.size();
			for (int i = 0; i < numGenes; ++i) {
				if(coveredSampleCount[i] > maxSample){
					maxSample = coveredSampleCount[i];
					maxGeneId = maxGeneIds[i];
					maxRef = maxRefs[i];
					sampleIdsOfMax = sampleIdsOfMaxs[i];
				}
			}
		}else{	//there is only one mutated gene that has maximum covering
			maxGeneId = maxGeneIds[0];
			maxRef = maxRefs[0];
			sampleIdsOfMax = sampleIdsOfMaxs[0];
		}

		//cout << "Driver = " << maxGeneId << " covered " << max << " (phenotype gene id, sample id)" << endl;
		DriverGene driverMax;
		driverMax.geneId = maxGeneId;
		for (unsigned i = 0; i < sampleIdsOfMax.size(); ++i) {
			driverMax.sampleIds.push_back(sampleIdsOfMax[i]);
		}
		driverGenes->push_back(driverMax);

		//deleted mutated gene from list
		mutatedGeneIdsList->erase(maxRef);
		// a list of covered genes to be deleted
		list<BipartitePhenotypeNode> coveredPhenotypeGens;
		coveredPhenotypeGens = bipartiteEdges->at(maxGeneId).phenotypeGeneIdsAndSampleIds;

		//deleted all edges connecting to the covered (phenotype gene id, sample id)s
		//for each mutated genes
		for (list<int>::iterator mutIt = mutatedGeneIdsList->begin(); mutIt != mutatedGeneIdsList->end(); mutIt++){
			int currentMutatedGeneId = *mutIt;
			list<BipartitePhenotypeNode>* currentNodes = &bipartiteEdges->at(currentMutatedGeneId).phenotypeGeneIdsAndSampleIds;

			//for each covered genes
			for (list<BipartitePhenotypeNode>::iterator coveredIt = coveredPhenotypeGens.begin(); coveredIt != coveredPhenotypeGens.end(); coveredIt++){
				int currentPhenotypeGeneIdUpDown = coveredIt->phenotypeGeneIdUpDown;
				int currentSampleId = coveredIt->sampleId;

				list<BipartitePhenotypeNode>::iterator nodeIt = currentNodes->begin();
				while(nodeIt != currentNodes->end()){
					if(nodeIt->sampleId == currentSampleId and nodeIt->phenotypeGeneIdUpDown == currentPhenotypeGeneIdUpDown){	//already covered, then delete
						nodeIt = currentNodes->erase(nodeIt); //erase the current element and return the pointer to next element
					}else{	//go to check next (phenotype gene id, sample id)
						nodeIt++;
					}
				}
			}
		}

		//check if all phenotypes gene in all samples have been covered (set coveredAll = true)
		//also delete mutated genes from list if all of their phenotype genes are already covered
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

