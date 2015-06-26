/*
 * explained_genes.cpp
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include "../header/explained_genes.h"

void getGeneExpressionFromSampleId(TDoubleMatrix* geneExpressionMatrix,
		vector<int>* genesEx, vector<double>* sampleExpression, int sampleId) {

	//initialization
	int totalGenes = sampleExpression->size();
	for (int i = 0; i < totalGenes; ++i) {
		sampleExpression->at(i) = 0;
	}
	
	int totalExGenes = genesEx->size();

	for (int i = 0; i < totalExGenes; ++i) {

		int currentGeneId = genesEx->at(i);

		sampleExpression->at(currentGeneId) =
				geneExpressionMatrix->at(i)[sampleId];
	}
}

void getMutatedGeneIdsFromSampleId(Mutations* mutations,
		vector<int>* mutatedGeneIds, int sampleId, vector<int>* genesMut) {
	// sampleId is column id
	int numGenes = genesMut->size();
	for (int i = 0; i < numGenes; ++i) {
		//if a current gene is mutated
		if (mutations->matrix->at(i)[sampleId] != 0) {
			mutatedGeneIds->push_back(genesMut->at(i));
		}
	}

}

void getExplainedGenes(vector<ExplainedGene>* explainedGenes,
		TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, double F) {

	int numMutatedGenes = mutatedGeneIds->size();
	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
		//BFS
		//cout << "For mutated gene: " << mutatedGeneIds->at(i) << endl;
		BFSforExplainedGenes(network, mutatedGeneIds->at(i), L, D, F,
				explainedGenes, sampleGeneExpression);
	}
}

void getExplainedGenesIdOnly(vector<int>* explainedGenesFrequency, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, double F){
	int numMutatedGenes = mutatedGeneIds->size();
	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
		//BFS
		//cout << "For mutated gene: " << mutatedGeneIds->at(i) << endl;
		BFSforExplainedGenesIdOnly(network, mutatedGeneIds->at(i), L, D, F,
				explainedGenesFrequency, sampleGeneExpression);
	}
}

//output: explained genes frequency
//up regulated gene frequency is in the first half [0,totalGenes - 1]
//down regulated gene frequency is in the second half [totalGenes, totalGenes * 2 - 1]
void getExplainedGenesIdOnlyUpDown(vector<bool>* isExplainedGeneUpDown, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, double F){
	int numMutatedGenes = mutatedGeneIds->size();
	int totalGenesUpDown = isExplainedGeneUpDown->size();
	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
		//BFS
//		cout << "For mutated gene: " << mutatedGeneIds->at(i) << endl;
		vector<bool> isExplaninedGeneUpDownForAMutatedGene(totalGenesUpDown);
		BFSforExplainedGenesIdOnlyUpDown(network, mutatedGeneIds->at(i), L, D, F,
				&isExplaninedGeneUpDownForAMutatedGene, sampleGeneExpression);

		for (int j = 0; j < totalGenesUpDown; ++j) {
			if(isExplaninedGeneUpDownForAMutatedGene[j]){
				isExplainedGeneUpDown->at(j) = true;
			}
		}


	}
}

void getExplainedGenesIdOnlyUpDownIncludingMutatedGene(vector<bool>* isExplainedGenesUpDown, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, double F, vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId){
	int numMutatedGenes = mutatedGeneIds->size();
	int totalGenesUpDown = isExplainedGenesUpDown->size();

	vector<bool> isExplaninedGeneUpDownForAMutatedGene(totalGenesUpDown, false);

	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
//		cout << "For mutated gene: " << mutatedGeneIds->at(i) << endl;

		//BFS
//		void BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(TIntAdjList* network, int mutatedGeneId, int L, int D,
//				double F, vector<bool>* isExplainedGenes, vector<double>* sampleGeneExpression, int currentSampleId, map<string, int>* geneSymbolToId);
		BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(network, mutatedGeneIds->at(i), L, D, F,
				&isExplaninedGeneUpDownForAMutatedGene, sampleGeneExpression, -1, geneIdToSymbol, geneSymbolToId);

		for (int j = 0; j < totalGenesUpDown; ++j) {
			if(isExplaninedGeneUpDownForAMutatedGene[j]){
				isExplainedGenesUpDown->at(j) = true;
			}
		}
	}
}

bool isExplainedInASample(int currentExplainedGeneId, TIntAdjList* network, vector<double>* sampleGeneExpression, int currentMutatedGeneId,
		int L, int D, double F, map<string, int>* geneSymbolToId){

	int totalGenes = network->size();

	// Mark all the vertices as not visited
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes); // store level of each nodes
	// Mark the current node as visited
	q.push(currentMutatedGeneId);
	visited[currentMutatedGeneId] = true;
	int currentGeneId;
	int currentLevel = 0;

	//consider all nodes that are in <= L distant (Note that the level value start from 0, so use < L)
	while (!q.empty() && currentLevel < L) {
		//read the root node
		currentGeneId = q.front();
		q.pop();
		currentLevel = levels[currentGeneId];

		if (currentLevel >= L) {
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGeneId][j]]) {

				int geneId = (*network)[currentGeneId][j];

				// check conditions before push to the queue
				// |fold change| >= F
//				cout << "gene id " << geneId << " is checking\n";

				if (fabs(sampleGeneExpression->at(geneId)) >= F) {

//					cout << "gene id " << geneId << " is explained\n";

					// is explained gene
					if(sampleGeneExpression->at(geneId) > 0.0){ 	// up regulated
						if(geneId == currentExplainedGeneId){
							return true;
						}
//						isExplainedGenesUpDown->at(geneId) = true;
//						countExpalinedGenes++;
					}else{											// down regulated
						if(geneId + totalGenes == currentExplainedGeneId){
							return true;
						}
//						isExplainedGenesUpDown->at(geneId + totalGenes) = true; //+9452
//						countExpalinedGenes++;
					}

					int degree = getNodeDegree(network, geneId);
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}

//				cout << "gene id " << geneId << " is ckecked\n";

				// mark visited
				visited[geneId] = true;

			}
		}
	}

	delete[] visited;

	return false;
}

//void getMutatedAndExplainedGenes(vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, TIntAdjList* network,
//		vector<double>* sampleGeneExpression, vector<int>* mutatedGeneIds, int L, int D, double F){
//	int numMutatedGenes = mutatedGeneIds->size();
//	int totalGenes = network->size();
//	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
//		int mutatedGeneId = mutatedGeneIds->at(i);
//		MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes->at(mutatedGeneId);
//		//BFS for explained genes of the current mutated gene
//		BFSforExplainedGenesIdOnly(network, mutatedGeneId, L, D, F,
//				meg->explainedGenesFreqency, sampleGeneExpression);
//	}
//
//}


void BFSforExplainedGenes(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<ExplainedGene>* explainedGenes, vector<double>* sampleGeneExpression) {
	// Mark all the vertices as not visited
	int totalGenes = network->size();
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes); // store level of each nodes
	// Mark the current node as visited
	q.push(geneId);
	visited[geneId] = true;
	int currentGene;
	int currentLevel = 0;

	//consider all nodes that are in <= L distant
	while (!q.empty() && currentLevel <= L) {
		//read the root node
		currentGene = q.front();
		q.pop();
		currentLevel = levels[currentGene];

		if (currentLevel > L) {
			break;
		}

		//cout << currentLevel << " : " << currentGene << "(" <<
		//		sampleGeneExpression->at(currentGene) << ")(" <<
		//		getNodeDegree(network, currentGene) << ")" << endl;

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGene];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGene][j]]) {
				// check conditions before push to the queue
				// |fold change| > F
				int geneId = (*network)[currentGene][j];
				if (fabs(sampleGeneExpression->at(geneId)) >= F) {	//in pl code used >=, but in paper said >
					// is explained gene
					explainedGenes->at(geneId).expression = sampleGeneExpression->at(geneId);
					explainedGenes->at(geneId).frequency++;
					int degree = getNodeDegree(network, geneId);
					explainedGenes->at(geneId).degree = degree;
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}
				// mark visited
				visited[geneId] = true;
			}
		}
	}

	delete[] visited;
}

void BFSforExplainedGenesIdOnly(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<int>* explainedGenesFrequency, vector<double>* sampleGeneExpression) {
	// Mark all the vertices as not visited
	int numNode = network->size();
	bool *visited = new bool[numNode];
	for (int j = 0; j < numNode; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(numNode); // store level of each nodes
	// Mark the current node as visited
	q.push(geneId);
	visited[geneId] = true;
	int currentGene;
	int currentLevel = 0;

	//consider all nodes that are in <= L distant
	while (!q.empty() && currentLevel <= L) {
		//read the root node
		currentGene = q.front();
		q.pop();
		currentLevel = levels[currentGene];

		if (currentLevel > L) {
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGene];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGene][j]]) {
				// check conditions before push to the queue
				// |fold change| > F
				int geneId = (*network)[currentGene][j];
				if (fabs(sampleGeneExpression->at(geneId)) >= F) {
					// is explained gene
					explainedGenesFrequency->at(geneId) += 1;
					int degree = getNodeDegree(network, geneId);
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}
				// mark visited
				visited[geneId] = true;
			}
		}
	}


	delete[] visited;
}

//up regulated gene frequency is in the first half [0,totalGenes - 1]
//down regulated gene frequency is in the second half [totalGenes, totalGenes * 2 - 1]
void BFSforExplainedGenesIdOnlyUpDown(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<bool>* isExplainedGenes, vector<double>* sampleGeneExpression) {

	//initialize explainedGenesFrequency
	int totalGenesUpDown = isExplainedGenes->size();
	for (int i = 0; i < totalGenesUpDown; ++i) {
		isExplainedGenes->at(i) = false;
	}

	int totalGenes = network->size();

	// Mark all the vertices as not visited
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes); // store level of each nodes
	// Mark the current node as visited
	q.push(geneId);
	visited[geneId] = true;
	int currentGeneId;
	int currentLevel = 0;

	//consider all nodes that are in <= L distant (Note that the level value start from 0, so use < L)
	while (!q.empty() && currentLevel < L) {
		//read the root node
		currentGeneId = q.front();
		q.pop();
		currentLevel = levels[currentGeneId];

		if (currentLevel >= L) {
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGeneId][j]]) {
				// check conditions before push to the queue
				// |fold change| > F
				int geneId = (*network)[currentGeneId][j];
				if (fabs(sampleGeneExpression->at(geneId)) >= F) {

					// is explained gene
					if(sampleGeneExpression->at(geneId) > 0.0){ 	// up regulated
						isExplainedGenes->at(geneId) = true;
					}else{											// down regulated
						isExplainedGenes->at(geneId + totalGenes) = true; //+9452
					}

					int degree = getNodeDegree(network, geneId);
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}
				// mark visited
				visited[geneId] = true;

			}
		}
	}


	delete[] visited;

}

//up regulated gene frequency is in the first half [0,totalGenes - 1]
//down regulated gene frequency is in the second half [totalGenes, totalGenes * 2 - 1]
//This also consider the deregulated mutated genes as an explained gene
void BFSforExplainedGenesIdOnlyUpDownIncludingMutatedGene(TIntAdjList* network, int mutatedGeneId, int L, int D,
		double F, vector<bool>* isExplainedGenesUpDown, vector<double>* sampleGeneExpression, int currentSampleId,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId) {

	int countExpalinedGenes = 0;

	//initialize explainedGenesFrequency
	int totalGenesUpDown = isExplainedGenesUpDown->size();
	for (int i = 0; i < totalGenesUpDown; ++i) {
		isExplainedGenesUpDown->at(i) = false;
	}

	int totalGenes = network->size();

	// Mark all the vertices as not visited
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes); // store level of each nodes
	// Mark the current node as visited
	q.push(mutatedGeneId);
	visited[mutatedGeneId] = true;
	int currentGeneId;
	int currentLevel = 0;

	//consider all nodes that are in <= L distant (Note that the level value start from 0, so use < L)
	while (!q.empty() && currentLevel < L) {
		//read the root node
		currentGeneId = q.front();
		q.pop();
		currentLevel = levels[currentGeneId];

		if (currentLevel >= L) {
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGeneId][j]]) {
				// check conditions before push to the queue
				// |fold change| > F
				int geneId = (*network)[currentGeneId][j];

//				cout << "gene id " << geneId << " is checking\n";

				if (fabs(sampleGeneExpression->at(geneId)) >= F) {

//					cout << "gene id " << geneId << " is explained\n";

					// is explained gene
					if(sampleGeneExpression->at(geneId) > 0.0){ 	// up regulated
						isExplainedGenesUpDown->at(geneId) = true;
						countExpalinedGenes++;
//						cout << geneIdToSymbol->at(geneId) << " is at level " << currentLevel << endl;
					}else{											// down regulated
						isExplainedGenesUpDown->at(geneId + totalGenes) = true; //+9452
						countExpalinedGenes++;
//						cout << geneIdToSymbol->at(geneId) << " is at level " << currentLevel << endl;
					}

					int degree = getNodeDegree(network, geneId);
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}

//				cout << "gene id " << geneId << " is ckecked\n";

				// mark visited
				visited[geneId] = true;

			}
		}
	}

	//because the permuted mutated genes label contain all the gene in the network,
	//it is possible that the current mutated gene is not in the gene expression matrix.. not a problem
//	cout << "before checking if the mutated gene is also deregulated\n";

	//add the mutated gene itself to the list of mutated genes
	if (countExpalinedGenes > 0 and fabs(sampleGeneExpression->at(mutatedGeneId)) >= F) {

		if(sampleGeneExpression->at(mutatedGeneId) > 0.0){ 	// up regulated
			isExplainedGenesUpDown->at(mutatedGeneId) = true;
		}else{											// down regulated
			isExplainedGenesUpDown->at(mutatedGeneId + totalGenes) = true; //+9452
		}
	}
//	cout << "after checking if the mutated gene is also deregulated\n";

	delete[] visited;
}

//use only one BFS for all L's values
void BFSforExplainedGenesUpDownIncludingMutatedGeneAllLengths(TIntAdjList* network, int mutatedGeneId, vector<int>* Ls, int D,
		double F, vector< vector<bool> >* isExplainedGenesUpDownAllLs, vector<double>* sampleGeneExpression, int currentSampleId,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId) {

//	cout << "for mutated gene " << mutatedGeneId << endl;

	int Lmax = Ls->at(Ls->size()-1);
	int numLs = Ls->size();
	vector<int> countExpalinedGenesForAllLs(numLs, 0);

	//isExplainedGenesUpDownAllL is initialized before sending to function
	//matrix: row = L value, col = gene id up down

	int totalGenes = network->size();

	// Mark all the vertices as not visited
	bool *visited = new bool[totalGenes];
	for (int j = 0; j < totalGenes; j++) {
		visited[j] = false;
	}

	// Create queue
	queue<int> q;
	vector<int> levels(totalGenes, 0); // store level of each nodes
	// Mark the current node as visited
	q.push(mutatedGeneId);
	visited[mutatedGeneId] = true;
	int currentGeneId;
	int currentLevel = 0;

	//consider all nodes that are in <= Lmax distant (Note that the level value start from 0, so use < Lmax)
	while (!q.empty() && currentLevel < Lmax) {
		//read the root node
		currentGeneId = q.front();
		q.pop();
		currentLevel = levels[currentGeneId];

		if (currentLevel >= Lmax) {
			break;
		}

		//explore all the connected nodes
		vector<int> adj = (*network)[currentGeneId];

		int numAdj = adj.size();
		for (int j = 0; j < numAdj; j++) {
			if (!visited[(*network)[currentGeneId][j]]) {

				int geneId = (*network)[currentGeneId][j];

				// check conditions before push to the queue
				// |fold change| >= F
				if (fabs(sampleGeneExpression->at(geneId)) >= F) {

					// is explained gene
					if(sampleGeneExpression->at(geneId) > 0.0){ 	// up regulated

						//for all L <= current level
						for (int li = 0; li < numLs; ++li) {
							if(currentLevel < Ls->at(li)){	//(use < L to get <= L)
								isExplainedGenesUpDownAllLs->at(li)[geneId] = true;
								countExpalinedGenesForAllLs[li]++;
							}
						}
					}else{											// down regulated
						//for all L <= current level
						for (int li = 0; li < numLs; ++li) {
							if(currentLevel < Ls->at(li)){	//(use < L to get <= L)
								isExplainedGenesUpDownAllLs->at(li)[geneId + totalGenes] = true;
								countExpalinedGenesForAllLs[li]++;
							}
						}
					}

					int degree = getNodeDegree(network, geneId);
					// check the degree
					if (degree <= D) {
						// push to queue
						q.push(geneId);
						// save level
						levels[geneId] = currentLevel + 1;
					}
				}

				// mark visited
				visited[geneId] = true;

			}
		}
	}

	//add the mutated gene itself to the list of mutated genes
	for (int li = 0; li < numLs; ++li) {
		if (countExpalinedGenesForAllLs[li] > 0 and fabs(sampleGeneExpression->at(mutatedGeneId)) >= F) {
			if(sampleGeneExpression->at(mutatedGeneId) > 0.0){ 	// up regulated
				isExplainedGenesUpDownAllLs->at(li)[mutatedGeneId] = true;
			}else{
				isExplainedGenesUpDownAllLs->at(li)[mutatedGeneId + totalGenes] = true;
			}
		}
	}

	delete[] visited;

}
