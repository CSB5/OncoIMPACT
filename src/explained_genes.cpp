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
#include "explained_genes.h"

void getGeneExpressionFromSampleId(TDoubleMatrix* geneExpressionMatrix,
		vector<int>* genesEx, vector<double>* sampleExpression, int sampleId, vector<string>* geneIdToSymbol) {
	int totalGene = geneExpressionMatrix->size();
	for (int i = 0; i < totalGene; ++i) {
		sampleExpression->at(genesEx->at(i)) =
				geneExpressionMatrix->at(i)[sampleId];
	}
}

void getMutatedGeneIdsFromSampleId(Mutations* mutations,
		vector<int>* mutatedGeneIds, int sampleId, vector<int>* genesMut) {
	// sampleId is column id
	int numGenes = genesMut->size();
	for (int i = 0; i < numGenes; ++i) {
		//if a current gene is mutated
		if (mutations->matrix->at(i)[sampleId] > 0) {
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

	//TO BE DELETED
	//get unique explained genes and delete mutated genes
	//delete all duplicate
//	set<int> uniqueGenes;
//	int size = explainGeneIds.size();
//	for (int i = 0; i < size; ++i) {
//		uniqueGenes.insert(explainGeneIds[i]);
//	}
//	//delete mutated genes
//	int numMut = mutatedGeneIds->size();
//	for (int i = 0; i < numMut; i++) {
//		uniqueGenes.erase(mutatedGeneIds->at(i));
//	}
//	//unique explained genes
//	set<int>::iterator it;
//	for (it = uniqueGenes.begin(); it != uniqueGenes.end(); it++) {
//		ExplainedGene eg;
//		eg.id = *it;
//		eg.expression = sampleGeneExpression->at(*it);
//		explainedGenes->push_back(eg);
//		//cout << eg.id << "\t" << eg.expression << "\t"
//		//		<< getNodeDegree(network, eg.id) << endl;
//	}
}

void getExplainedGenesOnlyId(vector<int>* explainedGeneIds, TIntAdjList* network, vector<double>* sampleGeneExpression,
		vector<int>* mutatedGeneIds, int L, int D, double F){
	int numMutatedGenes = mutatedGeneIds->size();
	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
		//BFS
		//cout << "For mutated gene: " << mutatedGeneIds->at(i) << endl;
		BFSforExplainedGenesIdOnly(network, mutatedGeneIds->at(i), L, D, F,
				explainedGeneIds, sampleGeneExpression);
	}
}

void getMutatedAndExplainedGenes(vector<MutatedAndExplianedGenes>* mutatedAndExplainedGenes, TIntAdjList* network,
		vector<double>* sampleGeneExpression, vector<int>* mutatedGeneIds, int L, int D, double F){
	int numMutatedGenes = mutatedGeneIds->size();
	int totalGenes = network->size();
	for (int i = 0; i < numMutatedGenes; ++i) {	// for each mutated genes
		int mutatedGeneId = mutatedGeneIds->at(i);
		MutatedAndExplianedGenes* meg = &mutatedAndExplainedGenes->at(mutatedGeneId);
		meg->explainedGenesFreqency = vector<int>(totalGenes);
		//BFS for explained genes of the current mutated gene
		BFSforExplainedGenesIdOnly(network, mutatedGeneId, L, D, F,
				&(meg->explainedGenesFreqency), sampleGeneExpression);
	}

//	cout << "in getMutatedAndExplainedGenes\n";
//	for (int i = 0; i < numMutatedGenes; ++i) {
//		for (int j = 0; j < totalGenes; ++j) {
//			if(mutatedAndExplainedGenes->at(mutatedGeneIds->at(i)).explainedGenesFreqency[j] > 0)
//				cout << mutatedAndExplainedGenes->at(mutatedGeneIds->at(i)).explainedGenesFreqency[j] << endl;
//		}
//	}
}


void BFSforExplainedGenes(TIntAdjList* network, int geneId, int L, int D,
		double F, vector<ExplainedGene>* explainedGenes, vector<double>* sampleGeneExpression) {
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
				if (fabs(sampleGeneExpression->at(geneId)) > F) {
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
				if (fabs(sampleGeneExpression->at(geneId)) > F) {
					// is explained gene
					explainedGenesFrequency->at(geneId) += 1;
//					cout << "gene " << geneId << " is explained" << endl;
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

//TO BE DELETED
//void getGeneFrequencyOfSamples(vector<int>* genesCount, int totalGenes,
//		vector<vector<ExplainedGene> >* explainedGenesList) {
//	int numSamples = explainedGenesList->size();
//
//	//consider up and down separately
//	vector<int> up(totalGenes);
//	vector<int> down(totalGenes);
//
//	for (int i = 0; i < numSamples; ++i) {
//		vector<ExplainedGene> explainedGenes = explainedGenesList->at(i);
//		int numGens = explainedGenes.size();
//		for (int j = 0; j < numGens; ++j) {
//			//consider up and down separately
//			ExplainedGene eg = explainedGenes[j];
//			if (eg.expression < 0) {
//				down[eg.id]++;
//			} else {
//				up[eg.id]++;
//			}
//		}
//	}
//
//	//get final frequency discarding gene id
//	for (int i = 0; i < totalGenes; ++i) {
//		if (up[i] > 0) {
//			genesCount->push_back(up[i]);
//		}
//		if (down[i] > 0) {
//			genesCount->push_back(down[i]);
//		}
//	}
//}
