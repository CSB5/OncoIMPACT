/*
 * driver_genes.h
 *
 *  Created on: 28 Mar, 2015
 *      Author: Nok
 */

#ifndef DRIVER_GENES_H_
#define DRIVER_GENES_H_

#include "utilities.h"
#include "explained_genes.h"
#include <list>

struct BipartitePhenotypeNode{
	int phenotypeGeneIdUpDown;
	int sampleId;
};

struct BipartiteEdge{
	//int mutatedGeneId;	//is known from the index
	list<BipartitePhenotypeNode> phenotypeGeneIdsAndSampleIds;
};

void getAllMutatedGenes(vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isMutatedGenes, list<int>* mutatedGeneIdsList);

void createBipartiteGraph(vector<vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal,
		vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenes,
		vector<BipartiteEdge>* bipartiteEdges, vector<string>* geneIdToSymbol);

void findDriverGenes(vector<BipartiteEdge>* bipartiteEdges, list<int>* mutatedGeneIdsList, vector<int>* driverGeneIds);
//Note mutatedGeneIdsListReal is a sample level mutated genes list, while mutatedGeneIdsList is a list of all mutated genes

#endif /* DRIVER_GENES_H_ */
