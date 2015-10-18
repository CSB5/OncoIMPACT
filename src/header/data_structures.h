/*
 * data_structures.cpp
 *
 *  Created on: May 19, 2015
 *      Author: nok
 */

#ifndef DATA_STRUCTURES_CPP_
#define DATA_STRUCTURES_CPP_

#include<vector>
#include<list>

using namespace std;

struct MutatedAndExplianedGenes{
	//int mutatedGeneId; this is equal to the index of the vector
	vector<bool>* isExplainedGenesUpDown;
};

struct BipartitePhenotypeNode{
	int phenotypeGeneIdUpDown;
	int sampleId;
};

struct BipartiteEdge{
	//int mutatedGeneId;	//is known from the index
	list<BipartitePhenotypeNode> phenotypeGeneIdsAndSampleIds;
};

struct DriverGene{
	int geneId;
	vector<int> sampleIds;
};


#endif /* DATA_STRUCTURES_CPP_ */
