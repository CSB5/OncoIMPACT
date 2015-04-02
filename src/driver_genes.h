/*
 * driver_genes.h
 *
 *  Created on: 28 Mar, 2015
 *      Author: Nok
 */

#ifndef DRIVER_GENES_H_
#define DRIVER_GENES_H_

#include "utilities.h"
#include <list>

struct BipartitePhenotypeNode{
	int phenotypeGeneId;
	int sampleId;
};

struct BipartiteEdge{
	//int mutatedGeneId;	//is known from the index
	list<BipartitePhenotypeNode> phenotypeGeneIdsAndSampleIds;
};



#endif /* DRIVER_GENES_H_ */
