/*
 * input.h
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#ifndef INPUT_H_
#define INPUT_H_

#include "utilities.h"

/*
 * Structure for importing data
 */

struct GeneExpression{
	TDoubleMatrix* matrix;
	vector<int>* genes;
	vector<int>* sample; //optional for saving sample id
};

struct Mutations{
	TIntegerMatrix* matrix; //0 if no mutation 1 if there are point mutations or CNVs
	vector<int>* genes;
	vector<int>* sample; //optional for saving sample id
};

struct PointMutations{
	TIntegerMatrix* matrix;
	vector<int>* genes;
	vector<int>* sample; //optional for saving sample id
};

struct CopyNumberVariation{
	TIntegerMatrix* matrix;
	vector<int>* genes;
	vector<int>* sample; //optional for saving sample id
};

/*
 * Functions for importing data
 */

void readGeneExpression(const char* filename, GeneExpression* geneExpression,
		char delim, map<string, int>* geneSymbolToId, vector<string>* sampleIdToName);
void readMutations(const char* filename, Mutations* pointMutations,
		char delim, map<string, int>* geneSymbolToId);
void readPointMutations(const char* filename, PointMutations* pointMutations,
		char delim, map<string, int>* geneSymbolToId);
void readCopyNumberVariation(const char* filename, CopyNumberVariation* copyNumberVariation,
		char delim, map<string, int>* geneSymbolToId);
void readGenesList(const char* filename, vector<int>* geneIds, map<string, int>* geneSymbolToId);
void readGenesListUpDown(const char* filename, vector<int>* geneIdsUpDown, map<string, int>* geneSymbolToId);

int findIndex(vector<int>* geneIds, int currentGeneId);

/*
 * Functions for calculate data statistics
 */


#endif /* INPUT_H_ */
