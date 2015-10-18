/*
 * input.h
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#ifndef INPUT_H_
#define INPUT_H_

#include "utilities.h"
#include "data_structures.h"
#include <set>

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
 * For discovery mode
 */

//GENE	DRIVER_FREQUENCY	DRIVER_SNV_FREQUENCY	DRIVER_DELTION_FREQUENCY	DRIVER_AMPLIFICATION_FREQUENCY	CANCER_CENSUS	PAN_CANCER	IMPACT	MUTATION_FREQUENCY	SNV_FREQUENCY	DELTION_FREQUENCY	AMPLIFICATION_FREQUENCY
//EGFR	0.393	0.262	0.000	0.241	Y	NA	120.1447552328	0.393	0.262	0.000	0.241
struct DriverGeneFromFile{
	//mutated gene id is the index of the vector<DriverGeneFromFile>
	double driverFreq;
	double driverSnpFreq;
	double driverDelFreq;
	double driverAmpFreq;
	bool isInCancerCensus;
	double impactScore;
	double mutFreq;
	double snpFreq;
	double delFreq;
	double ampFreq;
};

struct MutatedGeneFromFile{
	//mutated gene id is the index of the vector<MutatedGeneFromFile>
	double mutationFreq = 0.0;
	double snpFreq = 0.0;
	double delFreq = 0.0;
	double ampFreq = 0.0;
	double snpAndCnvFreq = 0.0;
	double driverFreq = 0.0;
	double driverSnpFreq = 0.0;
	double driverDelFreq = 0.0;
	double driverAmpFreq = 0.0;
	double driverSnpAndCnvFreq = 0.0;
	double impactScore = 0.0;
};

void readPhenotypeGenesFromFile(const char* filename, vector<int>* phenotypeGeneIdsUpDown, map<string, int>* geneSymbolToId);
void readDriverGenesFromFile(const char* filename, vector<DriverGeneFromFile>* driverGenesFromFile, map<string, int>* geneSymbolToId);
void readMutatedGenesFromFile(const char* filename, vector<MutatedGeneFromFile>* mutatedGenesFromFile, map<string, int>* geneSymbolToId);

void readDrugsListFromFile(const char* filename, map<int, string>* drugIdToName, map<string, int>* drugNameToId);
void readGenesDrugsAssocFromFile(const char* filename, vector< vector<int> >* geneDrugsAssocList, map<string, int>* geneSymbolToId);

void readModulesFromFile(string* moduleFileName, vector<string>* sampleIdToName, map<string, int>* sampleNameToId,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId,
		vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< vector<int> >* mutatedGeneIdsListReal);
void readModulesFromFileOnlySaveBipartiteGraph(string* moduleFileName, vector<string>* sampleIdToName, map<string, int>* sampleNameToId,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId,
		vector<BipartiteEdge>* bipartiteGraphDatabase, vector< vector<int> >* mutatedGeneIdsListReal, vector<bool>* isPhenotypeGenesUpDown);

int readCancerTypeList(const char* filename, map<string, string>* cancerTypeDataTypeMap);
bool replaceStr(std::string& str, const std::string& from, const std::string& to);


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
void readBenchmarkGeneList(string benchmarkGeneListFilename, vector<int>* cancerBenchmarkGenes, map<string, int>* geneSymbolToId, vector<string>* cancerBenchmarkGeneNames);

int findIndex(vector<int>* geneIds, int currentGeneId);


/*
 * Functions for calculate data statistics
 */


#endif /* INPUT_H_ */
