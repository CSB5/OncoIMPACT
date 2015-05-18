//============================================================================
// Name        : utilities.h
// Author      : Nok C Suphavilai
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include<string>
#include<vector>
#include<map>

using namespace std;

/*
 * For matrix
 */
// data types for matrix
typedef vector< vector<double> > TDoubleMatrix;
typedef vector< vector<int> > TIntegerMatrix;
typedef vector< vector<bool> > TBooleanMatrix;
// methods for matrix
void traceMatrixDouble(TDoubleMatrix* matrix);
void traceMatrixInt(TIntegerMatrix* matrix);
void readDoubleMatrix(TDoubleMatrix* matrix, const char* filename, char delim);
void readIntegerMatrix(TIntegerMatrix* matrix, const char* filename, char delim);

/*
 * For graph
 */
// data types for graph
typedef vector< vector<int> > TIntAdjList;
typedef vector< vector<string> > TStrEdge;
// methods for graph
int readNetwork(const char* filename, TIntAdjList* network,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId, char delim);
void printAdjacencyList(TIntAdjList* network);
void DFS(TIntAdjList* network, int geneId);
void DFSUtil(TIntAdjList* network, int geneId, bool visited[]);
void BFS(TIntAdjList* network, int geneId);
void BFSUtil(TIntAdjList* network, int geneId, bool visited[]);
int getNodeDegree(TIntAdjList* network, int nodeId);

/*
 * Output
 */
void writeStrVector(const char* filename, vector<string>* output);
void saveGeneSymbols(const char* filename, vector<int>* geneIds, vector<string>* geneIdToSymbol);
void printGeneSymbols(vector<int>* geneIds, vector<string>* geneIdToSymbol);
void writeToLogFile(ofstream* outLogStream, string outStr);

/*
 * Permutation
 */

// input: a list of int containing n zeroes
// output: a list of permuted rank containing [0. n)
void createPermutation(vector<int>* rank);
// input: originalGeneLabels (for example, a list of genes in gene expression / mutation matrix)
// output: a permuted gene labels
void permuteGeneLabels(vector<int>* originalGeneLabels, vector<int>* permutedGeneLabels);
void permutedGeneLabelsUsingAllGeneInNetwork(vector<int>* originalGeneLabels,
		vector<int>* permutedGeneLabels, int totalGenes);

/*
 * ID mapping
 */
void mapGeneSymbolToId(vector<string>* symbols, vector<int>* ids, map<string,int>* geneSymbolToId);
void mapGeneIdToSymbol(vector<int>* ids, vector<string>* symbols, vector<string>* geneIdToSymbol);

/*
 * Miscellaneous
 */

string intToStr(int i);
string doubleToStr(double i, int prec);
bool trimStr(std::string& str, const std::string& from);

string getCurrentDate();
string getCurrentDateAndTime();
string getCurrentTimestamp();

#endif /* UTILITIES_H_ */
