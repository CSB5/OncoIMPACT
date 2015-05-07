//============================================================================
// Name        : utilities.cpp
// Author      : Nok C Suphavilai
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <queue>
#include <algorithm>
#include <cstdlib>
#include "../header/utilities.h"

/*
 * For matrix
 */
void traceMatrixDouble(TDoubleMatrix* matrix) {
	int rowNum = matrix->size();
	int colNum = matrix->at(0).size();

	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			cout << matrix->at(i)[j] << " ";
		}
		cout << endl;
	}
}

void traceMatrixInt(TIntegerMatrix* matrix) {
	int rowNum = matrix->size();
	int colNum = matrix->at(0).size();

	cout << rowNum << "\t" << colNum;

	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			cout << matrix->at(i)[j] << " ";
		}
		cout << endl;
	}
}

void readDoubleMatrix(TDoubleMatrix* matrix, const char* filename, char delim) {
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	if (inFile.is_open()) {
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<double> row;

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				row.push_back(atof(s.c_str()));
			}

			matrix->push_back(row);
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readIntegerMatrix(TDoubleMatrix* matrix, const char* filename,
		char delim) {
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	if (inFile.is_open()) {
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<double> row;

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				row.push_back(atoi(s.c_str()));
			}

			matrix->push_back(row);
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

/*
 * For graph
 */

void readNetwork(const char* filename, TIntAdjList* network,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId,
		char delim) {
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	//read data
	vector<string> genesTemp;
	TStrEdge edgesTemp;

	if (inFile.is_open()) {
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<string> edge;
			vector<string> edgeBackward;

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				genesTemp.push_back(s);
				edge.push_back(s);
			}

			//add the other edge (b-a) in addition to (a-b)
			edgeBackward.push_back(edge[1]);
			edgeBackward.push_back(edge[0]);

			edgesTemp.push_back(edge);
			edgesTemp.push_back(edgeBackward);
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}

	cout << "Creating gene mapping ...\n";

	//delete all duplicate
	set<string> uniqueGenes;
	unsigned size = genesTemp.size();
	for (unsigned i = 0; i < size; ++i) {
		uniqueGenes.insert(genesTemp[i]);
	}

	//map id to symbol
	geneIdToSymbol->assign(uniqueGenes.begin(), uniqueGenes.end());

	int numNode = geneIdToSymbol->size();

	//map symbol to id
	for (int i = 0; i < numNode; ++i) {
		geneSymbolToId->insert(pair<string, int>(geneIdToSymbol->at(i), i));
		//cout << geneIdToSymbol->at(i) << " => " << i << endl;
	}

	cout << "Mapped " << numNode << " genes\n";
	cout << "Creating adjacency list of gene interaction networks ...\n";

	int numEdge = edgesTemp.size();
	network->resize(numNode);

//	cout << "total number of edges = " << numEdge << endl;

	int source = -1;
	int target = -1;
	for (int i = 0; i < numEdge; ++i) {
		source = geneSymbolToId->find(edgesTemp[i][0])->second;
		target = geneSymbolToId->find(edgesTemp[i][1])->second;
		(*network)[source].push_back(target);
	}

	//delete the duplicated (target) nodes in the adjacency list
	//(in case that the network file contains some backward edge)
	for (int i = 0; i < numNode; ++i) {
		//for a current gene i
		int size = network->at(i).size();

		//delete all duplicate
		set<int> uniqueIds;
		for (int j = 0; j < size; ++j) {
			uniqueIds.insert(network->at(i)[j]);
		}
		network->at(i).clear();
		network->at(i).assign(uniqueIds.begin(), uniqueIds.end());
	}
}

void printAdjacencyList(TIntAdjList* network) {
	int numNode = network->size();
	int numAdj = 0;
	cout << numNode << endl;
	for (int i = 0; i < numNode; ++i) {
		cout << i << " : ";
		numAdj = (*network)[i].size();
		cout << numAdj << " => ";
		for (int j = 0; j < numAdj; ++j) {
			cout << (*network)[i][j] << " ";
		}
		cout << endl;
	}
}

void DFS(TIntAdjList* network, int geneId) {
	// Mark all the vertices as not visited
	int numNode = network->size();
	bool *visited = new bool[numNode];
	for (int i = 0; i < numNode; i++)
		visited[i] = false;

	// Call the recursive helper function to print DFS traversal
	DFSUtil(network, geneId, visited);
	delete visited;
}

void DFSUtil(TIntAdjList* network, int geneId, bool visited[]) {
	// Mark the current node as visited and print it
	visited[geneId] = true;
	cout << "Current " << geneId << endl;

	// Recur for all the vertices adjacent to this vertex
	vector<int> adj = (*network)[geneId];
	int numAdj = adj.size();
	//cout << "Number of neighbors = " << numAdj << endl;
	for (int i = 0; i < numAdj; ++i) {
		if (!visited[(*network)[geneId][i]]) {
			DFSUtil(network, (*network)[geneId][i], visited);
		}
	}
}

void BFS(TIntAdjList* network, int geneId) {
	// Mark all the vertices as not visited
	int numNode = network->size();
	bool *visited = new bool[numNode];
	for (int i = 0; i < numNode; i++)
		visited[i] = false;

	// Call the helper function to print BFS traversal
	BFSUtil(network, geneId, visited);
	delete visited;
}

void BFSUtil(TIntAdjList* network, int geneId, bool visited[]) {
	// Create queue
	queue<int> q;
	// Mark the current node as visited and print it
	q.push(geneId);
	visited[geneId] = true;
	int current;

	while (!q.empty()) {
		current = q.front();
		q.pop();
		//cout << "Current " << current << endl;

		//explore all the connected nodes
		vector<int> adj = (*network)[current];
		int numAdj = adj.size();
		for (int i = 0; i < numAdj; i++) {
			if (!visited[(*network)[current][i]]) {
				q.push((*network)[current][i]);
				visited[(*network)[current][i]] = true;
			}
		}
	}
}

int getNodeDegree(TIntAdjList* network, int nodeId) {
	return network->at(nodeId).size();
}

/*
 * For output
 */

void writeStrVector(const char* filename, vector<string>* output) {
	ofstream outFile;
	outFile.open(filename);

	int size = output->size();
	for (int i = 0; i < size; i++) {
		outFile << (*output)[i] << endl;
	}
	outFile.close();
}

void printGeneSymbols(vector<int>* geneIds, vector<string>* geneIdToSymbol) {
	vector<string> geneSymbols;
	mapGeneIdToSymbol(geneIds, &geneSymbols, geneIdToSymbol);
	int size = geneSymbols.size();
	for (int i = 0; i < size; ++i) {
		cout << geneSymbols[i] << endl;
	}
}

void saveGeneSymbols(const char* filename, vector<int>* geneIds, vector<string>* geneIdToSymbol) {
	vector<string> geneSymbols;
	mapGeneIdToSymbol(geneIds, &geneSymbols, geneIdToSymbol);

	writeStrVector(filename, &geneSymbols);
}

/*
 * For permutation
 */

void createPermutation(vector<int>* rank) {
	int size = rank->size();
	for (int i = 0; i < size; i++) {
		rank->at(i) = i;
	}
	random_shuffle(rank->begin(), rank->end());
}

void permuteGeneLabels(vector<int>* originalGeneLabels,
		vector<int>* permutedGeneLabels) {
	int numGenes = originalGeneLabels->size();
	//put all genes into a new vector
	for (int i = 0; i < numGenes; ++i) {
		permutedGeneLabels->push_back(originalGeneLabels->at(i));
	}
	//permute the labels
	random_shuffle(permutedGeneLabels->begin(), permutedGeneLabels->end());
}

/*
 * For ID mapping
 */

void mapGeneSymbolToId(vector<string>* symbols, vector<int>* ids,
		map<string, int>* geneSymbolToId) {
	int size = symbols->size();
	for (int i = 0; i < size; ++i) {
		ids->push_back(geneSymbolToId->find(symbols->at(i))->second);
	}
}

void mapGeneIdToSymbol(vector<int>* ids, vector<string>* symbols,
		vector<string>* geneIdToSymbol) {
	int size = ids->size();
	for (int i = 0; i < size; ++i) {
		symbols->push_back(geneIdToSymbol->at(ids->at(i)));
	}

}

/*
 *
 */

string intToStr(int i){
	stringstream ss;
	ss << i;
	return ss.str();
}

string doubleToStr(double i, int prec){
	stringstream ss;
	ss.precision(prec);
	ss << fixed << i;
	return ss.str();
}

bool trimStr(std::string& str, const std::string& from) {
	size_t start_pos = str.find(from);
	if(start_pos == std::string::npos)
		return false;
	str.replace(start_pos, str.length(), "");
	return true;
}

