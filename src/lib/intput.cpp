/*
 * intput.cpp
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <cmath>
#include <map>
#include <cstdlib>
#include "../header/input.h"

void readGeneExpression(const char* filename, GeneExpression* geneExpression,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = geneExpression->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		while (inFile.good()) { //for each gene
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<double> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// is a currrent gene found in the network

			while (ss) {	//for each column
				string s;

				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						//cout << "\tfound " << s;
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atof(s.c_str()));
				}

				if(!found)	//not found in the network, so skip this gene (row)
					break;

				i++;
			}

			if(found){
				geneExpression->matrix->push_back(row);
				//cout << " and added" << endl;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readMutations(const char* filename, Mutations* mutations,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = mutations->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atoi(s.c_str()));
				}

				if(!found)
					break;

				i++;
			}

			if(found){
				mutations->matrix->push_back(row);
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readPointMutations(const char* filename, PointMutations* pointMutations,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = pointMutations->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atoi(s.c_str()));
				}

				if(!found)
					break;

				i++;
			}

			if(found){
				pointMutations->matrix->push_back(row);
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readCopyNumberVariation(const char* filename, CopyNumberVariation* copyNumberVariation,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = copyNumberVariation->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atoi(s.c_str()));
				}

				if(!found)
					break;

				i++;
			}

			if(found){
				copyNumberVariation->matrix->push_back(row);
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

int findIndex(vector<int>* geneIds, int currentGeneId){
	int index = -1;	//not found
	int size = geneIds->size();

	for (int i = 0; i < size; ++i) {
		if(geneIds->at(i) == currentGeneId){
			return i;
		}
	}

	return index;
}




