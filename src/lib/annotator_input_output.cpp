/*
 * input.cpp
 *
 *  Created on: Apr 27, 2015
 *      Author: nok
 */

#include "../header/annotator_input_output.h"
#include "../header/annotator_calculator.h"
#include "../annotator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

void readGeneList(string* filename, vector<string>* geneList){
	ifstream inFile;

	inFile.open(filename->c_str(), std::ifstream::in);

	int totalGenes = 0;

	if (inFile.is_open()) {

		while (inFile.good()) {
			string rowStr;	//for each gene set
			if (!getline(inFile, rowStr))
				break;
			geneList->push_back(rowStr);
			totalGenes++;

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}

//	cout << "total number of genes is " << totalGenes << endl;

}


void readGeneSets(vector<string>* filenames, vector<string>* geneSetNames, vector<string>* geneSetUrls,
		vector< vector<string> >* geneSetMembers, vector<string>* geneList){
	ifstream inFile;
	char delim = '\t';

//	set<string> added;
//	set<string> deleted;

	int numFiles = filenames->size();
	for (int f = 0; f < numFiles; ++f) {
		string filename = filenames->at(f);

		inFile.open(filename.c_str(), std::ifstream::in);

		if (inFile.is_open()) {

			while (inFile.good()) {
				string rowStr;	//for each gene set
				if (!getline(inFile, rowStr))
					break;

//				cout << rowStr << endl;

				istringstream ss(rowStr);

				int i = 0;	// for read first two columns
				vector<string> members;

				while (ss) {

					string s;
					if (!getline(ss, s, delim))
						break;

					if(i == 0){ //read gene set name
//						geneSetNames->push_back(intToStr(f) + "_" + s);
						geneSetNames->push_back(s);
					}else if(i == 1){
						geneSetUrls->push_back(s);
					}else{		//read gene symbol
						if(find(geneList->begin(), geneList->end(), s) != geneList->end()){
							members.push_back(s);
//							added.insert(s);
//						}else{
//							deleted.insert(s);
						}
					}

					i++;
				}

				geneSetMembers->push_back(members);

			}

			inFile.close();

		} else {
			cerr << "Error opening file\n";
		}
	}

//	cout << "added " << added.size() << " and deleted " << deleted.size() << endl;

}

void readModules(string* moduleFileName, vector<string>* moduleNames, vector< vector<string> >* moduleMembers,
		vector< vector<string> >* moduleDrivers, set<string>* driversList, set<string>* samplesList){
	ifstream inFile;
	char delim = '\t';


	inFile.open(moduleFileName->c_str(), std::ifstream::in);

	if (inFile.is_open()) {

		while (inFile.good()) {
			string rowStr;
			if (!getline(inFile, rowStr))
				break;

//				cout << rowStr << endl;
			vector<string> members;
			vector<string> drivers;

			istringstream rowStream(rowStr);

			int i = 0;	// for read gene in the first column;

			while (rowStream) {

				string colStr;
				if (!getline(rowStream, colStr, delim))
					break;

				if(i == 0){ 				//read module name
					moduleNames->push_back(colStr);
					trimStr(colStr, ".");
					string sampleName = colStr;
					samplesList->insert(sampleName);
				}else if(i == 1){			//read driver

					istringstream geneList(colStr);
					while(geneList){

						string gene;
						if (!getline(geneList, gene, ','))
							break;

						members.push_back(gene);
						drivers.push_back(gene);
						driversList->insert(gene);
					}
				}else if(i == 2 or i == 3){						//read phenotype or explained gene

						istringstream geneList(colStr);

						while(geneList){

							string gene;
							if (!getline(geneList, gene, ','))
								break;

							if(gene.compare("-") != 0){

								replaceStr( gene, "_UP", "");
								replaceStr( gene, "_DOWN", "");
								members.push_back(gene);
//								cout << gene << endl;
							}

						}
				}

				i++;
			}

			moduleMembers->push_back(members);
			moduleDrivers->push_back(drivers);
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}

}

//bool replaceStr(std::string& str, const std::string& from, const std::string& to) {
//	size_t start_pos = str.find(from);
//	if(start_pos == std::string::npos)
//		return false;
//	str.replace(start_pos, from.length(), to);
//	return true;
//}
//
//bool trimStr(std::string& str, const std::string& from) {
//	size_t start_pos = str.find(from);
//	if(start_pos == std::string::npos)
//		return false;
//	str.replace(start_pos, str.length(), "");
//	return true;
//}

//void writeStrVector(const char* filename, vector<string>* output) {
//	ofstream outFile;
//	outFile.open(filename);
//
//	int size = output->size();
//	for (int i = 0; i < size; i++) {
//		outFile << (*output)[i] << endl;
//	}
//	outFile.close();
//}

int findIndexStr(vector<string>* names, string name){
	int index = -1;	//not found
	int size = names->size();

	for (int i = 0; i < size; ++i) {
		if(names->at(i).compare(name) == 0){
//			cout << "matched " << names->at(i) << " and " << name << endl;
			return i;
		}
	}
	return index;
}

bool sortByCount(const GeneSetCount& first, const GeneSetCount& second) {
	if (first.count > second.count) {
		return true;
	} else {
		return false;
	}
}
