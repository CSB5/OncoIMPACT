/*
 * results.cpp
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#include <iostream>
#include <string>
#include "results.h"
#include "utilities.h"

void saveModules(vector< list<Module> > * modulesListOfAllSamples, string filename, vector<string>* geneIdToSymbol){
	int totalSamples = modulesListOfAllSamples->size();
	vector<string> outputStr;
	outputStr.push_back("GENE\tTYPE\tSAMPLE_ID\tMODULE");

	//get a list of nodes (gene symbol, type)
	vector<int> geneIds;
	for (int i = 0; i < totalSamples; ++i) {

		list<Module> modulesList = modulesListOfAllSamples->at(i);

		//for each module
		int j = 0;
		for(list<Module>::iterator it = modulesList.begin(); it != modulesList.end(); it++){

			//get driver genes
			for(list<int>::iterator g = it->driverGeneIds.begin(); g != it->driverGeneIds.end(); g++){
				string str = "" + geneIdToSymbol->at(*g) + "\tDRIVER\t" + intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			//get phenotype genes
			for(list<int>::iterator g = it->phenotypeGeneIds.begin(); g != it->phenotypeGeneIds.end(); g++){
				string str = "" + geneIdToSymbol->at(*g) + "\tPHENOTYPE\t" + intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			//get explained genes
			for(list<int>::iterator g = it->explainedGeneIds.begin(); g != it->explainedGeneIds.end(); g++){
				string str = "" + geneIdToSymbol->at(*g) + "\tEXPLAINED\t" + intToStr(i) + "\t" + intToStr(j);
				outputStr.push_back(str);
			}

			j++;
		}
	}

	//get a list of edges (gene symbol, gene symbol)
	//do not need now because the original network file can be used

	writeStrVector(filename.c_str(), &outputStr);
}

string geneIdsToNodesList(const vector<string>& geneSymbols) {
	string str;
	vector<string>::const_iterator it;
	for( it = geneSymbols.begin(); it != geneSymbols.end(); ++it){
		str += *it;
		str += "\n";
	}
	return str;
}

//st = myString.substr(0, myString.size()-1);



