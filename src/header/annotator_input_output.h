/*
 * input.h
 *
 *  Created on: Apr 27, 2015
 *      Author: nok
 */

#ifndef ANNOTATOR_INPUT_H_
#define ANNOTATOR_INPUT_H_

#include <vector>
#include <string>
#include <set>

using namespace std;

struct GeneSetCount{
	string name;
	int count;
};

void readGeneList(string* filename, vector<string>* geneList);
void readGeneSets(vector<string>* filenames, vector<string>* geneSetNames, vector<string>* geneSetUrls,
		vector< vector<string> >* geneSetMembers, vector<string>* geneList);
void readModules(string* moduleFileName, vector<string>* moduleNames, vector< vector<string> >* moduleMembers,
		vector< vector<string> >* moduleDrivers, set<string>* driversList, set<string>* samplesList);

//bool replaceStr(std::string& str, const std::string& from, const std::string& to);
//bool trimStr(std::string& str, const std::string& from);

//void writeStrVector(const char* filename, vector<string>* output);
int findIndexStr(vector<string>* names, string name);

bool sortByCount(const GeneSetCount& first, const GeneSetCount& second);

#endif /* INPUT_H_ */
