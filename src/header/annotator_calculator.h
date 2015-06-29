/*
 * calculator.h
 *
 *  Created on: Apr 27, 2015
 *      Author: nok
 */

#ifndef ANNOATATOR_CALCULATOR_H_
#define ANNOATATOR_CALCULATOR_H_

#include <vector>
#include <string>

using namespace std;

struct GeneSetPair{
	string geneSetName;
	string geneListName;
	string driverGeneList;
	int sizeOfGeneSet;
	int sizeOfGeneList;
	int numOverlap;
	long double pValue;
};


int countOverlap(vector<string>* currentModuleMembers, vector<string>* currentGeneSetMembers);
void calculateHypergeometric(GeneSetPair* pair, int N, int numGeneSets);
long double hypergeom(int n, int m, int N, int i);
long double logfact(int x);
long double gammln(int xx);

bool sortByPValues(const GeneSetPair& first, const GeneSetPair& second);


//string intToStr(int i);
//string doubleToStr(double i, int prec);

#endif /* CALCULATOR_H_ */
