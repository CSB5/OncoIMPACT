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
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "../header/input.h"

void readGeneExpression(const char* filename, GeneExpression* geneExpression,
		char delim, map<string, int>* geneSymbolToId, vector<string>* sampleIdToName){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = geneExpression->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // for a gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samplesList;
		getline(inFile, samplesList);
		istringstream ss(samplesList);
		string s;
		getline(ss, s, delim);	//discard the first column (GENE header)
		while (ss) {	//for each column (sample)
			if (!getline(ss, s, delim))
				break;
			sampleIdToName->push_back(s);
		}

		while (inFile.good()) { //for each row (gene)
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<double> row;	//to save expression values of each gene

			int i = 0;	// for read gene symbol in the first column
			bool found = false;	// is a current gene found in the network

			while (ss) {	//for each column (sample)
				string s;

				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);	//add the gene id to geneEx
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atof(s.c_str()));
				}

				if(!found)	//not found in the network, so skip this row (gene)
					break;

				i++;	//go to the next column (sample)
			}

			if(found){
				geneExpression->matrix->push_back(row);
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
			vector<int> row;	//to save mutation values of each gene

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
				}else{		//read mutation values
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
				}else{		//read mutation values
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

void readGenesList(const char* filename, vector<int>* geneIds, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			map<string, int>::iterator it = geneSymbolToId->find(s);
			geneIds->push_back(it->second);

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}


void readGenesListUpDown(const char* filename, vector<int>* geneIdsUpDown, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	int totalGenes = geneSymbolToId->size();

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
		  	size_t foundUp = s.find("UP");
		  	size_t foundDown = s.find("DOWN");

		  	trimStr(s, "_");
		  	string currentGeneName = s;

			map<string, int>::iterator it = geneSymbolToId->find(currentGeneName);

		  	if (foundUp!=std::string::npos){
		  		//the gene is upregulated
		  		geneIdsUpDown->push_back(it->second);
		  	}

		  	if (foundDown!=std::string::npos){
		  		//the gene is downregulated
		  		geneIdsUpDown->push_back(it->second + totalGenes);
		  	}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}


void readBenchmarkGeneList(string benchmarkGeneListFilename, vector<int>* cancerBenchmarkGenes,
		map<string, int>* geneSymbolToId, vector<string>* cancerBenchmarkGeneNames){
	ifstream inFile;
	inFile.open(benchmarkGeneListFilename.c_str(), std::ifstream::in);

	map<string, int>::iterator end = geneSymbolToId->end();

	if (inFile.is_open()) {

		//for each row
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream rowStr(s);

			int i = 0;
			bool found = false;
			while (rowStr) {	//for each column
				string s;

				if (!getline(rowStr, s, '\t'))
					break;
				if(i == 0){ //read gene symbols at the first column

					//trim the whitespace at the end of the string
					s.erase(s.find_last_not_of(" \n\r\t")+1);

					map<string, int>::iterator it = geneSymbolToId->find(s);
					cancerBenchmarkGeneNames->push_back(s);
					if(it != end){	//found in the network
						found = true;
						cancerBenchmarkGenes->push_back(it->second);	//add the gene id to geneEx
					}else{
						found = false;
					}
				}
				else if(i == 7){	//read tumor type (somatic mutation)
				}

				if(!found)	//not found in the network, so skip this row (gene)
					break;

				i++;	//go to the next column (sample)
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}

}

int readCancerTypeList(const char* filename, map<string, string>* cancerTypeDataTypeMap){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	if (inFile.is_open()) {

		//for each row
		while (inFile.good()) {

			string cancerType = "";
			string dataType = "";

			string s;
			if (!getline(inFile, s))
				break;

			istringstream rowStr(s);

			int i = 0;
			bool found = false;
			while (rowStr) {	//for each column
				string s;

				if (!getline(rowStr, s, '\t'))
					break;
				if(i == 0){ 	//cancer type
					cancerType = s;
				}else if(i == 1){	//data type
					dataType = s;
				}

				i++;	//go to the next column (sample)
			}

			cancerTypeDataTypeMap->insert(pair<string, string>(cancerType, dataType));

		}
		inFile.close();
		return 0;
	} else {
		cerr << "Error opening cancer_type_list.dat file\n";
		return 1;
	}
}

void readPhenotypeGenesFromFile(const char* filename, vector<int>* phenotypeGeneIdsUpDown, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	int totalGenes = geneSymbolToId->size();

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;

			istringstream lineStream(line);
			int i = 0;
			string geneSym;
			bool isPhenotype = false;

			while (lineStream) {	//for each column
				string token;
				if (!getline(lineStream, token, '\t'))
					break;
				if(i == 0){			//gene symbol
					geneSym = token;
				}else if(i == 9){	//phenotype status
					if(token.compare("0") == 0){
						isPhenotype = true;
					}
				}
				if(isPhenotype){
					size_t foundUp = geneSym.find("_UP");
					size_t foundDown = geneSym.find("_DOWN");
				  	trimStr(geneSym, "_");
					map<string, int>::iterator it = geneSymbolToId->find(geneSym);
					if (foundUp!=std::string::npos){
						//the gene is upregulated
						phenotypeGeneIdsUpDown->push_back(it->second);
					}
					if (foundDown!=std::string::npos){
						//the gene is downregulated
						phenotypeGeneIdsUpDown->push_back(it->second + totalGenes);
					}
				}
				i++;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file of explained and phenotype gene list (PHENO.dat) \n";
	}
}

void readDrugsListFromFile(const char* filename, map<int, string>* drugIdToName, map<string, int>* drugNameToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;

			istringstream lineStream(line);
			int i = 0;
			string drugName;
			int drugId;

			while (lineStream) {	//for each column
				string token;
				if (!getline(lineStream, token, '\t'))
					break;
				if(i == 0){			//drug id
					drugId = atoi(token.c_str());
				}else if(i == 1){	//drug name
					drugName = token;
					break;
				}
				i++;
			}

			//save drug id and name
			drugIdToName->insert(pair<int, string>(drugId, drugName));
			drugNameToId->insert(pair<string, int>(drugName, drugId));
//			cout << drugName << endl;
		}
		inFile.close();
	} else {
		cerr << "Error opening file (GDSC_drug_list.txt) \n";
	}
}

void readGenesDrugsAssocFromFile(const char* filename, vector< vector<int> >* geneDrugsAssocList, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	//read all genes
	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network
	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;

			istringstream lineStream(line);
			int i = 0;
			string geneName;
			int drugId;
			bool found = false;	// found in the network

			while (lineStream) {	//for each column
				string token;
				if (!getline(lineStream, token, '\t'))
					break;
				if(i == 0){			//drug id
					drugId = atoi(token.c_str());
				}else if(i == 1){	//gene name
					geneName = token;
					map<string, int>::iterator it = geneSymbolToId->find(geneName);
					if(it != end){	//found in the network
						found = true;
						int geneIds = it->second;
						geneDrugsAssocList->at(geneIds).push_back(drugId);
					}else{
						found = false;
					}
					break;
				}
				i++;
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file (GDSC_gene_drugs_assoc.txt) \n";
	}

}

void readDriverGenesFromFile(const char* filename, vector<DriverGeneFromFile>* driverGenesFromFile, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	int totalGenes = geneSymbolToId->size();
	int countDriver = 0;

	if (inFile.is_open()) {

		string header;
		if(inFile.good()){
			getline(inFile, header); //skip the header
		}

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;
			istringstream lineStream(line);
			int i = 0;
			string geneSym;
			bool isPhenotype = false;

			int currentDriverGeneId = -1;
			while (lineStream) {	//for each column (parameter)
				string token;
				if (!getline(lineStream, token, '\t'))
					break;

//				GENE	DRIVER_FREQUENCY	DRIVER_SNV_FREQUENCY	DRIVER_DELTION_FREQUENCY	DRIVER_AMPLIFICATION_FREQUENCY	CANCER_CENSUS	PAN_CANCER	IMPACT	MUTATION_FREQUENCY	SNV_FREQUENCY	DELTION_FREQUENCY	AMPLIFICATION_FREQUENCY
//				EGFR	0.393	0.262	0.000	0.241	Y	NA	120.1447552328	0.393	0.262	0.000	0.241
				if(i == 0){	//gene symbol
					map<string, int>::iterator it = geneSymbolToId->find(token);
					currentDriverGeneId = it->second;
					countDriver++;
				}else if(i == 1){
					driverGenesFromFile->at(currentDriverGeneId).driverFreq = atof(token.c_str());
				}else if(i == 2){
					driverGenesFromFile->at(currentDriverGeneId).driverSnpFreq = atof(token.c_str());
				}else if(i == 3){
					driverGenesFromFile->at(currentDriverGeneId).driverDelFreq = atof(token.c_str());
				}else if(i == 4){
					driverGenesFromFile->at(currentDriverGeneId).driverAmpFreq = atof(token.c_str());
				}else if(i == 5){
					if(token.compare("Y") == 0){
						driverGenesFromFile->at(currentDriverGeneId).isInCancerCensus = true;
					}else{
						driverGenesFromFile->at(currentDriverGeneId).isInCancerCensus = false;
					}
				}else if(i == 7){
					driverGenesFromFile->at(currentDriverGeneId).impactScore = atof(token.c_str());
				}else if(i == 8){
					driverGenesFromFile->at(currentDriverGeneId).mutFreq = atof(token.c_str());
				}else if(i == 9){
					driverGenesFromFile->at(currentDriverGeneId).snpFreq = atof(token.c_str());
				}else if(i == 10){
					driverGenesFromFile->at(currentDriverGeneId).delFreq = atof(token.c_str());
				}else if(i == 11){
					driverGenesFromFile->at(currentDriverGeneId).ampFreq = atof(token.c_str());
				}
				i++;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file of explained and phenotype gene list \n";
	}

}

void readMutatedGenesFromFile(const char* filename, vector<MutatedGeneFromFile>* mutatedGenesFromFile, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	int totalGenes = geneSymbolToId->size();
	int countDriver = 0;

	if (inFile.is_open()) {

//		string header;
//		if(inFile.good()){
//			getline(inFile, header); //skip the header
//		}

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;
			istringstream lineStream(line);
			int i = 0;
			string geneSym;
			bool isPhenotype = false;

			int currentDriverGeneId = -1;
			while (lineStream) {	//for each column (parameter)
				string token;
				if (!getline(lineStream, token, '\t'))
					break;

//gene, mutation freq, snp freq, del, amp, snp+cnv, driver freq, driver snp, driver del, driver amp, driver snp+cnv, cancer census, unused score, impact score
				if(i == 0){	//gene symbol
					map<string, int>::iterator it = geneSymbolToId->find(token);
					currentDriverGeneId = it->second;
					countDriver++;
				}else if(i == 1){
					mutatedGenesFromFile->at(currentDriverGeneId).mutationFreq = atof(token.c_str());
				}else if(i == 2){
					mutatedGenesFromFile->at(currentDriverGeneId).snpFreq = atof(token.c_str());
				}else if(i == 3){
					mutatedGenesFromFile->at(currentDriverGeneId).delFreq = atof(token.c_str());
				}else if(i == 4){
					mutatedGenesFromFile->at(currentDriverGeneId).ampFreq = atof(token.c_str());
				}else if(i == 5){
					mutatedGenesFromFile->at(currentDriverGeneId).snpAndCnvFreq = atof(token.c_str());
				}else if(i == 6){
					mutatedGenesFromFile->at(currentDriverGeneId).driverFreq = atof(token.c_str());
				}else if(i == 7){
					mutatedGenesFromFile->at(currentDriverGeneId).driverSnpFreq = atof(token.c_str());
				}else if(i == 8){
					mutatedGenesFromFile->at(currentDriverGeneId).driverDelFreq = atof(token.c_str());
				}else if(i == 9){
					mutatedGenesFromFile->at(currentDriverGeneId).driverAmpFreq = atof(token.c_str());
				}else if(i == 10){
					mutatedGenesFromFile->at(currentDriverGeneId).driverSnpAndCnvFreq = atof(token.c_str());
//				}else if(i == 11){
//					if(token.compare("CC") == 0){
//						mutatedGenesFromFile->at(currentDriverGeneId).isInCancerCensus = true;
//					}else{
//						mutatedGenesFromFile->at(currentDriverGeneId).isInCancerCensus = false;
//					}
				}else if(i == 13){
					mutatedGenesFromFile->at(currentDriverGeneId).impactScore = atof(token.c_str());
				}
				i++;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file of explained and phenotype gene list \n";
	}

}

void readModulesFromFile(string* moduleFileName, vector<string>* sampleIdToName, map<string, int>* sampleNameToId,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId,
		vector< vector<MutatedAndExplianedGenes> >* mutatedAndExplainedGenesListReal, vector< vector<int> >* mutatedGeneIdsListReal){
	ifstream inFile;
	char delim = '\t';

	int numInputSamples = sampleIdToName->size();
	int currentSampleId = numInputSamples;	//pool input sample modules and database modules together

	//read samples and their ids
	inFile.open(moduleFileName->c_str(), std::ifstream::in);

	if (inFile.is_open()) {

		while (inFile.good()) {
			string rowStr;
			if (!getline(inFile, rowStr))
				break;

			istringstream rowStream(rowStr);

			int i = 0;		// for read sample name in the first column;

			while (rowStream) {

				string colStr;
				if (!getline(rowStream, colStr, delim))
					break;

				if(i == 0){		//read sample name

					//check if sample's name already added to the module list
					vector<string>::iterator it = find(sampleIdToName->begin(), sampleIdToName->end(), colStr);
					if( it == sampleIdToName->end()){	//not found
						//add the sample to the list
						sampleIdToName->push_back(colStr);
						sampleNameToId->insert(pair<string, int>(colStr, currentSampleId));
//						cout << "sample " << currentSampleId << " = " << colStr << " is added\n";
						currentSampleId++;
					}else{	//found in the module list, then skip
					}

					break;
				}
				i++;
			}

		}
		inFile.close();

	} else {
		cerr << "Error opening file\n";
	}

	//read modules into mutatedAndExplainedGenesListReal

	inFile.open(moduleFileName->c_str(), std::ifstream::in);

	currentSampleId = -1;
	int currentMutatedGeneId = -1;
	int numDbSamples = sampleIdToName->size() - numInputSamples;
	int totalGenes = geneIdToSymbol->size();
	int totalGenesUpDown = totalGenes * 2;

	//initialize vector for samples in database
	for (int si = 0; si < numDbSamples; ++si) {
		mutatedGeneIdsListReal->push_back(vector<int>());
		vector<MutatedAndExplianedGenes> mutatedAndExplainedGenes(totalGenes);
		for (int gi = 0; gi < totalGenes; ++gi) {
			mutatedAndExplainedGenes[gi].isExplainedGenesUpDown = new vector<bool>(totalGenes * 2, false);
		}
		mutatedAndExplainedGenesListReal->push_back(mutatedAndExplainedGenes);
	}

	if (inFile.is_open()) {

		while (inFile.good()) {
			string rowStr;
			if (!getline(inFile, rowStr))
				break;

			vector<string> members;
			vector<string> drivers;

			istringstream rowStream(rowStr);

			int i = 0;	// for read gene in the first column;

			while (rowStream) {

				string colStr;
				if (!getline(rowStream, colStr, delim))
					break;

				if(i == 0){					//sample name
//					cout << colStr << endl;
					map<string, int>::iterator it = sampleNameToId->find(colStr);
					if(it == sampleNameToId->end()){	//not found, so this sample is one of the input samples
						currentSampleId = -1;
						break;	//skip to the next module (line)
					}else{
						currentSampleId = it->second;
					}

				}else if(i == 1){			//read mutated gene
//					cout << currentSampleId << " " << currentMutatedGeneId << endl;
					map<string, int>::iterator it = geneSymbolToId->find(colStr);
					currentMutatedGeneId = it->second;
					mutatedGeneIdsListReal->at(currentSampleId).push_back(currentMutatedGeneId);

				}else if(i == 2){			//read explained genes
					istringstream geneList(colStr);
					while(geneList){

						string gene;
						if (!getline(geneList, gene, ';'))
							break;

						size_t foundUp = gene.find("_UP");
						size_t foundDown = gene.find("_DOWN");
					  	trimStr(gene, "_");
						map<string, int>::iterator it = geneSymbolToId->find(gene);
						if (foundUp!=std::string::npos){
							//the gene is upregulated
							int currentExplainedGeneId = it->second;
							mutatedAndExplainedGenesListReal->at(currentSampleId)[currentMutatedGeneId].isExplainedGenesUpDown->at(currentExplainedGeneId) = true;
//							cout << gene << "_UP" << "\t";
						}
						if (foundDown!=std::string::npos){
							//the gene is downregulated
							int currentExplainedGeneId = it->second + totalGenes;
							mutatedAndExplainedGenesListReal->at(currentSampleId)[currentMutatedGeneId].isExplainedGenesUpDown->at(currentExplainedGeneId) = true;
//							cout << gene << "_DOWN" << "\t";
						}
					}
//					cout << endl;

					break;	//go to the next line
				}

				i++;	//go to the next column
			}
		}
		inFile.close();


	} else {
		cerr << "Error opening file\n";
	}

}

bool replaceStr(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if(start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}

