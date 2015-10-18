//============================================================================
// Name        : oncoIMPACT.cpp
// Author      : Nok C Suphavilai
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "discovery.h"
#include "database.h"
#include "annotator.h"
#include "header/utilities.h"

using namespace std;

int main( int argc, char *argv[] ) {

	/*
	 * START UP
	 */

	//TO RUN: oncoIMPACT --database|--discovery userConfig.cfg

	string modeStr;
	string configFilename = "userConfig.cfg";

	int oncoIMPACTMode = -1; //0 = database, 1 = discovery

	if ( argc != 3 ){ // argc should be 3 for correct execution
		cout<<"usage: "<< argv[0] <<" --discovery userConfig.cfg\n";
		return 1;
	}else{
		modeStr = argv[1];
		configFilename = argv[2];

		// check if the mode is valid
		string temp = modeStr.substr(2, modeStr.length()-2);
		size_t found = temp.find("discovery");
		if (found != string::npos){
			cout << "Start onocIMPACT Discovery Mode\n";
			oncoIMPACTMode = 1;
		}else{
			cout << "invalid mode\n";
			cout<<"usage: "<< argv[0] <<" --discovery <config_file>\n";
			return 1;
		}
	}


	/*
	 * DISCOVERY
	 */

	// output at current directory
	string outDir = "";

	// database directory name
	string dbPath = "database";

	// list of input files
	string expFilename = "EXPR.txt";
	string snpFilename = "SNP.txt";
	string cnvFilename = "CNV.txt";
	string networkFilename = dbPath + "/network_FIsInGene_041709.txt";
	string benchmarkGeneListFilename = dbPath + "/Census_all_04_06_2015.tsv";

	// Driver discovery modes can be re-considered in the future
	// 0 = sensitive, 1 = stringent
	// Now only the stringent mode is in used

	// Gene expression data type
	//0 = ARRAY, 1 = RNA_SEQ
	int dataType = -1;
	int inputDataType = -1;
	bool isDataTypeMatch = false;
	bool noFoldchangeCutoff = false;

	// Type of gene expression in the database
	string cancerType;
	map<string, string> cancerTypeDataTypeMap;

	// Read the cancer type list
	string cancerListFilename = dbPath + "/cancer_type_list.dat";
	if(readCancerTypeList(cancerListFilename.c_str(), &cancerTypeDataTypeMap) == 0){
		cout << "read cancer type list" << endl;
	}else{
		cerr << "read cancer type list" << endl;
		return 1;
	}

	// Open the config file
	ifstream inConfigFile;
	inConfigFile.open(configFilename.c_str(), std::ifstream::in);

	int countArvInConfigFile = 0;
	if (inConfigFile.is_open()) {

		// for each line
		while (inConfigFile.good()) {
			string line;
			if (!getline(inConfigFile, line))
				break;

			// Skip comment or empty lines
			if(line[0] == '#' || line.size() == 0){
				continue;
			}

			istringstream lineStrem(line);
			string configName;

			int i = 0;
			while (lineStrem) {	//for each column
				string token;

				if (!getline(lineStrem, token, '='))
					break;

				if(i == 0){ //read config name
					configName = token;
				}else if(i == 1){	//read config content

					// type of cancer
					if(configName.compare("cancerType") == 0){
						countArvInConfigFile++;
						cancerType = token;
						cout << configName << " = " << cancerType << endl;

						// check the data type of the input cancer type
						map<string, string>::iterator it = cancerTypeDataTypeMap.find(cancerType);
						if(it == cancerTypeDataTypeMap.end()){
							cerr << "oncoIMPACT does not have the input cancer type in database" << endl;
							return 1;
						}else{
							// type of gene expression in the database (for the input cancer type)
							string dataTypeStr = it->second;
							if(dataTypeStr.compare("ARRAY") == 0){	//use sensitive mode
								dataType = 0;
							}else if(dataTypeStr.compare("RNA_SEQ") == 0){
								dataType = 1;
							}

						}

					// type of gene expression data
					}else if(configName.compare("dataType") == 0){
						countArvInConfigFile++;
						if(token.compare("ARRAY") == 0){	//use sensitive mode
							inputDataType = 0;
						}else if(token.compare("RNA_SEQ") == 0){
							inputDataType = 1;
						}
						cout << configName << " = " << token << endl;

						//Check input data type
						if(dataType == inputDataType and dataType != -1){
							isDataTypeMatch = true;
							cout << "The data type of gene expression in the database matches with the input" << endl;
						}else{
							isDataTypeMatch = false;
							noFoldchangeCutoff = true;
							cout << "The data type of the input gene expression DOES NOT MATCH with the database" << endl;
						}
					}
				}

				i++;
			}

		}
		inConfigFile.close();
	} else {
		cerr << "Error opening configuration file\n";
		return 1;
	}

	if(countArvInConfigFile != 2){
		cout << "Incorrect configuration file\n";
		return 1;
	}else{

		/* Run oncoIMPACT discovery mode */

		// Check matching of data type
		if(!isDataTypeMatch){
			if(noFoldchangeCutoff){
				cout << "oncoIMPACT WILL NOT use any cutoff value for differential expression" << endl;
			}
		}

		cout << "Start constructing oncoIMPACT database ..." << endl;
		if(!discovery(outDir, networkFilename, expFilename, snpFilename, cnvFilename,
				benchmarkGeneListFilename, dbPath, cancerType, noFoldchangeCutoff)){
			return 1;
		}

		cout << "Start annotate the results ..." << endl;
		string mSigDbPath = dbPath + "/MSigDB/";
		string outputPrefix = outDir;
		string moduleFileName = "FINAL_MODULE.dat";
		double cutoff = 0.05;	// p-value for hypergeometric test
		int top = 10;	// number of gene sets most enriched for the modules
		string geneListFileName = dbPath + "/driver_net_background_gene_list.dat";
		annotator(mSigDbPath, moduleFileName, outputPrefix, cutoff, top, geneListFileName);

	}


	return 0;
}
