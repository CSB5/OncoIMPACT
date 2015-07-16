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
#include "database.h"
#include "discovery.h"
#include "annotator.h"
#include "header/utilities.h"

using namespace std;

int main( int argc, char *argv[] ) {

	/*
	 * START UP
	 */

	//oncoIMPACT --database|--discovery <config_file>

	string modeStr;
	string configFilename = "userConfig.cfg";

	int oncoIMPACTMode = -1; //0 = database, 1 = discovery

	if ( argc != 3 ){ // argc should be 3 for correct execution
		cout<<"usage: "<< argv[0] <<" --discovery <config_file>\n";
		return 1;
	}else{
		modeStr = argv[1];
		configFilename = argv[2];

		string temp = modeStr.substr (2, modeStr.length()-2);

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

	string outDir = "";

	string networkFilename;
	string expFilename = "EXPR.txt";
	string snpFilename = "SNP.txt";
	string cnvFilename = "CNV.txt";
	string benchmarkGeneListFilename;

	string dbPath = "database";

	int numThreads = 1;

	//0 = sensitive, 1 = stringent
	//int mode = 1; //default mode = stringent //Now only the stringent mode is in used

	//0 = ARRAY, 1 = RNA_SEQ
	int dataType = -1;
	int inputDataType = -1;
	bool isDataTypeMatch = false;
	bool noFoldchangeCutoff = false;
	string cancerType;
	map<string, string> cancerTypeDataTypeMap;

	//Read the cancer type list
	string cancerListFilename = dbPath + "/cancer_type_list.dat";
	if(readCancerTypeList(cancerListFilename.c_str(), &cancerTypeDataTypeMap) == 0){
		cout << "read cancer type list" << endl;
	}else{
		cerr << "read cancer type list" << endl;
		return 1;
	}

	ifstream inConfigFile;
	inConfigFile.open(configFilename.c_str(), std::ifstream::in);

	int countArvInConfigFile = 0;
	if (inConfigFile.is_open()) {

		//read the output from the
		while (inConfigFile.good()) {
			string line;
			if (!getline(inConfigFile, line))
				break;

			//comment or empty lines
			if(line[0] == '#' || line.size() == 0){
				continue;
			}

			istringstream lineStrem(line);
			string configName;

			int i = 0;
			while (lineStrem) {	//for each column (sample)
				string token;

				if (!getline(lineStrem, token, '='))
					break;

				if(i == 0){ //read config name
					configName = token;
				}else if(i == 1){	//read config content
//					if(configName.compare("outDir") == 0){
//						outDir = token;
//						countArvInConfigFile++;
//						cout << configName << " = " << outDir << endl;
//					}else if(configName.compare("numThreads") == 0){
//						countArvInConfigFile++;
//						numThreads = atoi(token.c_str());
//						cout << configName << " = " << numThreads << endl;
//					}else if(configName.compare("exp") == 0){
//						countArvInConfigFile++;
//						expFilename = token;
//						cout << configName << " = " << expFilename << endl;
//					}else if(configName.compare("snp") == 0){
//						countArvInConfigFile++;
//						snpFilename = token;
//						cout << configName << " = " << snpFilename << endl;
//					}else if(configName.compare("cnv") == 0){
//						countArvInConfigFile++;
//						cnvFilename = token;
//						cout << configName << " = " << cnvFilename << endl;
//					}else if(configName.compare("dbPath") == 0){
//						countArvInConfigFile++;
//						dbPath = token;
//						cout << configName << " = " << dbPath << endl;

//						//Read the cancer type list
//						string cancerListFilename = dbPath + "/cancer_type_list.dat";
//						if(readCancerTypeList(cancerListFilename.c_str(), &cancerTypeDataTypeMap) == 0){
//							cout << "read cancer type list" << endl;
//						}else{
//							cerr << "read cancer type list" << endl;
//							return 1;
//						}


					if(configName.compare("cancerType") == 0){
						countArvInConfigFile++;
						cancerType = token;
						cout << configName << " = " << cancerType << endl;

						//check the data type of the input cancer type
						map<string, string>::iterator it = cancerTypeDataTypeMap.find(cancerType);
						if(it == cancerTypeDataTypeMap.end()){
							cerr << "oncoIMPACT does not have the input cancer type in database" << endl;
							return 1;
						}else{
							string dataTypeStr = it->second;
							if(dataTypeStr.compare("ARRAY") == 0){	//use sensitive mode
								dataType = 0;
							}else if(dataTypeStr.compare("RNA_SEQ") == 0){
								dataType = 1;
							}

						}

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
							cout << "The data type of the database matches with the database" << endl;
						}else{
							isDataTypeMatch = false;
							cout << "The data type of the input DOES NOT MATCH with the database" << endl;
							cout << "Please uncomment 'noFoldChangeCutOff=yes' in the configuration file to allow oncoIMPACT run (ignore this if you already uncommented the line)" << endl;
						}
//					}else if(configName.compare("noFoldChangeCutOff") == 0){
//						if(token.compare("yes") == 0){
//							noFoldchangeCutoff = true;
//						}else{
//							noFoldchangeCutoff = false;
//						}
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

		//check matching of data type
		if(!isDataTypeMatch){
			if(noFoldchangeCutoff){
				cout << "oncoIMPACT WILL NOT use any cutoff value for differential expression" << endl;
			}else{
				//cout << "Please uncomment 'noFoldChangeCutOff=yes' in the configuration file to allow oncoIMPACT run" << endl;
				cout << "No fold change cutoff is used when the data type of the input DOES NOT MATCH with the database" << endl;
				return 1;
			}
		}


		networkFilename = dbPath + "/network_FIsInGene_041709.txt";
		benchmarkGeneListFilename = dbPath + "/Census_all_04_06_2015.tsv";

		cout << "Start constructing oncoIMPACT database ..." << endl;
		discovery(outDir, networkFilename, expFilename, snpFilename, cnvFilename,
				benchmarkGeneListFilename, dbPath, numThreads, cancerType, noFoldchangeCutoff);

		cout << "Start annotate the results ..." << endl;

		string mSigDbPath = dbPath + "/MSigDB/";
		//string outputPrefix = outDir + "/";
		string outputPrefix = outDir;
		//string moduleFileName = outDir + "/FINAL_MODULE.dat";
		string moduleFileName = "FINAL_MODULE.dat";
		double cutoff = 0.05;
		int top = 10;
		string geneListFileName = dbPath + "/driver_net_background_gene_list.dat";
		annotator(mSigDbPath, moduleFileName, outputPrefix, cutoff, top, geneListFileName);

	}


	return 0;
}
