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
#include "header/utilities.h"

using namespace std;

int main( int argc, char *argv[] ) {

	/*
	 * START UP
	 */

	//oncoIMPACT --database|--discovery <config_file>

	string modeStr;
	string configFilename;

	int oncoIMPACTMode = -1; //0 = database, 1 = discovery

	if ( argc != 3 ){ // argc should be 3 for correct execution
		cout<<"usage: "<< argv[0] <<" --discovery <config_file>\n";
		return 0;
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
			return 0;
		}
	}

//	/*
//	 * DATABASE
//	 */
//
//	if(oncoIMPACTMode == 0){
//
//		//variables
//		string outDir;
//		string networkFilename = "";
//		string expFilename;
//		string snpFilename;
//		string cnvFilename;
//		string benchmarkGeneListFilename;
//		int numThreads = 1;
//		//0 = sensitive, 1 = stringent
//		int mode = 0; //default mode = sensitive
//
//		//read config file
//		ifstream inConfigFile;
//		inConfigFile.open(configFilename.c_str(), std::ifstream::in);
//		int countArvInConfigFile = 0;
//		if (inConfigFile.is_open()) {
//
//			//read the output from the
//			while (inConfigFile.good()) {
//				string line;
//				if (!getline(inConfigFile, line))
//					break;
//
//				//comment or empty lines
//				if(line[0] == '#' || line.size() == 0){
//					continue;
//				}
//
//				istringstream lineStrem(line);
//				string configName;
//
//				int i = 0;
//				while (lineStrem) {	//for each column (sample)
//					string token;
//
//					if (!getline(lineStrem, token, '='))
//						break;
//
//					if(i == 0){ //read config name
//						configName = token;
//					}else if(i == 1){	//read config content
//						if(configName.compare("outDir") == 0){
//							outDir = token;
//							countArvInConfigFile++;
//							cout << configName << " = " << outDir << endl;
//						}else if(configName.compare("network") == 0){
//							countArvInConfigFile++;
//							networkFilename = token;
//							cout << configName << " = " << networkFilename << endl;
//						}else if(configName.compare("numThreads") == 0){
//							countArvInConfigFile++;
//							numThreads = atoi(token.c_str());
//							cout << configName << " = " << numThreads << endl;
//						}else if(configName.compare("exp") == 0){
//							countArvInConfigFile++;
//							expFilename = token;
//							cout << configName << " = " << expFilename << endl;
//						}else if(configName.compare("snp") == 0){
//							countArvInConfigFile++;
//							snpFilename = token;
//							cout << configName << " = " << snpFilename << endl;
//						}else if(configName.compare("cnv") == 0){
//							countArvInConfigFile++;
//							cnvFilename = token;
//							cout << configName << " = " << cnvFilename << endl;
//						}else if(configName.compare("benchmarkGeneList") == 0){
//							countArvInConfigFile++;
//							benchmarkGeneListFilename = token;
//							cout << configName << " = " << benchmarkGeneListFilename << endl;
//						}else if(configName.compare("dataType") == 0){
//							if(token.compare("ARRAY") == 0){	//use sensitive mode
//								mode = 0;
//							}else if(token.compare("RNA_SEQ") == 0){	//use stringent mode
//								mode = 1;
//							}
//							cout << configName << " = " << token << endl;
//						}
//					}
//
//					i++;
//				}
//
//			}
//			inConfigFile.close();
//		} else {
//			cerr << "Error opening configuration file\n";
//			return 1;
//		}
//
//		if(countArvInConfigFile != 7){
//			cout << "Incorrect configuration file\n";
//			return 1;
//		}else{
//			cout << "Start constructing oncoIMPACT database ..." << endl;
//			database(outDir, networkFilename, expFilename, snpFilename, cnvFilename,
//					benchmarkGeneListFilename, numThreads, mode);
//		}
//	}else{

	/*
	 * DISCOVERY
	 */

		string outDir;

		string networkFilename;
		string expFilename;
		string snpFilename;
		string cnvFilename;
		string benchmarkGeneListFilename;

		string dbPath = "";

		string cancerType;

		int numThreads = 1;
		//0 = sensitive, 1 = stringent
		int mode = 1; //default mode = stringent

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
						if(configName.compare("outDir") == 0){
							outDir = token;
							countArvInConfigFile++;
							cout << configName << " = " << outDir << endl;
//						}else if(configName.compare("network") == 0){
//							countArvInConfigFile++;
//							networkFilename = token;
//							cout << configName << " = " << networkFilename << endl;
						}else if(configName.compare("numThreads") == 0){
							countArvInConfigFile++;
							numThreads = atoi(token.c_str());
							cout << configName << " = " << numThreads << endl;
						}else if(configName.compare("exp") == 0){
							countArvInConfigFile++;
							expFilename = token;
							cout << configName << " = " << expFilename << endl;
						}else if(configName.compare("snp") == 0){
							countArvInConfigFile++;
							snpFilename = token;
							cout << configName << " = " << snpFilename << endl;
						}else if(configName.compare("cnv") == 0){
							countArvInConfigFile++;
							cnvFilename = token;
							cout << configName << " = " << cnvFilename << endl;
//						}else if(configName.compare("benchmarkGeneList") == 0){
//							countArvInConfigFile++;
//							benchmarkGeneListFilename = token;
//							cout << configName << " = " << benchmarkGeneListFilename << endl;
						}else if(configName.compare("dbPath") == 0){
							countArvInConfigFile++;
							dbPath = token;
							cout << configName << " = " << dbPath << endl;
						}else if(configName.compare("cancerType") == 0){
							countArvInConfigFile++;
							cancerType = token;
							cout << configName << " = " << cancerType << endl;
						}else if(configName.compare("dataType") == 0){
							if(token.compare("ARRAY") == 0){	//use sensitive mode
								mode = 0;
							}else if(token.compare("RNA_SEQ") == 0){
								mode = 1;
							}
							cout << configName << " = " << endl;
						}
					}

					i++;
				}

			}
			inConfigFile.close();
		} else {
			cerr << "Error opening configuration file\n";
		}

		if(countArvInConfigFile != 7){
			cout << "Incorrect configuration file\n";
		}else{

			networkFilename = dbPath + "/network_FIsInGene_041709.txt";
			benchmarkGeneListFilename = dbPath + "/Census_all_04_06_2015.tsv";

			cout << "Start constructing oncoIMPACT database ..." << endl;
			discovery(outDir, networkFilename, expFilename, snpFilename, cnvFilename,
					benchmarkGeneListFilename, dbPath, numThreads, cancerType);
		}


//	}	//end if-else for discovery mode

	return 0;
}
