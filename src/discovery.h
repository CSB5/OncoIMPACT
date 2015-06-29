/*
 * discovery.h
 *
 *  Created on: May 18, 2015
 *      Author: nok
 */

#ifndef DISCOVERY_H_
#define DISCOVERY_H_

#include <string>
#include "header/utilities.h"
#include "database.h"

int discovery(string outDir, string networkFilename, string expFilename, string snpFilename, string cnvFilename,
		string benchmarkGeneListFilename, string dbPath, int numThreads, string cancerType);


#endif /* DISCOVERY_H_ */
