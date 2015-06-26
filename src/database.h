/*
 * database.h
 *
 *  Created on: May 18, 2015
 *      Author: nok
 */

#ifndef DATABASE_H_
#define DATABASE_H_

#include <string>
#include "header/utilities.h"
#include "database.h"

int database(string outDir, string networkFilename, string expFilename, string snpFilename, string cnvFilename,
		string benchmarkGeneListFilename, int numThreads, int mode);

#endif /* DATABASE_H_ */
