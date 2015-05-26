/*
 * parameters.cpp
 *
 *  Created on: 3 Apr, 2015
 *      Author: Nok
 */

#include "../header/parameters.h"
#include "../header/sampling.h"
#include <cmath>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "omp.h"

void findParameters(vector<JSDivergence>* jsDivergences, vector<int>* Ls,
		vector<int>* Ds, vector<double>* Fs, int totalGenes,
		GeneExpression* geneExpression, Mutations* mutations,
		TIntAdjList* network, int numSamples, int numDatasets,
		vector<string>* geneIdToSymbol, map<string, int>* geneSymbolToId, int numThreads) {

	vector<int>* genesEx = geneExpression->genes;
	vector<int>* genesMut = mutations->genes;

	TDoubleMatrix* originalGeneExpressionMatrix = geneExpression->matrix;
	TIntegerMatrix* originalMutationMatrix = mutations->matrix;
	int totalSamples = originalMutationMatrix->at(0).size();

	int numLs = Ls->size();
	int numDs = Ds->size();
	int numFs = Fs->size();

	int totalGenesUpDown = totalGenes * 2;

	//ignore the choice (of F) in which the median number of deregulated genes is more than half of the gene in the network or <300
	double halfNumberOfGenesInNetwork = totalGenes / 2;
	vector<double> medianNumberOfDeregulatedGenes(Fs->size());
	for (unsigned i = 0; i < Fs->size(); ++i) {
		double F = Fs->at(i);
		medianNumberOfDeregulatedGenes[i] = getMedianNumberOfDeregulatedGenes(
				originalGeneExpressionMatrix, F);
	}

	/*
	 * Resampling numSamples samples for every round
	 */

	clock_t begin_time = clock();

	//gene expression submatrix
	vector<TDoubleMatrix>* subGeneExpressionMatrix = new vector<TDoubleMatrix>(numDatasets);
	vector<GeneExpression>* subGeneExpression = new vector<GeneExpression>(numDatasets);
	//mutation submatrix
	vector<TIntegerMatrix>* subMutationMatrix = new vector<TIntegerMatrix>(numDatasets);
	vector<Mutations>* subMutations = new vector<Mutations>(numDatasets);

	//for finding deregulated genes of all permutations and all possible F
	vector< vector< vector<bool> > > isDeregulatedGensAll(numFs, vector< vector <bool> >(numDatasets, vector<bool>(totalGenesUpDown, false)));

	//SUB
	//create 100 datasets containing 50 samples
	vector< vector<int> > sampleIdsOfDatasets(numDatasets, vector<int>(numSamples, -1));
	#pragma omp parallel for
	for (int i = 0; i < numDatasets; ++i) {

		//list of samples id to be used for tuning the parameters
		vector<int> rrank(totalSamples);
		createPermutation(&rrank);	//return a permutation of [0, totalSamples-1]

		//save samples ids for each dataset
		for (int ri = 0; ri < numSamples; ++ri) {
			sampleIdsOfDatasets[i][ri] = rrank[ri];
		}

		//create the gene expression matrix for sub-sample
		subGeneExpression->at(i).genes = genesEx;	// the same set of genes as the original gene expression matrix
		subGeneExpression->at(i).matrix = &subGeneExpressionMatrix->at(i);	//subset of samples
		randomlyChooseSamplesDouble(originalGeneExpressionMatrix,
				&subGeneExpressionMatrix->at(i), &rrank, numSamples);

		//create the gene expression matrix for sub-sample
		subMutations->at(i).genes = genesMut;	// the same set of genes as the combined mutation matrix
		subMutations->at(i).matrix = &subMutationMatrix->at(i);
		randomlyChooseSamplesInt(originalMutationMatrix, &subMutationMatrix->at(i),
				&rrank, numSamples);

		for (int fi = 0; fi < numFs; ++fi) {
			countNumberOfDeregulatedGenes(&subGeneExpressionMatrix->at(i), Fs->at(fi), &isDeregulatedGensAll[fi][i], genesEx);
		}
	}

//	//FULL
//	for (int i = 0; i < numDatasets; ++i) {
//		for (int fi = 0; fi < numFs; ++fi) {
//			countNumberOfDeregulatedGenes(originalGeneExpressionMatrix, Fs->at(fi), &isDeregulatedGensAll[fi][i], genesEx);
//		}
//	}

	cout << "\tDONE generating 100 datasets (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)	<< " sec)\n";

	int count = 0;	//count number of combinations
	//for each combination of parameters
	for (int di = 0; di < numDs; ++di) {
		for (int fi = 0; fi < numFs; ++fi) {

			int D = Ds->at(di);
			double F = Fs->at(fi);

			for (int li = 0; li < numLs; ++li) {
				cout << "\tcurrent parameters (L, D, F) is " << Ls->at(li) << ", " << D << ", " << F << "... \n";
			}

			begin_time = clock();

			//save values
			//for multiple L
			for (int li = 0; li < numLs; ++li) {
				jsDivergences->at(count).L = Ls->at(li);
				jsDivergences->at(count).D = D;
				jsDivergences->at(count).F = F;
				count++;
			}

			count -= numLs;	//reset the index of js divergence vector

			//Ignore the choice (of F) in which the median number of deregulated genes is more than half of the gene in the network or <300
			if (medianNumberOfDeregulatedGenes[fi] > halfNumberOfGenesInNetwork
					or medianNumberOfDeregulatedGenes[fi] < 300) {
				for (int li = 0; li < numLs; ++li) {
					jsDivergences->at(count).divergence = -1;
				}
				count -= numLs;	//reset the index of js divergence vector
				break;
			}


			/*
			 * Find the explained genes for all REAL samples
			 */

			//for multiple Ls
			vector< vector< vector<bool> > > isExplainedGenesInRealUpDownAllLs(numLs, vector< vector<bool> >(totalSamples, vector<bool>(totalGenesUpDown, false)));

			for (int sampleId = 0; sampleId < totalSamples; sampleId++) {

				vector<double> sampleGeneExpression(totalGenes); //expression of all genes in the network
				getGeneExpressionFromSampleId(originalGeneExpressionMatrix,
						genesEx, &sampleGeneExpression, sampleId);

				vector<int> mutatedGeneIds; // to store gene id of mutated genes
				getMutatedGeneIdsFromSampleId(mutations,
						&mutatedGeneIds, sampleId, genesMut);


				//For multiple L values
				vector< vector<bool> > isExplainedGenesUpDownForASampleAllLs(numLs, vector<bool>(totalGenesUpDown, false));	//collect all explained gene of this sample
				vector< vector<bool> > isExplaninedGeneUpDownForAMutatedGeneAllLs(numLs, vector<bool>(totalGenesUpDown, false));
				int numMutatedGenes = mutatedGeneIds.size();
				for (int mi = 0; mi < numMutatedGenes; ++mi) {	// for each mutated genes
					BFSforExplainedGenesUpDownIncludingMutatedGeneAllLengths(network, mutatedGeneIds[mi], Ls, D, F,
							&isExplaninedGeneUpDownForAMutatedGeneAllLs, &sampleGeneExpression, -1, geneIdToSymbol, geneSymbolToId);
					for (int li = 0; li < numLs; ++li) {
						for (int ei = 0; ei < totalGenesUpDown; ++ei) {
							if(isExplaninedGeneUpDownForAMutatedGeneAllLs[li][ei]){
								isExplainedGenesUpDownForASampleAllLs[li][ei] = true;
							}
						}
					}
				}

				//save the explained genes
				//for multiple Ls
				for (int li = 0; li < numLs; ++li) {
					for (int ei = 0; ei < totalGenesUpDown; ++ei) {//differentiate the up and down regulated genes
						isExplainedGenesInRealUpDownAllLs[li][sampleId][ei] = isExplainedGenesUpDownForASampleAllLs[li][ei];
					}
				}
			}	//end for loop of all real samples

			/*
			 * Crating distribution for both REAL and RANDOM
			 */

			//for multiple Ls
			vector< vector< vector<int> > > realDistributionAllLs(numLs, vector< vector<int> >(numDatasets, vector<int>(totalGenesUpDown, 0)));
			vector< vector< vector<int> > > randomDistributionAllLs(numLs, vector< vector<int> >(numDatasets, vector<int>(totalGenesUpDown, 0)));

			//100 iterations to generate the frequency distribution and compute JS divergence
			#pragma omp parallel for
			for (int i = 0; i < numDatasets; ++i) {

				/*
				 * find explained genes for REAL samples (without gene label permutation)
				 */

				//SUB
				//count the number of samples that each gene (up and down) is explained
				//for multiple Ls
				vector< vector<int> > realDistributionLs(numLs, vector<int>(totalGenesUpDown, 0)); //differentiate the up and down regulated genes

				for (int li = 0; li < numLs; ++li) {

					for (int si = 0; si < numSamples; si++) {

						int sampleId = sampleIdsOfDatasets[i][si];

						//update real distribution
						for (int ei = 0; ei < totalGenesUpDown; ++ei) {//differentiate the up and down regulated genes
							if (isExplainedGenesInRealUpDownAllLs[li][si][ei]) {
								realDistributionLs[li][ei] = realDistributionLs[li][ei] + 1;
							}
						}
					}	//end for loop of samples

					//save the real distribution of current dataset i
					for (int ei = 0; ei < totalGenesUpDown; ++ei) {
						realDistributionAllLs[li][i][ei] = realDistributionLs[li][ei];
					}
				}


//				//FULL
//				//count the number of samples that each gene (up and down) is explained
//				//for multiple Ls
//				vector< vector<int> > realDistributionLs(numLs, vector<int>(totalGenesUpDown, 0)); //differentiate the up and down regulated genes
//
//				for (int li = 0; li < numLs; ++li) {
//
//					for (int si = 0; si < totalSamples; si++) {
//
//						//update real distribution
//						for (int ei = 0; ei < totalGenesUpDown; ++ei) {//differentiate the up and down regulated genes
//							if (isExplainedGenesInRealUpDownAllLs[li][si][ei]) {
//								realDistributionLs[li][ei] = realDistributionLs[li][ei] + 1;
//							}
//						}
//					}	//end for loop of samples
//
//					//save the real distribution
//					for (int ei = 0; ei < totalGenesUpDown; ++ei) {
//						realDistributionAllLs[li][i][ei] = realDistributionLs[li][ei];
//					}
//				}

				/*
				 * find explained genes for RANDOM sub-sample (with gene label permutation)
				 */

				//TODO The permutation will be the same for all Ls, but different for diferent D and F
				//repermute the gene label of random a dataset
				//Create gene label permutation for both gene expression and mutation matrix
				//1. gene expression
				vector<int> permutedGeneLabelsEx;
				permuteGeneLabels(genesEx, &permutedGeneLabelsEx);

				//2. mutation
				vector<int> permutedGeneLabelsMut;
				//permuteGeneLabels(genesMut, &permutedGeneLabelsMut); [Use only gene in the mutation matrix]
				//use all the gene is the network
				permutedGeneLabelsUsingAllGeneInNetwork(genesMut, &permutedGeneLabelsMut, totalGenes);

				//count the number of samples that each gene (up and down) is explained
				//for multiple Ls
				vector< vector<int> > randomDistributionLs(numLs, vector<int>(totalGenesUpDown, 0)); //differentiate the up and down regulated genes

				//SUB
				//for multiple Ls
				for (int si = 0; si < numSamples; si++) {

					vector<double> sampleGeneExpression(totalGenes); // of all genes
					getGeneExpressionFromSampleId(&subGeneExpressionMatrix->at(i),
							&permutedGeneLabelsEx, &sampleGeneExpression,
							si);

					vector<int> mutatedGeneIds; // to store gene id of mutated genes
					getMutatedGeneIdsFromSampleId(&subMutations->at(i),
							&mutatedGeneIds, si,
							&permutedGeneLabelsMut);

					//for multiple L
					vector< vector<bool> > isExplainedGenesUpDownForASampleAllLs(numLs, vector<bool>(totalGenesUpDown, false));	//collect all explained gene of this sample
					vector< vector<bool> > isExplaninedGeneUpDownForAMutatedGeneAllLs(numLs, vector<bool>(totalGenesUpDown, false));
					//for multiple Ls
					int numMutatedGenes = mutatedGeneIds.size();
					for (int mi = 0; mi < numMutatedGenes; ++mi) {// for each mutated genes
						BFSforExplainedGenesUpDownIncludingMutatedGeneAllLengths(
								network, mutatedGeneIds[mi], Ls, D, F,
								&isExplaninedGeneUpDownForAMutatedGeneAllLs,
								&sampleGeneExpression, -1, geneIdToSymbol,
								geneSymbolToId);

						for (int li = 0; li < numLs; ++li) {
							for (int ei = 0; ei < totalGenesUpDown; ++ei) {
								if(isExplaninedGeneUpDownForAMutatedGeneAllLs[li][ei]){
									isExplainedGenesUpDownForASampleAllLs[li][ei] = true;
								}
							}
						}

					}

					//update random distribution
					for (int li = 0; li < numLs; ++li) {
						for (int j = 0; j < totalGenesUpDown; ++j) {//differentiate the up and down regulated genes
							if (isExplainedGenesUpDownForASampleAllLs[li][j]) {
								randomDistributionLs[li][j]++;
							}
						}
					}

				}	//end for loop of samples

//				//FULL
//				//for multiple Ls
//				for (int si = 0; si < totalSamples; si++) {
//
//					vector<double> sampleGeneExpression(totalGenes); // of all genes
//					getGeneExpressionFromSampleId(originalGeneExpressionMatrix,
//							&permutedGeneLabelsEx, &sampleGeneExpression,
//							si);
//
//					vector<int> mutatedGeneIds; // to store gene id of mutated genes
//					getMutatedGeneIdsFromSampleId(mutations,
//							&mutatedGeneIds, si,
//							&permutedGeneLabelsMut);
//
//					//for multiple L
//					vector< vector<bool> > isExplainedGenesUpDownForASampleAllLs(numLs, vector<bool>(totalGenesUpDown, false));	//collect all explained gene of this sample
//					vector< vector<bool> > isExplaninedGeneUpDownForAMutatedGeneAllLs(numLs, vector<bool>(totalGenesUpDown, false));
//					//for multiple Ls
//					int numMutatedGenes = mutatedGeneIds.size();
//					for (int mi = 0; mi < numMutatedGenes; ++mi) {// for each mutated genes
//						BFSforExplainedGenesUpDownIncludingMutatedGeneAllLengths(
//								network, mutatedGeneIds[mi], Ls, D, F,
//								&isExplaninedGeneUpDownForAMutatedGeneAllLs,
//								&sampleGeneExpression, -1, geneIdToSymbol,
//								geneSymbolToId);
//
//						for (int li = 0; li < numLs; ++li) {
//							for (int ei = 0; ei < totalGenesUpDown; ++ei) {
//								if(isExplaninedGeneUpDownForAMutatedGeneAllLs[li][ei]){
//									isExplainedGenesUpDownForASampleAllLs[li][ei] = true;
//								}
//							}
//						}
//
//					}
//
//					//update random distribution
//					for (int li = 0; li < numLs; ++li) {
//						for (int j = 0; j < totalGenesUpDown; ++j) {//differentiate the up and down regulated genes
//							if (isExplainedGenesUpDownForASampleAllLs[li][j]) {
//								randomDistributionLs[li][j]++;
//							}
//						}
//					}
//
//				}	//end for loop of samples


				//save the random distribution of current dataset i
				for (int li = 0; li < numLs; ++li) {
					for (int ei = 0; ei < totalGenesUpDown; ++ei) {
						randomDistributionAllLs[li][i][ei] = randomDistributionLs[li][ei];
					}
				}

			}	//end for loop of 100 round

			cout << "\t\tDONE finding explained gene for all random samples (" << (float(clock() - begin_time) / CLOCKS_PER_SEC)	<< " sec)\n";

			/*
			 * calculated JS divergence
			 */

			//SUB
			//Note: sometimes the divergence is not calculated because of the constraint of the number of deregulated genes
			//for multiple Ls
			for (int li = 0; li < numLs; ++li) {
				jsDivergences->at(count).divergence = calculateJSDivergence(&realDistributionAllLs[li], &randomDistributionAllLs[li],
						numSamples, &isDeregulatedGensAll[fi], Ls->at(li), D, F);

				cout << "\t\tdivergence = " << jsDivergences->at(count).divergence << endl;

				count++;
			}

//			//FULL
//			//for multiple Ls
//			for (int li = 0; li < numLs; ++li) {
//				jsDivergences->at(count).divergence = calculateJSDivergence(&realDistributionAllLs[li], &randomDistributionAllLs[li],
//						totalSamples, &isDeregulatedGensAll[fi], Ls->at(li), D, F);
//
//				cout << "\t\tdivergence = " << jsDivergences->at(count).divergence << endl;
//
//				count++;
//			}

		} //end for loop for Fs
	} //end for loop for Ds

	//gene expression submatrix
	delete subGeneExpressionMatrix;
	delete subGeneExpression;
	//mutation submatrix
	delete subMutationMatrix;
	delete subMutations;


}

double calculateJSDivergence(const vector<vector<int> >* realDistributionAll,
		const vector<vector<int> >* randomDistributionAll, int numSamples, //int sumOfNumDeregulatedGenes,
		vector< vector<bool> >* isDeregulatedGensAll, int L, int D, double F) {
	int round = randomDistributionAll->size();
	int totalGenesUpDown = randomDistributionAll->at(0).size();

	//create frequency distribution (x: #number of samples y: frequency) for both real and random samples
	vector<int> randomFrequencyDistribution(numSamples + 1, 0);	//+ 1 for count genes that are explained in all samples
	vector<int> realFrequencyDistribution(numSamples + 1, 0);

	int sumOfNumDeregulatedGenes = 0;

	// for each round
	for (int i = 0; i < round; ++i) {
		//for each genes, get the frequency
		for (int j = 0; j < totalGenesUpDown; ++j) {
			//only count the frequency for the deregulated genes
			if(isDeregulatedGensAll->at(i)[j]){
				int frequencyOfAGeneRandom = randomDistributionAll->at(i)[j];
				int frequencyOfAGeneReal = realDistributionAll->at(i)[j];
				//add the frequency to distribution
				randomFrequencyDistribution[frequencyOfAGeneRandom]++;
				realFrequencyDistribution[frequencyOfAGeneReal]++;
				sumOfNumDeregulatedGenes++;
			}
		}
	}

	//2. compute divergence
	//from perl code
	//$js += $P_v * log( 2 * $P_v / ( $P_v + $Q_v ) ) / $log_2 if($P_v != 0);   # P Log ( P /(P+Q)/2 )
	//$js += $Q_v * log( 2 * $Q_v / ( $P_v + $Q_v ) ) / $log_2 if($Q_v != 0);   # Q Log ( Q /(P+Q)/2 )
	double log2 = log(2);

	//TODO remove this
//	double sumPropP = 0.0;
//	double sumPropQ = 0.0;
//	string filename = "output/test_param/" + intToStr(L) + "_" + intToStr(D) + "_" + doubleToStr(F,1) + ".dat";
//	vector<string> outStr;

	double jsDivergence = 0;
	for (int i = 0; i < numSamples + 1; ++i) {	//+ 1 for count genes that are explained in all samples
		double pi = 1.0 * realFrequencyDistribution[i] / sumOfNumDeregulatedGenes;	// *2 because now we are considering up and down
//		sumPropP += pi;
		double qi = 1.0 * randomFrequencyDistribution[i] / sumOfNumDeregulatedGenes;
//		sumPropQ += qi;
		if (qi > 0 or pi > 0) {
			if (pi != 0) {	//real frequency is zero
				jsDivergence += pi * log(2 * pi / (pi + qi) ) / log2;
			}
			if (qi != 0) {	//random frequency is zero
				jsDivergence += qi * log(2 * qi / (pi + qi) ) / log2;
			}
		}
//		outStr.push_back(intToStr(realFrequencyDistribution[i]) + "\t" + intToStr(randomFrequencyDistribution[i]));
	}

	//TODO remove this
//	outStr.push_back(intToStr(sumOfNumDeregulatedGenes));	//add the sum of number of deregulated genes to the end of the file
//	writeStrVector(filename.c_str(), &outStr);

	//TODO remove this
//	cout << "\t\tsum pi = " << sumPropP << " qi = " << sumPropQ << endl;
//	cout << "\t\t# deregulated genes for 100 rounds " << sumOfNumDeregulatedGenes << endl;

	return jsDivergence / 2;
}

double getMedianNumberOfDeregulatedGenes(TDoubleMatrix* geneExpressionMatrix,
		double F) {
	int numSamples = geneExpressionMatrix->at(0).size();
	int numGenes = geneExpressionMatrix->size();

	//for each sample i
	vector<int> counts;
	for (int i = 0; i < numSamples; ++i) {
		//count deregulated genes
		int count = 0;
		//for each gene j
		for (int j = 0; j < numGenes; ++j) {
			if (fabs(geneExpressionMatrix->at(j)[i]) >= F) {
				count++;
			}
		}
		counts.push_back(count);
	}

	//find the median number of deregulated genes across samples
	sort(counts.begin(), counts.end());
	double median;
	if (counts.size() % 2 == 0) {
		median = (counts[counts.size() / 2 - 1] + counts[counts.size() / 2])
				/ 2.0;
	} else {
		median = counts[counts.size() / 2];
	}

	return median;
}

void findMaximumJsDivergence(vector<JSDivergence>* jsDivergences,
		JSDivergence* maxJs) {
	int total = jsDivergences->size();
	double max = 0;
	for (int i = 0; i < total; ++i) {
		if (jsDivergences->at(i).divergence > max) {
			maxJs->D = jsDivergences->at(i).D;
			maxJs->L = jsDivergences->at(i).L;
			maxJs->F = jsDivergences->at(i).F;
			maxJs->divergence = jsDivergences->at(i).divergence;
			max = maxJs->divergence;
		}
	}
}

int countNumberOfDeregulatedGenes(TDoubleMatrix* geneExpressionMatrix, double F, vector<bool>* isDeregulatedGens, vector<int>* geneEx){
	int count = 0;
	int totalExGenes = geneExpressionMatrix->size();
	int numSamples = geneExpressionMatrix->at(0).size();
	int totalGenes = isDeregulatedGens->size() / 2;

	for (int i = 0; i < totalExGenes; ++i) {
		bool isUpRegurated = false;
		bool isDownRegurated = false;
		for (int j = 0; j < numSamples; ++j) {
			if(geneExpressionMatrix->at(i)[j] >= F){
				isUpRegurated = true;
			}else if(geneExpressionMatrix->at(i)[j] <= -F){
				isDownRegurated = true;
			}
		}
		if(isUpRegurated){
			isDeregulatedGens->at(geneEx->at(i)) = true;
			count++;
		}
		if(isDownRegurated){
			isDeregulatedGens->at(geneEx->at(i) + totalGenes) = true;
			count++;
		}
	}
	return count;
}

