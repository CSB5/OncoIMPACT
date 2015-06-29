/*
 * calculator.cpp
 *
 *  Created on: Apr 27, 2015
 *      Author: nok
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include "../header/annotator_calculator.h"

int countOverlap(vector<string>* currentModuleMembers, vector<string>* currentGeneSetMembers){
	vector<string> overlap;

	sort(currentModuleMembers->begin(), currentModuleMembers->end());
	sort(currentGeneSetMembers->begin(), currentGeneSetMembers->end());
	set_intersection(currentGeneSetMembers->begin(), currentGeneSetMembers->end(),
			currentModuleMembers->begin(), currentModuleMembers->end(), back_inserter(overlap));

//	for(vector<string>::iterator it = overlap.begin(); it != overlap.end(); it++){
//		cout << *it << endl;
//	}

	return overlap.size();
}

void calculateHypergeometric(GeneSetPair* pair, int N, int numGeneSets){
//	[from pl script]
//	$a = min($m,$n);
//	$pvalue=0;
//	#From the overlap k to the minum a of(gene_set, gene_list)
//	#n gene list
//	#m gene set
//	#N backgroud
//	for ($j=$k; $j<=$a; $j++) {
//		$pvalue += hypergeom($m,($N-$m),$n,$j);
//	}

	int a = pair->sizeOfGeneList < pair->sizeOfGeneSet ? pair->sizeOfGeneList : pair->sizeOfGeneSet;
	int n = pair->sizeOfGeneList;
	int m = pair->sizeOfGeneSet;

	pair->pValue = 0;
	for (int j = pair->numOverlap; j <= a; ++j) {
		pair->pValue += hypergeom(m, (N-m), n, j);
	}

//	cout << "p-value = " << pair->pValue << " after correction = " << pair->pValue * 1.0 * numGeneSets << endl;
	pair->pValue *= 1.0 * numGeneSets;


}

long double hypergeom(int n, int m, int N, int i){
//	sub hypergeom {
//	   # There are m "bad" and n "good" balls in an urn.
//	   # Pick N of them. The probability of i or more successful selections
//
//	   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
//	   my ($n, $m, $N, $i) = @_;
//
//	   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
//	   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact(+$N-$i)+logfact($m+$n);
//	   return exp($loghyp1 - $loghyp2);
//	}

	long double loghyp1 = logfact(m)+logfact(n)+logfact(N)+logfact(m+n-N);
	long double loghyp2 = logfact(i)+logfact(n-i)+logfact(m+i-N)+logfact(N-i)+logfact(m+n);

	return exp(loghyp1 - loghyp2);
}

long double logfact(int x){
//	sub logfact {
//	   return gammln(shift(@_) + 1.0);
//	}

	return gammln(x + 1.0);
}

long double cof[] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,0.12086509738661e-2, -0.5395239384953e-5};

long double gammln(int xx){
//	sub gammln {
//	  my $xx = shift;
//	  my @cof = (76.18009172947146, -86.50532032941677,
//	             24.01409824083091, -1.231739572450155,
//	             0.12086509738661e-2, -0.5395239384953e-5);
//	  my $y = my $x = $xx;
//	  my $tmp = $x + 5.5;
//	  $tmp -= ($x + .5) * log($tmp);
//	  my $ser = 1.000000000190015;
//	  for my $j (0..5) {
//	     $ser += $cof[$j]/++$y;
//	  }
//	  -$tmp + log(2.5066282746310005*$ser/$x);
//	}

	long double y = xx;
	long double x = xx;
	long double tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	long double ser = 1.000000000190015;
	for (int j = 0; j <= 5; ++j) {
		ser += cof[j] / ++y;
	}

//	if(log(2.5066282746310005 * ser / x) > numeric_limits<long double>::max()){
//		cerr << "overflow > MAX\n";
//	}
//	if(log(2.5066282746310005 * ser / x) < numeric_limits<long double>::lowest()){
//		cerr << "overflow < MIN " << 2.5066282746310005 * ser / x << endl;cin.get();
//	}

	return -tmp + log(2.5066282746310005 * ser / x);
}

bool sortByPValues(const GeneSetPair& first, const GeneSetPair& second) {
	if (first.pValue < second.pValue) {
		return true;
	} else {
		return false;
	}
}

//string intToStr(int i){
//	stringstream ss;
//	ss << i;
//	return ss.str();
//}
//
//string doubleToStr(double i, int prec){
//	stringstream ss;
//	ss.precision(prec);
//	ss << fixed << i;
//	return ss.str();
//}


