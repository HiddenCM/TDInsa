/*
 * td2.cpp
 * 
 * camille marchet <camille.marchet@irisa.fr> INSA Rennes 2016
 *
 */


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "td2_etudiants.h"

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 3){
		runGetReads(argv);
	} else {
		runCorrectorHashTable(argv); // 1st part of TD2
		//~ runCorrectorBloomFilter(argv); // 2nd part of TD2
	}
	return 0;
}
