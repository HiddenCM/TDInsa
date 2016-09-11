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

using namespace std;

void readDataset(ifstream& readFile, vector<string>& readSet);
void writeCleanReads(ofstream& outFile, const vector<string>& cleanReadSet);

void pushKmerInMap(const vector<string>& readSet, unordered_map <string, uint>& occurrenceKmers, uint k);
void selectSolidKmers(const unordered_map <string, uint>& occurrenceKmers, unordered_set <string>& solidKmers);
void removeSpuriousReads(const vector<string>& readSet, vector<string>& cleanReadSet, const unordered_set <string>& solidKmers, uint k);


bool checkKmerInBloomFilter(const string& kmer, const vector <bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF);
void pushKmerInBloomFilter(const vector<string>& readSet, vector <bool>& BloomFilterSolid, vector <bool>& BloomFilterCheck, uint k, int nbHashFunctions, uint64_t sizeBF,  uint64_t sizeBF2);
void removeSpuriousReadsBloomFilter(const vector<string>& readSet, vector<string>& cleanReadSet, uint k, const vector<bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF);


void runCorrectorHashTable(char** argv);
void runCorrectorBloomFilter(char** argv);


void runGetReads(char** argv);

