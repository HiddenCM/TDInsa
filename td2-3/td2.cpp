/*
 * td2.cpp
 * 
 * camille marchet <camille.marchet@irisa.fr> INSA Rennes 2016
 *
 */


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

using namespace std;



void readDataset(ifstream& readFile, vector<string>& readSet){
	string sequence;
	while (not readFile.eof()){
        getline(readFile, sequence);
		getline(readFile, sequence);
		if (not sequence.empty()){
			readSet.push_back(sequence);
		}
	}
}


void pushKmerInMap(const vector<string>& readSet, unordered_map <string, uint>& occurrenceKmers, uint k){
	for (uint seq(0); seq < readSet.size(); ++seq){
		for (uint position(0); position < readSet[seq].size() - k + 1; ++position){
			// td: hash table
				string kmer(readSet[seq].substr(position, k));
				auto got = occurrenceKmers.find(kmer);
				if (got == occurrenceKmers.end()){
					occurrenceKmers.insert({kmer, 1});
				} else {
					got->second += 1;
				}
			// end td
		}
	}
}

// td : garder les kmers solides, avec table de hash
void selectSolidKmers(const unordered_map <string, uint>& occurrenceKmers, unordered_set <string>& solidKmers){
	for (auto iter(occurrenceKmers.begin()); iter != occurrenceKmers.end(); ++iter){
		if (iter->second > 1){
			solidKmers.insert(iter->first);
		}
	}
}


// td : retirer les reads erronés, avec table de hash
void removeSpuriousReads(const vector<string>& readSet, vector<string>& cleanReadSet, const unordered_set <string>& solidKmers, uint k){
	for (uint seq(0); seq < readSet.size(); ++seq){
		uint spuriousKmers(0);
		for (uint position(0); position < readSet[seq].size() - k + 1; ++position){
			string kmer(readSet[seq].substr(position, k));
			auto got = solidKmers.find(kmer);
			if (got == solidKmers.end()){
				++spuriousKmers;
			}
		}
		if (not spuriousKmers){
			cleanReadSet.push_back(readSet[seq]);
		}
	}
}

void writeCleanReads(ofstream& outFile, const vector<string>& cleanReadSet){
	uint count(0);
	for (uint seq(0); seq < cleanReadSet.size(); ++seq){
		outFile << "> read_" << count << endl;
		outFile << cleanReadSet[seq] << endl;
		++count;
	}
}

// td3 : muter un nucléotide
char randomNucleotide(){
	switch (rand() % 4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}


vector<char> switchNucleotide(char nucleotide){
	
	if (nucleotide == 'A'){
		return {'C','G','T'};
	}
	else if (nucleotide == 'C'){
		return {'A','G','T'};
	} else if (nucleotide == 'G'){
		return {'A','C','T'};
	} else {
		return {'A','C','G'};
	}
}


// td3 : muter un nucleotide du kmer
vector<string> getMutatedKmer(const string& kmer, uint k, uint positionOnKmer){
	char nucl(kmer[positionOnKmer]);
	vector<char> nuclVect = switchNucleotide(nucl);
	vector<string> resultVect;
	for (uint i(0); i < nuclVect.size(); ++i){
		string newKmer(kmer);
		newKmer[positionOnKmer] = nuclVect[i];
		resultVect.push_back(newKmer);
	}
	return resultVect;
}


//td 3: remplacer un  read par sa version corrigée
string replaceRead(const string& read, string& kmer, uint k, uint position){
	string firstPart(read.substr(0, position));
	string lastPart(read.substr(position + k));
	return firstPart + kmer + lastPart;
}


// td 3 : muter les reads erronés pour les corriger, avec table de hash, reprendre la fonction removeSpuriousReads
void correctSpuriousReads(const vector<string>& readSet, vector<string>& cleanReadSet, const unordered_set <string>& solidKmers, uint k){
	for (uint seq(0); seq < readSet.size(); ++seq){
		uint spuriousKmers(0);
		string newRead(readSet[seq]);
		uint position(0);
		while (position < readSet[seq].size() - k + 1){
			string kmer(readSet[seq].substr(position, k));
			auto got = solidKmers.find(kmer);
			if (got == solidKmers.end()){
				bool corrected(false);
				uint count(0);
				while(not corrected and count < kmer.size()){
					for (uint pos(0); pos < kmer.size(); ++pos){
						vector<string> mutatedKmers(getMutatedKmer(kmer, k, pos));
						for (uint i(0); i < mutatedKmers.size(); ++i){
							auto found = solidKmers.find(mutatedKmers[i]);
							if (found != solidKmers.end()){
								corrected = true;
								newRead = replaceRead(readSet[seq], mutatedKmers[i], k, position); //TODO il faut un while ! en corrigeant un kmer on corriger k positions
								position += k - 1;
								break;
							}
						}
						++count;
						if (corrected){
							break;
						}
					}
				}
				if (not corrected){
					++spuriousKmers;
				}
			}
			++position;
		}
		if (not spuriousKmers){
			cleanReadSet.push_back(newRead);
		}
	}
}

void runCorrector(char** argv){
	string fileName = argv[1];
	uint k = stoi(argv[2]);
	string outFileName = argv[3];
	ifstream readFile(fileName);
	ofstream outFile(outFileName);
	vector<string> readSet;
	// td2 : avec tables de hash
	//~ readDataset(readFile, readSet);
	//~ unordered_map<string, uint> occurrenceKmers;
	//~ pushKmerInMap(readSet, occurrenceKmers, k);
	//~ unordered_set<string> solidKmers;
	//~ selectSolidKmers(occurrenceKmers, solidKmers);
	//~ vector<string> cleanReadSet;
	//~ removeSpuriousReads(readSet, cleanReadSet, solidKmers, k);
	//~ writeCleanReads(outFile, cleanReadSet);
	// end td2
    // td3 : avec tables de hash
	readDataset(readFile, readSet);
	unordered_map<string, uint> occurrenceKmers;
	pushKmerInMap(readSet, occurrenceKmers, k);
	unordered_set<string> solidKmers;
	selectSolidKmers(occurrenceKmers, solidKmers);
	vector<string> cleanReadSet;
	correctSpuriousReads(readSet, cleanReadSet, solidKmers, k);
	writeCleanReads(outFile, cleanReadSet);
	// end td2
}



int main(int argc, char **argv)
{
	runCorrector(argv);
	return 0;
}

