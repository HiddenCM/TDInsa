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


// 512477 21-mers in lambdavirus_10X_err.fasta (400000 reads)
// 86100 solid 21-mers 


// correction avec hash:
//~ real	0m13.867s
//~ user	0m12.928s
//~ sys	0m0.939s
// mem: 134204 k (/usr/bin/time)


// correction avec BF:
//~ real	0m11.847s
//~ user	0m10.980s
//~ sys	0m0.863s
// mem: 18167k

// faire une fonction qui mute des reads parfaits
// faire rmq que c'est dur à valider !
// on peut aussi perdre un evt rare

// xorshift
uint64_t xorshift1(uint64_t x){
	x ^= x >> 12;
	x ^= x << 25;
	x ^= x >> 27;
	return x * UINT64_C(2685821657736338717);
}


// xorshift 2
uint64_t xorshift2(uint64_t x){
	x ^= (x << 21);
	x ^= (x >> 35);
	x ^= (x << 4);
	return x * UINT64_C(2685821657736338717);
}


// hash the integer u with the hash function number n
uint64_t multiHash(uint64_t u, int n){
	return (n + 1) * xorshift1(u) + xorshift2(u); 
}



uint64_t nuc2int(char c){
	switch(c){
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return 0;
}

uint64_t seq2int(const string& seq){
	uint64_t res(0);
	for(uint64_t i(0); i < seq.size(); ++i){
		res <<=  2;
		res += nuc2int(seq[i]);
	}
	return res;
}

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

// td : kmer in hash map
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


//TODO: estimer avant le nb de kmers pour avoir la taille du filtre de bloom en dur + le nb de fonctions de hash

//td : check kmer existence in bf
bool checkKmerInBloomFilter(const string& kmer, const vector <bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF){
	bool result(true);
	uint64_t integer(seq2int(kmer));  // get 64b integer corresponding to the string
	for (int n(0); n < nbHashFunctions; ++n){ // gatou : en mettre une seule 
		uint64_t hash(multiHash(integer, n) % sizeBF);
		if (BloomFilter[hash] == 1){
			result *= true;
		} else {
			result *= false;
		}
	}
	return result;
}



// td : kmer in bf
void pushKmerInBloomFilter(const vector<string>& readSet, vector <bool>& BloomFilterSolid, vector <bool>& BloomFilterCheck, uint k, int nbHashFunctions, uint64_t sizeBF,  uint64_t sizeBF2){
	for (uint seq(0); seq < readSet.size(); ++seq){
		for (uint position(0); position < readSet[seq].size() - k + 1; ++position){
			string kmer(readSet[seq].substr(position, k));
			uint64_t integer(seq2int(kmer));  // get 64b integer corresponding to the string
			bool solid(checkKmerInBloomFilter(kmer, BloomFilterSolid, nbHashFunctions, sizeBF)); // check if the kmer already exists in the first bf =>  solid
			
			for (int n(0); n < nbHashFunctions; ++n){ // gatou : en mettre une seule 
					if (solid){
						uint64_t hash(multiHash(integer, n) % sizeBF2);
						BloomFilterCheck[hash] = true;
					} else {
						uint64_t hash(multiHash(integer, n) % sizeBF);
						BloomFilterSolid[hash] = true;
				}
			}
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


// td : retirer les reads erronés, avec bf
void removeSpuriousReadsBloomFilter(const vector<string>& readSet, vector<string>& cleanReadSet, uint k, const vector<bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF){
	for (uint seq(0); seq < readSet.size(); ++seq){
		uint spuriousKmers(0);
		for (uint position(0); position < readSet[seq].size() - k + 1; ++position){
			string kmer(readSet[seq].substr(position, k));
			bool inBF(checkKmerInBloomFilter(kmer, BloomFilter, nbHashFunctions, sizeBF));
			if (not inBF){
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
		outFile << ">read_" << count + 1 << endl;
		outFile << cleanReadSet[seq] << endl;
		++count;
	}
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
								newRead = replaceRead(readSet[seq], mutatedKmers[i], k, position); 
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


// td 3 : BF VERSION muter les reads erronés pour les corriger, avec table de hash, reprendre la fonction removeSpuriousReads
void correctSpuriousReadsBloomFilter(const vector<string>& readSet, vector<string>& cleanReadSet, const vector<bool>& bloomFilter, uint k, int nbHashFunctions, uint64_t sizeBF){
	for (uint seq(0); seq < readSet.size(); ++seq){
		uint spuriousKmers(0);
		string newRead(readSet[seq]);
		uint position(0);
		while (position < readSet[seq].size() - k + 1){
			string kmer(readSet[seq].substr(position, k));
			//~ auto got = solidKmers.find(kmer);
			bool present(checkKmerInBF(kmer, bloomFilter, nbHashFunctions, sizeBF));
			if (not present){
				bool corrected(false);
				uint count(0);
				while(not corrected and count < kmer.size()){
					for (uint pos(0); pos < kmer.size(); ++pos){
						vector<string> mutatedKmers(getMutatedKmer(kmer, k, pos));
						for (uint i(0); i < mutatedKmers.size(); ++i){
							//~ auto found = solidKmers.find(mutatedKmers[i]);
							bool found(checkKmerInBloomFilter(mutatedKmers[i], bloomFilter, nbHashFunctions, sizeBF));
							if (found){
								corrected = true;
								newRead = replaceRead(readSet[seq], mutatedKmers[i], k, position);
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
	/* td2 : avec tables de hash */
	//~ readDataset(readFile, readSet);
	//~ unordered_map<string, uint> occurrenceKmers;
	//~ pushKmerInMap(readSet, occurrenceKmers, k);
	
	//~ unordered_set<string> solidKmers;
	//~ selectSolidKmers(occurrenceKmers, solidKmers);
	//~ cout << solidKmers.size() << endl;
	//~ vector<string> cleanReadSet;
	//~ removeSpuriousReads(readSet, cleanReadSet, solidKmers, k);
	//~ writeCleanReads(outFile, cleanReadSet);
	/* end td2 */
	/* td 2 with BF */
	//~ readDataset(readFile, readSet);
	//~ uint64_t sizeBF(500);
	//~ int nbHashFunctions(12);
	//~ vector<bool> bloomFilterSolid(sizeBF, false);
	//~ vector<bool> bloomFilterCheck(sizeBF, false);
	//~ pushKmerInBloomFilter(readSet, bloomFilterSolid, bloomFilterCheck, k, nbHashFunctions, sizeBF);
	//~ vector<string> cleanReadSet;
	//~ removeSpuriousReadsBloomFilter(readSet, cleanReadSet, k, bloomFilterCheck, nbHashFunctions, sizeBF);
	//~ writeCleanReads(outFile, cleanReadSet);
	/* end td2 with BF */
    /* td3 : avec tables de hash */
	//~ readDataset(readFile, readSet);
	//~ unordered_map<string, uint> occurrenceKmers;
	//~ pushKmerInMap(readSet, occurrenceKmers, k);
	//~ unordered_set<string> solidKmers;
	//~ selectSolidKmers(occurrenceKmers, solidKmers);
	//~ vector<string> cleanReadSet;
	//~ correctSpuriousReads(readSet, cleanReadSet, solidKmers, k);
	//~ writeCleanReads(outFile, cleanReadSet);
	/* end td3 hash */
	/* td3 with BF */
	readDataset(readFile, readSet);
	uint64_t sizeBF(5000000);
	uint64_t sizeBF2(860000);
	int nbHashFunctions(12);
	vector<bool> bloomFilterSolid(sizeBF, false);
	vector<bool> bloomFilterCheck(sizeBF2, false);
	pushKmerInBloomFilter(readSet, bloomFilterSolid, bloomFilterCheck, k, nbHashFunctions, sizeBF, sizeBF2);
	vector<string> cleanReadSet;
	correctSpuriousReadsBloomFilter(readSet, cleanReadSet, bloomFilterCheck, k, nbHashFunctions, sizeBF2);
	writeCleanReads(outFile, cleanReadSet);
	/* end td3 with BF */
}







/* READS */

void perfectsReadsFromRef(ifstream& ref, uint length, uint nbRead){
	ofstream out("perfect_reads.fa", ios::trunc);
	string line;
	getline(ref, line);
	getline(ref, line);
	for(uint i(0); i < nbRead; ++i){
		out<<">read_"<< i + 1 << endl;
		out << (line.substr(rand()%(line.size()-length), length)) << endl;
	}
}


char randNuc(){
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



string mutate(string& read){
	//~ for(uint i(0); i < n; ++i){
	int position(rand()%read.size());
	read[position] = randNuc();
		//~ }
	return read;
}


void mutateReads(ifstream& reads, uint mutRate){
	ofstream out("reads_error.fa", ios::trunc);
	string line;
	uint i(0);
	while (not reads.eof()){
		getline(reads, line);
		getline(reads, line);
		out<<">read_"<< i + 1 << endl;
		uint dice(rand() % 10000);
		if(dice <= mutRate){
			out << mutate(line) << endl;
		} else {
			out << line << endl;
		}
		++i;
	}
}

void runGetReads(char** argv){
	string refName = argv[1];
	ifstream refFile(refName);
	perfectsReadsFromRef(refFile, 100, 400000);
	ifstream reads("perfect_reads.fa");
	mutateReads(reads, 2);
}

int main(int argc, char **argv)
{
	if (argc < 3){
		runGetReads(argv);
	} else {
		runCorrector(argv);
	}
	return 0;
}

