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
#include "utils_etudiants.h"

using namespace std;


// IO functions

// lecture du fichier de séquences et stockage des séquences dans un vecteur de string
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

// ecriture des reads sans erreur dans un fichier de sortie
void writeCleanReads(ofstream& outFile, const vector<string>& cleanReadSet){
	uint count(0);
	for (uint seq(0); seq < cleanReadSet.size(); ++seq){
		outFile << ">read_" << count + 1 << endl;
		outFile << cleanReadSet[seq] << endl;
		++count;
	}
}

/* TD : with hash table */

//compléter cette fonction : mettre chaque k-mer dans une table de hash
void pushKmerInMap(const vector<string>& readSet, unordered_map <string, uint>& occurrenceKmers, uint k){
	for (uint seq(0); seq < readSet.size(); ++seq){
		for (uint position(0); position < readSet[seq].size() - k + 1; ++position){
			// récupérer le k-mer (utiliser substr())
			// l'ajouter à la table de hash occurenceKmers ainsi que le nombre de fois qu'il a été vu
		}
	}
}

// garder les kmers solides (d'occurrence > 1 dans occurrenceKmers) dans la table solidKmers
void selectSolidKmers(const unordered_map <string, uint>& occurrenceKmers, unordered_set <string>& solidKmers){
}

// remplir le vecteur cleanReadSet avec les reads ne contenant que des k-mers non erronés
void removeSpuriousReads(const vector<string>& readSet, vector<string>& cleanReadSet, const unordered_set <string>& solidKmers, uint k){
	// pour chaque read, lire tous ses k-mers. Si le read contient au moins un k-mer n'existe pas dans la table solidKmers, ne pas l'écrire dans le résultat final
}


/* TD : with Bloom Filter */


// vérifier si un k-mer appartient au filtre de bloom
bool checkKmerInBloomFilter(const string& kmer, const vector <bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF){
	// récupérer l'entier 64b qui correspond au kmer grace à une fonction de utils_etudiants.cpp
	// hash ce nombre nbHashFunctions fois grace à la fonction multiHash de utils_etudiants.cpp
	// vérifier que le bloom filter est à true pour chaque case correspondant aux hash, si oui renvoyer que l'on a trouvé le kmer, sinon renvoyer qu'il n'existe pas dans le filtre
}

// mettre les k-mers solides dans le filtre de bloom BloomFilterCheck en utilisant un filtre intermédiaire BloomFilterSolid
void pushKmerInBloomFilter(const vector<string>& readSet, vector <bool>& BloomFilterSolid, vector <bool>& BloomFilterCheck, uint k, int nbHashFunctions, uint64_t sizeBF,  uint64_t sizeBF2){
	// pour tous les reads de  readSet :
	// récupérer chaque kmer, et vérifier avec checkKmerInBloomFilter s'ils sont dans le filtre BloomFilterSolid
	// si oui, alors le k-mer est déjà apparu une fois au moins,, il est donc solide: l'ajouter au filtre BloomFilterCheck
	// sinon, l'ajouter à BloomFilterSolid
		// pour ajouter un k-mer au filtre : le hasher nbHashFunctions de fois et mettre les cases du filtre correspondant au hash à true
}

// remplir le vecteur cleanReadSet avec les reads ne contenant que des k-mers non erronés
void removeSpuriousReadsBloomFilter(const vector<string>& readSet, vector<string>& cleanReadSet, uint k, const vector<bool>& BloomFilter, int nbHashFunctions, uint64_t sizeBF){
}








// TD 2 : partie 1
void runCorrectorHashTable(char** argv){
	string fileName = argv[1];
	uint k = stoi(argv[2]); // size of k-mers (words of size k)
	string outFileName = argv[3];
	ifstream readFile(fileName); // input file with reads
	ofstream outFile(outFileName); // ouput file where we write clean reads
	vector<string> readSet;
	readDataset(readFile, readSet); // reads are stored in the vector of string readSet
	unordered_map<string, uint> occurrenceKmers;
	pushKmerInMap(readSet, occurrenceKmers, k); // 1st function to complete
}


// TD 2 : partie 2
void runCorrectorBloomFilter(char** argv){
}



// simulate reads - usage ./readCorrection referenceFile.fa
void runGetReads(char** argv){
	string refName = argv[1];
	ifstream refFile(refName);
	perfectsReadsFromRef(refFile, 100, 400000);
	ifstream reads("perfect_reads.fa");
	mutateReads(reads, 2);
}
