/*
 * 
 * td2_student.cpp
 *
 */


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;


//xorshift
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


// nucleotide to integer
uint64_t nuc2int(char c){
	switch(c){
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return 0;
}


// sequence to integer
uint64_t seq2int(const string& seq){
	uint64_t res(0);
	for(uint64_t i(0); i < seq.size(); ++i){
		res <<=  2;
		res += nuc2int(seq[i]);
	}
	return res;
}




// td3 : mutate a nucleotide
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
	int position(rand()%read.size());
	read[position] = randNuc();
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

