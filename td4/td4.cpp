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


struct readStruct{
  uint index;
  string sequence;
};

struct edge{
  uint index;
  string sequence;
};

struct compareEdge{
    bool operator()(const edge& seqL, const edge& seqR){
        return seqL.sequence <seqR.sequence;
    }
};

edge nPrefix(uint n, const edge& seq){
	return {seq.index, seq.sequence.substr(0, n)};
}


edge nSuffix(uint n, const edge& seq){
	return {seq.index, seq.sequence.substr(seq.sequence.size()-n, n)};
}



vector<edge> removeNotSingle(const vector<edge>& vect, uint k){
	vector<edge> vectResult;
	uint i(0);
	while (i < vect.size()){
		if (i == 0){
			if (vect[i].sequence != vect[i+1].sequence){
				vectResult.push_back(vect[i]);
			}
		} else if (i == vect.size()-1){
			if (vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			}
		} else {
			if (vect[i].sequence != vect[i+1].sequence and vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			}
		}
		++i;
	}
	return vectResult;
}


string compaction(const edge& seq1, const edge& seq2, uint k){	
	edge beg1(nPrefix(k, seq1));
	edge beg2(nPrefix(k, seq2));
	edge end1(nSuffix(k, seq1));
	edge end2 = nSuffix(k, seq2);
	
	if (end1.sequence == beg2.sequence){
		return seq1.sequence + seq2.sequence.substr(k);
	} else { //if (end2 == beg1)
		return seq2.sequence + seq1.sequence.substr(k);
	}
}


void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k){
	if (not vec[indexreadStruct1].sequence.empty()){
		if (not vec[indexreadStruct2].sequence.empty()){
			string c = compaction({vec[indexreadStruct1].index, vec[indexreadStruct1].sequence}, {vec[indexreadStruct2].index, vec[indexreadStruct2].sequence}, k);
			if (not c.empty()){
				vec[indexreadStruct1] = {vec[indexreadStruct1].index, c};
				vec[indexreadStruct2].index = vec[indexreadStruct1].index;
				vec[indexreadStruct2].sequence = "";
			}
		} else {
			compactInVector(vec, indexreadStruct1, vec[indexreadStruct2].index, k); //  each time a sequence is empty, the index leads to the sequence it's been compacted in-> recursive call until we find the sequence
		}
	} else {
		if (not vec[indexreadStruct2].sequence.empty()){
			compactInVector(vec, vec[indexreadStruct1].index, indexreadStruct2, k);
		} else {
			compactInVector(vec, vec[indexreadStruct1].index, vec[indexreadStruct2].index, k);
		}
	}
}


void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k){
	sort(left.begin(), left.end(), compareEdge());
	sort(right.begin(), right.end(), compareEdge());
	vector<edge> leftSingles;
	vector<edge> rightSingles;
	if (left.size()>1){
		leftSingles = removeNotSingle(left, k);
	} else {
		leftSingles = left;
	}
	if (right.size()>1){
		rightSingles = removeNotSingle(right, k);
	} else {
		rightSingles = right;
	}
	//~ cout << leftSingles.size() << " " << rightSingles.size() << endl;
	uint indexL(0),indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			if (leftSingles[indexL].index != rightSingles[indexR].index){
				compactInVector(readStructsVec, leftSingles[indexL].index, rightSingles[indexR].index, k);
			}
			++indexL;
			++indexR;
		} else {
			if (leftSingles[indexL].sequence < rightSingles[indexR].sequence){
				++indexL;
			} else {
				++indexR;
			}
		}
	}
}


/* fill vectors of prefixes k-mers coming from prefixes of readStructs */
void fillPrefVector(vector <edge>& vecLeft, const vector<readStruct>& vec, uint k){
	for (uint i(0); i < vec.size(); ++i){
		edge prefix = nPrefix(k, {vec[i].index, vec[i].sequence});
		vecLeft.push_back(prefix);
	}
}



void fillSuffVector(vector <edge>& vecRight, const vector<readStruct>& vec, uint k){
	for (uint i(0); i < vec.size(); ++i){
		edge suffix = nSuffix(k, {vec[i].index, vec[i].sequence});
		vecRight.push_back(suffix);
	}
}

/* remove duplicates in reads*/
void cleanDuplicatesInreadStructs(const vector <readStruct>& vec, vector <readStruct>& newVec){
	uint i(0);
	string previousSeq("");
	while(i<vec.size()){
		string temp(vec[i].sequence);
		if (temp != previousSeq){
			newVec.push_back(vec[i]);
		}
		++i;
	}
}


void setreadStructsIndex(vector <readStruct>& vec){
	for (uint i(0); i<vec.size(); ++i){
		if (not vec[i].sequence.empty()){
			vec[i].index = i;
		}
	}
}


void readDataset(ifstream& readFile, vector <readStruct>& vec){
	string sequence;
	while (not readFile.eof()){
        getline(readFile, sequence);
		getline(readFile, sequence);
		if (not sequence.empty()){
			vec.push_back({0, sequence});
		}
	}
}

void writeOutput(ofstream& outFile, const vector<readStruct>& vec){
	uint count(0);
	for (uint seq(0); seq < vec.size(); ++seq){
		if (not vec[seq].sequence.empty()){
			outFile << "> read_" << count << endl;
			outFile << vec[seq].sequence << endl;
			++count;
		}
	}
}


void runUnitiger(char** argv){
	string fileName = argv[1];
	uint k = stoi(argv[2]);
	string outFileName = argv[3];
	ifstream readFile(fileName);
	ofstream outFile(outFileName);
	vector <readStruct> tempReadSet;
	vector <readStruct> readSet;
	readDataset(readFile, tempReadSet);
	cleanDuplicatesInreadStructs(tempReadSet, readSet);
	setreadStructsIndex(readSet);
	vector <edge> left;
	vector <edge> right;
	fillPrefVector(left, readSet, k);
	fillSuffVector(right, readSet, k);
	parseVector(left, right, readSet, k);
	writeOutput(outFile, readSet);
}

int main(int argc, char **argv)
{
	runUnitiger(argv);
	return 0;
}
