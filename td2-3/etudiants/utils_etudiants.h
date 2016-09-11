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


// bloom filter hashing functions
uint64_t xorshift1(uint64_t x);
uint64_t xorshift2(uint64_t x);
uint64_t multiHash(uint64_t u, int n);

// conversion functions
uint64_t nuc2int(char c);
uint64_t seq2int(const string& seq);

// reads generation
void perfectsReadsFromRef(ifstream& ref, uint length, uint nbRead);
char randNuc();
string mutate(string& read);
void mutateReads(ifstream& reads, uint mutRate);
