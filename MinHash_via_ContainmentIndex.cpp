#include<iostream>
#include<vector>
#include<string>
#include<limits.h>
#include<unordered_map>
#include<fstream>
#include <time.h>  
#include "MurmurHash3.h"
#include "MinHash.h"
#include "ContainmentHash.h"
#define KMER_LENGTH 3
#define SEED 42
#define FALSE_POSITIVE_PROBABILITY 0.00001
#define LONG_READ_FILE_NAME "long_read.txt"
#define HASH_FUNCTIONS 100
using namespace std;

string get_long_read(){
	string line;
	string output = "";
	ifstream myfile (LONG_READ_FILE_NAME);
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{	
			output += line;
		}
		myfile.close();
	}

	else 
		cout << "Unable to open file"; 

	return line;
}

void print_mbits(std::vector<int> mbits){
	for(int i=0;i<mbits.size();i++)
		cout<<mbits[i]<<" ";
	cout<<endl;
}

int main(int argc, char const *argv[])
{
	MinHash *minhash_obj = new MinHash();
	ContainmentHash *containmenthash_obj = new ContainmentHash();
	string input1 = get_long_read();

	uint32_t seed = SEED;
	int kmer_length = KMER_LENGTH;
	int hash_functions = HASH_FUNCTIONS;
	//cout<<"hash_functions"<<"    "<<"actual_jaccardian"<<"    "<<"minhash_jaccardian_distance"<<"    "<<"	diff"<<endl;

	vector<string> input_kmers_1 = minhash_obj->get_kmers(input1,kmer_length);

	vector<vector<uint32_t> > characteristic_matrix = minhash_obj->get_characteristic_matrix(hash_functions,input_kmers_1,seed,kmer_length);
	vector<uint32_t> sketch_vector = minhash_obj->get_sketch_vector(characteristic_matrix);

	minhash_obj->save_object(minhash_obj);

	containmenthash_obj->populate_bloom_filter(input_kmers_1,hash_functions);
	//containmenthash_obj->save_object(containmenthash_obj);
	ofstream file_obj;
	file_obj.open("./ContainmentHash_Obj_Save.txt", ios::binary);
	file_obj.write((char*)containmenthash_obj, sizeof(ContainmentHash));
	//file_obj.close();
	//print_mbits(containmenthash_obj.get_mbits());

	return 0;
}