#include<iostream>
#include<vector>
#include<string>
#include<limits.h>
#include<unordered_map>
#include<fstream>
#include<stdlib.h>
#include "MinHash.h"
#include "MurmurHash3.h"
#include "ContainmentHash.h"
//#define KMER_LENGTH 11
//#define SEED 42
//#define FALSE_POSITIVE_PROBABILITY 0.00001
#define SHORT_READ_FILE_NAME "short_read.txt"
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


vector<string>  get_short_reads(){

	string line;
	string output = "";
	vector<string> short_reads;
	ifstream myfile (SHORT_READ_FILE_NAME);
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{	
			short_reads.push_back(line);
		}
		myfile.close();
	}

	else 
		cout << "Unable to open file"; 

	return short_reads;
}

double get_min_hash_value(MinHash *new_obj,vector<uint32_t> sketch_vector,string input2,int kmer_length,int hash_functions,int seed){
	//MinHash minhash_obj;
	//MinHash new_obj = minhash_obj.restore_object();
	vector<string> input_kmers_2 = new_obj->get_kmers(input2,kmer_length);
	
	vector<vector<uint32_t> > characteristic_matrix_2 = new_obj->get_characteristic_matrix(hash_functions,input_kmers_2,seed,kmer_length);
	//double actual_jaccardian = new_obj.calculate_jaccard_similarity(input_kmers_1,input_kmers_2);
	vector<uint32_t> sketch_vector_2 = new_obj->get_sketch_vector(characteristic_matrix_2);

	int intersect_count = 0;
	for(int i=0;i<sketch_vector.size();i++)
		for (int j = 0; j < sketch_vector_2.size(); ++j)
			if(sketch_vector[i] == sketch_vector_2[j]){
				intersect_count = intersect_count + 1;
			}

	int union_total = sketch_vector.size()+sketch_vector_2.size() - intersect_count;
	double minhash_jaccardian_distance = intersect_count/(union_total*1.0);

	return minhash_jaccardian_distance;
	//double diff = actual_jaccardian - minhash_jaccardian_distance;
	//cout<<hash_functions<<"    			"<<actual_jaccardian<<"   	 	"<<minhash_jaccardian_distance<<"    		"<<diff<<endl;
}

void print_mbits(std::vector<int> mbits){
	for(int i=0;i<mbits.size();i++)
		cout<<mbits[i]<<" ";
	cout<<endl;
}


ContainmentHash* restore_containment_hash_object(){
	ContainmentHash *new_obj = (ContainmentHash *)malloc(sizeof(ContainmentHash));
	ifstream file_obj;
	file_obj.open("./ContainmentHash_Obj_Save.txt", ios::in);
	file_obj.read((char*)new_obj, sizeof(ContainmentHash));
	//print_mbits(new_obj->get_mbits());
	//file_obj.close();
	return  new_obj;
}

double get_containment_hash_value(ContainmentHash *new_obj,string input1,string input2,int hash_functions){

	//ContainmentHash *new_obj = restore_containment_hash_object();
	
	int kmer_length = KMER_LENGTH;
	vector<string> long_read_kmers = new_obj->get_kmers(input1,kmer_length);
	vector<string> short_read_kmers = new_obj->get_kmers(input2,kmer_length);

	int size_A = short_read_kmers.size();
	int size_B = long_read_kmers.size();

	double containment_index = new_obj->compute_containment_index(short_read_kmers,hash_functions);
	double jaccard_index = new_obj->get_jaccard_index(long_read_kmers, short_read_kmers);
	double jaccard_estimate = (size_A*containment_index)/(size_A+size_B-(size_A*containment_index));
	double relative_error = ((jaccard_estimate - jaccard_index) / jaccard_index)*100;

	return jaccard_estimate;
	//cout<<hash_functions<<"				"<<containment_index<<"				"<<jaccard_estimate<<"			"<<jaccard_index<<"			"<<relative_error<<endl;

}

int main(int argc, char const *argv[])
{

	string input1 = get_long_read();
	vector<string> short_reads = get_short_reads();
	

	int hash_functions = HASH_FUNCTIONS;
	uint32_t seed = SEED;
	int kmer_length = KMER_LENGTH;

	
	ContainmentHash *new_obj = new ContainmentHash();
	MinHash *minhash_obj = new MinHash();
	vector<string> long_read_kmers = new_obj->get_kmers(input1,kmer_length);
	vector<vector<uint32_t> > characteristic_matrix = minhash_obj->get_characteristic_matrix(hash_functions,long_read_kmers,seed,kmer_length);
	vector<uint32_t> sketch_vector = minhash_obj->get_sketch_vector(characteristic_matrix);
	

	new_obj->populate_bloom_filter(long_read_kmers,hash_functions);
	//cout<<"hash_functions"<<"    "<<"actual_jaccardian"<<"    "<<"minhash_jaccardian_distance"<<"    "<<"	diff"<<endl;
	
	cout<<"actual_jaccardian"<<" "<<"containment_hash_jaccard_estimate"<<" "<<"minhash_jaccardian_distance"<<endl;
	for(int i=0;i<short_reads.size();i++){

		string input2 = short_reads[i];
		//cout<<input1<<endl;
		//cout<<input2<<endl;

		double containment_hash_jaccard_estimate = get_containment_hash_value(new_obj,input1,input2,hash_functions);
		double minhash_jaccardian_distance = get_min_hash_value(minhash_obj,sketch_vector,input2,kmer_length,hash_functions,seed);

		vector<string> short_read_kmers = new_obj->get_kmers(input2,kmer_length);
		double actual_jaccardian = minhash_obj->calculate_jaccard_similarity(long_read_kmers,short_read_kmers);
		cout<<actual_jaccardian<<" "<<containment_hash_jaccard_estimate<<" "<<minhash_jaccardian_distance<<endl;
	}
	return 0;
}