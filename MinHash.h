#ifndef _MINHASH_H_
#define _MINHASH_H_
#include<iostream>
#include<vector>
#include<string>
#include<limits.h>
#include<unordered_map>
#include<fstream>
#include "MurmurHash3.h"
#include "MinHash.h"
#define KMER_LENGTH 3
#define SEED 42
#define FALSE_POSITIVE_PROBABILITY 0.00001
#define LONG_READ_FILE_NAME "long_read.txt"
using namespace std;

class MinHash{
	vector<uint32_t> act_sketch_vector;
	public:
	vector<string> get_kmers(string input_string,int kmer_length){
		int end = input_string.size() - kmer_length + 1;
		vector<string> input_kmers;
		for(int i=0;i<end;i++){
			string kmer = input_string.substr(i,kmer_length);
			input_kmers.push_back(kmer);
		}
		return input_kmers;
	}

	void print_vector(vector<string> input_kmers){
		for(int i=0;i<input_kmers.size();i++){
			cout<<input_kmers[i]<<" ";
		}
		cout<<endl;
	}


	void print_sketch_vector(vector<uint32_t> sketch_vector){
		for(int i=0;i<sketch_vector.size();i++){
			cout<<sketch_vector[i]<<" ";
		}
		cout<<endl;
	}

	uint32_t xorshift32(uint32_t state)
	{
		uint32_t x = state;
		x ^= x << 13;
		x ^= x >> 17;
		x ^= x << 5;
		//state[0] = x;
		return x;
	}


	uint32_t xorshift_min(uint32_t a,uint32_t b)
	{
		if(a<b)
			return a;
		else
			return b;
	}

	void print_characteristic_matrix(vector<vector<uint32_t> >& characteristic_matrix){
		for (int i = 0; i < characteristic_matrix.size(); ++i){
			for (int j = 0; j < characteristic_matrix[i].size(); ++j){
				cout<<characteristic_matrix[i][j]<<" ";
			}
			cout<<endl;
		}
		cout<<endl<<endl;
	}

	vector<vector<uint32_t> > get_characteristic_matrix(int hash_functions,vector<string>& input_kmers,uint32_t seed,int kmer_length){
		vector<vector<uint32_t> > characteristic_matrix;
		uint32_t mod = 999999937;
		for(int i=0;i<input_kmers.size();i++){
			string kmer = input_kmers[i];
			vector<uint32_t> hashed_vector;
			uint32_t hash[4];
			char *s = &kmer[0];
			MurmurHash3_x86_128(s, kmer_length, seed, hash);	
			hashed_vector.push_back((hash[0]+hash[1]+hash[2]+hash[3])%mod);
			for(int i=0;i<hash_functions-1;i++){
				hashed_vector.push_back(xorshift32(hashed_vector[i]));
			}
			characteristic_matrix.push_back(hashed_vector);
		}
		return characteristic_matrix;
	}

	vector<uint32_t> get_sketch_vector(vector<vector<uint32_t> >& characteristic_matrix){
		vector<uint32_t> sketch_vector(characteristic_matrix[0].size(), UINT_MAX);
		//cout<<"UINT_MAX "<<UINT_MAX<<endl;
		//print_sketch_vector(sketch_vector);
		for (int i = 0; i < characteristic_matrix.size(); ++i){
			for (int j = 0; j < characteristic_matrix[i].size(); ++j){
				sketch_vector[j] = xorshift_min(sketch_vector[j],characteristic_matrix[i][j]);
			}
			//print_sketch_vector(sketch_vector);
		}
		act_sketch_vector = sketch_vector;
		return sketch_vector;
	}

	vector<uint32_t> get_stored_sketch_vector(){
		return act_sketch_vector;
	}

	unordered_map<string,int> getMap(vector<string> input_kmers_1){
		unordered_map<string,int> mymap_1;
	 	 for(int i=0;i<input_kmers_1.size();i++){
	 	 	string kmer = input_kmers_1[i];
	 	 	//unordered_map<string,int>::const_iterator got = mymap_1.find (kmer);
	 	 	//if(got == mymap_1.end()){
	 	 		mymap_1[kmer] = 1;
	 	 	//}else{
	 	 	//	mymap_1[kmer] = got->second+1;
	 	 	//}
	 	 }
	 	 return mymap_1;
	}


	double calculate_jaccard_similarity(vector<string>& input_kmers_1,vector<string>& input_kmers_2)
	{
		 unordered_map<string,int> mymap_1 = getMap(input_kmers_1);
		 unordered_map<string,int> mymap_2 = getMap(input_kmers_2);
		 double val = 0.0;
		 int union_tot = 0;
		 int intersection = 0;
		 for (auto it : mymap_1){ 
		 	string kmer_1 = it.first;
		 	bool found = false;
	    	for(auto it2 : mymap_2){
	    		string kmer_2 = it2.first;
	    		if(kmer_1.compare(kmer_2) == 0){
	    			found = true;
	    			intersection += 1;
	    			break;
	    		}
	    	}

	     }
	     int size1 = mymap_1.size();
	     int size2 = mymap_2.size();
		 union_tot = size1 + size2 - intersection;     
	     //cout<<"intersection "<<intersection<<endl;
	     //cout<<"union_tot "<<union_tot<<endl;
	     val = (intersection*1.0)/(union_tot*1.0);
	 	 return val;
	}

	void save_object(MinHash *obj){
		ofstream file_obj;
 
		// Opening file in append mode
		file_obj.open("./MinHash_Obj_Save.txt", ios::binary);
		file_obj.write((char*)obj, sizeof(MinHash));
	}


};
#endif
