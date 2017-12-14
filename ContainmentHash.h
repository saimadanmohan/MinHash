#ifndef _CONTAINMENTHASH_H_
#define _CONTAINMENTHASH_H_

#include<iostream>
#include<vector>
#include<string>
#include<limits.h>
#include<unordered_map>
#include<fstream>
#include "MurmurHash3.h"
#define KMER_LENGTH 3
#define SEED 42
#define FALSE_POSITIVE_PROBABILITY 0.00001
#define LONG_READ_FILE_NAME "long_read.txt"
using namespace std;



class ContainmentHash{
	public:
	int optimal_length_of_bloom_filter_m;
	int* m_ptr;
	//vector<int>  m_bitz;
	int *m_bitz;
	int size ;


	void print_mbits(int* mbitz){
		for(int i=0;i<size;i++)
			cout<<mbitz[i]<<" ";
		cout<<endl;
	}

	std::vector<int> get_mbits(){
		cout<<"inside_mbits"<<size<<endl;
		print_mbits(m_bitz);
		vector<int> vec;	
		return vec;
	}
	uint32_t nthHash(int n, uint32_t hashA, uint32_t hashB, int filterSize) {
		//double hashing technique
	    return (hashA + n * hashB) % filterSize;
	}

	vector<int> add_elements_to_bloom_filter(vector<string> string_set,int num_hash_functions_k)
	{
		vector<int> m_vec;
		int m_bits[optimal_length_of_bloom_filter_m];

		for(int i=0;i<optimal_length_of_bloom_filter_m;i++){
			m_bits[i] = 0;
		}

		for(int i=0; i<string_set.size(); i++){
			uint32_t hash[1];
			string string_i = string_set[i];
			char* data = &string_i[0];
			MurmurHash3_x86_32 (data, 3 ,0, hash);
			for(int i=0; i< num_hash_functions_k; i++)
			{
				m_bits[nthHash(i, hash[0], hash[0], optimal_length_of_bloom_filter_m)] = 1;
			}
		}

		for(int i=0;i<optimal_length_of_bloom_filter_m;i++){
			m_vec.push_back(m_bits[i]);
		}
		//cout<<endl;	

		return m_vec;
	}

	int test_for_membership(string str, int* m_bits,int num_hash_functions_k){
		uint32_t hash[1];
		MurmurHash3_x86_32(&str[0], str.length(),0,hash);

		for(int i=0; i<num_hash_functions_k;i++){
			int index = nthHash(i,hash[0], hash[0], optimal_length_of_bloom_filter_m);
			int val = m_bits[index];
			if(!val){
				return 0;
			}
		}

		return 1;
	}

	double compute_containment_index(vector<string> string_set,int num_hash_functions_k){
		int y_values[num_hash_functions_k];
		double containment_index;
		int yk=0;
		vector<vector<uint32_t> > matrix;

		for(int i=0; i<string_set.size(); i++){
			vector<uint32_t> hash_vector;
			uint32_t hash[1];
			string string_i = string_set[i];
			char* data = &string_i[0];
			MurmurHash3_x86_32 (data, 3 ,0, hash); 
			for(int j=0;j<num_hash_functions_k; j++){
				hash_vector.push_back(nthHash(j, hash[0], hash[0], optimal_length_of_bloom_filter_m));
			}
			matrix.push_back(hash_vector);
		}

		/*for(int i =0; i<string_set.size(); i++){
			for(int j=0; j<num_hash_functions_k; j++){
				//cout<<matrix[i][j]<<"    ";
			}
			//cout<<endl;
		}*/

		for(int i=0;i<num_hash_functions_k;i++){
			uint32_t min = matrix[0][i];
			string min_kmer=string_set[0];
			for(int j=0;j<string_set.size();j++){
				if(matrix[j][i] < min){
					min = matrix[j][i];
					min_kmer = string_set[j];
				}

			}
			y_values[i] = test_for_membership(min_kmer, m_bitz,num_hash_functions_k);
			yk += y_values[i];
		}

		//cout<<"yk:"<<yk<<endl;
		//cout<<"num_hash_functions_k:"<<num_hash_functions_k<<endl;
		containment_index = (float(yk)/(float)num_hash_functions_k) - FALSE_POSITIVE_PROBABILITY;

		return containment_index;
	}


	/*
	Cest=(Yk/k)-p; 
	where, Yk = E(i=1 to k)Yi
	where Yi = 1 if min(hi(A)) is in bloom filter; 0 otherwise
	*/
	void populate_bloom_filter(vector<string> string_set_A,int hash_functions){

		int size_B = string_set_A.size();
		//cout<<"sizeA"<<size_A<<endl;
		//cout<<"size_B"<<size_B<<endl;

		int cardinality =  string_set_A.size();
		//int cardinality = (size_A<=size_B) ? string_set_B.size() : string_set_A.size();

		optimal_length_of_bloom_filter_m = (-cardinality) * (log(FALSE_POSITIVE_PROBABILITY)/pow(log(2), 2));
		//num_hash_functions_k = (optimal_length_of_bloom_filter_m/cardinality) * log(2);
		
		//cout << "optimal_length_of_bloom_filter_m=" << optimal_length_of_bloom_filter_m << endl ;
		//cout << "num_hash_functions_k=" << num_hash_functions_k << endl ;

	 	//vector<int>  m_bits = (size_A<=size_B) ? add_elements_to_bloom_filter(string_set_B) : add_elements_to_bloom_filter(string_set_A);
		
		vector<int> m_bits = add_elements_to_bloom_filter(string_set_A,hash_functions);
		m_bitz = (int *)malloc(m_bits.size()*sizeof(int));
		for(int i=0;i<m_bits.size();i++)
			m_bitz[i] = m_bits[i];
		size = m_bits.size();

		//return compute_containment_index(string_set_B, m_bits);

	}


	double get_jaccard_index(vector<string>string_set_A, vector<string>string_set_B){
		int size_A = string_set_A.size();
		int size_B = string_set_B.size();
		vector<string> union_set;
		vector<string> intersection_set;

		for(int i=0; i<size_A; i++){
			union_set.push_back(string_set_A[i]);
		}

		int flag;
		for(int i=0; i<size_B; i++){
			flag=0;
			for(int j=0; j<size_A; j++){
				if(string_set_B[i] == string_set_A[j])
					flag = 1;
			}
			if(!flag)
				union_set.push_back(string_set_B[i]);
		}

		for(int i=0; i<size_B; i++){
			flag =0;
			for(int j=0; j<size_A && !flag; j++){
				if(string_set_B[i] == string_set_A[j]){
					flag = 1;
				}
			}
			if(flag)
				intersection_set.push_back(string_set_B[i]);
		}


		//cout<<"intersect_count="<<intersection_set.size()<<endl;
		//cout<<"union_total="<<union_set.size()<<endl;
		
		return (double)intersection_set.size()/(double)union_set.size();
	}


	vector<string> get_kmers(string input_string,int kmer_len){
		int end = input_string.size() - kmer_len + 1;
		vector<string> input_kmers;
		for(int i=0;i<end;i++){
			string kmer = input_string.substr(i,kmer_len);
			input_kmers.push_back(kmer);
		}
		return input_kmers;
	}	

	void save_object(ContainmentHash *obj){
		ofstream file_obj;
 
		// Opening file in append mode
		file_obj.open("./ContainmentHash_Obj_Save.txt", ios::binary);
		file_obj.write((char*)obj, sizeof(ContainmentHash));
		//file_obj.close();
	}

	
};



#endif
