#include<iostream>
#include<vector>
#include<string>
#include<limits.h>
#include<unordered_map>
#include<fstream>
#include<math.h>
#include<time.h>
#include "MurmurHash3.h"

#define KMER_LENGTH 3
#define LONG_READ_FILE_NAME "long_read.txt"
#define FALSE_POSITIVE_PROBABILITY 0.00001
using namespace std;

class MinHash{
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
		return sketch_vector;
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

	void save_object(MinHash obj){
		ofstream file_obj;
 
		// Opening file in append mode
		file_obj.open("./MinHash_Obj_Save.txt", ios::app);
		file_obj.write((char*)&obj, sizeof(obj));
	}

	MinHash restore_object(){
		MinHash new_obj;
		ifstream file_obj;
		file_obj.open("Input.txt", ios::in);
		file_obj.read((char*)&new_obj, sizeof(new_obj));
		return  new_obj;
	}

};


class ContainmentHash{
	int num_hash_functions_k;
	int optimal_length_of_bloom_filter_m;
	int* m_ptr;

	uint32_t nthHash(int n, uint32_t hashA, uint32_t hashB, int filterSize) {
		//double hashing technique
	    return (hashA + n * hashB) % filterSize;
	}

	vector<int> add_elements_to_bloom_filter(vector<string> string_set)
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

	int test_for_membership(string str, vector<int>  m_bits){
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

	double compute_containment_index(vector<string> string_set, vector<int>  m_bits){
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

		for(int i =0; i<string_set.size(); i++){
			for(int j=0; j<num_hash_functions_k; j++){
				//cout<<matrix[i][j]<<"    ";
			}
			//cout<<endl;
		}

		for(int i=0;i<num_hash_functions_k;i++){
			uint32_t min = matrix[0][i];
			string min_kmer=string_set[0];
			for(int j=0;j<string_set.size();j++){
				if(matrix[j][i] < min){
					min = matrix[j][i];
					min_kmer = string_set[j];
				}

			}
			y_values[i] = test_for_membership(min_kmer, m_bits);
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
	double get_containment_index(vector<string> string_set_A, vector<string> string_set_B, int num_hash_functions_k){

		int size_A = string_set_A.size();
		int size_B = string_set_B.size();
		//cout<<"sizeA"<<size_A<<endl;
		//cout<<"size_B"<<size_B<<endl;

		int cardinality = (size_A<=size_B) ? string_set_B.size() : string_set_A.size();

		optimal_length_of_bloom_filter_m = (-cardinality) * (log(FALSE_POSITIVE_PROBABILITY)/pow(log(2), 2));
		//num_hash_functions_k = (optimal_length_of_bloom_filter_m/cardinality) * log(2);
		
		//cout << "optimal_length_of_bloom_filter_m=" << optimal_length_of_bloom_filter_m << endl ;
		//cout << "num_hash_functions_k=" << num_hash_functions_k << endl ;

	 	vector<int>  m_bits = (size_A<=size_B) ? add_elements_to_bloom_filter(string_set_B) : add_elements_to_bloom_filter(string_set_A);

		return (size_A<=size_B) ? compute_containment_index(string_set_A, m_bits) : compute_containment_index(string_set_B, m_bits);

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


	vector<string> get_kmers(string input_string){
		int end = input_string.size() - KMER_LENGTH + 1;
		vector<string> input_kmers;
		for(int i=0;i<end;i++){
			string kmer = input_string.substr(i,KMER_LENGTH);
			input_kmers.push_back(kmer);
		}
		return input_kmers;
	}	


	
};


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

int main(int argc, char const *argv[])
{
	MinHash minhash_obj;

	string input1 = get_long_read();
	string input2 = "ATTAACGTCTATCCGCCTCACTAGGACGACAGGATTTTTAATATCCGCCG";
	//TGGGTCGCGGCCGCGCCTTACTCGCAGCTTAGTTGTAAGAGTATGACACGGATCCACGTAGAACTAGCTGCGCCCCCGCGTACGACCAGACAGTCAGAGAATAAATCGCGATCCAATAGCGGTCATTTCTGGTATGAACTACATGACAAACTGCACTGACTTTTGAGCCTGGGTTTTTTTAAGCCTGATGGTTCAAATGCATGCACCGGGCAGCCGCCACTGCACACCGGAAGACGGGATCCAGAGACAGGTGCTTTAGAGATTCCCCATTGTGTGTCGGGCGATTAGTGACAATCGAGATGGATCGGGCACTAGTTTCGGGATCGTCCTTTCCAACTGCTCGGAAAAACGGGATTATCCTCCCGTCAGCGCCTACATGCCTTGTTTAGTTTTATCTTATCCTGTTTTTAGCCAGTCATCGGAGTATATTATGG
	int hash_functions = 100;
	uint32_t seed = 42;
	int kmer_length = KMER_LENGTH;

	cout<<"hash_functions"<<"    "<<"actual_jaccardian"<<"    "<<"minhash_jaccardian_distance"<<"    "<<"	diff"<<endl;
	for(int k=0;k<50;k++){
		time_t t1 = time(0);
		vector<string> input_kmers_1 = minhash_obj.get_kmers(input1,kmer_length);
		vector<string> input_kmers_2 = minhash_obj.get_kmers(input2,kmer_length);

		double actual_jaccardian = minhash_obj.calculate_jaccard_similarity(input_kmers_1,input_kmers_2);
		vector<vector<uint32_t> > characteristic_matrix = minhash_obj.get_characteristic_matrix(hash_functions,input_kmers_1,seed,kmer_length);
		minhash_obj.save_object(minhash_obj);

		MinHash new_obj = minhash_obj.restore_object();

		vector<vector<uint32_t> > characteristic_matrix_2 = new_obj.get_characteristic_matrix(hash_functions,input_kmers_2,seed,kmer_length);

		//print_characteristic_matrix(characteristic_matrix);
		//print_characteristic_matrix(characteristic_matrix_2);
		vector<uint32_t> sketch_vector = new_obj.get_sketch_vector(characteristic_matrix);
		vector<uint32_t> sketch_vector_2 = new_obj.get_sketch_vector(characteristic_matrix_2);

		//print_sketch_vector(sketch_vector);
		//print_sketch_vector(sketch_vector_2);
		int intersect_count = 0;
		for(int i=0;i<sketch_vector.size();i++)
			for (int j = 0; j < sketch_vector_2.size(); ++j)
				if(sketch_vector[i] == sketch_vector_2[j]){
					//cout<<sketch_vector[i]<<" "<<sketch_vector_2[j]<<endl;
					intersect_count = intersect_count + 1;
				}

		//cout<<"sketch_vector size"<<sketch_vector.size()<<endl;
		int union_total = sketch_vector.size()+sketch_vector_2.size() - intersect_count;
		//cout<<"intersect_count "<<intersect_count<<endl;
		//cout<<"union_total "<<union_total<<endl;
		double minhash_jaccardian_distance = intersect_count/(union_total*1.0);
		//cout<<"jaccardian_distance "<<jaccardian_distance<<endl;

		double diff = abs(actual_jaccardian - minhash_jaccardian_distance)/(actual_jaccardian*1.0);
		//cout<<hash_functions<<"    		"<<diff<<endl;
		time_t t2 = time(0);
		double diff2 = difftime(t2, t1);
		cout<<hash_functions<<" "<<diff2<<endl;
		hash_functions += 100;
		
	}
	return 0;
}