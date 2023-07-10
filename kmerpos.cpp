#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <bitset>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cstring>
#include "thread_pool.hpp"
#include "kmerpos.hpp"
#include "stringtools.hpp"

using namespace std;

mutex m;
int k_size;          // k-mer len
unordered_map<uint64_t, string> kmer;

std::unordered_map<char, uint64_t> ctoi = {
	{'A', 0ull},
	{'T', 3ull},
	{'C', 1ull},
	{'G', 2ull}
};

// 字符串转换成二进制
uint64_t tobin(string* s){
    uint64_t k_value = 0;
    int i = 0;
    for (auto c = s->begin(); c < s->end(); c++){
        uint64_t v = ctoi[*c];
        k_value = k_value << 2 | v;
    }
    return k_value;
}


int read_kmer(const char *k_path, int k)
{
    ifstream k_mer(k_path);
    uint64_t v;
    uint64_t k_value = 0;
    // printf("%s",s);
	string line;
    while(std::getline(k_mer, line)){
		string tmp = line.substr(0, k);	
		v = tobin(&tmp);
		kmer[v] = "";
	}
	printf("Total kmer number: %ld\n", kmer.size());
    return 1;
}

/* store the sbin(bin of kmer),rbin(reverse bin of kmer), pos at contig and contig name*/
class KMER{
public:
    uint64_t sbin, rbin;
    uint32_t pos;
    string nameid;
    int is_ukmer(){
        if (kmer.find(sbin) != kmer.end()){
		   if (kmer[sbin] != "")
				cout << "+KMER不是唯一" << endl; 
		   kmer[sbin] = nameid;	
		}
	    if (kmer.find(rbin) != kmer.end()){
			if (kmer[rbin] != "")
				cout << "-KMER不是唯一" << endl; 
            kmer[rbin] = nameid;
        }
        return 0;
    }
};

int search_kmer(string line, string name, int n){
    KMER k;
	k.nameid = name;
    uint64_t base;
    uint32_t len = line.length();
    string qkmer = line.substr(0, k_size);
    k.sbin = tobin(&qkmer);
    reverse(qkmer.begin(),qkmer.end()); // 第一个kmer的反向kmer
    k.rbin = ~tobin(&qkmer); // 第一个kmer的反向互补
	// 遍历该序列寻找unique kmer
    for (uint32_t i = 21; i < len; i++){
        k.pos = i - 21;
        k.is_ukmer();
        base = ctoi[line[i]];
        k.sbin = (k.sbin << 2 | base) & 0x3ffffffffffull;
        k.rbin = (k.rbin >> 2) | ((3ull-base) << 40);
    }
    return 1;
}

int build_pos(const char *fasta_file, int k, int t){
    cout << "Number of threads: " << t << "\n"
         << "Reading the kmer to dicts\n";
    ThreadPool fpool(t);
    fpool.init();
    ifstream fa(fasta_file);
    string line, name;
	k_size = k;
	cout << "kmer size:"<<k_size<<"\n";
    // auto size=sizeof(KMER);
    cout << "Searching kmer in fasta file!\n";
    while (getline(fa, line)){
        if (line[0] == '>')
            name = line.substr(1, line.find(" ")-1);
        else
            fpool.submit(search_kmer, line, name, t);
	}
    fpool.shutdown();
    return 1;
}



int read_belong(string seq, string readid, vector<pair<string, string>> *readhash){
	unordered_map<string, int> read2contig;
	uint32_t len = seq.length();
	uint64_t base;
    string qkmer = seq.substr(0, 21);
    uint64_t sbin = tobin(&qkmer);
    reverse(qkmer.begin(),qkmer.end()); // 第一个kmer的反向kmer
    uint64_t rbin = ~tobin(&qkmer); // 第一个kmer的反向互补
	//遍历reads，寻找具有共同uniquekmer的contig并记录各个contig出现的次数。
    for (uint32_t i = 21; i < len; i++){
       	if (kmer.find(sbin) != kmer.end())
			if (read2contig.find(kmer[sbin]) != read2contig.end())
				read2contig [kmer[sbin]]++;
			else
				read2contig [kmer[sbin]]=1;
		if (kmer.find(rbin) != kmer.end())
			if (read2contig.find(kmer[rbin]) != read2contig.end())
				read2contig [kmer[rbin]]++;
			else
        		read2contig [kmer[rbin]]=1;
		base = ctoi[seq[i]];
        sbin = (sbin << 2 | base) & 0x3ffffffffffull;
        rbin = (rbin >> 2) | ((3ull-base) << 40);
	} //for
	if (read2contig.size() > 3 || read2contig.size() == 0){
		m.lock();
		(*readhash).push_back(make_pair(readid, ""));
		m.unlock();
		return 0;
	}
	auto maxValue = read2contig.begin()->second; // 初始化最大值为第一个值
    auto maxContig = read2contig.begin()->first;
	for (auto p : read2contig) 
        if (p.second > maxValue){ 
            maxValue = p.second;
			maxContig = p.first;
		}
    // add contig of maxValue biger than 2 to read2contig
	if (maxValue >= 2){
		m.lock();
		(*readhash).push_back(make_pair(readid, maxContig));
		m.unlock();
	}
	else{
		m.lock();
		(*readhash).push_back(make_pair(readid, ""));
		m.unlock();
	}
	return 0;
}


int position_reads(const char *filepath, const char *output, vector<pair<string, string>> *read_ownership, int t){
	
	ifstream inputFile(filepath);
	if (!inputFile) {
		std::cerr << "ERROR opening file." << endl;
		exit(0);
		return 1;
	}
	ThreadPool fpool(t);
    fpool.init();
	std::string line;
	int lineCount = 0;
	string name;
	while (getline(inputFile, line)) {
		lineCount++;
		if (lineCount % 4 == 1)
			name = line.substr(1, line.find(" ")-1);
		else if(lineCount %4 == 2) 
			fpool.submit(read_belong, line, name, read_ownership);
	}
	fpool.shutdown();
	FILE* outfile = fopen(output, "w");
	if (outfile == NULL)
		std::cerr << "ERROR: output file not exist" << std::endl;
	cout << "read  to contig size()\t"<<(*read_ownership).size()<<endl;
	std::sort((*read_ownership).begin(), (*read_ownership).end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
	for (auto p:(*read_ownership))
		fprintf(outfile, "%s\t%s\n", p.first.c_str(), p.second.c_str());
	fclose(outfile);
	return 1;
}

int contribute_hic_map(const char *read1, const char *read2, const char *output, int t){
	/*
	 * This function can construct the hic map of contigs and hic reads.
	 * 1. Obtain the belonging of every read.
	 * 2. Connect the contig of paired reads 
	 * 3. 
	 * */
	
	vector<pair<string, string>> r1contig; 
	vector<pair<string, string>> r2contig;
	char r1output[100];
	char r2output[100];
	char hicmapout[100];
	strcpy(hicmapout, output);
	strcat(hicmapout, "/raw_hicmap.tab");
	
	strcpy(r1output, output);
	strcat(r1output,"/read1_of_contig.tab");
	strcpy(r2output, output);
	strcat(r2output,"/read2_of_contig.tab");

	cout << hicmapout << "\t" << r1output << " " << r2output << endl;
	
	position_reads(read1, r1output, &r1contig, t);
	position_reads(read2, r2output, &r2contig, t);
	map<pair<string, string>, int> hicmap;
	cout << "Hi-C map constructing start" << endl; 
	for (int i = 0;i < r1contig.size();i++){
		if (r1contig[i].first != r2contig[i].first){
			cout << "ERROR: Read name does not correspond!" << endl;
			exit(0);
		}
		if (r1contig[i].second == "" || r2contig[i].second == "")
			continue;
		pair<string, string> key = make_pair(r1contig[i].second, r2contig[i].second);
		pair<string, string> key2 = make_pair(r2contig[i].second, r1contig[i].second);
		auto it = hicmap.find(key);
		if (it != hicmap.end())
			it->second++;
		else if (hicmap.find(key2) != hicmap.end())
			hicmap[key2]++;
		else
			hicmap[key] = 1;
	}
	cout<< "Hi-C map have been constructed\n"<<"Preparing to write Hi-C map to file"<<endl;
	FILE *outfile = fopen(hicmapout, "w");
	for (auto item : hicmap){
		fprintf(outfile, "%s\t%s\t%d\n", item.first.first.c_str(), item.first.second.c_str(), item.second);
	}
	fclose(outfile);
	cout<< "Hi-C map already have been written.\n"<<endl;
	return 0;
}

