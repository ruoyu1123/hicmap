# include "kmerpos.hpp"
# include <iostream>
# include "clipp.h"
# include <vector>
using namespace clipp; using std::cout; using std::string; using std::vector;


int main(int argc, char *argv[]){
    string kmerfile;        // jellyfish result dump
    string reference;
    string output = "./";          //output path
    string fq1;
    string fq2;
    int t = 8;
    int k_size = 21;        //k-mer length

    auto cli = (
        option("-t", "--threads").doc("Number of threads [8]") & value("threads", t),
        option("-o").doc("Output path [./]") & value("outpath", output),
        value("kmer_dump", kmerfile).doc("Dump file of jellyfish result"),
        value("reference", reference).doc("fasta file of reference"),
        value("read1.fq", fq1).doc("Hi-C reads1"),
        value("read2.fq", fq2).doc("Hi-C reads2")
    );

    auto fmt = doc_formatting{}
		.first_column(8)                           //left border column for text body
		.doc_column(30)                            //column where parameter docstring starts
		.last_column(100);

    if(!parse(argc, const_cast<char **>(argv), cli)) {
		cout << "Usage:\n" << usage_lines(cli, "HiCmap")
     << "\nParameter:\n" << documentation(cli,fmt) << "\nERROR: Required parameter missing\n";
		// throw "Division by zero condition!";
		exit(0);
	}
    read_kmer(kmerfile.c_str(), k_size);
    build_pos(reference.c_str(), k_size, t);
    contribute_hic_map(fq1.c_str(), fq2.c_str(), output.c_str(), t);
    return 0; 
}
