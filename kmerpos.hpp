#include <string>
using namespace std;

// add the unique kmer to a hash table
int read_kmer(const char *k_path, int k);

// search the unique kmer for reference; Using 64 bit binary int to store the k-mer 
int build_pos(const char  *fasta_file, int k, int t);

// read hic reads and constructet hic map depend on the max number of the unique kmer.
int contribute_hic_map(const char *read1, const char *read2, const char *output, int t);
