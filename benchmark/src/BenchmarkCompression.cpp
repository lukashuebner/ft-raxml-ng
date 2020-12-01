#include <iostream>
#include <string>
#include "../../src/MSACompression.hpp"

using namespace std;

size_t constexpr RETURN_SUCCESS = 0;
size_t constexpr RETURN_INVALID_ARGUMENTS = 1;
size_t constexpr RETURN_CONSISTENCY_CHECK_FAILED = 2;

void load_msa_and_build_parsimony_tree(string filename, PartitionedMSA& part_msa, Tree& tree, unsigned int random_seed) {
    // Load MSa file into memory
    MSA load_msa = msa_load_from_file(filename, FileFormat::autodetect);

    // Build the parsimony tree
    unsigned int score;

    unsigned int attrs = sysutil_simd_autodetect();
    //attrs |= PLL_ATTRIB_PATTERN_TIP;

    part_msa.emplace_part_info("p1", DataType::dna, "GTR");
    part_msa.full_msa(std::move(load_msa));
    part_msa.split_msa();
    
    tree = Tree::buildParsimony(part_msa, random_seed, attrs, &score);
}

double proportion_of_gaps(const MSA& msa) {
    uint64_t nGaps = 0;
    for (auto sequence = msa.begin(); sequence != msa.end(); sequence++) {
        for (auto site: *sequence) {
            nGaps += (site == '-');
        }
    }

    return (double)nGaps / (msa.num_sites() * msa.size());
}

double identity_score(const string& seq1, const string& seq2) {
    uint64_t identical_sites = 0;
    assert(seq1.length() == seq2.length());

    for (size_t idx = 0; idx < seq1.length(); idx++) {
        identical_sites += (seq1[idx] == seq2[idx]);
    }
    
    return (double)identical_sites / seq1.length();
}

double average_pairwise_identity_score(const MSA& msa) {
    double identity_score_sum = 0;
    uint64_t nPairs = 0;

    for (auto seq1 = msa.begin(); seq1 != msa.end(); seq1++) {
        for (auto seq2 = msa.begin(); seq2 != msa.end(); seq2++) {
            if (seq1 == seq2) {
                continue;
            }
            identity_score_sum += identity_score(*seq1, *seq2);
            nPairs++;
        }
    }

    return identity_score_sum / nPairs;
}

void print_msa_stats(const MSA& msa, const string msa_name) {
    cout << "msa,proportionOfGaps,nPatterns,nSites,avgPairwiseIdentityScore" << endl;

    // Name of the msa
    cout << msa_name << ",";

    // Proportion of gaps
    cout << proportion_of_gaps(msa) << ",";

    // Number of patterns
    cout << msa.num_patterns() << ",";

    // Sequence length
    cout << msa.num_sites() << ",";

    // Average pairwise identity score
    cout << average_pairwise_identity_score(msa) << endl;
}

int benchmark_compression(string msa_file, unsigned int random_seed) {
    PartitionedMSA part_msa;
    Tree tree;

    string msa_name = msa_file.substr(msa_file.find_last_of("/\\") + 1);
    
    // cout << "Loading file and building parsimony tree" << endl;
    load_msa_and_build_parsimony_tree(msa_file, part_msa, tree, random_seed);
    //part_msa.compress_patterns();
    
    const MSA& msa = part_msa.full_msa();

    print_msa_stats(msa, msa_name);
    //return RETURN_SUCCESS;

    // Compress
    // cout << "Compressing" << endl;
    auto start_compression = chrono::high_resolution_clock::now();
    MSATreeCompression compressed_msa(tree);
    size_t compressed_size_in_bytes = compressed_msa.compress(msa);
    auto end_compression = chrono::high_resolution_clock::now();
    auto time_for_compression = (end_compression - start_compression);

    // Decompress and check if decompression ° compression = identity 
    // cout << "Decompressing" << endl;
    auto start_decompression = chrono::high_resolution_clock::now();
    RangeList all_sites = { SiteRange(0, msa.num_sites()) };
    auto decompressed_msa = compressed_msa.decompress(all_sites);
    auto end_decompression = chrono::high_resolution_clock::now();
    auto time_for_decompression = (end_decompression - start_decompression);

    // Consistency checks
    // cout << "Consistency checks" << endl;
    if (msa.size() == 0) {
        cout << "The MSA has zero taxons." << endl;
        return RETURN_CONSISTENCY_CHECK_FAILED;
    } else if (decompressed_msa->size() != msa.size()) {
        cout << "The input and decompressed MSA have different number of taxons." << endl;
        return RETURN_CONSISTENCY_CHECK_FAILED;
    } else {
        for (string taxon: msa.labels()) {
            string no_question_marks_sequence = msa.at(taxon);
            replace(no_question_marks_sequence.begin(), no_question_marks_sequence.end(), '?', '-');
            if (decompressed_msa->at(taxon) != no_question_marks_sequence) {
                cout << "The decompressed sequence did not match the original sequence" << endl;
                return RETURN_CONSISTENCY_CHECK_FAILED;
            }
        }
    }

    // Output stats
    uint8_t constexpr MSA_BYTES_PER_NUCLEOTIDE_STATE = 1;
    size_t uncompressed_size_in_bytes = msa.size() * msa.num_sites() * MSA_BYTES_PER_NUCLEOTIDE_STATE;
    cout << "msa,compressed_size_in_bytes,uncompressed_size_in_bytes,time_for_compression_in_ms,time_for_decompression_in_ms" << endl;
    cout << msa_name << "," << compressed_size_in_bytes << "," << uncompressed_size_in_bytes << "," << time_for_compression / chrono::milliseconds(1) << "," << time_for_decompression / chrono::milliseconds(1) << endl;
    return RETURN_SUCCESS;
}

int main(int argc, char** argv) {
    if (argc != 3 || string(argv[1]) != "msa-compression") {
        cout << "Usage " << argv[0] << " msa-compression <msa-file>" << endl;
        return RETURN_INVALID_ARGUMENTS;
   }

    string msa_file(argv[2]);
    for (uint64_t random_seed = 0; random_seed < 1; random_seed++) {
        int ret =  benchmark_compression(msa_file, random_seed);
        if (ret != RETURN_SUCCESS) {
            return ret;
        }
    } 
}

