#ifndef __COMMON_H__
#define __COMMON_H__


#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <list>
#include <string>
#include <numeric>
#include <unordered_map>
#include <vector>


#include <armadillo>


// Struct to hold OTU counts and associated information
struct OtuTable {
    // Samples names and OTU names
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_names;

    // Dimensions of the table
    unsigned int sample_number = 0;
    unsigned int otu_number = 0;

    // OTU counts
    arma::Mat<double> otu_counts;

    // OTU total read counts
    arma::Row<double> otu_sum_totals;
};


struct FastaRecord {
    std::string description;
    std::string sequence;
};


// Struct to contain a merged set of OTUs
struct MergeOtu {
    FastaRecord *fasta;
    double abundance;
};


// Read OTU counts from file
OtuTable read_otu_table_from_file(std::string &otu_count_fp);


// Read FASTA from file
std::unordered_map<std::string,FastaRecord> read_fasta_from_file(std::string &fasta_fp);


// Function to convert optarg into double
double double_from_optarg(const char *optarg);


#endif
