#ifndef __COMMON_H__
#define __COMMON_H__


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <list>
#include <random>
#include <string>
#include <numeric>
#include <unordered_map>
#include <vector>


#include <armadillo>
#include <gsl/gsl_cdf.h>


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
    long long unsigned int count_index;
    std::vector<long long unsigned int> member_count_indices;
    FastaRecord *fasta;
    arma::Col<double> otu_counts;
    double abundance;
};


// Struct to hold all OTU-related data
struct OtuData {
    OtuTable *table;
    std::unordered_map<std::string,FastaRecord> *fasta;
    arma::Col<long long unsigned int> *otu_indices_abundance;
    std::vector<MergeOtu> *merged_otus;
};


// Read OTU counts from file
OtuTable read_otu_table_from_file(std::string &otu_count_fp);


// Read FASTA from file
std::unordered_map<std::string,FastaRecord> read_fasta_from_file(std::string &fasta_fp);


// Write merged OTU counts to file
void write_otu_table_to_file(std::vector<MergeOtu> &merged_otus, OtuData &otu_data, std::string &output_fp);


// Write merged OTU members to file
void write_merged_otu_members_to_file(std::vector<MergeOtu> &merged_otus, OtuData &otu_data, std::string &output_fp);


// Function to convert optarg into double
double double_from_optarg(const char *optarg);


#endif
