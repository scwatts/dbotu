#ifndef __DBOTU_H__
#define __DBOTU_H__


#include "common.h"
#include "levenshtein.h"


struct OtuData {
    OtuTable *table;
    std::unordered_map<std::string,FastaRecord> *fasta;
    arma::Col<long long unsigned int> *otu_indices_abundance;
    std::vector<MergeOtu> *merged_otus;
};


bool merge_otu(long long unsigned int &otu_index, double &otu_min_abundance, OtuData &otu_data);


#endif
