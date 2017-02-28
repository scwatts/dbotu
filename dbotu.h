#ifndef __DBOTU_H__
#define __DBOTU_H__


#include "common.h"
#include "levenshtein.h"


// Core algorithm for merging OTUs
bool merge_otu(long long unsigned int &otu_index, double &otu_min_abundance, double &min_distance, double &max_pvalue, OtuData &otu_data);

// Slow pair wise vector addition
inline std::vector<double> element_wise_sum(std::vector<double> &vec_a, std::vector<double> &vec_b);

// Function used in calculation of the p value
inline double d_helper(std::vector<double> &counts);


#endif
