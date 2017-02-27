#ifndef __DBOTU_H__
#define __DBOTU_H__


#include "common.h"
#include "levenshtein.h"


// Core algorithm for merging OTUs
bool merge_otu(long long unsigned int &otu_index, double &otu_min_abundance, double &min_distance, double &max_pvalue, OtuData &otu_data);


// Function used in calculation of the p value
inline double d_helper(arma::Col<double> counts);


#endif
