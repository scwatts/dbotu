#ifndef __DBOTU_H__
#define __DBOTU_H__


#include "common.h"


bool merge_otu(long long unsigned int otu_index, OtuTable &otu_table, std::unordered_map<std::string,FastaRecord> &fasta_records);

#endif
