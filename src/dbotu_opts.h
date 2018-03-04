#ifndef __DBOTU_OPTS_H__
#define __DBOTU_OPTS_H__


#include <getopt.h>


#include "common.h"
#include "version.h"


struct DbotuOptions {
    // Set some default parameters
    double min_distance = 0.1f;
    double min_abundance = 10.0f;
    double max_pvalue = 0.0005f;

    // Declare some important variables
    std::string input_otu_counts_fp;
    std::string input_fasta_fp;
    std::string output_otu_counts_fp;
    std::string output_membership_fp;

    // Number of threads
    unsigned int threads = 1;
};


// Help and version text
void print_help();
void print_version();


DbotuOptions get_commandline_arguments(int argc, char **argv);


#endif
