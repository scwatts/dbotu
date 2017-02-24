#include "dbotu.h"
#include "dbotu_opts.h"


int main(int argc, char **argv) {
    // Get and parse commandline arguments
    DbotuOptions dbotu_options = get_commandline_arguments(argc, argv);


    // Load OTU count table from file and sum OTU counts
    OtuTable otu_table = read_otu_table_from_file(dbotu_options.input_otu_counts_fp);
    otu_table.otu_sum_totals = arma::sum(otu_table.otu_counts, 0);

    // Load FASTA records
    std::unordered_map<std::string,FastaRecord> fasta_records = read_fasta_from_file(dbotu_options.input_fasta_fp);


    // Collect indexes in order of most abundant OTU
    arma::Col<long long unsigned int> otu_indices_abundance = arma::sort_index(otu_table.otu_sum_totals, "descend");

    // Iterate from the most abundant OTU to the least and assign OTUs
    std::vector<MergeOtu> merged_otus;

    for (auto &otu_index : otu_indices_abundance) {
        // TODO: check that it really makes sense to be accummulating total OTU counts here

        // If not OTU merged (function call with boolean retval), then create new OTU group
        if (!merge_otu(otu_index, otu_table, fasta_records)) {
            // Create new OTU
        }

        break;
    }

    return 0;
}


bool merge_otu(long long unsigned int otu_index, OtuTable &otu_table, std::unordered_map<std::string,FastaRecord> &fasta_records) {
    // Collect name of OTU sequence
    std::string otu_name = otu_table.otu_names[otu_index];

    // Get FASTA record using otu_name to lookup
    FastaRecord otu_fasta = fasta_records[otu_name];

    // Abundance test

    // Distance test

    // Statistical significance test

    // Merge if all passed

    printf("%s\n", otu_fasta.description.c_str());
    printf("%s\n", otu_fasta.sequence.c_str());

    return true;
}
