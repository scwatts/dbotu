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


    // TODO: check that all OTUs have a matching FASTA record


    // Collect indexes in order of most abundant OTU
    arma::Col<long long unsigned int> otu_indices_abundance = arma::sort_index(otu_table.otu_sum_totals, "descend");


    // Place a pointer to related data into a struct
    std::vector<MergeOtu> merged_otus;
    OtuData otu_data { &otu_table, &fasta_records, &otu_indices_abundance, &merged_otus };


    // Iterate from the most abundant OTU to the least and assign OTUs
    for (auto &otu_index : otu_indices_abundance) {
        // TODO: check that it really makes sense to be accummulating total OTU counts here


        // Determine minimum abundance to merge current OTU
        double otu_abundance = static_cast<double>(otu_data.otu_indices_abundance->at(otu_index));
        double otu_min_abundance = otu_abundance * dbotu_options.min_abundance;

        // If not OTU merged (function call with boolean retval), then create new OTU group
        if (!merge_otu(otu_index, otu_min_abundance, otu_data)) {
            // Create new OTU
            // TODO: refactor to reduce code re-use
            std::string *otu_name = &otu_data.table->otu_names[otu_index];
            otu_data.merged_otus->push_back(MergeOtu { &otu_data.fasta->at(*otu_name), otu_abundance });
        }

    }

    return 0;
}


bool merge_otu(long long unsigned int &otu_index, double &otu_min_abundance, OtuData &otu_data) {
    // Get a pointer to OTU name, FASTA, and counts
    std::string *otu_name = &otu_data.table->otu_names[otu_index];
    FastaRecord *otu_fasta = &otu_data.fasta->at(*otu_name);

    // Running list of candidate OTUs to merge into
    // TODO: Use a list of pointers rather than copies of the MergeOtu struct
    std::list<MergeOtu> merged_otu_candidates;

    // Discover all OTUs which pass the abundance test
    for (auto &merged_otu : *otu_data.merged_otus) {
        if (merged_otu.abundance > otu_min_abundance) {
            merged_otu_candidates.push_back(merged_otu);
        }
    }

    // If not candidates are found, return
    if (merged_otu_candidates.empty()) {
        return false;
    }

    // Order candidate merge OTUs by increasing genetic distance
    std::vector<double> merged_otu_candidate_distances;
    for (auto &merged_otu_candidate : merged_otu_candidates) {
        // Calculate levenshtein distance
        size_t lev_distance = lev_edit_distance(merged_otu_candidate.fasta->sequence.length(),
                                                merged_otu_candidate.fasta->sequence.c_str(),
                                                otu_fasta->sequence.length(),
                                                otu_fasta->sequence.c_str(), 0);

        // Length-adjusted Levenshtein distance
        size_t sequence_length_sum = merged_otu_candidate.fasta->sequence.length() + otu_fasta->sequence.length();
        double distance = static_cast<double>(lev_distance) / (0.5 * static_cast<double>(sequence_length_sum));

        // Add distance to vector
        merged_otu_candidate_distances.push_back(distance);
    }

    // Get ordered indicies
    std::vector<unsigned int> distance_indices(merged_otu_candidate_distances.size());
    std::iota(distance_indices.begin(), distance_indices.end(), 0);
    std::sort(distance_indices.begin(), distance_indices.end(), [&merged_otu_candidate_distances](size_t i1, size_t i2) { return merged_otu_candidate_distances[i1] < merged_otu_candidate_distances[i2];});


    // Statistical significance test

    // Merge if all passed

    printf("%llu\n", otu_index);
    printf("%s\n", otu_fasta->description.c_str());
    printf("%s\n", otu_fasta->sequence.c_str());

    return true;
}
