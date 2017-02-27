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
        // Determine minimum abundance to merge current OTU
        double otu_abundance = static_cast<double>(otu_data.table->otu_sum_totals.at(otu_index));
        double otu_min_abundance = otu_abundance * dbotu_options.min_abundance;

        // If not OTU merged (function call with boolean retval), then create new OTU group
        if (!merge_otu(otu_index, otu_min_abundance, dbotu_options.min_distance, dbotu_options.max_pvalue, otu_data)) {
            // Create new OTU
            // TODO: refactor to reduce code re-use
            std::string *otu_name = &otu_data.table->otu_names[otu_index];
            MergeOtu merged_otu = MergeOtu { otu_index, std::vector<long long unsigned int> {otu_index},
                                             &otu_data.fasta->at(*otu_name), otu_data.table->otu_counts.col(otu_index),
                                             otu_abundance };
            otu_data.merged_otus->push_back(merged_otu);
        }

    }


    // Write out results
    write_otu_table_to_file(merged_otus, otu_data, dbotu_options.output_otu_counts_fp);
    write_merged_otu_members_to_file(merged_otus, otu_data, dbotu_options.output_membership_fp);


    return 0;
}


bool merge_otu(long long unsigned int &otu_index, double &otu_min_abundance, double &min_distance, double &max_pvalue, OtuData &otu_data) {
    // Get a pointer to OTU name, FASTA, and counts
    std::string *otu_name = &otu_data.table->otu_names[otu_index];
    FastaRecord *otu_fasta = &otu_data.fasta->at(*otu_name);

    // Running list of candidate OTUs to merge into
    // TODO: Use a list of pointers rather than copies of the MergeOtu struct
    std::vector<MergeOtu*> merged_otu_candidates;

    // Discover all OTUs which pass the abundance test
    for (auto &merged_otu : *otu_data.merged_otus) {
        if (merged_otu.abundance > otu_min_abundance) {
            merged_otu_candidates.push_back(&merged_otu);
        }
    }

    // If not candidates are found, return
    if (merged_otu_candidates.empty()) {
        return false;
    }

    // Order candidate merge OTUs by increasing genetic distance
    std::vector<double> merged_otu_candidate_distances;
    for (auto merged_otu_candidate : merged_otu_candidates) {
        // Calculate levenshtein distance
        size_t lev_distance = lev_edit_distance(merged_otu_candidate->fasta->sequence.length(),
                                                merged_otu_candidate->fasta->sequence.c_str(),
                                                otu_fasta->sequence.length(),
                                                otu_fasta->sequence.c_str(), 0);

        // Length-adjusted Levenshtein distance
        size_t sequence_length_sum = merged_otu_candidate->fasta->sequence.length() + otu_fasta->sequence.length();
        double distance = static_cast<double>(lev_distance) / (0.5 * static_cast<double>(sequence_length_sum));

        // Add distance to vector
        merged_otu_candidate_distances.push_back(distance);
    }

    // Get the indices of merged OTUs in order of increasing genetic distance
    std::vector<unsigned int> merged_otu_indices(merged_otu_candidate_distances.size());
    std::iota(merged_otu_indices.begin(), merged_otu_indices.end(), 0);
    std::sort(merged_otu_indices.begin(), merged_otu_indices.end(), [&merged_otu_candidate_distances](size_t i1, size_t i2) { return merged_otu_candidate_distances[i1] < merged_otu_candidate_distances[i2];});


    // Statistical significance test; iterate through distances from smallest to largest distance
    for (auto &merged_otu_index : merged_otu_indices) {
        // Test distance; if we fail this then all following will fail, exit early
        double *distance = &merged_otu_candidate_distances[merged_otu_index];
        if (*distance > min_distance) {
            break;
        }

        // Significance test; if pass return true
        // TODO: use pointer here
        MergeOtu *merged_otu = merged_otu_candidates[merged_otu_index];
        arma::Col<double> merged_otu_counts = merged_otu->otu_counts;
        arma::Col<double> otu_counts = otu_data.table->otu_counts.col(otu_index);

        // Determine degress of freedom for chi^2 distribution and likelyhood ratio
        unsigned int df = static_cast<unsigned int>(otu_counts.n_elem) - 1;
        double lr = -2.0 * (d_helper(otu_counts + merged_otu_counts) - d_helper(otu_counts) - d_helper(merged_otu_counts));

        // Calculate p-value
        double p_value = 1 - gsl_cdf_chisq_P(lr, df);

        // If p value passes threshold hold, merge and then return true, else test next genetically
        // closest OTU
        if (p_value > max_pvalue) {
            // Merge; add counts, abundance, and name to successful candidate
            // TODO: do this in a nicer way (accessed in calling function)
            merged_otu->otu_counts += otu_counts;
            merged_otu->abundance += static_cast<double>(otu_data.otu_indices_abundance->at(otu_index));
            merged_otu->member_count_indices.push_back(otu_index);

            // Return true
            return true;
        }
    }

    // If we reach this point then return false
    return false;
}


inline double d_helper(arma::Col<double> counts) {
    // Get non-zero counts and calculate
    arma::Col<double> nz_counts = counts.elem(arma::find(counts > 0.0));
    return arma::sum(nz_counts % arma::log(nz_counts)) - (arma::sum(nz_counts) * std::log(arma::sum(nz_counts)));
}
