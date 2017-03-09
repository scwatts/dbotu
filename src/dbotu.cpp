#include "dbotu.h"
#include "dbotu_opts.h"


int main(int argc, char **argv) {
    // Get and parse commandline arguments
    DbotuOptions dbotu_options = get_commandline_arguments(argc, argv);


    // Set number of threads to use
    omp_set_num_threads(dbotu_options.threads);


    // Load OTU count table from file and sum OTU counts
    OtuTable otu_table = read_otu_table_from_file(dbotu_options.input_otu_counts_fp);
    otu_table.otu_sum_totals.reserve(otu_table.otu_counts.size());
    for (auto &row_counts : otu_table.otu_counts) {
        // Sum row
        double otu_sum = 0;
        for (auto &count : row_counts) {
            otu_sum += count;
        }

        // Add row sum to otu_sum_totals vector
        otu_table.otu_sum_totals.push_back(otu_sum);
    }


    // Load FASTA records
    std::unordered_map<std::string,FastaRecord> fasta_records = read_fasta_from_file(dbotu_options.input_fasta_fp);


    // TODO: check that all OTUs have a matching FASTA record


    // Collect indexes in order of most abundant OTU; Initialise output vector and fill with 0..n where n is input_vector.size() - 1
    std::vector<long long unsigned int> otu_indices_abundance(otu_table.otu_sum_totals.size());
    std::iota(otu_indices_abundance.begin(), otu_indices_abundance.end(), 0);

    // Perform sort
    std::sort(otu_indices_abundance.begin(), otu_indices_abundance.end(), [&otu_table](size_t i1, size_t i2) { return otu_table.otu_sum_totals[i1] > otu_table.otu_sum_totals[i2]; });


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
            MergeOtu merged_otu = MergeOtu { otu_index, otu_data.table->otu_names.at(otu_index),
                                             std::vector<long long unsigned int> {otu_index},
                                             &otu_data.fasta->at(*otu_name),
                                             otu_data.table->otu_counts.at(otu_index),
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


    // Initialise value and index vectors for distance
    std::vector<MergeOtuDistancePair> merged_otu_distance_pairs;

    // Calculate distances to merge candidates
#pragma omp parallel for
    for (long long unsigned int i = 0; i < merged_otu_candidates.size(); ++i) {
        // Get pointer to merge OTU candidate
        MergeOtu *merged_otu_candidate = merged_otu_candidates[i];

        // Calculate levenshtein distance
        size_t lev_distance = lev_edit_distance(merged_otu_candidate->fasta->sequence.length(),
                                                merged_otu_candidate->fasta->sequence.c_str(),
                                                otu_fasta->sequence.length(),
                                                otu_fasta->sequence.c_str(), 0);

        // Length-adjusted Levenshtein distance
        size_t sequence_length_sum = merged_otu_candidate->fasta->sequence.length() + otu_fasta->sequence.length();
        double distance = static_cast<double>(lev_distance) / (0.5 * static_cast<double>(sequence_length_sum));

        // Add distance and index to vector
#pragma omp critical
        merged_otu_distance_pairs.push_back(MergeOtuDistancePair {distance, merged_otu_candidate});
    }

    // Get the indices of merged OTUs in order of increasing genetic distance
    sort_merged_otu_distance_pair(merged_otu_distance_pairs);


    // Statistical significance test; iterate through distances from smallest to largest distance
    for (auto &merged_otu_distance_pair : merged_otu_distance_pairs) {
        // Test distance; if we fail this then all following will fail, exit early
        double *distance = &merged_otu_distance_pair.distance;
        if (*distance > min_distance) {
            break;
        }

        // Significance test; if pass return true
        // TODO: use pointer here
        MergeOtu *merged_otu = merged_otu_distance_pair.merged_otu;
        std::vector<double> *merged_otu_counts = &merged_otu->otu_counts;
        std::vector<double> *otu_counts = &otu_data.table->otu_counts.at(otu_index);

        // Determine degress of freedom for chi^2 distribution and likelyhood ratio
        unsigned int df = static_cast<unsigned int>(otu_counts->size()) - 1;
        std::vector<double> count_vector_sum = element_wise_sum(*otu_counts, *merged_otu_counts);
        double lr = -2.0 * (d_helper(count_vector_sum) - d_helper(*otu_counts) - d_helper(*merged_otu_counts));

        // Calculate p-value
        double p_value = 1 - gsl_cdf_chisq_P(lr, df);

        // If p value passes threshold hold, merge and then return true, else test next genetically
        // closest OTU
        if (p_value > max_pvalue) {
            // Merge; add counts, abundance, and name to successful candidate
            // TODO: do this in a nicer way (accessed in calling function)
            // TODO: are we wasting cycles doing an assignment here (can we sum inplace?)
            merged_otu->otu_counts = element_wise_sum(merged_otu->otu_counts, *otu_counts);
            merged_otu->abundance += static_cast<double>(otu_data.otu_indices_abundance->at(otu_index));
            merged_otu->member_count_indices.push_back(otu_index);

            // Return true
            return true;
        }
    }

    // If we reach this point then return false
    return false;
}


inline std::vector<double> element_wise_sum(std::vector<double> &vec_a, std::vector<double> &vec_b) {
    // Ensure incoming vectors are of the same length
    assert(vec_a.size() == vec_b.size());

    // Return variable
    std::vector<double> out_vec;
    out_vec.reserve(vec_a.size());

    // Sum
    for (unsigned int i = 0; i < vec_a.size(); ++i) {
        out_vec.push_back(vec_a[i] + vec_b[i]);
    }

    // Return
    return out_vec;
}

inline double d_helper(std::vector<double> &counts) {
    // Total sum of non-zero counts and total sum of log(nz_count) * nz_count)
    double total = 0;
    double log_product_total = 0;

    // Get values
    for (auto &count : counts) {
        if (count > 0.0) {
            // Log and multiply
            log_product_total += (std::log(count) * count);

            // Add count to total
            total += count;
        }
    }

    // Calculate and return
    return log_product_total - (total * std::log(total));
}
