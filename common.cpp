#include "common.h"


// Convert optarg type to double
double double_from_optarg(const char *local_optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(local_optarg);
    std::string string_double = optstring.substr(0, 8);

    for (std::string::iterator it = string_double.begin(); it != string_double.end(); ++it) {
        if (!isdigit(*it) && (*it) != '.') {
            fprintf(stderr, "This doesn't look like a usable float: %s\n", local_optarg);
            exit(1);
        }
    }

    return std::atof(string_double.c_str());
}


// Load an OTU table from file
OtuTable read_otu_table_from_file(std::string &otu_count_fp) {
    //Return variable
    OtuTable otu_table;

    // Variables for storing lines and tokens
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    std::vector<double> counts_vector;
    bool id;

    // Open file stream
    std::ifstream otu_count_fh;
    otu_count_fh.open(otu_count_fp);

    // Process header
    std::getline(otu_count_fh, line);
    line_stream.str(line);

    // Iterate header columns
    while(std::getline(line_stream, ele, '\t')) {
        // Skip the #OTU ID column (first column)
        if (ele == "#OTU ID") {
            continue;
        }

        // Store samples and increment count
        otu_table.sample_names.push_back(ele);
        ++otu_table.sample_number;
    }

    // Process sample counts, need to get OTU names first
    while(std::getline(otu_count_fh, line)) {
        // (Re)sets variables for loop
        id = true;
        line_stream.clear();

        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Grab the #OTU ID
            if (id) {
                otu_table.otu_names.push_back(ele);
                id = false;
                continue;
            }

            // Add current element to OTU count after converting to double; some OTUs may be corrected and therefore a double
            counts_vector.push_back(std::stod(ele));
        }

        ++otu_table.otu_number;
    }

    // Finally construct the OTU observation matrix directly from the vector memory
    // TODO: check if operator= assigning via copy or is it being optimised out
    otu_table.otu_counts = arma::Mat<double>(&counts_vector[0], otu_table.sample_number, otu_table.otu_number, false, true);

    return otu_table;
}


// Read FASTA sequence from file returning an unordered_map with keys as title and value as sequence
// Ported from Biopython; un-optimised
std::unordered_map<std::string,FastaRecord> read_fasta_from_file(std::string &fasta_fp) {
    // Return variable
    std::unordered_map<std::string,FastaRecord> fasta_records;

    // Variables for storing lines
    std::string line;
    std::stringstream line_stream;

    // Open file stream
    std::ifstream fasta_fh;
    fasta_fh.open(fasta_fp);

    // Skip any text before the first record (e.g. blank lines, comments)
    while (true) {
        std::getline(fasta_fh, line);
        if (line.empty()) {
            // Premature end of file, or just empty?
            fprintf(stderr, "First line of FASTA file was empty, exiting early\n");
            exit(1);
        }
        if (line[0] == '>') {
            break;
        }
    }

    // Begin parsing
    while (true) {
        // Make sure FASTA description starts with > character
        if (line[0] != '>') {
            fprintf(stderr, "Records in Fasta files should start with '>' character\n");
            exit(1);
        }

        // Parse description
        std::string description = line.substr(1, line.size());
        std::string sequence;

        // Collect all FASTA lines in current record
        while (std::getline(fasta_fh, line)) {
            if (line[0] == '>') {
                break;
            }

            // Collecting FASTA
            sequence += line;
        }

        // Create new FastaRecord and add to unordered_map
        FastaRecord fasta_record { description, sequence };
        fasta_records.emplace(description, fasta_record);

        // If line empty, then it's _likely_ EOF and we're done here
        if (line.empty()) {
            // Finish iteration and return
            return fasta_records;
        }
    }
}
