#include "dbotu_opts.h"


void print_help() {
    fprintf(stderr, "Program: dbOTU3 (c++ implementation)\n");
    fprintf(stderr, "Version 0.1a\n");
    fprintf(stderr, "C++ Implementation: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
    fprintf(stderr, "Original Author: Scott Olesen (swo@alum.mit.edu)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  dbotu3 [options] --input_otu_table <if> --output_otu_table <of> --output_membership <mf>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t <if>, --input_otu_table <if>\n");
    fprintf(stderr, "                Input OTU count table (absolute values in BIOM TSV format)\n");
    fprintf(stderr, "  -f <ff>, --input_fasta <if>\n");
    fprintf(stderr, "                Input FASTA file of representative OTU sequences\n");
    fprintf(stderr, "  -o <of>, --output_otu_table <of>\n");
    fprintf(stderr, "                Output merged OTU counts file\n");
    fprintf(stderr, "  -m <mf>, --output_membership <mf>\n");
    fprintf(stderr, "                Output membership file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -d <float>, --distance <float>\n");
    fprintf(stderr, "                Minimum distance for merging (default 0.1)\n");
    fprintf(stderr, "  -a <float>, --abundance <float>\n");
    fprintf(stderr, "                Minimum fold abundance for merging (default: 10.0)\n");
    fprintf(stderr, "  -p <float>, --pvalue <float>\n");
    fprintf(stderr, "                Minimum p value for merging (default: 0.0005)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -h        --help\n");
    fprintf(stderr, "                Display this help and exit\n");
    fprintf(stderr, "  -v        --version\n");
    fprintf(stderr, "                Display version information and exit\n");
}


void print_version() {
    fprintf(stderr, "Program: dbOTU3 (c++ implementation)\n");
    fprintf(stderr, "Version 0.1a\n");
    fprintf(stderr, "C++ Implementation: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
    fprintf(stderr, "Original dbOTU3 Author: Scott Olesen (swo@alum.mit.edu)\n");
}


DbotuOptions get_commandline_arguments(int argc, char **argv) {
    // Get instance of DbotuOptionsOptions
    DbotuOptions dbotu_options;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"input_otu_table", required_argument, NULL, 't'},
            {"input_fasta", required_argument, NULL, 'f'},
            {"output_otu_table", required_argument, NULL, 'o'},
            {"output_membership", required_argument, NULL, 'm'},
            {"distance", required_argument, NULL, 'd'},
            {"abundance", required_argument, NULL, 'a'},
            {"pvalue", required_argument, NULL, 'p'},
            {"version", no_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Parse commandline arguments
    while (1) {
        // Parser variables
        int option_index = 0;
        int c;

        // Parser
        c = getopt_long(argc, argv, "hvt:f:o:m:d:a:p:", long_options, &option_index);

        // If no more arguments to parse, break
        if (c == -1) {
            break;
        }

        // Process current arguments
        switch(c) {
            case 't':
                dbotu_options.input_otu_counts_fp = optarg;
                break;
            case 'f':
                dbotu_options.input_fasta_fp = optarg;
                break;
            case 'o':
                dbotu_options.output_otu_counts_fp = optarg;
                break;
            case 'm':
                dbotu_options.output_membership_fp = optarg;
                break;
            case 'd':
                dbotu_options.min_distance = double_from_optarg(optarg);
                break;
            case 'a':
                dbotu_options.min_abundance = double_from_optarg(optarg);
                break;
            case 'p':
                dbotu_options.min_pvalue = double_from_optarg(optarg);
                break;
            case 'v':
                print_version();
                exit(0);
            case 'h':
                print_help();
                exit(0);
            default:
                exit(1);
        }
    }


    // Check if have an attempt at arguments
    if (argc < 9) {
        print_help();
        fprintf(stderr,"\n%s: error: option -t/--input_otu_table, -f/--input_fasta, -o/--output_otu_table, and -m/--output_membership are required\n", argv[0]);
        exit(1);
    }

    // Abort execution if given unknown arguments
    if (optind < argc){
        print_help();
        fprintf(stderr, "\n%s: invalid argument: %s\n", argv[0], argv[optind++]);
    }


    // Make sure we have filenames; Input OTU table
    if (dbotu_options.input_otu_counts_fp.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -t/--input_otu_table is required\n", argv[0]);
        exit(1);
    }

    // Input FASTA file
    if (dbotu_options.input_fasta_fp.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -f/--input_fasta is required\n", argv[0]);
        exit(1);
    }


    // Check that the input OTU count file exists; Input OTU table
    std::ifstream check_fh;
    check_fh.open(dbotu_options.input_otu_counts_fp);
    if (!check_fh.good()) {
        print_help();
        fprintf(stderr, "\n%s: error: input OTU count table %s does not exist\n", argv[0], dbotu_options.input_otu_counts_fp.c_str());
        exit(1);
    }
    check_fh.close();

    // Input FASTA file
    check_fh.open(dbotu_options.input_fasta_fp);
    if (!check_fh.good()) {
        print_help();
        fprintf(stderr, "\n%s: error: input FASTA file %s does not exist\n", argv[0], dbotu_options.input_fasta_fp.c_str());
        exit(1);
    }


    // Float comparison done below should be okay for our purposes
    // Ensure distance threshold is greater or equal to zero
    if (dbotu_options.min_distance < 0.0f) {
        print_help();
        fprintf(stderr,"\n%s: error: -d/--distance must be greater or equal to 0\n", argv[0]);
        exit(1);
    }

    // Ensure abundance threshold is greater or equal to zero
    if (dbotu_options.min_abundance < 0.0f) {
        print_help();
        fprintf(stderr,"\n%s: error: -a/--abundance must be greater or equal to 0\n", argv[0]);
        exit(1);
    }

    // Ensure p value threshold is between [0.0, 1.0]
    if (dbotu_options.min_pvalue < 0.0f || dbotu_options.min_pvalue > 1.0f) {
        print_help();
        fprintf(stderr,"\n%s: error: -p/--pvalue must be between 0.0 and 1.0\n", argv[0]);
        exit(1);
    }


    return dbotu_options;
}
