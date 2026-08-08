#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <getopt.h>
#include <limits.h>

static void usage(const char* progname)
{
    fprintf(stderr,
            "====== RIsearch1 ver 1.2 ======\n= RNA-RNA interaction search =\n");
    fprintf(stderr,
            "=   Contact: wenzel@rth.dk   =\n==============================\n\n");
    fprintf(stderr, "Usage: \t%s [ARGUMENTS]\n", progname);
    fprintf(stderr, "\n   [INPUT]\n");
    fprintf(stderr, "\t-q <file> Fasta file containing query sequence(s)\n");
    fprintf(stderr,
            "\t-t <file> Fasta file containing target sequence(s) or specify '-t -' to pass the fasta formatted input to STDIN\n");
    fprintf(stderr, "\t-Q <str>  Query sequence only, direct as string\n");
    fprintf(stderr, "\t-T <str>  Target sequence as string\n");
    fprintf(stderr,
            "   If q is given, Q will be ignored; same for t over T.\n");
    fprintf(stderr, "\n   [OPTIONS]\n");
    fprintf(stderr,
            "\t-d <int>  per-nucleotide extension penalty given in dacal/mol (recommended: 30, default: 0)\n");
    fprintf(stderr,
            "\t-s <int>  threshold for suboptimal duplexes (minimum score, NOT energy)\n");
    fprintf(stderr,
            "\t-n <int>  neighborhood (only backtrack from best position within this range - to omit many overlapping results), default is 0, backtrack all\n");
    fprintf(stderr,
            "\t-f <int> force the interaction to start and end, respectively, at the 3'and 5' end of target and end at the query's 3' end. Takes input int >0.\n");
    fprintf(stderr,
            "\t            Use a high value, e.g. 200*max(length(query),length(target))\n");
    fprintf(stderr,
            "\t            It is currently not possible to compute weighted interactions without a fixed start and end.\n");
    fprintf(stderr, "\t            This option requires -w.\n");
    fprintf(stderr, "\t-l <int>  max trace back length (default: 40)\n");
    fprintf(stderr,
            "\t-m <str>  matrix to use, t99 or t04(def) or su95 or su95_noGU or slh04_noGU \n");
    fprintf(stderr,
            "\t-w <str>  weights vector to use, CRISPR_20nt_5p_3p or noweights. \n");
    fprintf(stderr,
            "\t            Weights length mush be >= than the length of the query -1.\n");
    fprintf(stderr,
            "\t            This option requires -f. Please set -f appropriately, do not try to use a low -f value to avoid the force start. \n");
    fprintf(stderr,
            "\t            Launched with this option, RIsearch does not perform step 1 (memory optimization) of the algorithm (see Wenzel et al. 2012).\n");
    fprintf(stderr, "\t-e <num>  energy threshold, checked after backtrack\n");
    fprintf(stderr,
            "\t            An interaction is only printed if the predicted hybridization energy is lower than (or equal to) this threshold.\n");
    fprintf(stderr,
            "\t            Also the 'best hit' per query/target pair might be filtered out, appears only once if at all (and not necessarily as first).\n");
    fprintf(stderr,
            "\t-p          switch for short output, for backwards compatibility, same as -p1\n");
    fprintf(stderr,
            "\t-R         Transpose the scoring matrix");
    fprintf(stderr, "\t-p[1-3]     different shorter output modes:\n");
    fprintf(stderr,
            "\t\t-p1     one line per hit, incl. interaction string, still header for each pair (query / target)\n");
    fprintf(stderr,
            "\t\t-p2     one line per hit, tab seperated 'Qname Qbeg Qend Tname Tbeg Tend score energy'; no header; 'best' hit only once, not first\n");
    fprintf(stderr,
            "\t\t-p3     one line per pair with number of hits that would have been printed, tab seperated 'Qname Tname hit-count'\n");
    fprintf(stderr,
            "\t\t        without p (and with p1), the 'best' interaction (per pair) is always shown first, and repeated in the list of results\n");
    fprintf(stderr, "\n\n");
    exit(1);
}

/* values to be overwritten by command line parameters */
typedef struct config
{
    int transpose_matrix_flag;
    int extension_penalty;
    int all_vs_all;
    const char* seq1_file_name;
    const char* seq2_file_name;
    const char* seq1_cli;
    const char* seq2_cli;
    const char* mat_name;
    const char* pos_weights;
    int min_score;
    int doSubopt;
    double max_energy;
    int filter_e;
    int weighted_positions;
    int vicinity;
    char printShort;
    int force_start_val;
    int tblen;
} config_st;

static const char* DEFAULT_MAT_NAME = "t04";
static const char* DEFAULT_POS_WEIGHTS = "CRISPR_20nt_3p_5p";

/*TODO possibly several print styles, Vienna-like (one line, but still IA) */

static void getArgs(int argc, char* argv[], config_st* config)
{
    config->transpose_matrix_flag = 0;
    config->extension_penalty = 0; /* extension penalty; used to compute dsm */
    config->seq1_file_name = NULL;
    config->seq2_file_name = NULL;
    config->seq1_cli = NULL;
    config->seq2_cli = NULL;
    config->all_vs_all = 1; /* run all queries vs all targets */
    config->mat_name = DEFAULT_MAT_NAME;
    config->min_score = INT_MAX; /* score cutoff; if not set, only print best */
    config->doSubopt = 0; /* flag : minScore in use */
    config->max_energy = INT_MAX; /*energy cutoff, not even print 'best' if it's not lower than that! */
    config->filter_e = 0; /* flag : filterE in use */
    config->weighted_positions = 0; /*true if each position of interaction has a certain weight. */
    config->pos_weights = DEFAULT_POS_WEIGHTS; /*name of the array being used to assign weights to positions */
    config->vicinity = 0; /* to omit neighboring hits (subalignments) */
    config->printShort = 0; /* switch p to print 1 line per IA, only pos&E, not IA itself */
    config->force_start_val = -1;
    /* values used to unitialize the first column of the M matrix. If sufficiently high, can force the interaction to start at position 0 of the DNA. */
    config->tblen = 40; /* trace-back length, that many nucleotides before 'maxHit' */

    char c;
    while ((c = getopt(argc, argv, "1q:t:Q:T:d:X:m:s:e:n:w:l:f:p:R::")) != -1)
        switch (c)
        {
        case '1':
            config->all_vs_all = 0;
            break;
        case 'q':
            config->seq1_file_name = optarg;
            break;
        case 'R':
            config->transpose_matrix_flag = 1;
            break;
        case 't':
            config->seq2_file_name = optarg;
            break;
        case 'Q':
            config->seq1_cli = optarg;
            break;
        case 'T':
            config->seq2_cli = optarg;
            break;
        case 'd':
            config->extension_penalty = atoi(optarg);
            break;
        case 'X':
            break; /*silent var */
        case 'm':
            config->mat_name = optarg;
            break;
        case 's':
            config->min_score = atoi(optarg);
            config->doSubopt = 1;
            break;
        case 'e':
            config->max_energy = atof(optarg);
            config->filter_e = 1;
            break;
        case 'n':
            config->vicinity = atoi(optarg);
            break;
        case 'w':
            config->weighted_positions = 1;
            config->pos_weights = optarg;
            break;
        case 'l':
            config->tblen = atoi(optarg);
            break;
        case 'f':
            config->force_start_val = atoi(optarg);
            break;
        case 'p':
            if (optarg)
                config->printShort = atoi(optarg);
            else
                config->printShort = 1;
            break;
        case '?':
            usage(argv[0]);
            break;
        default:
            abort();
        }

    if (!((config->seq1_file_name || config->seq1_cli)
        && (config->seq2_file_name || config->seq2_cli)))
    {
        fprintf(stderr,
                "\nYou need to provide a query (see -Q or -q option) and a target (-T/-t)\n\n");
        usage(argv[0]);
    }
    if (config->seq1_file_name && (!strcmp(config->seq1_file_name, "-")))
    {
        fprintf(stderr,
                "\nQuery can currently not be read from STDIN, only target can!\n\n");
        usage(argv[0]);
    }
}
