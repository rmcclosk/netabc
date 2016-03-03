#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>

#include "util.h"
#include "tree.h"
#include "mmpp.h"
#include "stats.h"

struct mmpp_options {
    FILE *tree_file;
    FILE *output;
    scaling scale_branches;
    char *cmaes_settings;
    int seed;
    int trace;
    int nrates;
    int bg_states;
    int trans_at_nodes;
    int use_tips;
    double lbound_branch;
    double ubound_branch;
    double lbound_trans;
    double ubound_trans;
    model_selector ms;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"scale-branches", required_argument, 0, 'b'},
    {"seed", required_argument, 0, 's'},
    {"cmaes-settings", required_argument, 0, 'c'},
    {"trace", no_argument, 0, 't'},
    {"rates", required_argument, 0, 'r'},
    {"model-selector", required_argument, 0, 'm'},
    {"bg-states", required_argument, 0, 'l'},
    {"output", required_argument, 0, 'o'},
    {"trans-at-nodes", no_argument, 0, 'n'},
    {"ignore-tips", no_argument, 0, 'i'},
    {"lbound-branch", required_argument, 0, '1'},
    {"ubound-branch", required_argument, 0, '2'},
    {"lbound-trans", required_argument, 0, '3'},
    {"ubound-trans", required_argument, 0, '4'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: pcbr [options] [tree]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -b, --scale-branches      how to scale branches (mean/median/max/none)\n");
    fprintf(stderr, "  -s, --seed                random seed\n");
    fprintf(stderr, "  -c, --cmaes-settings      file containing settings for CMA-ES\n");
    fprintf(stderr, "  -t, --trace               show trace of tested parameter values\n");
    fprintf(stderr, "  -r, --rates               number of rates to fit\n");
    fprintf(stderr, "  -m, --model-selector      test used to select number of rates\n");
    fprintf(stderr, "                            (lrt/aic/bic)\n");
    fprintf(stderr, "  -l, --bg-states           number of lowest states to use for background\n");
    fprintf(stderr, "  -o, --output              write clustering results and rates here\n");
    fprintf(stderr, "  -n, --trans-at-nodes      assume transitions happen at nodes, not along edges\n");
    fprintf(stderr, "  -i, --ignore-tips         do not factor in terminal branches\n");
    fprintf(stderr, "  -1, --lbound-branch       lower bound for branching rates\n");
    fprintf(stderr, "  -2, --ubound-branch       upper bound for branching rates\n");
    fprintf(stderr, "  -3, --lbound-trans        lower bound for transition rates\n");
    fprintf(stderr, "  -4, --ubound-trans        upper bound for transition rates\n");
}

void display_results(int nrates, double *theta, double branch_scale)
{
    int i, j;
    int *rate_order = malloc(nrates * sizeof(int));

    order(theta, rate_order, sizeof(double), nrates, compare_doubles);

    fprintf(stderr, "rates: ");
    for (i = 0; i < nrates; ++i)
	    fprintf(stderr, "%f ", theta[rate_order[i]] * branch_scale);
    fprintf(stderr, "\n");

    fprintf(stderr, "Q: ");
    for (i = 0; i < nrates; ++i) {
        fprintf(stderr, "[ ");
        for (j = 0; j < nrates; ++j) {
            if (i == j)
                fprintf(stderr, "   *   ");
            else if (i > j)
	            fprintf(stderr, "%f", theta[nrates + rate_order[i]*(nrates-1) + rate_order[j]] * branch_scale);
            else
	            fprintf(stderr, "%f", theta[nrates + rate_order[i]*(nrates-1) + rate_order[j] - 1] * branch_scale);
        }
        if (i < nrates - 1)
            fprintf(stderr, " ]\n   ");
        else
            fprintf(stderr, " ]\n");
    }
    free(rate_order);
}

struct mmpp_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct mmpp_options opts = {
        .scale_branches = MAX,
        .tree_file = stdin,
        .output = stdout,
        .cmaes_settings = NULL,
        .seed = -1,
        .trace = 0,
        .nrates = 0,
        .bg_states = 1,
        .trans_at_nodes = 0,
        .use_tips = 1,
        .lbound_branch = 1e-10,
        .ubound_branch = 1e10,
        .lbound_trans = 1e-10,
        .ubound_trans = 1e10,
        .ms = LRT
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hs:b:c:itr:m:nl:o:1:2:3:4:", long_options, &i);

        switch (c)
        {
            case 0:
            case -1:
                break;
            case '1':
                opts.lbound_branch = atof(optarg);
                break;
            case '2':
                opts.ubound_branch = atof(optarg);
                break;
            case '3':
                opts.lbound_trans = atof(optarg);
                break;
            case '4':
                opts.ubound_trans = atof(optarg);
                break;
            case 'b':
                if (strcmp(optarg, "mean") == 0) {
                    opts.scale_branches = MEAN;
                }
                else if (strcmp(optarg, "median") == 0) {
                    opts.scale_branches = MEDIAN;
                }
                else if (strcmp(optarg, "none") == 0) {
                    opts.scale_branches = NONE;
                }
                break;
            case 'c':
                opts.cmaes_settings = optarg;
                break;
            case 'i':
                opts.use_tips = 0;
                break;
            case 'l':
                opts.bg_states = atoi(optarg);
                break;
            case 'm':
                if (strcmp(optarg, "aic") == 0) {
                    opts.ms = AIC;
                }
                else if (strcmp(optarg, "bic") == 0) {
                    opts.ms = BIC;
                }
            case 'n':
                opts.trans_at_nodes = 1;
                break;
            case 'o':
                opts.output = fopen(optarg, "w");
                break;
            case 's':
                opts.seed = atoi(optarg);
                break;
            case 't':
                opts.trace = 1;
                break;
            case 'r':
                opts.nrates = atoi(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case '?':
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    if (optind < argc)
        opts.tree_file = fopen(argv[optind++], "r");

    return opts;
}

int main (int argc, char **argv)
{
    struct mmpp_options opts = get_options(argc, argv);
    double branch_scale;
    double *theta = malloc(opts.nrates * opts.nrates * sizeof(double));
    double bounds[4] = {opts.lbound_branch, opts.ubound_branch, opts.lbound_trans, opts.ubound_trans};
    int i, error, *states, *clusters;
    igraph_t *tree;

    set_seed(opts.seed < 0 ? time(NULL) : opts.seed);
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    tree = parse_newick(opts.tree_file);
    ladderize(tree);
    branch_scale = scale_branches(tree, opts.scale_branches);
    states = malloc(igraph_vcount(tree) * sizeof(int));
    for (i = 0; i < 4; ++i)
        bounds[i] /= branch_scale;

    error = fit_mmpp(tree, &opts.nrates, &theta, opts.trace, opts.cmaes_settings,
            states, opts.ms, opts.use_tips, opts.trans_at_nodes, bounds);
    display_results(opts.nrates, theta, branch_scale);

    clusters = malloc(igraph_vcount(tree) * sizeof(int));
    get_clusters(tree, states, clusters, opts.bg_states);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        fprintf(opts.output, "%s\t%e\t%d\n", VAS(tree, "id", i),
                theta[states[i]] * branch_scale, clusters[i]);
    }

    igraph_destroy(tree);
    free(theta);
    free(states);
    free(clusters);
    if (opts.tree_file != stdin)
        fclose(opts.tree_file);
    if (opts.output != stdin)
        fclose(opts.output);
    return error;
}
