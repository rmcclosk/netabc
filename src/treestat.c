#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>

#include "tree.h"
#include "treestats.h"

typedef enum {
    TREESTAT_HEIGHT,
    TREESTAT_NTIP,
    TREESTAT_SACKIN,
    TREESTAT_COLLESS,
    TREESTAT_COPHENETIC,
    TREESTAT_LADDER_LENGTH,
    TREESTAT_IL_NODES,
    TREESTAT_WIDTH,
    TREESTAT_BMI,
    TREESTAT_MAX_DELTA_WIDTH,
    TREESTAT_CHERRIES,
    TREESTAT_PROP_UNBALANCED,
    TREESTAT_AVG_UNBALANCE,
    TREESTAT_GAMMA
} tree_statistic;

struct treestat_options {
    tree_statistic stat;
    int ladderize;
    int use_branch_lengths;
    int yule;
    int normalize;
    scaling scale_branches;
    FILE *tree_file;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"statistic", required_argument, 0, 's'},
    {"ignore-branch-lengths", no_argument, 0, 'i'},
    {"ladderize", no_argument, 0, 'd'},
    {"scale-branches", required_argument, 0, 'b'},
    {"yule", no_argument, 0, 'y'},
    {"normalize", no_argument, 0, 'n'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: treestat [options] [tree1] [tree2]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -s, --statistic           statistic to compute (see options below)\n");
    fprintf(stderr, "  -d, --ladderize           ladderize the tree before computing the statistic\n");
    fprintf(stderr, "  -b, --scale-branches      type of branch scaling to apply (mean/median/max/none, default none)\n");
    fprintf(stderr, "  -y, --yule                normalize by expected value under Yule model\n");
    fprintf(stderr, "                            for sackin, colless, and cophenetic\n");
    fprintf(stderr, "  -n, --normalize           divide by number of nodes or tips (depending on statistic)\n");
    fprintf(stderr, "  -i, --ignore-branch-lengths  treat all branches as if they have unit length\n\n");

    fprintf(stderr, "Available statistics:\n");
    fprintf(stderr, "  height                    tree height\n");
    fprintf(stderr, "  ntip                      number of tips\n");
    fprintf(stderr, "  sackin                    Sackin's index\n");
    fprintf(stderr, "  colless                   Colless' index\n");
    fprintf(stderr, "  cophenetic                total cophenetic index\n");
    fprintf(stderr, "  ladder                    maximum ladder length\n");
    fprintf(stderr, "  il                        number of 'IL' nodes\n");
    fprintf(stderr, "  width                     tree width\n");
    fprintf(stderr, "  bmi                       width / height\n");
    fprintf(stderr, "  max-delta-width           maximum width difference\n");
    fprintf(stderr, "  cherries                  number of cherries\n");
    fprintf(stderr, "  prop-unbalanced           proportion of unbalanced subtrees\n");
    fprintf(stderr, "  unbalance                 average unbalance\n");
    fprintf(stderr, "  gamma                     Pybus & Harvey's gamma statistic\n");
}

struct treestat_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct treestat_options opts = {
        .stat = TREESTAT_NTIP,
        .ladderize = 0,
        .use_branch_lengths = 1,
        .yule = 0,
        .normalize = 0,
        .scale_branches = NONE,
        .tree_file = stdin
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hs:k:nidb:y", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 's':
                if (strcmp(optarg, "height") == 0) {
                    opts.stat = TREESTAT_HEIGHT;
                }
                else if (strcmp(optarg, "sackin") == 0) {
                    opts.stat = TREESTAT_SACKIN;
                }
                else if (strcmp(optarg, "colless") == 0) {
                    opts.stat = TREESTAT_COLLESS;
                }
                else if (strcmp(optarg, "ntip") == 0) {
                    opts.stat = TREESTAT_NTIP;
                }
                else if (strcmp(optarg, "cophenetic") == 0) {
                    opts.stat = TREESTAT_COPHENETIC;
                }
                else if (strcmp(optarg, "ladder") == 0) {
                    opts.stat = TREESTAT_LADDER_LENGTH;
                }
                else if (strcmp(optarg, "il") == 0) {
                    opts.stat = TREESTAT_IL_NODES;
                }
                else if (strcmp(optarg, "width") == 0) {
                    opts.stat = TREESTAT_WIDTH;
                }
                else if (strcmp(optarg, "bmi") == 0) {
                    opts.stat = TREESTAT_BMI;
                }
                else if (strcmp(optarg, "max-delta-width") == 0) {
                    opts.stat = TREESTAT_MAX_DELTA_WIDTH;
                }
                else if (strcmp(optarg, "cherries") == 0) {
                    opts.stat = TREESTAT_CHERRIES;
                }
                else if (strcmp(optarg, "prop-unbalanced") == 0) {
                    opts.stat = TREESTAT_PROP_UNBALANCED;
                }
                else if (strcmp(optarg, "unbalance") == 0) {
                    opts.stat = TREESTAT_AVG_UNBALANCE;
                }
                else if (strcmp(optarg, "gamma") == 0) {
                    opts.stat = TREESTAT_GAMMA;
                }
                else {
                    fprintf(stderr, "Unrecognized tree statistic \"%s\"\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'i':
                opts.use_branch_lengths = 0;
                break;
            case 'b':
                if (strcmp(optarg, "mean") == 0) {
                    opts.scale_branches = MEAN;
                }
                else if (strcmp(optarg, "median") == 0) {
                    opts.scale_branches = MEDIAN;
                }
                else if (strcmp(optarg, "max") == 0) {
                    opts.scale_branches = MAX;
                }
                break;
            case 'n':
                opts.normalize = 1;
                break;
            case 'y':
                opts.yule = 1;
                break;
            case '?':
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // TODO: safety
    if (optind < argc) {
        opts.tree_file = fopen(argv[optind++], "r");
    }
    return opts;
}

int main (int argc, char **argv)
{
    int i;
    double s = 0;
    struct treestat_options opts = get_options(argc, argv);

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_t *t = parse_newick(opts.tree_file);

    if (opts.ladderize) {
        ladderize(t);
    }

    if (!opts.use_branch_lengths) {
        for (i = 0; i < igraph_ecount(t); ++i) {
            SETEAN(t, "length", i, 1.0);
        }
    }
    else {
        scale_branches(t, opts.scale_branches);
    }

    switch (opts.stat) 
    {
        case TREESTAT_NTIP:
            s = (igraph_vcount(t) + 1) / 2;
            break;
        case TREESTAT_HEIGHT:
            s = height(t);
            if (opts.normalize) {
                s /= NTIP(t);
            }
            break;
        case TREESTAT_SACKIN:
            s = sackin(t, 1);
            if (opts.yule) {
                s = (s - SACKIN_YULE(NTIP(t))) / NTIP(t);
            }
            break;
        case TREESTAT_COLLESS:
            s = colless(t);
            if (opts.yule) {
                s = (s - COLLESS_YULE(NTIP(t))) / NTIP(t);
            }
            break;
        case TREESTAT_COPHENETIC:
            s = cophenetic(t, 1);
            if (opts.yule) {
                s = (s - COPHENETIC_YULE(NTIP(t))) / NTIP(t);
            }
            break;
        case TREESTAT_LADDER_LENGTH:
            s = ladder_length(t);
            if (opts.normalize) {
                s /= NTIP(t);
            }
            break;
        case TREESTAT_IL_NODES:
            s = il_nodes(t);
            if (opts.normalize) {
                s /= NTIP(t) - 1;
            }
            break;
        case TREESTAT_WIDTH:
            s = width(t);
            if (opts.normalize) {
                s /= NTIP(t);
            }
            break;
        case TREESTAT_BMI:
            if (!igraph_ecount(t)) {
                s = 1;
            }
            else {
                s = width(t) / height(t);
            }
            break;
        case TREESTAT_MAX_DELTA_WIDTH:
            s = max_delta_width(t);
            if (opts.normalize) {
                s /= NTIP(t);
            }
            break;
        case TREESTAT_CHERRIES:
            s = cherries(t);
            if (opts.normalize) {
                s /= NTIP(t);
            }
            break;
        case TREESTAT_PROP_UNBALANCED:
            s = prop_unbalanced(t);
            break;
        case TREESTAT_AVG_UNBALANCE:
            s = avg_unbalance(t);
            break;
        case TREESTAT_GAMMA:
            s = pybus_gamma(t);
            break;
        default:
            fprintf(stderr, "Unrecognized tree statistic\n");
            return 1;
    }

    if (s == (int) s) {
        printf("%d\n", (int) s);
    }
    else {
        printf("%f\n", s);
    }

    igraph_destroy(t);

    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }
    return EXIT_SUCCESS;
}
