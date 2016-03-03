#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "simulate.h"
#include "util.h"
#include "tree.h"

struct nettree_options {
    double sim_time;
    int sim_nodes;
    double tree_height;
    int extant_only;
    int ntip;
    double sample_baseline;
    double sample_peer;
    FILE *net_file;
    FILE *tree_file;
    int seed;

    int nsample;
    double *sample_prop;
    double *sample_time;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"sim-time", required_argument, 0, 's'},
    {"sim-nodes", required_argument, 0, 'n'},
    {"tree-height", required_argument, 0, 'r'},
    {"tree-tips", required_argument, 0, 't'},
    {"sample-baseline", required_argument, 0, 'b'},
    {"sample-peer", required_argument, 0, 'x'},
    {"extant-only", no_argument, 0, 'e'},
    {"seed", required_argument, 0, 'd'},
    {"sample-time", required_argument, 0, 'm'},
    {"sample-prop", required_argument, 0, 'p'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: nettree [options] [graph]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -s, --sim-time            amount of epidemic time to simulate (default infinity)\n");
    fprintf(stderr, "  -n, --sim-nodes           stop simulation after this many nodes are infected (default infinity)\n");
    fprintf(stderr, "  -r, --tree-height         cut the tree at this height (default infinity)\n");
    fprintf(stderr, "  -t, --tree-tips           sample this many tips from the tree (default infinity)\n");
    fprintf(stderr, "  -b, --sample-baseline     baseline sampling probability (default 1)\n");
    fprintf(stderr, "  -x, --sample-peer         additional peer-driven sampling probability (default 0)\n");
    fprintf(stderr, "  -e, --extant-only         sample only extant tips (default false)\n");
    fprintf(stderr, "  -d, --seed                random seed\n");
    fprintf(stderr, "  -m, --sample-time         sample some tips at this time\n");
    fprintf(stderr, "  -p, --sample-prop         sample this proportion of tips at some time point\n");
}

struct nettree_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct nettree_options opts = {
        .sim_time = 0, 
        .sim_nodes = 0,
        .tree_height = DBL_MAX,
        .extant_only = 0,
        .ntip = 0,
        .sample_baseline = 1,
        .sample_peer = 0,
        .net_file = stdin,
        .tree_file = stdout,
        .seed = -1,
        .nsample = 0,
        .sample_prop = NULL,
        .sample_time = NULL
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "b:hs:m:p:n:r:t:ed:x:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'b':
                opts.sample_baseline = atof(optarg);
                break;
            case 's':
                opts.sim_time = atof(optarg);
                break;
            case 'm':
                opts.nsample++;
                opts.sample_prop = safe_realloc(opts.sample_prop, opts.nsample * sizeof(double));
                opts.sample_time = safe_realloc(opts.sample_time, opts.nsample * sizeof(double));
                opts.sample_prop[opts.nsample - 1] = atof(optarg);
                break;
            case 'p':
                opts.sample_time[opts.nsample - 1] = atof(optarg);
                break;
            case 'n':
                opts.sim_nodes = atoi(optarg);
                break;
            case 'r':
                opts.tree_height = atof(optarg);
                break;
            case 't':
                opts.ntip = atoi(optarg);
                break;
            case 'e':
                opts.extant_only = 1;
                break;
            case 'd':
                opts.seed = atoi(optarg);
                break;
            case 'x':
                opts.sample_peer = atof(optarg);
                break;
            case '?':
            case 0:
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // TODO: safety
    if (optind < argc)
    {
        opts.net_file = fopen(argv[optind++], "r");
        if (optind < argc)
        {
            opts.tree_file = fopen(argv[optind++], "w");
        }
    }

    return opts;
}

int net_ok(igraph_t *net, int min_comp_size)
{
    int ok = 1;
    igraph_vector_t csize;
    igraph_vector_init(&csize, igraph_vcount(net));

    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "remove"))
    {
        fprintf(stderr, "Vertices must have the 'remove' attribute\n");
        ok = 0;
    }
    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "id"))
    {
        fprintf(stderr, "Vertices must have the 'id' attribute\n");
        ok = 0;
    }
    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_EDGE, "transmit"))
    {
        fprintf(stderr, "Edges must have the 'transmit' attribute\n");
        ok = 0;
    }

    igraph_clusters(net, NULL, &csize, NULL, IGRAPH_WEAK);
    if (igraph_vector_max(&csize) < min_comp_size)
    {
        fprintf(stderr, "Largest connected component has %d nodes, but need %d\n",
                (int) igraph_vector_max(&csize), min_comp_size);
        ok = 0;
    }
    igraph_vector_destroy(&csize);

    return ok;
}

int main (int argc, char **argv)
{
    int i, ok, numeric_ids = 0, min_comp_size = 0;
    char buf[128];
    struct nettree_options opts = get_options(argc, argv);
    gsl_rng *rng = set_seed(opts.seed);
    igraph_t net, *tree = malloc(sizeof(igraph_t));
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_t gtypes, vtypes, etypes;

    // validate command line options
    if (opts.sim_nodes > 0 && opts.ntip > opts.sim_nodes) {
        fprintf(stderr, "Cannot simulate a tree with %d tips from an epidemic on %d nodes\n",
                opts.ntip, opts.sim_nodes);
        return EXIT_FAILURE;
    }

    // parse the graph
    // NB: this will segfault on an invalid graph even if I use error handling
    // it's a bug in igraph
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_gml(&net, opts.net_file);

    // find all the attributes attached to the graph
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);
    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
    igraph_cattribute_list(&net, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);

    // check if the ids are numeric or strings
    // the accessors are different for string vs. numeric attributes, so we
    // have to know which type to use
    for (i = 0; i < igraph_strvector_size(&vnames); ++i) {
        if (strcmp(STR(vnames, i), "id") == 0 && VECTOR(vtypes)[i] == 1) {
            numeric_ids = 1;
            break;
        }
    }

    // calculate the minimum connected component size we need
    if (opts.sim_nodes > 0) {
        min_comp_size = opts.sim_nodes;
    }
    else if (opts.ntip > 0) {
        min_comp_size = opts.ntip;
    }
    else {
        min_comp_size = 1;
    }

    ok = net_ok(&net, min_comp_size);
    if (ok)
    {
        // reject until we get a tree with enough nodes
        // we already checked for a sufficiently large connected component,
        // so this should succeed eventually
        simulate_phylogeny(tree, &net, rng, opts.sim_time, opts.sim_nodes, numeric_ids);
        while (igraph_vcount(tree) < opts.ntip) {
            igraph_destroy(tree);
            simulate_phylogeny(tree, &net, rng, opts.sim_time, opts.sim_nodes, numeric_ids);
        }

        // post-process the tree
        cut_at_time(tree, opts.tree_height, opts.extant_only);
        if (opts.sample_peer == 0) {
            subsample_tips(tree, opts.ntip, rng);
        }
        else {
            subsample_tips_peerdriven(tree, &net, opts.sample_baseline,
                    opts.sample_peer, opts.ntip, rng);
        }

        if (opts.nsample > 0) {
            subsample(tree, opts.nsample, opts.sample_prop, opts.sample_time, rng);
        }

        ladderize(tree);
        write_tree_newick(tree, opts.tree_file);
        igraph_destroy(tree);
        free(tree);
    }

    // clean up
    if (opts.net_file != stdin) {
        fclose(opts.net_file);
    }
    if (opts.tree_file != stdout) {
        fclose(opts.tree_file);
    }
    gsl_rng_free(rng);
    igraph_destroy(&net);
    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);
    igraph_vector_destroy(&gtypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&etypes);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
