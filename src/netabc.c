/** \file netabc.c
 * \brief Main program for the netabc binary.
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <yaml.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>

#include "tree.h"
#include "treestats.h"
#include "simulate.h"
#include "smc.h"
#include "util.h"
#include "stats.h"

#define MAX_PARAMS 7
#define MAX_REJECTIONS 100
#ifndef INFINITY
#define INFINITY DBL_MAX
#endif

/** Available contact network models. */
typedef enum {
    NET_TYPE_PA = 0,
    NET_TYPE_GNP = 1,
    NET_TYPE_SMALLWORLD = 2
} net_type;

/** 
 * Number of parameters each network model has. Must be kept in sync with
 * net_type enum.
 */
static const int NUM_PARAMS[3] = {7, 6, 7};

/** 
 * All the network parameters. Parameters for individual networks come after
 * the UNIVERSAL ones and must be in sequential order.
 */
typedef enum {
    UNIVERSAL_N = 0,
    UNIVERSAL_I = 1,
    UNIVERSAL_TIME = 2,
    UNIVERSAL_TRANSMIT_RATE = 3,
    UNIVERSAL_REMOVE_RATE = 4,

    // for preferential attachment networks
    PA_M = 5,
    PA_ALPHA = 6,

    // for random networks
    GNP_P = 5,

    // for smallworld networks
    SMALLWORLD_NEI = 5,
    SMALLWORLD_P = 6,
} net_parameter;

/**
 * Names for the network parameters. Must be kept in sync with the
 * net_parameter and net_type enums.
 */
static char *PARAM_NAMES[3][7] = {
    {"N", "I", "time", "transmit_rate", "remove_rate", "m", "alpha"},
    {"N", "I", "time", "transmit_rate", "remove_rate", "p"},
    {"N", "I", "time", "transmit_rate", "remove_rate", "nei", "p"}
};

static const char ZEROES[BUFSIZ] = {0};

/*******************************************************************************
 * Command line options.                                                       *
 ******************************************************************************/

struct netabc_options {
    FILE *tree_file;
    FILE *yaml_file;
    FILE *trace_file;
    int nthread;
    int seed;
    int nparticle;
    int nsample;
    int nltt;
    net_type net;
    double decay_factor;
    double rbf_variance;
    double quality;
    double final_epsilon;
    double final_accept_rate;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"num-threads", required_argument, 0, 't'},
    {"seed", required_argument, 0, 's'},
    {"decay-factor", required_argument, 0, 'l'},
    {"rbf-variance", required_argument, 0, 'g'},
    {"nltt", no_argument, 0, 'c'},
    {"num-particles", required_argument, 0, 'n'},
    {"num-samples", required_argument, 0, 'p'},
    {"quality", required_argument, 0, 'q'},
    {"trace", required_argument, 0, 'd'},
    {"net-type", required_argument, 0, 'm'},
    {"final-epsilon", required_argument, 0, 'e'},
    {"final-accept-rate", required_argument, 0, 'a'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: netabc [options] [tree] [yaml_file]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -t, --num-threads         number of threads\n");
    fprintf(stderr, "  -s, --seed                random seed\n");
    fprintf(stderr, "  -l, --decay-factor        decay factor for tree kernel\n");
    fprintf(stderr, "  -g, --rbf-variance        variance for tree kernel radial basis function\n");
    fprintf(stderr, "  -c, --nltt                multiply tree kernel by nLTT statistic\n");
    fprintf(stderr, "  -n, --num-particles       number of particles for SMC\n");
    fprintf(stderr, "  -p, --num-samples         number of sampled datasets per particle\n");
    fprintf(stderr, "  -q, --quality             tradeoff between speed and accuracy (0.9=fast, 0.99=accurate)\n");
    fprintf(stderr, "  -d, --trace               write population and weights at each iteration to this file\n");
    fprintf(stderr, "  -m, --net-type            type of network (pa/gnp/smallworld)\n");
    fprintf(stderr, "  -e, --final-epsilon       last epsilon for SMC\n");
    fprintf(stderr, "  -a, --final-accept-rate   stop when particle acceptance rate drops below this value\n");
}

struct netabc_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct netabc_options opts = {
        .tree_file = stdin,
        .yaml_file = NULL,
        .trace_file = NULL,
        .nthread = 1,
        .seed = -1,
        .nparticle = 1000,
        .nsample = 5,
        .net = NET_TYPE_PA,
        .decay_factor = 0.3,
        .rbf_variance = 4,
        .nltt = 0,
        .quality = 0.95,
        .final_epsilon = 0.01,
        .final_accept_rate = 0.015
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "ha:cd:e:g:l:m:n:p:q:s:t:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'a':
                opts.final_accept_rate = atof(optarg);
                break;
            case 'c':
                opts.nltt = 1;
                break;
            case 'd':
                opts.trace_file = fopen(optarg, "a");
                break;
            case 'e':
                opts.final_epsilon = atof(optarg);
                break;
            case 'g':
                opts.rbf_variance = atof(optarg);
                break;
            case 'l':
                opts.decay_factor = atof(optarg);
                break;
            case 'n':
                opts.nparticle = atoi(optarg);
                break;
            case 'm':
                if (strcmp(optarg, "pa") == 0) {
                    opts.net = NET_TYPE_PA;
                }
                else if (strcmp(optarg, "gnp") == 0) {
                    opts.net = NET_TYPE_GNP;
                }
                else if (strcmp(optarg, "smallworld") == 0) {
                    opts.net = NET_TYPE_SMALLWORLD;
                }
                else {
                    fprintf(stderr, "Error: unrecognized network type \"%s\"\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'p':
                opts.nsample = atoi(optarg);
                break;
            case 'q':
                opts.quality = atof(optarg);
                break;
            case 's':
                opts.seed = atoi(optarg);
                break;
            case 't':
                opts.nthread = atoi(optarg);
                break;
            case '?':
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    if (optind < argc) {
        opts.tree_file = fopen(argv[optind++], "r");
    }
    if (optind < argc) {
        opts.yaml_file = fopen(argv[optind++], "r");
    }

    return opts;
}

/*******************************************************************************
 * Priors.                                                                     *
 ******************************************************************************/

struct prior_data {
    int nparam;
    distribution *dist;
    double **params;
};

struct prior_data *read_priors(FILE *f, net_type net)
{
    int key = 1, seq = 0, seq_cur = 0, param;
    yaml_parser_t parser;
    yaml_event_t event;

    struct prior_data *pdata = malloc(sizeof(struct prior_data));

    pdata->dist = malloc(NUM_PARAMS[net] * sizeof(distribution));
    pdata->nparam = NUM_PARAMS[net];
    pdata->params = malloc(NUM_PARAMS[net] * sizeof(double *));

    for (param = 0; param < pdata->nparam; ++param) {
        pdata->dist[param] = 0;
    }

    if (f == NULL) {
        return;
    }

    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, f);
    do
    {
        yaml_parser_parse(&parser, &event);
        switch (event.type)
        {
            case YAML_SEQUENCE_START_EVENT:
                seq = 1;
                key = 0;
                seq_cur = 0;
                break;
            case YAML_SEQUENCE_END_EVENT:
                seq = 0;
                key = 1;
                break;
            case YAML_SCALAR_EVENT:
                if (key) {
                    if (strcmp(event.data.scalar.value, "N") == 0) {
                        param = UNIVERSAL_N;
                    }
                    else if (strcmp(event.data.scalar.value, "I") == 0) {
                        param = UNIVERSAL_I;
                    }
                    else if (strcmp(event.data.scalar.value, "time") == 0) {
                        param = UNIVERSAL_TIME;
                    }
                    else if (strcmp(event.data.scalar.value, "transmit_rate") == 0) {
                        param = UNIVERSAL_TRANSMIT_RATE;
                    }
                    else if (strcmp(event.data.scalar.value, "remove_rate") == 0) {
                        param = UNIVERSAL_REMOVE_RATE;
                    }
                    else if (strcmp(event.data.scalar.value, "p") == 0) {
                        param = net == NET_TYPE_GNP ? GNP_P : SMALLWORLD_P;
                    }
                    else if (strcmp(event.data.scalar.value, "m") == 0) {
                        param = PA_M;
                    }
                    else if (strcmp(event.data.scalar.value, "alpha") == 0) {
                        param = PA_ALPHA;
                    }
                    else if (strcmp(event.data.scalar.value, "nei") == 0) {
                        param = SMALLWORLD_NEI;
                    }
                    else {
                        fprintf(stderr, "Error: unrecognized parameter \"%s\" in YAML file\n",
                                event.data.scalar.value);
                        exit(EXIT_FAILURE);
                    }
                }
                else if (seq) {
                    if (seq_cur == 0) {
                        if (strcmp(event.data.scalar.value, "uniform") == 0) {
                            pdata->dist[param] = UNIFORM;
                            pdata->params[param] = malloc(2 * sizeof(double));
                        }
                        else if (strcmp(event.data.scalar.value, "gaussian") == 0) {
                            pdata->dist[param] = GAUSSIAN;
                            pdata->params[param] = malloc(2 * sizeof(double));
                        }
                        else if (strcmp(event.data.scalar.value, "delta") == 0) {
                            pdata->dist[param] = DELTA;
                            pdata->params[param] = malloc(sizeof(double));
                        }
                        else {
                            fprintf(stderr, "Error: unrecognized distribution \"%s\"\n",
                                    event.data.scalar.value);
                            exit(EXIT_FAILURE);
                        }
                    }
                    else {
                        pdata->params[param][seq_cur-1] = atof(event.data.scalar.value);
                    } 
                    ++seq_cur;
                }
                else {
                    pdata->dist[param] = DELTA;
                    pdata->params[param] = malloc(sizeof(double));
                    pdata->params[param][0] = atof(event.data.scalar.value);
                }

                if (!seq) {
                    key = !key;
                }
                break;
            default:
                break;
        }
        if (event.type != YAML_STREAM_END_EVENT) {
            yaml_event_delete(&event);
        }
    }
    while (event.type != YAML_STREAM_END_EVENT);

    yaml_event_delete(&event);
    yaml_parser_delete(&parser);

    for (param = 0; param < pdata->nparam; ++param) {
        if (pdata->dist[param] == 0) {
            fprintf(stderr, "ERROR: values for some priors are missing\n");
            exit(EXIT_FAILURE);
        }
    }
    return pdata;
}

void print_priors(struct prior_data *pdata, net_type net) 
{
    int i;

    for (i = 0; i < pdata->nparam; ++i) {
        fprintf(stderr, "%-13s ", PARAM_NAMES[net][i]);
        switch(pdata->dist[i]) {
            case UNIFORM:
                fprintf(stderr, "Uniform(%.1f, %.1f)\n", pdata->params[i][0], pdata->params[i][1]);
                break;
            case GAUSSIAN:
                fprintf(stderr, "Gaussian(%.1f, %.1f)\n", pdata->params[i][0], pdata->params[i][1]);
                break;
            case DELTA:
                fprintf(stderr, "Delta(%.1f)\n", pdata->params[i][0]);
                break;
            default:
                fprintf(stderr, "ERROR: unknown distribution type %d\n", pdata->dist[i]);
                exit(EXIT_FAILURE);
        }
    }
}

/*******************************************************************************
 * User-supplied functions for SMC.                                            *
 ******************************************************************************/

struct kernel_data {
    int ntip;
    double decay_factor;
    double rbf_variance;
    int nltt;
    net_type type;
};

void propose(gsl_rng *rng, double *theta, const void *feedback, const void *arg)
{
    int i;
    double *var = (double *) feedback;
    int nparam = * ((int *) arg);
    for (i = 0; i < nparam; ++i) {
        theta[i] += gsl_ran_gaussian(rng, sqrt(2*var[i]));
    }
}

double proposal_density(const double *from, const double *to, const void *feedback, const void *arg)
{
    int i;
    double *var = (double *) feedback;
    int nparam = * ((int *) arg);
    double p = 1;

    if (to[UNIVERSAL_I] > to[UNIVERSAL_N]) {
        return 0;
    }

    for (i = 0; i < nparam; ++i) {
        if (fabs(var[i]) > FLT_EPSILON || fabs(to[i] - from[i]) > FLT_EPSILON) {
            p *= gsl_ran_gaussian_pdf(to[i] - from[i], sqrt(2*var[i]));
        }
    }
    return p;
}

void sample_dataset(gsl_rng *rng, const double *theta, const void *arg, void *X)
{
    int i, failed = 0;
    igraph_t net, *tree = (igraph_t *) X;
    igraph_vector_t v;
    igraph_rng_t igraph_rng;
    unsigned long int igraph_seed = gsl_rng_get(rng);

    struct kernel_data *karg = (struct kernel_data *) arg;
    int ntip = karg->ntip;
    double decay_factor = karg->decay_factor;
    double rbf_variance = karg->rbf_variance;
    int nltt = karg->nltt;

    igraph_rng_init(&igraph_rng, &igraph_rngtype_mt19937);
    igraph_rng_seed(&igraph_rng, igraph_seed);

    igraph_vector_init(&v, (int) theta[UNIVERSAL_N]);

    switch (karg->type) {
        case NET_TYPE_PA:
            igraph_barabasi_game(&net, (int) theta[UNIVERSAL_N], theta[PA_ALPHA],
                    (int) theta[PA_M], NULL, 0, 1, 0,
                    IGRAPH_BARABASI_PSUMTREE, NULL, &igraph_rng);
            break;
        case NET_TYPE_GNP:
            igraph_erdos_renyi_game(&net, IGRAPH_ERDOS_RENYI_GNP, (int)
                    theta[UNIVERSAL_N], theta[GNP_P], 0, 0);
            break;
        case NET_TYPE_SMALLWORLD:
            igraph_watts_strogatz_game(&net, 1, (int) theta[UNIVERSAL_N], 
                    (int) theta[SMALLWORLD_NEI], theta[SMALLWORLD_P], 0, 0);
            break;
        default:
            fprintf(stderr, "BUG: unknown network type %d\n", karg->type);
            memset(&net, 0, sizeof(igraph_t));
            break;
    }
    igraph_to_directed(&net, IGRAPH_TO_DIRECTED_MUTUAL);
    
    igraph_vector_fill(&v, theta[UNIVERSAL_REMOVE_RATE]);
    SETVANV(&net, "remove", &v);
    for (i = 0; i < (int) theta[UNIVERSAL_N]; ++i) {
        VECTOR(v)[i] = i;
    }
    SETVANV(&net, "id", &v);

    igraph_vector_resize(&v, igraph_ecount(&net));
    igraph_vector_fill(&v, theta[UNIVERSAL_TRANSMIT_RATE]);
    SETEANV(&net, "transmit", &v);
    
    simulate_phylogeny(tree, &net, rng, theta[UNIVERSAL_TIME], theta[UNIVERSAL_I], 1);
    i = 0;
    while (igraph_vcount(tree) < (ntip - 1) / 2) {
        if (i == 20) {
            fprintf(stderr, "Too many tries to simulate a tree\n");
            failed = 1;
            break;
        }
        igraph_destroy(tree);
        simulate_phylogeny(tree, &net, rng, theta[UNIVERSAL_TIME], theta[UNIVERSAL_I], 1);
        ++i;
    }

    if (failed) {
        memset(tree, 0, sizeof(igraph_t));
    }
    else {
        subsample_tips(tree, ntip, rng);
        ladderize(tree);
        scale_branches(tree, MEAN);
        SETGAN(tree, "kernel", kernel(tree, tree, decay_factor, rbf_variance, 1));
    }
    igraph_destroy(&net);

    igraph_rng_destroy(&igraph_rng);
    igraph_vector_destroy(&v);
}

double distance(const void *x, const void *data, const void *arg)
{
    igraph_t *gx = (igraph_t *) x;
    igraph_t *gy = (igraph_t *) data;
    struct kernel_data *kdata = (struct kernel_data *) arg;
    double k, kx, ky, dist;

    if (memcmp(x, ZEROES, sizeof(igraph_t)) == 0 ||
        memcmp(data, ZEROES, sizeof(igraph_t)) == 0) {
        dist = INFINITY;
    }
    else {
        ky = GAN(gy, "kernel");
        kx = GAN(gx, "kernel");
        k = kernel(gx, gy, kdata->decay_factor, kdata->rbf_variance, 1);
        if (kdata->nltt) {
            k *= (1.0 - nLTT(gx, gy));
        }
        dist = 1.0 - k / sqrt(kx) / sqrt(ky);
    }

    return dist;
}

void feedback(const double *theta, int nparticle, void *params, const void *arg)
{
    int i;
    double *var = (double *) params;
    int nparam = *((int *) arg);
    for (i = 0; i < nparam; ++i) {
        var[i] = gsl_stats_variance(&theta[i], nparam, nparticle);
    }
}

void destroy_dataset(void *z)
{
    igraph_destroy((igraph_t *) z);
}

void sample_from_prior(gsl_rng *rng, double *theta, const void *arg)
{
    struct prior_data *parg = (struct prior_data *) arg;
    int tries = 0;

    do {
        sample_distribution(parg->nparam, theta, rng, parg->dist, parg->params);
        ++tries;
    } while (theta[UNIVERSAL_N] < theta[UNIVERSAL_I] && tries < MAX_REJECTIONS);

    if (tries == MAX_REJECTIONS) {
        fprintf(stderr, "ERROR: not enough prior density on N > I\n");
        exit(1);
    }
}

double prior_density(double *theta, const void *arg)
{
    struct prior_data *parg = (struct prior_data *) arg;
    if (theta[UNIVERSAL_N] < theta[UNIVERSAL_I]) {
        return 0;
    }
    return density_distribution(parg->nparam, theta, parg->dist, parg->params);
}

smc_functions functions = {
    .propose = propose,
    .proposal_density = proposal_density,
    .sample_dataset = sample_dataset,
    .distance = distance,
    .feedback = feedback,
    .destroy_dataset = destroy_dataset,
    .sample_from_prior = sample_from_prior,
    .prior_density = prior_density
};

smc_config config = {
    .step_tolerance = 1e-5,
    .dataset_size = sizeof(igraph_t)
};

/*******************************************************************************
 * Main.                                                                       *
 ******************************************************************************/

int main (int argc, char **argv)
{
    int i;
    struct netabc_options opts = get_options(argc, argv);
    igraph_t *tree;
    smc_result *result;

    double *darg = malloc(3 * sizeof(double));
    struct kernel_data kdata;
    struct prior_data *pdata;

#if IGRAPH_THREAD_SAFE == 0
    if (opts.nthread > 1)
    {
        fprintf(stderr, "Warning: igraph is not thread-safe\n");
        fprintf(stderr, "Disabling multithreading\n");
        fprintf(stderr, "To fix this, install a recent igraph compiled with thread-local storage\n");
        opts.nthread = 1;
    }
#endif

    // pars the priors
    pdata = read_priors(opts.yaml_file, opts.net);
    fclose(opts.yaml_file);
    fprintf(stderr, "Priors\n======\n");
    print_priors(pdata, opts.net);
    fprintf(stderr, "\n");

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    tree = parse_newick(opts.tree_file);
    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }

    ladderize(tree); 
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", kernel(tree, tree, opts.decay_factor, opts.rbf_variance, 1));

    // SMC configuration
    config.nparticle = opts.nparticle;
    config.ess_tolerance = opts.nparticle / 2;
    config.nsample = opts.nsample;
    config.quality = opts.quality;
    config.final_epsilon = opts.final_epsilon;
    config.final_accept_rate = opts.final_accept_rate;
    config.nparam = NUM_PARAMS[opts.net];
    config.feedback_size = config.nparam * sizeof(double);

    // fill and attach arguments for the user-supplied functions
    kdata.ntip = NTIP(tree);
    kdata.decay_factor = opts.decay_factor;
    kdata.rbf_variance = opts.rbf_variance;
    kdata.nltt = opts.nltt;
    kdata.type = opts.net;

    config.propose_arg = &config.nparam;
    config.proposal_density_arg = &config.nparam;
    config.sample_dataset_arg = &kdata;
    config.distance_arg = &kdata;
    config.feedback_arg = &config.nparam;
    config.sample_from_prior_arg = pdata;
    config.prior_density_arg = pdata;

    if (opts.trace_file != NULL) {
        fprintf(opts.trace_file, "iter\tweight\tN\tI\ttime\ttransmit\tremove\t");
        switch (opts.net) {
            case NET_TYPE_PA:
                fprintf(opts.trace_file, "m\talpha");
                break;
            case NET_TYPE_SMALLWORLD:
                fprintf(opts.trace_file, "nei\tp");
                break;
            case NET_TYPE_GNP:
                fprintf(opts.trace_file, "p");
                break;
            default:
                fprintf(stderr, "BUG: unrecognized network type %d\n", opts.net);
                exit(EXIT_FAILURE);
        }

        for (i = 0; i < opts.nsample; ++i) {
            fprintf(opts.trace_file, "\tX%d", i);
        }
        fprintf(opts.trace_file, "\n");
    }
    else {
        fprintf(stderr, "WARNING: no trace file specified, output will not be recorded\n");
    }

    result = abc_smc(config, functions, opts.seed, opts.nthread, tree, opts.trace_file);

    if (opts.trace_file != NULL) {
        fclose(opts.trace_file);
    }

    igraph_destroy(tree);
    smc_result_free(result);
    return EXIT_SUCCESS;
}
