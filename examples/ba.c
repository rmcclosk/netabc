/**
 *
 * This file demonstrates how to use ABC-SMC it the Barabasi-Albert network
 * model to a phylogeny. 
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <igraph/igraph.h>
#include <netabc.h>

// model parameters
#define NPARAM 4
#define ALPHA 0
#define I 1
#define M 2
#define N 3

// kernel parameters
#define LAMBDA 0.3
#define SIGMA 4

void propose(gsl_rng *rng, double *theta, const void *feedback, const void *arg)
{
    // feedback is the empirical variance of the particles
    double *var = (double *) feedback;

    // continuous parameters have Gaussian proposals
    theta[ALPHA] += gsl_ran_gaussian(rng, sqrt(2*var[ALPHA]));
    theta[I] += gsl_ran_gaussian(rng, sqrt(2*var[I]));
    theta[N] += gsl_ran_gaussian(rng, sqrt(2*var[N]));

    // use a Poisson proposal for discrete parameters
    // use gsl_ran_flat instead of rand for thread safety
    if (gsl_ran_flat(rng, 0, 1) < 0.5) {
        theta[M] += gsl_ran_poisson(rng, sqrt(2*var[M]));
    }
    else {
        theta[M] -= gsl_ran_poisson(rng, sqrt(2*var[M]));
    }
}

double proposal_density(const double *from, const double *to, const void *feedback, const void *arg)
{
    double *var = (double *) feedback, p;

    p  = gsl_ran_gaussian_pdf(to[ALPHA] - from[ALPHA], sqrt(2*var[ALPHA]));
    p *= gsl_ran_gaussian_pdf(to[I] - from[I], sqrt(2*var[I]));
    p *= gsl_ran_poisson_pdf((int) fabs(to[M] - from[M]), sqrt(2*var[M])) / 2;
    p *= gsl_ran_gaussian_pdf(to[N] - from[N], sqrt(2*var[N]));

    return p;
}

typedef struct {
    int ntip;
    igraph_rng_t rng;
} sample_dataset_arg;

void sample_dataset(gsl_rng *rng, const double *theta, const void *arg, void *X)
{
    int i;

    // we passed in the number of tips and the RNG as the extra argument
    sample_dataset_arg *sarg = (sample_dataset_arg *) arg;

    igraph_t net, *tree = (igraph_t *) X;
    igraph_vector_t v;

    igraph_vector_init(&v, (int) theta[N]);

    // make the network
    igraph_barabasi_game(&net, (int) theta[N], theta[ALPHA], (int) theta[M],
                         NULL, 0, 1, 0, IGRAPH_BARABASI_PSUMTREE, NULL,
                         &sarg->rng);
    igraph_to_directed(&net, IGRAPH_TO_DIRECTED_MUTUAL);
    
    // set removal rate to zero
    igraph_vector_fill(&v, 0);
    SETVANV(&net, "remove", &v);

    // give each node an ID (will need this later to record which node in the
    // tree came from which node in the network)
    for (i = 0; i < (int) theta[N]; ++i) {
        VECTOR(v)[i] = i;
    }
    SETVANV(&net, "id", &v);

    // set the transmission rate to 1
    igraph_vector_resize(&v, igraph_ecount(&net));
    igraph_vector_fill(&v, 1);
    SETEANV(&net, "transmit", &v);
    
    // simulate a tree
    simulate_phylogeny(tree, &net, rng, 0, (int) theta[I], 1);
    subsample_tips(tree, sarg->ntip, rng);
    ladderize(tree);
    scale_branches(tree, MEAN);
    
    // clean up
    igraph_destroy(&net);
    igraph_vector_destroy(&v);
}

double distance(const void *x, const void *data, const void *arg)
{
    double k, kx, ky;
    igraph_t *gx = (igraph_t *) x;
    igraph_t *gy = (igraph_t *) data;

    kx = kernel(gx, gx, LAMBDA, SIGMA, 1);
    ky = kernel(gy, gy, LAMBDA, SIGMA, 1);
    k = kernel(gx, gy, LAMBDA, SIGMA, 1);
    return 1.0 - k / sqrt(kx) / sqrt(ky);
}

void feedback(const double *theta, int nparticle, void *fdbk, const void *arg)
{
    int i;
    double *var = (double *) fdbk;

    for (i = 0; i < NPARAM; ++i) {
        var[i] = gsl_stats_variance(&theta[i], NPARAM, nparticle);
    }
}

void destroy_dataset(void *z)
{
    igraph_destroy((igraph_t *) z);
}

void sample_from_prior(gsl_rng *rng, double *theta, const void *arg)
{
    do {
        theta[ALPHA] = gsl_ran_flat(rng, 0, 2);
        theta[I] = gsl_ran_flat(rng, 0, 10000);
        theta[M] = (int) gsl_ran_flat(rng, 1, 6);
        theta[N] = gsl_ran_flat(rng, 0, 10000);
    } while (theta[I] > theta[N]);
}

double prior_density(double *theta, const void *arg)
{
    double p = 1;

    // don't allow I > N
    if (theta[N] < theta[I]) {
        return 0;
    }

    p *= gsl_ran_flat_pdf(theta[ALPHA], 0, 2);
    p *= gsl_ran_flat_pdf(theta[I], 0, 10000);
    p *= gsl_ran_flat_pdf(theta[M], 1, 6);
    p *= gsl_ran_flat_pdf(theta[N], 0, 10000);
    return p;
}

int main (void) {
    igraph_rng_t rng;
    igraph_rng_init(&rng, &igraph_rngtype_mt19937);

    FILE *f = fopen("test.nwk", "r");
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_t *tree = parse_newick(f);
    fclose(f);

    sample_dataset_arg sarg = {
        .ntip = igraph_vcount(tree),
        .rng = rng
    };

    smc_config config = {
        .nparam = NPARAM,
        .nparticle = 100,
        .nsample = 1,
        .ess_tolerance = 50,

        .final_epsilon = 0,
        .final_accept_rate = 0.015,
        .quality = 0.95,
        .step_tolerance = 1e-4,

        .dataset_size = sizeof(igraph_t),
        .feedback_size = NPARAM * sizeof(double),

        .propose = propose,
        .proposal_density = proposal_density,
        .sample_dataset = sample_dataset,
        .distance = distance,
        .feedback = feedback,
        .destroy_dataset = destroy_dataset,
        .sample_from_prior = sample_from_prior,
        .prior_density= prior_density,

        .propose_arg = NULL,
        .proposal_density_arg = NULL,
        .sample_dataset_arg = &sarg,
        .distance_arg = NULL,
        .feedback_arg = NULL,
        .sample_from_prior_arg = NULL,
        .prior_density_arg = NULL
    };

    FILE *t = fopen("ba_trace.tsv", "w");
    abc_smc(config, 0, 1, tree, t);
    fclose(t);

    igraph_rng_destroy(&rng);
    return 0;
}
