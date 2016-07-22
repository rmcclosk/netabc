/**
 * This file demonstrates how to use ABC-SMC to estimate the toy distribution 
 * used for demonstration by Sisson et al. 2007 and del Moral et al. 2012.
 */

#include <stdio.h>
#include <math.h>
#include <netabc.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

void propose(gsl_rng *rng, double *theta, const void *feedback, const void *arg)
{
    // feedback is the empirical variance of the particles
    double var = *((double *) feedback);
    *theta += gsl_ran_gaussian(rng, sqrt(2*var));
}

double proposal_density(const double *from, const double *to, const void *feedback, const void *arg)
{
    double var = * ((double *) feedback);
    return gsl_ran_gaussian_pdf(*to - *from, sqrt(2*var));
}

void sample_dataset(gsl_rng *rng, const double *theta, const void *arg, void *X)
{
    double *x = (double *) X;
    if (gsl_rng_get(rng) % 2) {
        *x = (*theta + gsl_ran_gaussian(rng, 1)); 
    }
    else {
        *x = (*theta + gsl_ran_gaussian(rng, 0.1));
    }
}

double distance(const void *x, const void *data, const void *arg)
{
    return fabs(*(double *) x - *(double *) data);
}

void feedback(const double *theta, int nparticle, void *fdbk, const void *arg)
{
    double *var = (double *) fdbk;
    *var = gsl_stats_variance(theta, 1, nparticle);
}

void destroy_dataset(void *z)
{
    return;
}

void sample_from_prior(gsl_rng *rng, double *theta, const void *arg)
{
    *theta = gsl_ran_flat(rng, -10, 10);
}

double prior_density(double *theta, const void *arg)
{
    return gsl_ran_flat_pdf(*theta, -10, 10);
}

int main (void) {
    smc_config config = {
        .nparam = 1,
        .nparticle = 1000,
        .nsample = 1,
        .ess_tolerance = 500,

        .final_epsilon = 0,
        .final_accept_rate = 0.015,
        .quality = 0.95,
        .step_tolerance = 1e-4,

        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double),

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
        .sample_dataset_arg = NULL,
        .distance_arg = NULL,
        .feedback_arg = NULL,
        .sample_from_prior_arg = NULL,
        .prior_density_arg = NULL
    };

    double y = 0;
    FILE *t = fopen("toy_trace.tsv", "w");
    abc_smc(config, 0, 1, &y, t);
    fclose(t);

    return 0;
}
