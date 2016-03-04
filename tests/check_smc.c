#include <check.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include "../src/smc.h"
#include "../src/util.h"

void toy_propose(gsl_rng *rng, double *theta, const void *feedback, const void *arg)
{
    double var = *((double *) feedback);
    *theta += gsl_ran_gaussian(rng, sqrt(2*var));
}

double toy_proposal_density(const double *from, const double *to, const void *feedback, const void *arg)
{
    double var = * ((double *) feedback);
    return gsl_ran_gaussian_pdf(*to - *from, sqrt(2*var));
}

void toy_sample_dataset(gsl_rng *rng, const double *theta, const void *data, void *X)
{
    double x;
    if (rand() % 2) {
        x = (*theta + gsl_ran_gaussian(rng, 1)); 
    }
    else {
        x = (*theta + gsl_ran_gaussian(rng, 0.1));
    }
    memcpy(X, &x, sizeof(double));
}

double toy_distance(const void *x, const void *y, const void *arg)
{
    return fabs(*(double *) x - *(double *) y);
}

void toy_feedback(const double *theta, int nparticle, void *feedback, const void *arg)
{
    double var = gsl_stats_variance(theta, 1, nparticle);
    memcpy(feedback, &var, sizeof(double));
}

void toy_destroy_dataset(void *z)
{
    return;
}

void toy_sample_from_prior(gsl_rng *rng, double *theta, const void *arg)
{
    *theta = gsl_ran_flat(rng, -10, 10);
}

double toy_prior_density(double *theta, const void *arg)
{
    gsl_ran_flat_pdf(*theta, -10, 10);
}

smc_functions toy_functions = {
    .propose = toy_propose,
    .proposal_density = toy_proposal_density,
    .sample_dataset = toy_sample_dataset,
    .distance = toy_distance,
    .feedback = toy_feedback,
    .destroy_dataset = toy_destroy_dataset,
    .sample_from_prior = toy_sample_from_prior,
    .prior_density = toy_prior_density
};

void write_tsv(char **hdr, double **data, int nrow, int ncol, const char *fn)
{
    int i, j;
    FILE *f = fopen(fn, "w");

    for (i = 0; i < ncol-1; ++i) {
        fprintf(f, "%s\t", hdr[i]);
    }
    fprintf(f, "%s\n", hdr[ncol-1]);

    for (i = 0; i < nrow; ++i) {
        for (j = 0; j < ncol-1; ++j) {
            fprintf(f, "%e\t", data[j][i]);
        } 
        fprintf(f, "%e\n", data[ncol-1][i]);
    }
    fclose(f);
}

Suite *smc_suite(void);

START_TEST (test_smc_toy)
{
    double y;
    int i;
    char *hdr[2] = {"theta", ""};
    double *data[2];
    FILE *trace = fopen("check_smc_trace.tsv", "w");

    smc_config config = {
        .nparam = 1,
        .nparticle = 10000,
        .nsample = 1,
        .ess_tolerance = 5000,
        .final_epsilon = 0.01,
        .quality = 0.95,
        .step_tolerance = 1e-9,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    fprintf(trace, "iter\tweight\ttheta\tX0\n");
    smc_result *res = abc_smc(config, toy_functions, 0, 1, (void *) &y, trace);
    fclose(trace);

    for (i = 0; i < config.nparticle; ++i) {
        ck_assert(res->theta[res->niter][i] > -10 && res->theta[res->niter][i] < 10);
    }
    data[0] = res->theta[res->niter];
    write_tsv(hdr, data, config.nparticle, 1, "check_smc_theta.tsv");

    // approx. same number of steps as in Del Moral 2012
    ck_assert(res->niter > 120 && res->niter < 140);
    data[0] = &res->epsilon[1];
    data[1] = &res->acceptance_rate[1];
    hdr[0] = "epsilon";
    hdr[1] = "acceptance_rate";
    write_tsv(hdr, data, res->niter-1, 2, "check_smc_iter.tsv");

    smc_result_free(res);
}
END_TEST

Suite *smc_suite(void)
{
    Suite *s;
    TCase *tc_io, *tc_smc;

    s = suite_create("smc");

    tc_smc = tcase_create("Core");
    tcase_add_test(tc_smc, test_smc_toy);
    tcase_set_timeout(tc_smc, 60);
    suite_add_tcase(s, tc_smc);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = smc_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
