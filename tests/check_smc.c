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

void plot(const double *x, const double *y, int nx, const char *pdf, 
          const char *plot_cmd)
{
    int i, fd;
    char fn[BUFSIZ], cmd[BUFSIZ];

    sprintf(fn, "%s/smcXXXXXX", P_tmpdir);
    fd = mkstemp(fn);
    for (i = 0; i < nx; ++i) {
        if (y) {
            dprintf(fd, "%e\t%e\n", x[i], y[i]);
        } 
        else {
            dprintf(fd, "%e\n", x[i]);
        }
    }
    close(fd);
    
    sprintf(cmd, "R --vanilla --silent -e 'd <- read.table(\"%s\")' ", fn);
    sprintf(cmd, "%s -e 'pdf(\"%s\")' ", cmd, pdf);
    sprintf(cmd, "%s -e '%s' ", cmd, plot_cmd);
    sprintf(cmd, "%s -e 'dev.off()'\n", cmd);
    i = system(cmd);

    unlink(fn);
}

Suite *smc_suite(void);

START_TEST (test_smc_toy)
{
    double y = 0;
    int i, fd;

    smc_config config = {
        .nparam = 1,
        .nparticle = 10000,
        .nsample = 5,
        .ess_tolerance = 5000,
        .final_epsilon = 0.01,
        .quality = 0.95,
        .step_tolerance = 1e-9,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    smc_result *res = abc_smc(config, toy_functions, 0, 1, (void *) &y, NULL);

    for (i = 0; i < config.nparticle; ++i) {
        ck_assert(res->theta[res->niter][i] > -10 && res->theta[res->niter][i] < 10);
    }

    plot(res->theta[res->niter], NULL, config.nparticle, "check_smc_hist.pdf",
         "plot(density(d[,1]), xlim=c(-3, 3), ylim=c(0, 2.5), xlab=\"theta\", ylab=\"density\", main=NA); polygon(density(d[,1]), col=\"gray\"); x <- seq(-3, 3, 0.01); lines(x, 0.5*dnorm(x, sd=1) + 0.5*dnorm(x, sd=0.1))");
    smc_result_free(res);
}
END_TEST

START_TEST (test_smc_toy_steps)
{
    int i;
    double y = 0;
    smc_config toy_config = {
        .nparam = 1,
        .nparticle = 10000,
        .nsample = 1,
        .ess_tolerance = 5000,
        .final_epsilon = 0.01,
        .quality = 0.95,
        .step_tolerance = 1e-5,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    smc_result *res = abc_smc(toy_config, toy_functions, 0, 8, (void *) &y, NULL);

    // same number of steps as in Del Moral 2012
    ck_assert(res->niter > 120 && res->niter < 140);
    plot(&res->epsilon[1], NULL, res->niter-1, "check_smc_epsilon.pdf", 
         "plot(d[,1], xlab=\"time index\", ylab=\"epsilon\", main=NA, type=\"l\")");
    plot(&res->acceptance_rate[1], NULL, res->niter-1, "check_smc_accept.pdf", 
         "plot(d[,1], xlab=\"time index\", ylab=\"acceptance rate\", main=NA, type=\"l\")");
    smc_result_free(res);
}
END_TEST

START_TEST (test_smc_distance)
{
    double y = 0;
    int i, fd;
    char fn[] = "/tmp/test_smc_distance", cmd[BUFSIZ];
    FILE *trace = fopen(fn, "w+");

    smc_config config = {
        .nparam = 1,
        .nparticle = 1000,
        .nsample = 2,
        .ess_tolerance = 500,
        .final_epsilon = 0.0,
        .final_accept_rate = 0.015,
        .quality = 0.95,
        .step_tolerance = 1e-9,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    fprintf(trace, "iter\tweight\ttheta\tX0\tX1\n");
    smc_result *res = abc_smc(config, toy_functions, 0, 1, (void *) &y, trace);
    fflush(trace);
    fseek(trace, 0L, SEEK_SET);

    sprintf(cmd, "R --vanilla --silent -e 'd <- read.table(\"%s\", header=TRUE)' ", fn);
    sprintf(cmd, "%s -e 'library(ggplot2) ' ", cmd);
    sprintf(cmd, "%s -e 'd <- subset(d, iter %%%% 11 == 0) ' ", cmd);
    sprintf(cmd, "%s -e 'pdf(\"check_smc_distance.pdf\")' ", cmd);
    sprintf(cmd, "%s -e 'ggplot(d, aes(x=abs(theta), y=X0)) + geom_point() + facet_wrap(~iter) + theme_bw()' ", cmd);
    sprintf(cmd, "%s -e 'dev.off()'\n", cmd);
    i = system(cmd);

    fclose(trace);
    unlink(fn);
}
END_TEST

Suite *smc_suite(void)
{
    Suite *s;
    TCase *tc_io, *tc_smc;

    s = suite_create("smc");

    tc_smc = tcase_create("Core");
    tcase_add_test(tc_smc, test_smc_toy);
    tcase_add_test(tc_smc, test_smc_toy_steps);
    tcase_add_test(tc_smc, test_smc_distance);
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
