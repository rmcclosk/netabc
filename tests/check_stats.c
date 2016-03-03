#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/stats.h"
#include "../src/util.h"

Suite *stats_suite(void);

START_TEST (test_sample_distribution)
{
    distribution dist[3] = {UNIFORM, DELTA, GAUSSIAN};
    double *params[3];
    double theta[3];
    gsl_rng *rng = set_seed(0);
    int i;

    for (i = 0; i < 3; ++i) {
        params[i] = malloc(2*sizeof(double));
    }
    params[0][0] = 0;
    params[0][1] = 1;
    params[1][0] = 2;
    params[2][0] = 3;
    params[2][1] = 1;

    for (i = 0; i < 10; ++i) {
        sample_distribution(3, theta, rng, dist, params);
        ck_assert(theta[0] >= 0 && theta[0] <= 1);
        ck_assert(theta[1] == 2);
        ck_assert(theta[2] >= 0 && theta[2] <= 6);
    }
}
END_TEST

START_TEST (test_density_distribution)
{
    distribution dist[3] = {UNIFORM, DELTA, GAUSSIAN};
    double *params[3];
    double theta[3] = {1, 2, 3};
    int i;

    for (i = 0; i < 3; ++i) {
        params[i] = malloc(2*sizeof(double));
    }
    params[0][0] = 0;
    params[0][1] = 2;
    params[1][0] = 2;
    params[2][0] = 3;
    params[2][1] = 1;

    ck_assert(fabs(density_distribution(3, theta, dist, (double **) params) - 0.1994711) < 1e-5);
    theta[0] = -1;
    ck_assert(density_distribution(3, theta, dist, (double **) params) == 0);
    theta[0] = 1; theta[1] = 0;
    ck_assert(density_distribution(3, theta, dist, (double **) params) == 0);
    theta[1] = 2; theta[2] = 2;
    ck_assert(fabs(density_distribution(3, theta, dist, (double **) params) - 0.1209854) < 1e-5);
}
END_TEST

Suite *stats_suite(void)
{
    Suite *s;
    TCase *tc_core;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = suite_create("stats");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_sample_distribution);
    tcase_add_test(tc_core, test_density_distribution);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = stats_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
