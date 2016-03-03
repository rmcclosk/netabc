#include <stdio.h>
#include <string.h>
#include <check.h>
#include <math.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "../src/tree.h"
#include "../src/mmpp.h"

Suite *mmpp_suite(void); 

// TODO: this is copypasta from check_tree
FILE *newick_file(const char *newick)
{
    FILE *f = tmpfile();
    fprintf(f, "%s", newick);
    fseek(f, 0, SEEK_SET);
    return f;
}

igraph_vector_t new_vector(void)
{
    igraph_vector_t v;
    igraph_vector_init(&v, 0);
    return v;
}

igraph_t *tree_from_newick(const char *newick)
{
    FILE *f = newick_file(newick);
    igraph_t *tree = parse_newick(f);
    fclose(f);
    return tree;
}

void _calculate_pi(int nrates, const double *theta, double *pi)
{
    int i, j, cur = nrates;
    double sum;

    gsl_matrix *Q = gsl_matrix_alloc(nrates, nrates);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(nrates);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(nrates, nrates);
    gsl_eigen_nonsymmv_workspace *ew = gsl_eigen_nonsymmv_alloc(nrates);

    for (i = 0; i < nrates; ++i) {
        sum = 0;
        for (j = 0; j < nrates; ++j) {
            if (i != j) {
                gsl_matrix_set(Q, i, j, theta[cur]);
                sum += theta[cur++];
            }
        }
        gsl_matrix_set(Q, i, i, -sum);
    }

    gsl_matrix_transpose(Q);
    gsl_eigen_nonsymmv(Q, eval, evec, ew);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    sum = 0;
    for (i = 0; i < nrates; ++i)
    {
        pi[i] = GSL_REAL(gsl_matrix_complex_get(evec, i, 0));
        sum += pi[i];
    }
    for (i = 0; i < nrates; ++i)
        pi[i] /= sum;

    gsl_eigen_nonsymmv_free(ew);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_free(Q);
}

double _calculate_P(int nrates, const double *theta, int pstate, int cstate,
        double t, int trans_at_nodes)
{
    gsl_matrix *Q = gsl_matrix_alloc(nrates, nrates);
    gsl_matrix *expQt = gsl_matrix_alloc(nrates, nrates);
    double sum, ans;
    int i, j, cur = nrates;

    for (i = 0; i < nrates; ++i) {
        sum = 0;
        for (j = 0; j < nrates; ++j) {
            if (i != j) {
                gsl_matrix_set(Q, i, j, t * theta[cur]);
                sum += t * theta[cur++];
            }
        }
        if (trans_at_nodes)
            gsl_matrix_set(Q, i, i, -sum);
        else
            gsl_matrix_set(Q, i, i, -sum-t*theta[i]);
    }
    gsl_linalg_exponential_ss(Q, expQt, GSL_PREC_DOUBLE);
    ans = gsl_matrix_get(expQt, pstate, cstate);
    gsl_matrix_free(Q);
    gsl_matrix_free(expQt);
    return ans;
}

double lik_states(const igraph_t *tree, int nrates, int *states, double *theta, 
        double *pi, igraph_adjlist_t *al, int trans_at_nodes)
{
    double lik = pi[states[root(tree)]], lt, rt;
    int i, lc, rc, lterm, rterm, edge;
    igraph_vector_int_t *children;

    for (i = 0; i < igraph_vcount(tree); ++i) {
        children = igraph_adjlist_get(al, i);
        if (igraph_vector_int_size(children) > 0) {
            lc = VECTOR(*children)[0];
            rc = VECTOR(*children)[1];

            igraph_get_eid(tree, &edge, i, lc, 1, 1);
            lt = EAN(tree, "length", edge);
            lik *= _calculate_P(nrates, theta, states[i], states[lc], lt, trans_at_nodes);

            igraph_get_eid(tree, &edge, i, rc, 1, 1);
            rt = EAN(tree, "length", edge);
            lik *= _calculate_P(nrates, theta, states[i], states[rc], rt, trans_at_nodes);

            children = igraph_adjlist_get(al, lc);
            lterm = igraph_vector_int_size(children) == 0;
            children = igraph_adjlist_get(al, rc);
            rterm = igraph_vector_int_size(children) == 0;

            if (trans_at_nodes) {
                lik *= exp(-theta[states[i]] * lt);
                lik *= exp(-theta[states[i]] * rt);
                if (!lterm)
                    lik *= theta[states[i]];
                if (!rterm)
                    lik *= theta[states[i]];
            }
            else {
                if (!lterm)
                    lik *= theta[states[lc]];
                if (!rterm)
                    lik *= theta[states[rc]];
            }
        }
    }
    return lik;
}

double reconstruct_long(const igraph_t *tree, int nrates, double *theta, 
        int *states, int trans_at_nodes)
{
    mmpp_workspace *w = mmpp_workspace_create(tree, nrates);
    igraph_adjlist_t al;

    int i, j, icopy;
    int *v = malloc(igraph_vcount(tree) * sizeof(int));
    double lik = 0, likmax = 0;
    double *pi = malloc(nrates * sizeof(double));

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    _calculate_pi(nrates, theta, pi);

    for (i = 0; i < (int) pow(nrates, igraph_vcount(tree)); ++i) {
        icopy = i;
        j = 0;
        memset(v, 0, igraph_vcount(tree)*sizeof(int));
        while (icopy > 0) {
            v[j++] = icopy % nrates;
            icopy /= nrates;
        }
        if (trans_at_nodes)
            lik = lik_states(tree, nrates, v, theta, pi, &al, trans_at_nodes);
        else
            lik = lik_states(tree, nrates, v, theta, pi, &al, trans_at_nodes);
        if (likmax < lik) {
            likmax = lik;
            memcpy(states, v, igraph_vcount(tree)*sizeof(int));
        }
    }

    igraph_adjlist_destroy(&al);
    mmpp_workspace_free(w);
    free(v);
}

START_TEST(test_guess_parameters)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    double theta[4];
    guess_parameters(tree, 2, theta);
    ck_assert(fabs(theta[0] - 1.0/0.88) < 1e-5);
    ck_assert(fabs(theta[1] - 1.0/0.46) < 1e-5);
    ck_assert(fabs(theta[2] - 0.007463) < 1e-5);
    ck_assert(fabs(theta[3] - 0.007463) < 1e-5);
}
END_TEST

START_TEST(test_likelihood_toy)
{
    igraph_t *tree = tree_from_newick("(t1:1,t2:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 0, 1, 0)) - 0.0883478) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_likelihood_toy3)
{
    igraph_t *tree = tree_from_newick("((t1:0.5,t2:1):0.75,t3:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 0, 1, 0)) - 0.02248663) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_likelihood_nodes_toy)
{
    igraph_t *tree = tree_from_newick("(t1:1,t2:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 1, 1, 0)) - 0.1060804) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_likelihood_nodes_toy3)
{
    igraph_t *tree = tree_from_newick("((t1:0.5,t2:1):0.75,t3:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 1, 1, 0)) - 0.02633081) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_get_clusters)
{
    igraph_t *tree = tree_from_newick("(t1:1,((t5:5,t2:2),(t4:4,t3:3)));");
    int states[9] = {0, 2, 2, 2, 1, 1, 1, 1, 0};
    int clusters[9];
    int i;

    // clustering on state 1
    get_clusters(tree, states, clusters, 1);
    ck_assert_int_eq(clusters[0], 0);
    ck_assert_int_eq(clusters[8], 0);
    for (i = 1; i < 8; ++i)
        ck_assert_int_eq(clusters[i], 1);

    // clustering on state 2
    get_clusters(tree, states, clusters, 2);
    ck_assert_int_eq(clusters[0], 0);
    for (i = 1; i < 4; ++i)
        ck_assert_int_eq(clusters[i], 1);
    for (i = 4; i < 9; ++i)
        ck_assert_int_eq(clusters[i], 0);

    igraph_destroy(tree);
}
END_TEST

START_TEST(test_reconstruct)
{
    igraph_t *tree = tree_from_newick("(t2:10.94,(t4:0.04,(t1:0.75,(t3:0.29,t5:0.8):0.39):0.21):0.8);");
    double theta[4] = {0.1, 1, 1, 2};
    int states1[9];
    int states2[9];
    int i;
    double pi[2];
    igraph_adjlist_t al;
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);

    _calculate_pi(2, theta, pi);
    igraph_adjlist_init(tree, &al, IGRAPH_OUT);

    reconstruct(tree, 2, theta, w, states1, 1, 0);
    reconstruct_long(tree, 2, theta, states2, 0);

    for (i = 0; i < 9; ++i) {
        ck_assert_int_eq(states1[i], states2[i]);
    }
    ck_assert(fabs(lik_states(tree, 2, states1, theta, pi, &al, 0)) -
                   pow(10, likelihood(tree, 2, theta, w, 0, 1, 1)) < 1e-5);
    igraph_adjlist_destroy(&al);
    mmpp_workspace_free(w);
}
END_TEST

START_TEST(test_node_assignment)
{
    igraph_t *tree = tree_from_newick("(t1:1,t2:1);");
    double theta[4] = {1, 2, 1, 3};
    int states[3] = {0, 0, 0};
    double pi[2];
    igraph_adjlist_t al;
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    _calculate_pi(2, theta, pi);
    ck_assert(fabs(lik_states(tree, 2, states, theta, pi, &al, 1) -
                   0.05779385) < 1e-5);

    states[2] = 1;
    fprintf(stderr, "%f\n", lik_states(tree, 2, states, theta, pi, &al, 1));
    ck_assert(fabs(lik_states(tree, 2, states, theta, pi, &al, 1) -
              0.002482152) < 1e-6);

    states[1] = 1;
    ck_assert(fabs(lik_states(tree, 2, states, theta, pi, &al, 1) -
                   0.0008891312) < 1e-7);

    ck_assert(fabs(pow(10, reconstruct(tree, 2, theta, w, states, 1, 1)) -
              0.05779385) < 1e-6);
    ck_assert_int_eq(states[0], 0);
    ck_assert_int_eq(states[1], 0);
    ck_assert_int_eq(states[2], 0);

    igraph_adjlist_destroy(&al);
    igraph_destroy(tree);
    mmpp_workspace_free(w);
}
END_TEST

START_TEST(test_reconstruct_nodes)
{
    igraph_t *tree = tree_from_newick("(t2:10.94,(t4:0.04,(t1:0.75,(t3:0.29,t5:0.8):0.39):0.21):0.8);");
    double theta[4] = {0.1, 1, 1, 2};
    int states1[9];
    int states2[9];
    int i;
    double pi[2];
    igraph_adjlist_t al;
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);

    _calculate_pi(2, theta, pi);
    igraph_adjlist_init(tree, &al, IGRAPH_OUT);

    pow(10, reconstruct(tree, 2, theta, w, states1, 1, 1));
    reconstruct_long(tree, 2, theta, states2, 1);

    for (i = 0; i < 9; ++i) {
        fprintf(stderr, "%d\t%d\n", states1[i], states2[i]);
    }
    ck_assert(fabs(lik_states(tree, 2, states1, theta, pi, &al, 1) -
                   pow(10, likelihood(tree, 2, theta, w, 1, 1, 1))) < 1e-8);
    mmpp_workspace_free(w);
}
END_TEST

Suite *mmpp_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("mmpp");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_likelihood_nodes_toy3);
    tcase_add_test(tc_core, test_likelihood_nodes_toy);
    tcase_add_test(tc_core, test_likelihood_toy);
    tcase_add_test(tc_core, test_guess_parameters);
    tcase_add_test(tc_core, test_likelihood_toy3);
    tcase_add_test(tc_core, test_get_clusters);
    tcase_add_test(tc_core, test_reconstruct);
    tcase_add_test(tc_core, test_node_assignment);
    tcase_add_test(tc_core, test_reconstruct_nodes);
    tcase_set_timeout(tc_core, 10);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = mmpp_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
