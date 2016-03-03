#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/util.h"

Suite *util_suite(void);

START_TEST (test_get_scale_nonzero)
{
    double x[2] = { 1e5, 1e-10 };
    int scale = get_scale(x, 2);
    ck_assert_int_eq(scale, -2);
}
END_TEST

START_TEST (test_get_scale_zero)
{
    double x[2] = { 1e5, 0 };
    int scale = get_scale((double *) &x, 2);
    ck_assert_int_eq(scale, 5);
    x[0] = 0; x[1] = 1e5;
    scale = get_scale((double *) &x, 2);
    ck_assert_int_eq(scale, 5);
}
END_TEST

START_TEST (test_get_scale_tiny)
{
    double x[2] = { 1e5, 1e-200 };
    int scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
    x[0] = 1e-200; x[1] = 1e5;
    scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
}
END_TEST

START_TEST (test_set_seed)
{
    gsl_rng *rng = set_seed(0);
    ck_assert(rand() == 1804289383);
    ck_assert(gsl_rng_get(rng) == 4293858116);
    ck_assert(igraph_rng_get_integer(igraph_rng_default(), 0, INT_MAX) == 2146929058);
}
END_TEST

START_TEST (test_order_ints)
{
    int x[5] = {4, 2, 1, 5, 3};
    int ord[5];
    order(x, ord, sizeof(int), 5, compare_ints);

    ck_assert_int_eq(ord[0], 2);
    ck_assert_int_eq(ord[1], 1);
    ck_assert_int_eq(ord[2], 4);
    ck_assert_int_eq(ord[3], 0);
    ck_assert_int_eq(ord[4], 3);
}
END_TEST

START_TEST (test_order_doubles)
{
    double x[5] = {4, 2, 1, 5, 3};
    int ord[5];
    order(x, ord, sizeof(double), 5, compare_doubles);

    ck_assert_int_eq(ord[0], 2);
    ck_assert_int_eq(ord[1], 1);
    ck_assert_int_eq(ord[2], 4);
    ck_assert_int_eq(ord[3], 0);
    ck_assert_int_eq(ord[4], 3);
}
END_TEST

START_TEST (test_rotl_int)
{
    int x[5] = {0, 1, 2, 3, 4};
    rotl(x, sizeof(x), 2 * sizeof(int));
    ck_assert_int_eq(x[0], 2);
    ck_assert_int_eq(x[1], 3);
    ck_assert_int_eq(x[2], 4);
    ck_assert_int_eq(x[3], 0);
    ck_assert_int_eq(x[4], 1);
}
END_TEST

START_TEST (test_rotl_double)
{
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    rotl(x, sizeof(x), 2 * sizeof(double));
    ck_assert(x[0] == 2.0);
    ck_assert(x[1] == 3.0);
    ck_assert(x[2] == 4.0);
    ck_assert(x[3] == 0.0);
    ck_assert(x[4] == 1.0);
}
END_TEST

START_TEST (test_permute_vector)
{
    int i;
    int perm[10] = {3, 7, 2, 5, 9, 1, 8, 6, 4, 0};
    int perm2[10] = {3, 1, 2, 5, 9, 0, 8, 6, 4, 7};
    igraph_vector_t v;
    igraph_vector_init(&v, 10);

    for (i = 0; i < 10; ++i) {
        VECTOR(v)[i] = i;
    }
    permute(&v, sizeof(igraph_real_t), 10, perm, get_igraph_vector_t, set_igraph_vector_t);
    for (i = 0; i < 10; ++i) {
        ck_assert(VECTOR(v)[perm[i]] == i);
    }

    for (i = 0; i < 10; ++i) {
        VECTOR(v)[i] = i;
    }
    permute(&v, sizeof(igraph_real_t), 10, perm2, get_igraph_vector_t, set_igraph_vector_t);
    for (i = 0; i < 10; ++i) {
        ck_assert(VECTOR(v)[perm2[i]] == i);
    }

    igraph_vector_destroy(&v);
}
END_TEST

START_TEST (test_permute_strvector)
{
    int i;
    int perm[10] = {3, 7, 2, 5, 9, 1, 8, 6, 4, 0};
    igraph_strvector_t v;
    igraph_strvector_init(&v, 10);
    char str[10];

    for (i = 0; i < 10; ++i) {
        sprintf(str, "%d", i);
        igraph_strvector_set(&v, i, str);
    }
    permute(&v, 2, 10, perm, get_igraph_strvector_t, set_igraph_strvector_t);
    for (i = 0; i < 10; ++i) {
        sprintf(str, "%d", i);
        ck_assert(strcmp(STR(v, perm[i]), str) == 0);
    }
    igraph_strvector_destroy(&v);
}
END_TEST

START_TEST (test_permute_vector_bool)
{
    int i;
    int perm[10] = {3, 7, 2, 5, 9, 1, 8, 6, 4, 0};
    igraph_vector_bool_t v;
    igraph_vector_bool_init(&v, 10);

    for (i = 0; i < 10; ++i) {
        VECTOR(v)[i] = i % 2;
    }
    permute(&v, 2, 10, perm, get_igraph_vector_bool_t, 
            set_igraph_vector_bool_t);
    for (i = 0; i < 10; ++i) {
        ck_assert(VECTOR(v)[perm[i]] % 2 == i % 2);
    }
    igraph_vector_bool_destroy(&v);
}
END_TEST

START_TEST (test_id_type_numeric)
{
    int i;
    igraph_t g;
    igraph_empty(&g, 5, 0);
    for (i = 0; i < igraph_vcount(&g); ++i) {
        SETVAN(&g, "id", i, i);
    }
    ck_assert_int_eq(get_igraph_id_type(&g), IGRAPH_ATTRIBUTE_NUMERIC);
    igraph_destroy(&g);
}
END_TEST

START_TEST (test_id_type_string)
{
    int i;
    igraph_t g;
    igraph_empty(&g, 5, 0);
    char buf[3];
    for (i = 0; i < igraph_vcount(&g); ++i) {
        sprintf(buf, "t%d", i);
        SETVAS(&g, "id", i, buf);
    }
    ck_assert_int_eq(get_igraph_id_type(&g), IGRAPH_ATTRIBUTE_STRING);
    igraph_destroy(&g);
}
END_TEST

START_TEST (test_id_type_none)
{
    int i;
    igraph_t g;
    igraph_empty(&g, 5, 0);
    ck_assert_int_eq(get_igraph_id_type(&g), IGRAPH_ATTRIBUTE_DEFAULT);
    igraph_destroy(&g);
}
END_TEST

START_TEST (test_match_vector)
{
    int i, index[3];
    igraph_vector_t x, tbl;

    igraph_vector_init(&x, 3);
    igraph_vector_init(&tbl, 3);

    for (i = 0; i < 3; ++i) {
        VECTOR(x)[i] = i;
        VECTOR(tbl)[i] = i+1;
    }

    match(&x, &tbl, index, sizeof(double), 3, 3, get_igraph_vector_t, compare_doubles);

    ck_assert_int_eq(index[0], -1);
    ck_assert_int_eq(index[1], 0);
    ck_assert_int_eq(index[2], 1);

    igraph_vector_destroy(&x);
    igraph_vector_destroy(&tbl);
}
END_TEST

START_TEST (test_match_strvector)
{
    int i, index[3];
    igraph_strvector_t x, tbl;

    igraph_strvector_init(&x, 0);
    igraph_strvector_init(&tbl, 0);

    igraph_strvector_add(&x, "e");
    igraph_strvector_add(&x, "a");
    igraph_strvector_add(&x, "c");

    igraph_strvector_add(&tbl, "a");
    igraph_strvector_add(&tbl, "b");
    igraph_strvector_add(&tbl, "c");
    igraph_strvector_add(&tbl, "d");
    igraph_strvector_add(&tbl, "e");

    match(&x, &tbl, index, 1, 3, 5, get_igraph_strvector_t, compare_strings);

    ck_assert_int_eq(index[0], 4);
    ck_assert_int_eq(index[1], 0);
    ck_assert_int_eq(index[2], 2);

    igraph_strvector_destroy(&x);
    igraph_strvector_destroy(&tbl);
}
END_TEST

START_TEST (test_sample)
{
    int x[5] = {0, 1, 2, 3, 4};
    double prob[5] = {0, 0, 1, 0, 0};
    int s[2];
    gsl_rng *r = set_seed(0);

    sample_weighted(x, s, 1, 5, sizeof(int), prob, 0, r);
    ck_assert_int_eq(s[0], 2);

    sample_weighted(x, s, 2, 5, sizeof(int), prob, 1, r);
    ck_assert_int_eq(s[0], 2);
    ck_assert_int_eq(s[1], 2);

    prob[0] = 1;
    sample_weighted(x, s, 2, 5, sizeof(int), prob, 0, r);
    ck_assert( (s[0] == 0 && s[1] == 2) || (s[0] == 2 && s[1] == 0) );

    gsl_rng_free(r);
}
END_TEST

Suite *util_suite(void)
{
    Suite *s;
    TCase *tc_core;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = suite_create("util");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_get_scale_zero);
    tcase_add_test(tc_core, test_get_scale_nonzero);
    tcase_add_test(tc_core, test_get_scale_tiny);
    tcase_add_test(tc_core, test_set_seed);
    tcase_add_test(tc_core, test_order_ints);
    tcase_add_test(tc_core, test_order_doubles);
    tcase_add_test(tc_core, test_rotl_int);
    tcase_add_test(tc_core, test_permute_vector);
    tcase_add_test(tc_core, test_permute_strvector);
    tcase_add_test(tc_core, test_permute_vector_bool);
    tcase_add_test(tc_core, test_id_type_numeric);
    tcase_add_test(tc_core, test_id_type_string);
    tcase_add_test(tc_core, test_id_type_none);
    tcase_add_test(tc_core, test_match_vector);
    tcase_add_test(tc_core, test_match_strvector);
    tcase_add_test(tc_core, test_sample);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = util_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
