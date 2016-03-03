#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/tree.h"
#include "../src/treestats.h"
#include "../src/util.h"

Suite *tree_suite(void);

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

START_TEST(test_kernel)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("(((1:0.25,2:0.25)5:0.5,3:0.25)6:0.5,4:0.25)7;");

    ck_assert(fabs(kernel(t1, t2, 0.5, 1, 1) - 1.125 * (1 + exp(-0.0625))) < 1e-5);
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_nLTT)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("((((1:0.4,2:0.4)6:0.3,3:0.4)7:0.2,4:0.4)8:0.1,5:0.4)9;");
    ck_assert(fabs(nLTT(t1, t2) - 0.111111) < 1e-5);
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_nLTT_identical)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    ck_assert(fabs(nLTT(t1, t1)) < 1e-5);
    igraph_destroy(t1);
}
END_TEST

START_TEST(test_nLTT_uniform_diffsizes)
{
    igraph_t *t1 = tree_from_newick("((((:1,:1):1,:1):1,:1):1,:1);");
    igraph_t *t2 = tree_from_newick("(((:1,:1):1,:1):1,:1);");
    ck_assert(fabs(nLTT(t1, t2)) < 1e-5);
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_sackin)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(sackin(tree, 0) == 13);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_sackin_yule)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    double s = (sackin(tree, 0) - SACKIN_YULE(NTIP(tree))) / NTIP(tree);
    ck_assert(fabs(s - 0.0333333) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_sackin_branch_lengths)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(fabs(sackin(tree, 1) - 5.81) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_colless)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert_int_eq(colless(t), 2);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_colless_yule)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    double c = (colless(t) - COLLESS_YULE(NTIP(t))) / NTIP(t);
    ck_assert(fabs(c + 0.09350639) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_cophenetic)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(cophenetic(t, 0) == 5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_cophenetic_branch_lengths)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(fabs(cophenetic(t, 1) - 3.04) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_ladder_length)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert_int_eq(ladder_length(t), 5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_il_nodes)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert_int_eq(il_nodes(t), 2);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_width)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert_int_eq(width(t), 6);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_max_delta_width)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert_int_eq(max_delta_width(t), 4);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_cherries)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert_int_eq(cherries(t), 3);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_unbalanced)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert(fabs(prop_unbalanced(t) - 2.0 / 7.0) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_avg_unbalance)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert(fabs(avg_unbalance(t) - 0.8333333) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_gamma)
{
    igraph_t *t = tree_from_newick("(((t4:0.29,(t6:0.03,t1:0.03):0.26):0.01,(t3:0.29,t5:0.29):0.02):0.51,((t8:0.19,t2:0.19):0.13,t7:0.32):0.5);");
    ck_assert(fabs(pybus_gamma(t) + 0.2562976) < 1e-5);
    igraph_destroy(t);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("tree");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_kernel);
    tcase_add_test(tc_core, test_nLTT);
    tcase_add_test(tc_core, test_nLTT_identical);
    tcase_add_test(tc_core, test_nLTT_uniform_diffsizes);
    tcase_add_test(tc_core, test_sackin);
    tcase_add_test(tc_core, test_sackin_branch_lengths);
    tcase_add_test(tc_core, test_sackin_yule);
    tcase_add_test(tc_core, test_colless);
    tcase_add_test(tc_core, test_colless_yule);
    tcase_add_test(tc_core, test_cophenetic);
    tcase_add_test(tc_core, test_cophenetic_branch_lengths);
    tcase_add_test(tc_core, test_ladder_length);
    tcase_add_test(tc_core, test_il_nodes);
    tcase_add_test(tc_core, test_width);
    tcase_add_test(tc_core, test_max_delta_width);
    tcase_add_test(tc_core, test_cherries);
    tcase_add_test(tc_core, test_unbalanced);
    tcase_add_test(tc_core, test_avg_unbalance);
    tcase_add_test(tc_core, test_gamma);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = tree_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
