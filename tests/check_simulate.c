#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/tree.h"
#include "../src/util.h"
#include "../src/simulate.h"

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

START_TEST (test_simulate_chain)
{
    igraph_t net, tree;
    char buf[2];
    int i, from, to;
    double diff;
    gsl_rng *rng = set_seed(0);

    igraph_empty(&net, 5, 1);
    for (i = 0; i < 5; ++i)
    {
        if (i < 4)
        {
            igraph_add_edge(&net, i, i+1);
            igraph_add_edge(&net, i+1, i);
        }
        sprintf(buf, "%d", i);
        SETVAS(&net, "id", i, buf);
        SETVAN(&net, "remove", i, 0);
    }
    for (i = 0; i < igraph_ecount(&net); ++i)
    {
        SETEAN(&net, "transmit", i, 1);
    }
    simulate_phylogeny(&tree, &net, rng, 100, 100, 0);

    ck_assert_int_eq(igraph_vcount(&tree), 9);
    ck_assert_int_eq(igraph_ecount(&tree), 8);
    for (i = 0; i < igraph_ecount(&tree); ++i)
    {
        igraph_edge(&tree, i, &from, &to);
        diff = fabs(atoi(VAS(&tree, "id", from)) - atoi(VAS(&tree, "id", to)));
        ck_assert(diff == 0.0 || diff == 1.0);
    }

    igraph_destroy(&net);
    igraph_destroy(&tree);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_io, *tc_tree;

    s = suite_create("simulate");

    tc_tree = tcase_create("Core");
    tcase_add_test(tc_tree, test_simulate_chain);
    suite_add_tcase(s, tc_tree);

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
