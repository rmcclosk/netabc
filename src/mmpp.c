#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> 
#include <igraph/igraph.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_cblas.h>
#include "../c-cmaes/cmaes_interface.h"
#include "../c-cmaes/boundary_transformation.h"
#include "mmpp.h"
#include "tree.h"
#include "util.h"
#include "stats.h"

#define CMAES_POP_SIZE 100
#define MAX_NRATES 6

struct ck_params {
    int nrates;
    double *Q;
};

struct mmpp_workspace {
    struct ck_params ckpar;
    double *y;
    double *P;
    double *L;
    double *Li;
    int *C;
    double *pi;
    double *branch_lengths;
    int *scale;
    int *bl_order;
    gsl_matrix *Q;
    gsl_vector_complex *eval;
    gsl_matrix_complex *evec;
    gsl_eigen_nonsymmv_workspace *ew;
    igraph_adjlist_t al;
};

int ckdiff(double t, const double y[], double f[], void *params);
void calculate_P(const igraph_t *tree, int nrates, const double *theta, 
        mmpp_workspace *w, int trans_at_nodes, int use_tips);
void calculate_pi(const igraph_t *tree, int nrates, const double *theta, 
        mmpp_workspace *w);
void mmpp_workspace_set_params(mmpp_workspace *w, const double *theta,
        int trans_at_nodes);
int _fit_mmpp(const igraph_t *tree, int nrates, double *theta, int trace,
             const char *cmaes_settings, int *states, double *loglik,
             int use_tips, int trans_at_nodes, double bounds[4]);

int fit_mmpp(const igraph_t *tree, int *nrates, double **theta, int trace,
             const char *cmaes_settings, int *states, model_selector sel,
             int use_tips, int trans_at_nodes, double bounds[4])
{
    int i, accept, error = 0;
    double loglik, prev_loglik, test_stat;
    double *prev_theta = malloc(MAX_NRATES * MAX_NRATES * sizeof(double));

    if (*nrates > 0)
    {
        error = _fit_mmpp(tree, *nrates, *theta, trace, cmaes_settings, states,
                          &loglik, use_tips, trans_at_nodes, bounds);
        fprintf(stderr, "log likelihood for %d state model is %f\n", *nrates, loglik);
        return error;
    }

    error = _fit_mmpp(tree, 1, prev_theta, trace, cmaes_settings, states,
                      &prev_loglik, use_tips, trans_at_nodes, bounds);
    fprintf(stderr, "log likelihood for 1 state model is %f\n", prev_loglik);
    *nrates = 1;

    for (i = 2; i <= MAX_NRATES && !error; ++i)
    {
        *theta = safe_realloc(*theta, i * i * sizeof(double));

        // if we failed to fit a model with more states, fall back to the previous one
        if (_fit_mmpp(tree, i, *theta, trace, cmaes_settings, states, &loglik,
                      use_tips, trans_at_nodes, bounds)) {
            fprintf(stderr, "Warning: parameter estimates for %d state model did not converge\n", i);
        }
        fprintf(stderr, "log likelihood for %d state model is %f\n", i, loglik);
        
        if (sel == LRT) {
            test_stat = lrt(prev_loglik, loglik, (i-1)*(i-1), i*i);
            fprintf(stderr, "P-value for %d state model is %f\n", i, test_stat);
            accept = test_stat < 0.05;
        }
        else if (sel == AIC) {
            test_stat = aic(loglik, i*i) - aic(prev_loglik, (i-1)*(i-1));
            fprintf(stderr, "delta AIC for %d state model is %f\n", i, test_stat);
            accept = test_stat <= -2;
        }
        else if (sel == BIC) {
            test_stat = bic(loglik, i*i, igraph_ecount(tree)) - bic(prev_loglik, (i-1)*(i-1), igraph_ecount(tree));
            fprintf(stderr, "delta BIC for %d state model is = %f\n", i, test_stat);
            accept = test_stat <= -2;
        }

        if (!accept)
        {
            fprintf(stderr, "%d state model is not supported\n", i);
            break;
        }

        memcpy(prev_theta, *theta, i * i * sizeof(double));
        prev_loglik = loglik;
        *nrates = i;
    }

    memcpy(*theta, prev_theta, *nrates * *nrates * sizeof(double));
    free(prev_theta);
    return error;
}

void get_clusters(const igraph_t *tree, const int *states, int *clusters,
        int cluster_state)
{
    igraph_adjlist_t al;
    igraph_vector_int_t *children;
    int i, lchild, rchild, ccount = 0;

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);

    memset(clusters, 0, igraph_vcount(tree) * sizeof(int));

    if (states[root(tree)] >= cluster_state) {
        clusters[root(tree)] = ++ccount;
    }

    for (i = igraph_vcount(tree)-1; i >= 0; --i) {
        children = igraph_adjlist_get(&al, i);
        if (igraph_vector_int_size(children) > 0)
        {
            lchild = VECTOR(*children)[0];
            rchild = VECTOR(*children)[1];

            if (states[lchild] < cluster_state)
                clusters[lchild] = 0;
            else if (states[i] < cluster_state)
                clusters[lchild] = ++ccount;
            else
                clusters[lchild] = clusters[i];

            if (states[rchild] < cluster_state)
                clusters[rchild] = 0;
            else if (states[i] < cluster_state)
                clusters[rchild] = ++ccount;
            else
                clusters[rchild] = clusters[i];
        }
    }

    igraph_adjlist_destroy(&al);
}

mmpp_workspace *mmpp_workspace_create(const igraph_t *tree, int nrates)
{
    struct mmpp_workspace *w = malloc(sizeof(struct mmpp_workspace));
    igraph_vector_t vec;

    w->ckpar.nrates = nrates;
    w->ckpar.Q = malloc(nrates * nrates * sizeof(double));
    w->y = malloc(nrates * nrates * sizeof(double));
    w->P = malloc(nrates * nrates * igraph_vcount(tree) * sizeof(double));
    w->L = malloc(nrates * igraph_vcount(tree) * sizeof(double));
    w->Li = malloc(nrates * nrates * sizeof(double));
    w->C = malloc(nrates * igraph_vcount(tree) * sizeof(int));
    w->scale = malloc(igraph_vcount(tree) * sizeof(int));
    w->pi = malloc(nrates * sizeof(double));
    w->branch_lengths = malloc(igraph_ecount(tree) * sizeof(double));
    w->bl_order = malloc(igraph_ecount(tree) * sizeof(int));
    igraph_adjlist_init(tree, &w->al, IGRAPH_OUT);

    w->Q = gsl_matrix_alloc(nrates, nrates);
    w->eval = gsl_vector_complex_alloc(nrates);
    w->evec = gsl_matrix_complex_alloc(nrates, nrates);
    w->ew = gsl_eigen_nonsymmv_alloc(nrates);

    // collect and order branch lengths
    igraph_vector_init(&vec, igraph_ecount(tree));
    EANV(tree, "length", &vec);
    order(VECTOR(vec), w->bl_order, sizeof(double), igraph_ecount(tree), compare_doubles);
    memcpy(w->branch_lengths, VECTOR(vec), igraph_ecount(tree) * sizeof(double));

    igraph_vector_destroy(&vec);
    return w;
}

void mmpp_workspace_free(mmpp_workspace *w)
{
    free(w->ckpar.Q);
    free(w->y);
    free(w->P);
    free(w->L);
    free(w->Li);
    free(w->C);
    free(w->scale);
    free(w->pi);
    free(w->bl_order);
    free(w->branch_lengths);
    igraph_adjlist_destroy(&w->al);
    gsl_matrix_free(w->Q);
    gsl_vector_complex_free(w->eval);
    gsl_matrix_complex_free(w->evec);
    gsl_eigen_nonsymmv_free(w->ew);
    free(w);
}

void guess_parameters(const igraph_t *tree, int nrates, double *theta)
{
    int i, nbranch = (igraph_vcount(tree) - 3)/2, cur = 0;
    int part_size = (int) ceil((float) nbranch / (float) nrates);
    igraph_inclist_t il;
    igraph_vector_t degree;
    igraph_vector_int_t *edge;
    double *int_edges = malloc(nbranch * sizeof(double));
    double intsum = 0;

    // find out-degree of all nodes
    igraph_vector_init(&degree, igraph_vcount(tree));
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    // collect the internal branch lengths
    igraph_inclist_init(tree, &il, IGRAPH_IN);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        edge = igraph_inclist_get(&il, i);
        if (igraph_vector_int_size(edge) > 0 && (int) VECTOR(degree)[i] > 0) {
            int_edges[cur++] = EAN(tree, "length", VECTOR(*edge)[0]);
            intsum += int_edges[cur-1];
        }
    }

    // for branching rates, partition internal branch lengths, and take means of partitions
    qsort(int_edges, nbranch, sizeof(double), compare_doubles);
    memset(theta, 0, nrates * sizeof(double));
    for (i = 0; i < nbranch; ++i)
        theta[nrates - (i / part_size) - 1] += int_edges[i] / part_size;

    // the guessed rates are the inverse of the means
    for (i = 0; i < nrates; ++i)
        theta[i] = 1.0 / fmax(theta[i], 1e-4);

    // choose transition rates all equal such that the expected number of
    // transitions is 1% of the number of branches
    for (i = 0; i < nrates * (nrates - 1); ++i)
        theta[nrates+i] = nbranch / nrates / intsum / 100.0;

    // clean up
    igraph_inclist_destroy(&il);
    free(int_edges);
}

double likelihood(const igraph_t *tree, int nrates, const double *theta,
                  mmpp_workspace *w, int trans_at_nodes, int use_tips,
                  int reconstruct)
{ 
    int i, j, rt = root(tree);
    double lik = 0;
    int lchild, rchild, pstate, cstate, new_scale;
    igraph_vector_int_t *children;

    mmpp_workspace_set_params(w, theta, trans_at_nodes);
    calculate_P(tree, nrates, theta, w, trans_at_nodes, use_tips);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        children = igraph_adjlist_get(&w->al, i);
        if (igraph_vector_int_size(children) == 0)
        {
            lchild = -1; rchild = -1;
            w->scale[i] = 0;
        }
        else
        {
            lchild = VECTOR(*children)[0];
            rchild = VECTOR(*children)[1];
            w->scale[i] = w->scale[lchild] + w->scale[rchild];
        }

        for (pstate = 0; pstate < nrates; ++pstate)
        {
            for (cstate = 0; cstate < nrates; ++cstate) {
                if (lchild != -1) {
                    w->Li[cstate] = w->L[lchild * nrates + cstate] *
                                    w->L[rchild * nrates + cstate] *
                                    w->P[i * nrates * nrates + pstate * nrates + cstate];
                }
                else {
                    w->Li[cstate] = w->P[i * nrates * nrates + pstate * nrates + cstate];
                }
            }

            if (reconstruct) {
                w->C[i * nrates + pstate] = which_max(w->Li, nrates);
                w->L[i * nrates + pstate] = w->Li[w->C[i * nrates + pstate]];
            }
            else {
                w->L[i * nrates + pstate] = sum_doubles(w->Li, nrates);
            }
        }
        new_scale = get_scale(&w->L[i * nrates], nrates);
        for (pstate = 0; pstate < nrates; ++pstate) {
            w->L[i * nrates + pstate] /= pow(10, new_scale);
        }
        w->scale[i] += new_scale;
    }

    lik = (reconstruct ? max_doubles : sum_doubles)(&w->L[rt * nrates], nrates);
    return log10(lik) + w->scale[rt];
}

double reconstruct(const igraph_t *tree, int nrates, const double *theta,
        mmpp_workspace *w, int *states, int use_tips, int trans_at_nodes)
{
    int i, rt = root(tree), lchild, rchild;
    double lik = likelihood(tree, nrates, theta, w, trans_at_nodes, use_tips, 1);
    igraph_vector_int_t *children;

    states[rt] = which_max(&w->L[rt * nrates], nrates);
    for (i = rt; i >= 0; --i)
    {
        children = igraph_adjlist_get(&w->al, i);
        if (igraph_vector_int_size(children) > 0)
        {
            lchild = VECTOR(*children)[0];
            rchild = VECTOR(*children)[1];
            states[lchild] = w->C[lchild * nrates + states[i]];
            states[rchild] = w->C[rchild * nrates + states[i]];
        }
    }
    return lik;
}

/* Private. */
int _fit_mmpp(const igraph_t *tree, int nrates, double *theta, int trace,
             const char *cmaes_settings, int *states, double *loglik, 
             int use_tips, int trans_at_nodes, double bounds[4])
{
    int i, j, dimension = nrates * nrates, error = 0, cur = nrates;
    int *state_order;
    double *lbound = malloc(dimension * sizeof(double));
    double *ubound = malloc(dimension * sizeof(double));
    double *init_sd = malloc(dimension * sizeof(double));
    double *funvals, *tmp, *const *pop;
    struct mmpp_workspace *w = mmpp_workspace_create(tree, nrates);
    cmaes_t evo;
    cmaes_boundary_transformation_t trbound;

    for (i = 0; i < nrates; ++i)
    {
        lbound[i] = log(bounds[0]);
        ubound[i] = log(bounds[1]);
        init_sd[i] = 1;
    }
    for (i = nrates; i < dimension; ++i)
    {
        lbound[i] = log(bounds[2]);
        ubound[i] = log(bounds[3]);
        init_sd[i] = 1;
    }

    cmaes_boundary_transformation_init(&trbound, lbound, ubound, dimension);

    guess_parameters(tree, nrates, theta);
    for (i = 0; i < dimension; ++i)
        theta[i] = log(theta[i]);

    funvals = cmaes_init(&evo, dimension, theta, init_sd, 0, CMAES_POP_SIZE, cmaes_settings);

	while (!cmaes_TestForTermination(&evo)) {

		pop = cmaes_SamplePopulation(&evo);
		for (i = 0; i < CMAES_POP_SIZE; ++i) {
            cmaes_boundary_transformation(&trbound, pop[i], theta, dimension);
            for (j = 0; j < dimension; ++j)
            {
                theta[j] = exp(theta[j]);
                if (trace)
                    fprintf(stderr, "%f\t", theta[j]);
            }
            funvals[i] = -likelihood(tree, nrates, theta, w, trans_at_nodes, use_tips, 0);
            if (funvals[i] != funvals[i])
                funvals[i] = FLT_MAX;
            if (trace)
                fprintf(stderr, "%f\n", -funvals[i]);
        }
		cmaes_UpdateDistribution(&evo, funvals);
    }

    if (strncmp(cmaes_TestForTermination(&evo), "TolFun", 6) != 0)
    {
        error = 1;
        fprintf(stderr, "%s", cmaes_TestForTermination(&evo));
    }

    cmaes_boundary_transformation(&trbound, 
        (double const *) cmaes_GetPtr(&evo, "xbestever"), theta, dimension);

    state_order = malloc(nrates * sizeof(int));
    tmp = malloc(dimension * sizeof(double));
    for (i = 0; i < dimension; ++i)
        tmp[i] = exp(theta[i]);
    order(tmp, state_order, sizeof(double), nrates, compare_doubles);
    for (i = 0; i < nrates; ++i)
    {
        theta[i] = tmp[state_order[i]];
        for (j = 0; j < nrates; ++j)
        {
            if (state_order[i] > state_order[j])
	            theta[cur++] = tmp[nrates + state_order[i]*(nrates-1) + state_order[j]];
            else if (state_order[i] < state_order[j])
	            theta[cur++] = tmp[nrates + state_order[i]*(nrates-1) + state_order[j]-1];
        }
    }
    loglik[0] = likelihood(tree, nrates, theta, w, trans_at_nodes, 1, 0);
    if (states != NULL)
        reconstruct(tree, nrates, theta, w, states, 1, trans_at_nodes);

    cmaes_exit(&evo);
    cmaes_boundary_transformation_exit(&trbound);
    free(lbound);
    free(ubound);
    free(init_sd);
    free(state_order);
    free(tmp);
    mmpp_workspace_free(w);
    return error;
}

void mmpp_workspace_set_params(mmpp_workspace *w, const double *theta, 
        int trans_at_nodes)
{
    double rowsum = 0;
    int i, j, nrates = w->ckpar.nrates, cur = nrates;

    for (i = 0; i < nrates; ++i)
    {
        rowsum = 0;
        for (j = 0; j < nrates; ++j)
        {
            if (j != i)
            {
                w->y[i * nrates + j] = 0;
                w->ckpar.Q[i * nrates + j] = theta[cur];
                rowsum += theta[cur++];
            }
        }
        if (trans_at_nodes)
            w->ckpar.Q[i * nrates + i] = -rowsum;
        else
            w->ckpar.Q[i * nrates + i] = -rowsum - theta[i];
        w->y[i * nrates + i] = 1;
    }
}

void calculate_P(const igraph_t *tree, int nrates, const double *theta,
        struct mmpp_workspace *w, int trans_at_nodes, int use_tips)
{
    int i, j, k, from, to, cur = nrates; 
    double t = 0.0;
    igraph_vector_int_t *children;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver *d;

    // set the values at the root of the tree to the equilibrium frequencies
    calculate_pi(tree, nrates, theta, w);
    memset(&w->P[root(tree) * nrates * nrates], 0, nrates * nrates * sizeof(double));
    for (i = 0; i < nrates; ++i) {
        w->P[root(tree) * nrates * nrates + i * nrates + i] = w->pi[i];
    }

    // set up ODE
    sys.function = ckdiff;
    sys.jacobian = NULL;
    sys.dimension = nrates * nrates;
    sys.params = &w->ckpar;
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    // calculate values for each branch
    for (i = 0; i < igraph_ecount(tree); ++i)
    {
        if (w->branch_lengths[w->bl_order[i]] != t)
            gsl_odeiv2_driver_apply(d, &t, w->branch_lengths[w->bl_order[i]], w->y);
        igraph_edge(tree, w->bl_order[i], &from, &to);
        memcpy(&w->P[to * nrates * nrates], w->y, nrates * nrates * sizeof(double));

        // for non-terminal nodes, multiply by branching rate
        children = igraph_adjlist_get(&w->al, to);
        for (j = 0; j < nrates; ++j) {
            for (k = 0; k < nrates; ++k) {
                if (igraph_vector_int_size(children) == 0 && !use_tips)
                    w->P[to * nrates * nrates + j * nrates + k] = (j == k);
                else if (igraph_vector_int_size(children) > 0 && trans_at_nodes)
                    w->P[to * nrates * nrates + j * nrates + k] *= theta[j] * exp(-theta[j] * t);
                else if (igraph_vector_int_size(children) == 0 && trans_at_nodes)
                    w->P[to * nrates * nrates + j * nrates + k] *= exp(-theta[j] * t);
                else if (igraph_vector_int_size(children) > 0 && !trans_at_nodes)
                    w->P[to * nrates * nrates + j * nrates + k] *= theta[k];
            }
        }
    }
    gsl_odeiv2_driver_free(d);
}

void calculate_pi(const igraph_t *tree, int nrates, const double *theta, mmpp_workspace *w)
{
    int i, j, cur = nrates;
    double sum;

    // calculate Q
    for (i = 0; i < nrates; ++i) {
        sum = 0;
        for (j = 0; j < nrates; ++j) {
            if (i != j) {
                gsl_matrix_set(w->Q, i, j, theta[cur]);
                sum += theta[cur++];
            }
        }
        gsl_matrix_set(w->Q, i, i, -sum);
    }

    // calculate pi
    gsl_matrix_transpose(w->Q);
    gsl_eigen_nonsymmv(w->Q, w->eval, w->evec, w->ew);
    gsl_eigen_nonsymmv_sort(w->eval, w->evec, GSL_EIGEN_SORT_ABS_ASC);

    sum = 0;
    for (i = 0; i < nrates; ++i)
    {
        w->pi[i] = GSL_REAL(gsl_matrix_complex_get(w->evec, i, 0));
        sum += w->pi[i];
    }
    for (i = 0; i < nrates; ++i)
        w->pi[i] /= sum;
}

int ckdiff(double t, const double y[], double f[], void *params)
{
    struct ck_params *ckpar = (struct ck_params *) params;
    int n = ckpar->nrates;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y, n, 
                ckpar->Q, n, 0.0, f, n);
    return GSL_SUCCESS;
}
