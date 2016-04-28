#include <math.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "stats.h"

double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt)
{
    double teststat = - 2 * log10lik_null / log10(exp(1)) + 2 * log10lik_alt / log10(exp(1));
    int df = nparam_alt - nparam_null;
    return 1 - gsl_cdf_chisq_P(teststat, df);
}

double bic(double log10lik, int nparam, int ndata)
{
    return -2 * log10lik / log10(exp(1)) + nparam * log(ndata);
}

double aic(double log10lik, int nparam)
{
    return -2 * log10lik / log10(exp(1)) + 2 * nparam;
}

void sample_distribution(int n, double *theta, gsl_rng *rng, 
                         const distribution *dist, double **params)
{
    int i;

    for (i = 0; i < n; ++i)
    {
        theta[i] = params[i][0];
        switch (dist[i])
        {
            case UNIFORM:
                theta[i] += gsl_ran_flat(rng, 0, params[i][1] - params[i][0]);
                break;
            case GAUSSIAN:
                theta[i] += gsl_ran_gaussian(rng, params[i][1]);
                break;
            case DELTA:
                break;
            case EXPONENTIAL:
                theta[i] += gsl_ran_exponential(rng, 1.0 / params[i][1]);
                break;
            case LAPLACE:
                theta[i] += gsl_ran_laplace(rng, params[i][1]);
                break;
            case EXPONENTIAL_POWER:
                theta[i] += gsl_ran_exppow(rng, params[i][1], params[i][2]);
                break;
            case CAUCHY:
                theta[i] += gsl_ran_cauchy(rng, params[i][1]);
                break;
            case RAYLEIGH:
                theta[i] += gsl_ran_rayleigh(rng, params[i][1]);
                break;
            case GAMMA:
                theta[i] = gsl_ran_gamma(rng, params[i][0], params[i][1]);
                break;
            case LOGNORMAL:
                theta[i] = gsl_ran_lognormal(rng, params[i][0], params[i][1]);
                break;
            case CHI_SQUARED:
                theta[i] += gsl_ran_chisq(rng, params[i][1]);
                break;
            case F:
                theta[i] += gsl_ran_fdist(rng, params[i][1], params[i][2]);
                break;
            case STUDENT_T:
                theta[i] += gsl_ran_tdist(rng, params[i][1]);
                break;
            case BETA:
                theta[i] += gsl_ran_beta(rng, params[i][1], params[i][2]) * params[i][3];
                break;
            case LOGISTIC:
                theta[i] += gsl_ran_logistic(rng, params[i][1]);
                break;
            case PARETO:
                theta[i] = gsl_ran_pareto(rng, params[i][1], params[i][0]);
                break;
            case WEIBULL:
                theta[i] += gsl_ran_weibull(rng, params[i][1], params[i][2]);
                break;
            case POISSON:
                theta[i] += gsl_ran_poisson(rng, params[i][1]);
                break;
            case DISCRETE_UNIFORM:
                theta[i] += floor(gsl_ran_flat(rng, 0, params[i][1] + 1 - params[i][0]));
                break;
            default:
                fprintf(stderr, "BUG: tried to sample from unknown distribution\n");
                theta[i] = 0;
                break;
        }
    }
}

double density_distribution(int n, const double *theta,
                            const distribution *dist, double **params)
{
    int i;
    double dens = 1;

    for (i = 0; i < n; ++i)
    {
        switch (dist[i])
        {
            case UNIFORM:
                dens *= gsl_ran_flat_pdf(theta[i], params[i][0], params[i][1]);
                break;
            case GAUSSIAN:
                dens *= gsl_ran_gaussian_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case DELTA:
                dens *= fabs(theta[i] - params[i][0]) < FLT_EPSILON;
                break;
            case EXPONENTIAL:
                dens *= gsl_ran_exponential_pdf(theta[i] - params[i][0], 1.0 / params[i][1]);
                break;
            case LAPLACE:
                dens *= gsl_ran_laplace_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case EXPONENTIAL_POWER:
                dens *= gsl_ran_exppow_pdf(theta[i] - params[i][0], params[i][1], params[i][2]);
                break;
            case CAUCHY:
                dens *= gsl_ran_cauchy_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case RAYLEIGH:
                dens *= gsl_ran_rayleigh_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case GAMMA:
                dens *= gsl_ran_gamma_pdf(theta[i], params[i][0], params[i][1]);
                break;
            case LOGNORMAL:
                dens *= gsl_ran_lognormal_pdf(theta[i], params[i][0], params[i][1]);
                break;
            case CHI_SQUARED:
                dens *= gsl_ran_chisq_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case F:
                dens *= gsl_ran_fdist_pdf(theta[i] - params[i][0], params[i][1], params[i][2]);
                break;
            case STUDENT_T:
                dens *= gsl_ran_tdist_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case BETA:
                dens *= gsl_ran_beta_pdf((theta[i] - params[i][0]) / params[i][3], params[i][1], params[i][2]);
                break;
            case LOGISTIC:
                dens *= gsl_ran_logistic_pdf(theta[i] - params[i][0], params[i][1]);
                break;
            case PARETO:
                dens *= gsl_ran_pareto_pdf(theta[i], params[i][1], params[i][0]);
                break;
            case WEIBULL:
                dens *= gsl_ran_weibull_pdf(theta[i] - params[i][0], params[i][1], params[i][2]);
                break;
            case POISSON:
                dens *= gsl_ran_poisson_pdf( (int) (theta[i] - params[i][0]), params[i][1] );
                break;
            case DISCRETE_UNIFORM:
                dens *= gsl_ran_flat_pdf(theta[i], params[i][0], params[i][1]+1);
                break;
            default:
                fprintf(stderr, "BUG: tried to calculate density for unknown distribution\n");
                break;
        }
    }
    return dens;
}

distribution parse_distribution(const char *s)
{
    if (strcmp(s, "uniform") == 0) {
        return UNIFORM;
    }
    else if (strcmp(s, "gaussian") == 0) {
        return GAUSSIAN;
    }
    else if (strcmp(s, "delta") == 0) {
        return DELTA;
    }
    else if (strcmp(s, "exponential") == 0) {
        return EXPONENTIAL;
    }
    else if (strcmp(s, "laplace") == 0) {
        return LAPLACE;
    }
    else if (strcmp(s, "exponential_power") == 0) {
        return EXPONENTIAL_POWER;
    }
    else if (strcmp(s, "cauchy") == 0) {
        return CAUCHY;
    }
    else if (strcmp(s, "rayleigh") == 0) {
        return RAYLEIGH;
    }
    else if (strcmp(s, "gamma") == 0) {
        return GAMMA;
    }
    else if (strcmp(s, "lognormal") == 0) {
        return LOGNORMAL;
    }
    else if (strcmp(s, "chi_squared") == 0) {
        return CHI_SQUARED;
    }
    else if (strcmp(s, "f") == 0) {
        return F;
    }
    else if (strcmp(s, "student_t") == 0) {
        return STUDENT_T;
    }
    else if (strcmp(s, "beta") == 0) {
        return BETA;
    }
    else if (strcmp(s, "logistic") == 0) {
        return LOGISTIC;
    }
    else if (strcmp(s, "pareto") == 0) {
        return PARETO;
    }
    else if (strcmp(s, "weibull") == 0) {
        return WEIBULL;
    }
    else if (strcmp(s, "poisson") == 0) {
        return POISSON;
    }
    else if (strcmp(s, "discrete_uniform") == 0) {
        return DISCRETE_UNIFORM;
    }
    return 0;
}
