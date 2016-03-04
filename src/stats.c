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
            default:
                fprintf(stderr, "BUG: tried to calculate density for unknown distribution\n");
                break;
        }
    }
    return dens;
}

