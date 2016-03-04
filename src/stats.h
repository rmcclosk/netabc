#ifndef STATS_H
#define STATS_H

#include <gsl/gsl_rng.h>

typedef enum {
    AIC,
    BIC,
    LRT
} model_selector;

typedef enum {
    UNIFORM = 1,
    GAUSSIAN = 2,
    DELTA = 3,
    EXPONENTIAL = 4,
    LAPLACE = 5,
    EXPONENTIAL_POWER = 6,
    CAUCHY = 7,
    RAYLEIGH = 8,
    GAMMA = 9,
    LOGNORMAL = 10,
    CHI_SQUARED = 11,
    F = 12,
    STUDENT_T = 13,
    BETA = 14,
    LOGISTIC = 15,
    PARETO = 16,
    WEIBULL = 17
} distribution;

/** Calculate a likelihood ratio test.
 *
 * \param[in] log10lik_null log10 likelihood of null model
 * \param[in] log10lik_alt log10 likelihood of alternative model
 * \param[in] nparam_null number of parameters of null model
 * \param[in] nparam_alt number of parameters of alternative model
 * \return p-value for support of the alternative model
 */
double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt);

/** Calculate the Bayesian information criterion for a model fit.
 * 
 * \param[in] log10lik log10 likelihood of the fitted model
 * \param[in] nparam number of parameters of the model
 * \param[in] ndata number of data points of the model
 * \return the BIC
 */
double bic(double log10lik, int nparam, int ndata);

/** Calculate the Akaike information criterion for a model fit.
 *
 * \param[in] log10lik log 10 likelihood of the fitted model
 * \param[in] nparam number of parameters of the model
 */
double aic(double log10lik, int nparam);

/** Sample from a joint distribution with independent components.
 *
 * This samples a point (x_1, ..., x_n) from an n-dimensional distribution,
 * where each component is fully described by its marginal. The parameters for
 * each distribution are:
 *
 *  - GAUSSIAN: mean, variance
 *  - UNIFORM: start, end
 *  - DELTA: value
 *
 * \param[in] n the number of dimensions in the distribution
 * \param[out] theta place to put the sampled values
 * \param[in] rng GSL random number generator object
 * \param[in] dist distributions for each component
 * \param[in] params parameters for the distributions in dist
 */
void sample_distribution(int n, double *theta, gsl_rng *rng, 
                         const distribution *dist, double **params);

/** Calculate the density of a point from its distribution.
 *
 * Given a point (x_1, ..., x_n), where each x_i was independently drawn from a
 * separate probability distribution, calculate the density at the point.
 *
 * \param[in] n the number of dimensions in the distribution
 * \param[in] theta point to calculate density at
 * \param[in] dist distributions for each component
 * \param[in] params parameters for the distributions in dist
 * \return the density at the point
 * \sa sample_distribution()
 */
double density_distribution(int n, const double *theta, 
                            const distribution *dist, double **params);

/** Parse a distribution from a string.
 *
 * This takes a string (like "exponential") and return the corresponding
 * element of the distribution enum, or 0 if the string does not match any
 * known distribution.
 *
 * \param[in] s string to parse
 * \return the distribution named by s
 */
distribution parse_distribution(const char *s);

#endif
