#ifndef SMC_H
#define SMC_H

#include <gsl/gsl_rng.h>
#include "stats.h"

#define MAX_DIST_PARAMS 2

typedef struct {
    int nparam; /**< number of parameters in the model */
    int nparticle; /**< number of particles used to approximate the posterior */
    int nsample; /**< number of sampled data points per particle */
    int ess_tolerance; /**< ESS below this value triggers resampling */

    double final_epsilon; /**< tolerance level to end at */
    double final_accept_rate; /**< MCMC acceptance rate to stop at */
    double quality; /**< between 0 and 1, where 0 is fast and coarse, 1 is slow and accurate */
    double step_tolerance; /**< tolerance for bisection solution of next epsilon */

    size_t dataset_size; /**< size of data objects */
    size_t feedback_size; /**< size of feedback objects */

    void *propose_arg; /**< extra argument for propose function */
    void *proposal_density_arg; /**< extra argument for proposal_density function */
    void *sample_dataset_arg; /**< extra argument for sample_dataset function */
    void *distance_arg; /**< extra argument for distance function */
    void *feedback_arg; /**< extra argument for feedback function */
    void *sample_from_prior_arg; /**< extra argument for sample_from_prior function */
    void *prior_density_arg; /**< extra argument for sample_from_prior function */
} smc_config;

typedef struct {
    void   (*propose)           (gsl_rng *rng, double *, const void *, const void *);
    double (*proposal_density)  (const double *, const double *, const void *, const void *);
    void   (*sample_dataset)    (gsl_rng *rng, const double *, const void *, void *);
    double (*distance)          (const void *, const void *, const void *);
    void   (*feedback)          (const double *, int, void *, const void *);
    void   (*destroy_dataset)   (void *);
    void   (*sample_from_prior) (gsl_rng *rng, double *, const void *);
    double (*prior_density)     (double *, const void *);
} smc_functions;

typedef struct {
    int niter;
    double *epsilon;
    double *acceptance_rate;
    double **theta;
    double **W;
} smc_result;

/** Perform ABC-SMC.
 *
 * This implements the adaptive tolerance algorithm from DelMoral et al. 2012.
 * The object returned must be passed to smc_result_free() when you are done
 * with it.
 *
 * \param[in] config control parameters for the algorithm
 * \param[in] functions functions to use in the algorithm
 * \param[in] seed random seed (if negative, use time)
 * \panam[in] nthread number of threads to use
 * \param[in] data true data
 * \param[in] trace_file file to record particle populations to at each iteration
 * \return an smc_result object with the posterior distribution (theta) and
 * other information about the run
 */
smc_result *abc_smc(const smc_config config, const smc_functions functions,
                    int seed, int nthread, const void *data, FILE *trace_file);

/** Free an smc_result object.
 *
 * \param[in] r object to free
 */
void smc_result_free(smc_result *r);

#endif
