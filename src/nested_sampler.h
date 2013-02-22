#ifndef MAXILLARIA_NESTED_SAMPLER_H
#define MAXILLARIA_NESTED_SAMPLER_H

#include"live_point.h"
#include"mt19937.h"

void run_maxillaria_nested_sampler(unsigned num_live,unsigned num_dim,unsigned num_par,
	double z_tol,unsigned update_interval,unsigned long seed,
	char* prefix_to_files,int feedback,int resume,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew),
	void (*dumper)(double* log_z));
/* just uniform sampling*/	
void explore_prior_space(live_point* livpnt,double llstar,unsigned num_dim,
	unsigned num_par,ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew));


/* mcmc exploration */
void explore_prior_space_with_mcmc(live_point* livpnt,double llstar,unsigned num_dim,
	unsigned num_par,ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew));


/*
 *This function implements an approximation to Mukherjee's Ellipsoidal methods.
 */
void explore_prior_space_with_mcmc_var_tuned(live_point* livpnt,double llstar,
	live_point* phys_live_points,unsigned num_phys_live,unsigned num_dim,unsigned num_par,	ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew));

double log_add(double x,double y);

void print_post_stats(live_point* psamps,unsigned long num_post_samps,double logz,unsigned num_dim);


#endif/*MAXILLARIA_NESTED_SAMPLER_H*/
