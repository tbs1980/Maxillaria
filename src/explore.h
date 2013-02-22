#ifndef MAXILLARIA_EXPLORE_H
#define MAXILLARIA_EXPLORE_H

#include"live_point.h"
#include"mt19937.h"

/*
 * Sample the parameter space uniformly. Very slow!
 */ 
void explore_prior_space(live_point* livpnt,double llstar,unsigned num_dim,
unsigned num_par,ellipsis_mt19937_rng* rng,
void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew));

/*
 * Explore the parameter space with MCMC.
 * See Sivia & Skilling 2006, page 193 for more details
 */
void explore_prior_space_with_mcmc(live_point* livpnt,double llstar,unsigned num_dim,
	unsigned num_par,ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew));



#endif/*MAXILLARIA_EXPLORE_H*/
