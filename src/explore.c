#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>

#include"explore.h"
#include"mt19937.h"
/*
 * Sample the parameter space uniformly. Very slow!
 */ 
void explore_prior_space(live_point* livpnt,double llstar,unsigned num_dim,
unsigned num_par,ellipsis_mt19937_rng* rng,
void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew))
{
	/*Evolve the object within the likelihood constraint*/
	
	unsigned i;
	double llik;
	
	live_point* trial_lvpnt=(live_point*)malloc(sizeof(live_point));
	init_live_point(trial_lvpnt,num_dim);
	
	for(;;)
	{
	
		for(i=0;i<num_dim;++i)
		{
			trial_lvpnt->x[i]=genrand_uniform(rng);
		}
			
		log_likelihood(trial_lvpnt->x,num_dim,num_par,&llik);
		trial_lvpnt->log_lik=llik;
		
		/*accpet if and only if within l > lstar */
		if(trial_lvpnt->log_lik > llstar)
		{
			copy_live_point(livpnt,trial_lvpnt);
			break;
		}
	}
	
	free(trial_lvpnt);
	
}


/*
 * Explore the parameter space with MCMC.
 * See Sivia & Skilling 2006, page 193 for more details
 */
void explore_prior_space_with_mcmc(live_point* livpnt,double llstar,unsigned num_dim,
	unsigned num_par,ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew))
{
	/*Evolve the object within the likelihood constraint*/
	
	unsigned num_samples=20;
	double step=0.1;
	unsigned accept=0;
	unsigned reject=0;
	unsigned i;
	double llik;
	unsigned samps=0;
	double* uvals;
	double* uvals_new;
	double* xvals;
	double* xvals_new;
	double llik_new;
	
	uvals=(double*)malloc(num_dim*sizeof(double));
	uvals_new=(double*)malloc(num_dim*sizeof(double));
	xvals=(double*)malloc(num_dim*sizeof(double));
	xvals_new=(double*)malloc(num_dim*sizeof(double));

	
	for(i=0;i<num_dim;++i)
	{
		uvals[i]=livpnt->u[i];
		xvals[i]=livpnt->x[i];
		llik_new=livpnt->log_lik;
	}
	
	
	for(;;)
	{
	
		for(i=0;i<num_dim;++i)
		{
			uvals_new[i]=uvals[i]+step*(2.*genrand_uniform(rng)-1.);
			uvals_new[i]-=floor(uvals_new[i]);
			xvals_new[i]=uvals_new[i];
		}
			
		log_likelihood(xvals_new,num_dim,num_par,&llik);
		
		/*accpet if and only if within l > lstar */
		if(llik > llstar)
		{
			for(i=0;i<num_dim;++i)
			{
				uvals[i]=uvals_new[i];
				xvals[i]=xvals_new[i];
			}
			llik_new=llik;
			++accept;
		}
		else
		{
			++reject;
		}
		
		/*refine step size to let acceptance ratio converge around 50% */
		if(accept>reject)
		{
			step*=exp(1./(double)accept);
		}
		if(accept<reject)
		{
			step/=exp(1./(double)reject);
		}
		
		
		if(accept>0 && samps>=num_samples)
		{
			for(i=0;i<num_dim;++i)
			{
				livpnt->u[i]=uvals[i];
				livpnt->x[i]=xvals[i];
			}
			livpnt->log_lik=llik_new;
			break;
		}
		++samps;
	}
	
	free(uvals);
	free(uvals_new);
	free(xvals);
	free(xvals_new);
}


/*
 *This function implements an approximation to Mukherjee's Ellipsoidal methods.
 * EXPERIMENTAL -> DO NOT USE
 */
void explore_prior_space_with_mcmc_var_tuned(live_point* livpnt,double llstar,
	live_point* phys_live_points,unsigned num_phys_live,unsigned num_dim,unsigned num_par,	ellipsis_mt19937_rng* rng,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew))
{
	/*Evolve the object within the likelihood constraint*/
	
	unsigned num_samples=20;
	unsigned i,j;
	double llik;
	unsigned samps=0;
	unsigned accept=0;
	unsigned reject=0;
	double step=2.;
	double* sigma;
	double* mean;
	double sum;
	double sum2;
	double enlarge_fact=2.;
	double* uvals;
	double* uvals_new;
	double* xvals;
	double* xvals_new;
	double llik_new;
	unsigned best=0;
	
	/*From the available live points, estimate the variance in each dimension*/
	sigma=(double*)malloc(num_dim*sizeof(double));
	mean=(double*)malloc(num_dim*sizeof(double));
	
	for(i=0;i<num_phys_live;++i)
	{
		if(phys_live_points[i].log_lik>phys_live_points[best].log_lik)
		{
			best=i;
		}
	}
	
	/*printf("\nbest= %d\tval=%f\n",best,phys_live_points[best].log_lik);*/
	
	for(i=0;i<num_dim;++i)
	{
		sum=0;
		sum2=0;
		for(j=0;j<num_phys_live;++j)
		{
			sum+=phys_live_points[j].u[i];
			sum2+=phys_live_points[j].u[i]*phys_live_points[j].u[i];
		}
		sigma[i]=sqrt( (sum2-sum*sum/(double)num_phys_live)/(double)num_phys_live );
		/*mean[i]=sum/(double)num_phys_live;*/
		mean[i]=(phys_live_points[best].u[i]+sum/(double)num_phys_live)*0.5;
		/*printf("%f\t",mean[i]);*/
	}
	/*printf("\n");*/
	
	uvals=(double*)malloc(num_dim*sizeof(double));
	uvals_new=(double*)malloc(num_dim*sizeof(double));
	xvals=(double*)malloc(num_dim*sizeof(double));
	xvals_new=(double*)malloc(num_dim*sizeof(double));
	
	/*printf("\nllstar= %f\n",llstar);*/
	

	for(i=0;i<num_dim;++i)
	{
		uvals[i]=livpnt->u[i];
		xvals[i]=livpnt->x[i];
		llik_new=livpnt->log_lik;
	}
	
	
	for(;;)
	{
	
		for(i=0;i<num_dim;++i)
		{
			uvals_new[i]=(uvals[i]+mean[i])*0.5+sigma[i]*genrand_uniform(rng);
			xvals_new[i]=uvals_new[i];
		}
			
		log_likelihood(xvals_new,num_dim,num_par,&llik);
		
		/*trial_lvpnt->log_lik=llik*/
		
		/*accpet if and only if within l > lstar */
		if(llik > llstar)
		{
			for(i=0;i<num_dim;++i)
			{
				uvals[i]=uvals_new[i];
				xvals[i]=xvals_new[i];
			}
			llik_new=llik;
			++accept;
		}
		else
		{
			++reject;
		}
		
		/*refine step size to let acceptance ratio converge around 50% */
		if(accept>reject)
		{
			step*=exp(1./(double)accept);
		}
		if(accept<reject)
		{
			step/=exp(1./(double)reject);
		}
		
		if(accept>0 && samps>=num_samples)
		{
			/*copy_live_point(livpnt,trial_lvpnt);*/
			for(i=0;i<num_dim;++i)
			{
				livpnt->u[i]=uvals[i];
				livpnt->x[i]=xvals[i];
			}
			livpnt->log_lik=llik_new;
			break;
		}
		++samps;
		
		/*printf("accept=%d\tsamps= %d\n",accept,samps);*/
	}
	
	free(uvals);
	free(uvals_new);
	free(xvals);
	free(xvals_new);
	
	free(sigma);
	free(mean);
	
}
