#include<stdio.h>
#include<stdlib.h>

#include"nested_sampler.h"
#include"priors.h"

#define NUM_DIM 2

int myid = 0;
double glog_z;

void log_lik(double* cube,unsigned num_dim, unsigned num_par, double *lnew);


void my_dumper(double* log_z)
{
	glog_z=*log_z;
}


void log_lik(double* cube,unsigned num_dim, unsigned num_par, double *lnew)
{
	/* 2-D Gaussian Distribution with mu=0 and sigma=1*/
	double priors[NUM_DIM][2];
	int prior_type[NUM_DIM];
	double para[NUM_DIM];
	double val;
	double loglik=0;
	unsigned i;
	
	prior_type[0] = 0;  priors[0][0] = -15.;	priors[0][1] = 15.;
	prior_type[1] = 0;  priors[1][0] = -15.;	priors[1][1] = 15.;

	for(i=0;i<num_dim;++i)
	{
		val=apply_prior(prior_type[i], cube[i], priors[i]);
		cube[i]=para[i]=val;
	}
	
	for(i=0;i<num_dim;++i)
	{
		loglik+=(0.-para[i])*(0.-para[i]);
	}
	*lnew=-0.5*loglik;	
}

int main(int argc, char *argv[])
{
	unsigned nlive=1000;
	unsigned ndim=NUM_DIM;
	unsigned npar=NUM_DIM;
	double ztol=1e-4;
	unsigned update_int=100;
	unsigned long seed=12457;
	char* root="chains/test";
	int feedback=0;
	int resume=0;
	
	run_maxillaria_nested_sampler(nlive,ndim,npar,ztol,update_int,seed,root,feedback,resume,&log_lik,&my_dumper);
	
	return 0;
	 
}
