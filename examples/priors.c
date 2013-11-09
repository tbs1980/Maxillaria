#include<math.h>
#include<assert.h>
#include"priors.h"

double flat_prior(double r,double x1,double x2)
{
	return x1+r*(x2-x1);
}

double log_prior(double r, double x1, double x2)
{
	double lx1, lx2;
	lx1=log(x1);
	lx2=log(x2);
	return exp(lx1+r*(lx2-lx1));
}

double apply_prior(int prior_type,double r,double prior_arg[])
{

	if(prior_type==0)
		return flat_prior(r,prior_arg[0],prior_arg[1]);
	else if(prior_type==1)
		return log_prior(r,prior_arg[0],prior_arg[1]);
	return 0;
}

