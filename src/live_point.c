#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>

#include"live_point.h"

void init_live_point(live_point* lp, unsigned num_dim)
{
	lp->x=(double*)malloc(num_dim*sizeof(double));
	lp->u=(double*)malloc(num_dim*sizeof(double));
	lp->log_lik=-1e90;
	lp->log_weight=0;
	lp->num_dim=num_dim;
}

void free_live_point(live_point* lp)
{
	if(lp->x != NULL)
		free(lp->x);
	if(lp->u != NULL)
		free(lp->u);
	lp->log_lik=-1e90;
	lp->log_weight=0;
	lp->num_dim=0;
}

void copy_live_point(live_point* destination,live_point* source)
{
	unsigned i;
	assert(destination->num_dim==source->num_dim);
	/*memcpy(destination->x,source->x,source->num_dim*sizeof(double));*/
	for(i=0;i<destination->num_dim;++i)
	{
		destination->u[i]=source->u[i];
		destination->x[i]=source->x[i];
	}
	destination->log_lik=source->log_lik;
	destination->log_weight=source->log_weight;
}


