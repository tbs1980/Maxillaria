#ifndef MAXILLARIA_LIVE_POINT_H
#define MAXILLARIA_LIVE_POINT_H

typedef struct
{
	double* u;
	double* x;
	double log_lik;
	double log_weight;
	unsigned num_dim;
} live_point;

void init_live_point(live_point* lp, unsigned num_dim);
void free_live_point(live_point* lp);
void copy_live_point(live_point* destination,live_point* source);




#endif/*MAXILLARIA_LIVE_POINT_H*/
