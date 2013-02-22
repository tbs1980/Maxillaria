#ifndef MAXILLARIA_POST_STATS_H
#define MAXILLARIA_POST_STATS_H

#include"live_point.h"

void read_post_samples(char* file_name,live_point* psamps,unsigned long num_post_samps,unsigned num_dim);

/*print posterior statistics*/
void print_post_stats(live_point* psamps,unsigned long num_post_samps,double logz,unsigned num_dim);

void print_max_like_estimate(live_point* lvpnts,unsigned num_live,unsigned num_dim);

/*
 * Write sampling and posterior statistics to a file 
 */
void write_stats_file(live_point* psamps,unsigned long num_post_samps,live_point* lvpnts,
	unsigned num_live,unsigned num_dim,double logz,double dlogz,double H,char* stats_file_name);
#endif/*MAXILLARIA_POST_STATS_H*/
