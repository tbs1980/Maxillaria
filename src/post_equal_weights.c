#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"utils.h"
#include"live_point.h"
#include"post_equal_weights.h"
#include"mt19937.h"

void calc_post_equal_weights_samples(live_point* psamps,unsigned num_post_samps,double logz,unsigned num_dim,char* file_name)
{
	unsigned i,k;
	unsigned num_eff_posts;
	ellipsis_mt19937_rng* rand;
	unsigned long rand_init[4]={0x123, 0x234, 0x345, 0x456};
	int rand_length=4;
	double uni_rand;
	double S_k;
	double eff_log_num_samples=0;
	double* weight=(double*)malloc(num_post_samps*sizeof(double));
	double sum_wt;
	live_point* post_equal_weight_samples;

	/* allocate memroy for random number generator*/
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));
	init_by_array(rand,rand_init, rand_length);
	
	/*calculate effective number of samples */
	for(i=0;i<num_post_samps;++i)
	{
		weight[i]=exp(psamps[i].log_weight-logz);
		if(weight[i]>0)
		{
			eff_log_num_samples-=weight[i]*log(weight[i]);
		}
	}	
	num_eff_posts=(unsigned)floor((double)exp(eff_log_num_samples));
	
	/*check if it is less than number of posterior samples*/
	if(num_post_samps<num_eff_posts)
		num_eff_posts=num_post_samps;
	
	/*allocate eqully weighted posterior samples */	
	post_equal_weight_samples=(live_point*)malloc(num_eff_posts*sizeof(live_point));
	for(i=0;i<num_eff_posts;++i)
	{
		init_live_point(&post_equal_weight_samples[i],num_dim);
	}
	
	/*staircase sampling */
	sum_wt=0;
	k=0;
	for(i=0;i<num_post_samps;++i)
	{
		uni_rand=genrand_uniform(rand);
		sum_wt+=weight[i];
		S_k=uni_rand+num_eff_posts*sum_wt;
		if(S_k-(double)k > 1.)
		{
			copy_live_point(&post_equal_weight_samples[k],&psamps[i]);
			++k;
			
		}
	}
	
	/* once again check if have the expected number of posterior samples */
	if(k<num_eff_posts)
		num_eff_posts=k;
	
	/*write to txt file*/
	write_live_points_to_txt_file(file_name,post_equal_weight_samples,num_eff_posts,num_dim);
	
	/*clean up*/
	free(rand);
	free(weight);
	free(post_equal_weight_samples);
	
}
