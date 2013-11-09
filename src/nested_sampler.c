#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "live_point.h"
#include "maxillaria.h"
#include "nested_sampler.h"
#include "nested_diagnostic.h"
#include "mt19937.h"
#include "utils.h"
#include "post_stats.h"
#include "post_equal_weights.h"
#include "explore.h"

double log_add(double x,double y)
{
	/*printf("x= %e\ty= %e\n",x,y);*/
	/*printf("log(1.+exp(x-y)= %e\n",log(1.+exp(x-y)));*/
	return (x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) );
}

void run_maxillaria_nested_sampler(unsigned num_live,unsigned num_dim,unsigned num_par,
	double z_tol,unsigned update_interval,unsigned long seed,
	char* prefix_to_files,int feedback,int resume,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew),
	void (*dumper)(double* log_z))
{
	unsigned long iter;
	char diag_file_name[128];
	char rand_file_name[128];
	char phys_live_file_name[128];
	char post_equal_weights_file_name[128];
	char extract_file_name[128];
	char stats_file_name[128];
	FILE* diag_out_file;
	FILE* rand_out_file;
	FILE* ext_out_file;
	live_point* live_points;
	live_point* post_samples;
	live_point* all_samples;
	unsigned i,j;
	ellipsis_mt19937_rng* rand;
	unsigned long rand_init[4]={0x123, 0x234, 0x345, 0x456};
	unsigned long rand_length=4;
	nested_diagnostic_data* nested_diag_data;
	int run_sampler=1;
	/*int perform_sampling=1;*/
	double log_lik_new;
	double log_width;
	double log_z_new;
	double log_z=-DBL_MAX;
	double delta_log_z;
	double log_lik_star;
	double info_h=0.;
	unsigned worst;
	unsigned copy;
	unsigned isamp;
	
	/*const unsigned long num_post_samps= (update_interval>0 ? update_interval : 1);*/
	
	assert(update_interval>0);
	
	
	printf("\n-------------------------------------------------\n");
	printf("             Sree's Nested Sampler\n");
	printf("                Version %d.%d\n",MAXILLARIA_VERSION_MAJOR,MAXILLARIA_VERSION_MINOR);
	printf("   Sreekumar Thaithara Balan (tbs1980@gmail.com)\n");
	printf("-------------------------------------------------\n");

	
	/* allocate memroy for random number generator*/
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));
	init_by_array(rand,rand_init, rand_length);

	/* assign the file names */
	strcpy(diag_file_name,prefix_to_files);
	strcat(diag_file_name,".diag.txt");
	strcpy(rand_file_name,prefix_to_files);
	strcat(rand_file_name,".rand.txt");
	strcpy(phys_live_file_name,prefix_to_files);
	strcat(phys_live_file_name,".phys_live_points.txt");
	strcpy(post_equal_weights_file_name,prefix_to_files);
	strcat(post_equal_weights_file_name,".post_equal_weights.txt");
	strcpy(extract_file_name,prefix_to_files);
	strcat(extract_file_name,".extract.txt");
	strcpy(stats_file_name,prefix_to_files);
	strcat(stats_file_name,".stats.txt");

	

	/* are we resuming from pervious state? */
	if(resume)
	{
		printf("\nResuming from previous state...\n");
		
		/* check if resume files exist files for output */
		diag_out_file=fopen(diag_file_name,"r");
		rand_out_file=fopen(rand_file_name,"r");
		
		if(diag_out_file == NULL || rand_out_file == NULL)
		{
			printf("No resume files foud...\n");
			printf("Starting sampling from scratch...\n");
		
			/* initialise the random number generator */
			init_genrand(rand,seed);				
		}
		else
		{
			/* close the opened files*/
			fclose(diag_out_file);
			fclose(rand_out_file);

			/* read the diagnostic data */
			printf("Reading diagnostic data...\t");
			read_diagnostic_data_from_file(nested_diag_data,diag_file_name);
			printf("Done!\n");

			/* has the algorithm converged already? */
			if(!nested_diag_data->keep_sampling)
			{
				printf("Algorithm has already converged!\n");
				run_sampler=0;
			}
			else
			{
		
				/* read the random number data */
				printf("Reading the random number state...\t");
				read_rand_state(rand,rand_file_name);
				printf("Done!\n");		
			}		
		}		
	}
	else
	{
		printf("Starting sampling from scratch...\n");
		init_genrand(rand,seed);
		run_sampler=1;
	}	

	/*---------------------------run sampler ----------------------------*/
	if(run_sampler)
	{
		ext_out_file=fopen(extract_file_name,"w");
		if(ext_out_file!=NULL)
		{
			fclose(ext_out_file);
		}
		else
		{
			printf("Cannot open extract file for writing. Please check the path\n");
			return;
		}
		
		/* open extract file */
		ext_out_file=fopen(extract_file_name,"a");
		
		/*allocate live points */
		live_points=(live_point*)malloc(num_live*sizeof(live_point));
		for(i=0;i<num_live;++i)
		{
			init_live_point(&live_points[i],num_dim);
		}

		/*allocate posterior samples*/
		post_samples=(live_point*)malloc(update_interval*sizeof(live_point));
		for(i=0;i<update_interval;++i)
		{
			init_live_point(&post_samples[i],num_dim);
		}
	
		/*set prior objects */
		for(i=0;i<num_live;++i)
		{
			for(j=0;j<num_dim;++j)
			{
				live_points[i].x[j]=genrand_uniform(rand);
				live_points[i].u[j]=live_points[i].x[j];
			}
			
			/*caclulate the log lik*/
			log_likelihood(live_points[i].x,num_dim,num_par,&log_lik_new);
			live_points[i].log_lik=log_lik_new;
			
		}
		
		/*outermost interval of prior mass*/
		log_width=log(1.-exp(-1./(double)num_live));
		
		/*set counter*/
		isamp=0;
		iter=0;
		
		/*=====================nested sampling loop==================*/
		for(;;)
		{
			/*worst object in collection, with weight = width*likelihood*/
			worst=0;
			for(i=0;i<num_live;++i)
			{
				if(live_points[i].log_lik<live_points[worst].log_lik)
				{
					worst=i;
				}
			}
			live_points[worst].log_weight=log_width+live_points[worst].log_lik;
			
			/*update Evidence Z and information H*/
			log_z_new=log_add(log_z,live_points[worst].log_weight);
			info_h=exp(live_points[worst].log_weight-log_z_new)*live_points[worst].log_lik
				+exp(log_z-log_z_new)*(info_h+log_z)-log_z_new;
				
			delta_log_z=fabs(log_z_new-log_z);

			/*update log z*/
			log_z=log_z_new;
			
			/*posterior samples*/
			copy_live_point(&post_samples[isamp],&live_points[worst]);
			++isamp;			
						
			/* feedback, update files */
			if(isamp==update_interval)
			{
				
				if(feedback)
				{
					printf("\nNumber of iterates     = %ld\n",iter+1);
					printf("Evidence,log(Z)        = %g +- %g\n",log_z,sqrt(info_h/(double)num_live));
					printf("Dlog(Z)                = %g\n",delta_log_z);
				}
				
				/*write the physical live points to a file same functionaliy as MultiNest (Feroz et al 2008)*/
				write_live_points_to_txt_file(phys_live_file_name,live_points,num_live,num_dim);

				/*write the extract same functionaliy as MultiNest (Feroz et al 2008)*/
				write_live_points_to_txt_file_app(ext_out_file,post_samples,isamp,num_dim);
				
				/*call dumper same functionaliy as MultiNest (Feroz et al 2008)*/
				dumper(&log_z);
				
				/* terminating crieterion */
				if(delta_log_z <z_tol)
				{
					break;
				}

				isamp=0;

			}

			/*kill the worst object in favor of copy of different survivor*/
			do
			{
				copy=(unsigned)(num_live*genrand_uniform(rand))%num_live;
			}
			while(copy==worst && num_live > 1);
			
			log_lik_star=live_points[worst].log_lik;

			copy_live_point(&live_points[worst],&live_points[copy]);
			
			/*evolve the copied object within constraint*/
			/*explore_prior_space(&live_points[worst],log_lik_star,num_dim,num_par,rand,log_likelihood);*/
			explore_prior_space_with_mcmc(&live_points[worst],log_lik_star,num_dim,num_par,rand,log_likelihood);
			/*explore_prior_space_with_mcmc_var_tuned(&live_points[worst],log_lik_star,live_points,num_live,num_dim,num_par,rand,log_likelihood);*/
			
			/*shrink interval*/
			log_width-=1./(double)num_live;
			
			++iter;
		}
		
		/*close the extract file*/
		if(ext_out_file!=NULL)
		{
			fclose(ext_out_file);
		}
		
		/*exit with eidence Z, information H and optional posterior samples*/
		printf("\n-------------------------------------------------\n");
		printf("Sampling finished!\n");
		printf("Iteraions          = %ld\n",iter+1);
		printf("Evidence:log(Z)    = %g +- %g\n",log_z,sqrt(info_h/(double)num_live));
		printf("Information: H     = %f nats(= %f bits)\n",info_h,info_h/log(2.));
		printf("-------------------------------------------------\n");
		
		/*allocate memory for all samples written*/
		all_samples=(live_point*)malloc((iter+1)*sizeof(live_point));
		for(i=0;i<(iter+1);++i)
		{
			init_live_point(&all_samples[i],num_dim);
		}
		
		/*read in all samples */
		read_post_samples(extract_file_name,all_samples,(iter+1),num_dim);

		/*print posterior stats*/
		/*print_post_stats(all_samples,(iter+1),log_z,num_dim);*/
		
		/*print the maximum likelihood parameters*/
		/*print_max_like_estimate(live_points,num_live,num_dim);*/
		
		write_stats_file(all_samples,(iter+1),live_points,num_live,num_dim,log_z,
			sqrt(info_h/(double)num_live),info_h,stats_file_name);
		
		/*write the equally weighted posterior samples by staircase sampling*/
		calc_post_equal_weights_samples(all_samples,(iter+1),log_z,num_dim,post_equal_weights_file_name);
		
		free(all_samples);
		
	}
	
	/*free live points*/
	free(live_points);
	free(post_samples);
	free(rand);
}



