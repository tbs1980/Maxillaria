#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"post_stats.h"

void read_post_samples(char* file_name,live_point* psamps,unsigned long num_post_samps,unsigned num_dim)
{
	FILE* infile;
	unsigned long i;
	unsigned j;
	double val;
	infile=fopen(file_name,"r");
	
	if(infile!=NULL)
	{
		for(i=0;i<num_post_samps;++i)
		{
			if(feof(infile))
				break;
			for(j=0;j<num_dim;++j)
			{
				fscanf(infile,"%le",&val);
				psamps[i].x[j]=val;			
			}
			fscanf(infile,"%le",&val);
			psamps[i].log_lik=val;
			fscanf(infile,"%le",&val);
			psamps[i].log_weight=val;
		
		}
		
		if(i!=num_post_samps)
		{
			printf("\t ERROR! Number of posterior samples in %s (%ld) \n\t!= Number of posterior samples expected (%ld) (aborting!)\n",file_name,i,num_post_samps);
			exit(-1);
		}
	
		fclose(infile);
	}
	else
	{
		printf("Unable to open file: %s\n",file_name);
		exit(-1);
	}
}

void print_post_stats(live_point* psamps,unsigned long num_post_samps,double logz,unsigned num_dim)
{
	double* mom1=(double*)malloc(num_dim*sizeof(double));
	double* mom2=(double*)malloc(num_dim*sizeof(double));
	double* weight=(double*)malloc(num_post_samps*sizeof(double));
	unsigned i,j;
	
	/*memset(mom1,0,num_dim*sizeof(double));*/
	/*memset(mom2,0,num_dim*sizeof(double));*/
	
	for(i=0;i<num_dim;++i)
	{
		mom1[i]=0.;
		mom2[i]=0.;
	}
	
	for(i=0;i<num_post_samps;++i)
	{
		weight[i]=exp(psamps[i].log_weight-logz);
		for(j=0;j<num_dim;++j)
		{
			mom1[j]+=weight[i]*psamps[i].x[j];
			mom2[j]+=weight[i]*psamps[i].x[j]*psamps[i].x[j];
		}
	}	
	
	printf("\n\nPosterior Statistics\n");
	/*printf("--------------------------------------------------------------------------\n");*/
	printf("Dimension           Mean                            Standard Deviation\n");
	printf("--------------------------------------------------------------------------\n");
	for(i=0;i<num_dim;++i)
	{
		printf("%4d                %12.10e                %12.10e\n",i,mom1[i],sqrt(mom2[i]-mom1[i]*mom1[i]));
	}
	/*printf("--------------------------------------------------------------------------\n");*/
	
	
	free(mom1);
	free(mom2);
	free(weight);
}

void print_max_like_estimate(live_point* lvpnts,unsigned num_live,unsigned num_dim)
{
	unsigned best=0;
	unsigned i;
	for(i=0;i<num_live;++i)
	{
		if(lvpnts[i].log_lik>lvpnts[best].log_lik)
		{
			best=i;
		}
	}
	printf("\nMax Like parameters\n");
	printf("Dimension           Value\n");
	printf("--------------------------------------------------------------------------\n");
	for(i=0;i<num_dim;++i)
	{
		printf("%4d                %12.10e\n",i,lvpnts[best].x[i]);
	}
	
}


/*
 * Write sampling and posterior statistics to a file 
 */
void write_stats_file(live_point* psamps,unsigned long num_post_samps,live_point* lvpnts,
	unsigned num_live,unsigned num_dim,double logz,double dlogz,double H,char* stats_file_name)
{
	double* mom1=(double*)malloc(num_dim*sizeof(double));
	double* mom2=(double*)malloc(num_dim*sizeof(double));
	double* weight=(double*)malloc(num_post_samps*sizeof(double));
	unsigned i,j;
	unsigned best=0;
	FILE* outfile;
	
	/*memset(mom1,0,num_dim*sizeof(double));*/
	/*memset(mom2,0,num_dim*sizeof(double));*/
	
	/*calculate the mean and stadnard deviation*/
	for(i=0;i<num_dim;++i)
	{
		mom1[i]=0.;
		mom2[i]=0.;
	}
	
	for(i=0;i<num_post_samps;++i)
	{
		weight[i]=exp(psamps[i].log_weight-logz);
		for(j=0;j<num_dim;++j)
		{
			mom1[j]+=weight[i]*psamps[i].x[j];
			mom2[j]+=weight[i]*psamps[i].x[j]*psamps[i].x[j];
		}
	}	
	
	/*find the maximum likelihood */
	for(i=0;i<num_live;++i)
	{
		if(lvpnts[i].log_lik>lvpnts[best].log_lik)
		{
			best=i;
		}
	}
	
	/*write statistics */
	outfile=fopen(stats_file_name,"w");	
	if(outfile!=NULL)
	{
		fprintf(outfile,"Number of iterates               = %d\n",num_post_samps);
		fprintf(outfile,"Evidence: log(Z)                 = %g +/- %g\n",logz,dlogz);
		fprintf(outfile,"Information: H                   = %g nats(= %g bits)\n",H,H/log(2.));
		fprintf(outfile,"\n\nPosterior Statistics\n");
		fprintf(outfile,"Dimension           Mean                            Standard Deviation\n");
		fprintf(outfile,"--------------------------------------------------------------------------\n");
		for(i=0;i<num_dim;++i)
		{
			fprintf(outfile,"%4d                %12.10e                %12.10e\n",i,mom1[i],sqrt(mom2[i]-mom1[i]*mom1[i]));
		}
		
		fprintf(outfile,"\nMaximum Likelihood Parameters\n");
		fprintf(outfile,"Dimension           Value\n");
		fprintf(outfile,"--------------------------------------------------------------------------\n");
		for(i=0;i<num_dim;++i)
		{
			fprintf(outfile,"%4d                %12.10e\n",i,lvpnts[best].x[i]);
		}
		
		fclose(outfile);
	}
	else
	{
		printf("Unable to open file: %s (aborting)\n",stats_file_name);
	}
}
