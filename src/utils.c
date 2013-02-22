#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<assert.h>

#include"utils.h"
#include"live_point.h"

void write_vector_to_txt_file(FILE* out_file,double* x,unsigned long num_dim)
{
	unsigned long i;
	for(i=0;i<num_dim;++i)
	{
		fprintf(out_file,"%#.12e\t",x[i]);
	}
	fprintf(out_file,"\n");
}

void write_live_points_to_txt_file(char* file_name,live_point* lvs,unsigned num_live,unsigned num_dim)
{

	
	FILE* out_file;
	unsigned i,j;
	
	out_file=fopen(file_name,"w");
	
	
	for(i=0;i<num_live;++i)
	{
		for(j=0;j<num_dim;++j)
		{
			fprintf(out_file,"%#.12e\t",lvs[i].x[j]);
		}
		fprintf(out_file,"%#.12e\t",lvs[i].log_lik);
		fprintf(out_file,"%#.12e\t",lvs[i].log_weight);
		fprintf(out_file,"\n");
	}

	if(out_file!=NULL)
		fclose(out_file);
}

/*Append files*/
void write_live_points_to_txt_file_app(FILE* out_file,live_point* lvs,unsigned num_live,unsigned num_dim)
{
	unsigned i,j;
	if(out_file != NULL)
	{
		
		for(i=0;i<num_live;++i)
		{
			for(j=0;j<num_dim;++j)
			{
				fprintf(out_file,"%#.12e\t",lvs[i].x[j]);
			}
			fprintf(out_file,"%#.12e\t",lvs[i].log_lik);
			fprintf(out_file,"%#.12e\t",lvs[i].log_weight);
			fprintf(out_file,"\n");
		}
	}
	else
	{
		printf("Unable to open file: (aborting!)\n");
	}
}

void sort_live_points(live_point* lvs,unsigned num_live,unsigned num_dim)
{
	double* llvals;
	unsigned long* inds;
	unsigned i;
	live_point* temp_lvs;

	/*allocate memory*/
	llvals=(double*)malloc(num_live*sizeof(double));
	inds=(unsigned long*)malloc(num_live*sizeof(unsigned long));
	temp_lvs=(live_point*)malloc(num_live*sizeof(live_point));
	
	/*init live points */
	for(i=0;i<num_live;++i)
	{
		init_live_point(&temp_lvs[i],num_dim);
	}
	
	/* store a temp copy for later use*/
	for(i=0;i<num_live;++i)
	{
		copy_live_point(&temp_lvs[i],&lvs[i]);
	}

	/*sort the lvs according to the log likelihood */
	for(i=0;i<num_live;++i)
	{
		llvals[i]=lvs[i].log_lik;
		inds[i]=(unsigned long)i;	
	}
	
	quick_sort(llvals,inds,num_live);
	
	/*now apply the sorted indices*/
	for(i=0;i<num_live;++i)
	{
		copy_live_point(&lvs[i],&temp_lvs[inds[i]]);
	}
	
	
	free(llvals);
	free(inds);
	free(temp_lvs);
	
}


void quick_sort(double *arr, unsigned long *brr, unsigned long n) 
{
	unsigned long i;
	unsigned long ir=n;
	unsigned long j;
	unsigned long k;
	unsigned long l=0;
	unsigned long M=7;
	unsigned long *istack;
	unsigned long jstack=0;
	unsigned long NSTACK=50;
	double a,dble_temp;
	unsigned long ul_temp,b;
	istack=(unsigned long*)malloc(NSTACK*sizeof(unsigned long));
	
	for (;;)
	{
		if (ir-l < M)
		{
			for (j=l+1;j<=ir;j++)
			{
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=1;i--)
				{
					if (arr[i] <= a)
					{
						break;
					}
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack)
			{
				/*free(istack);*/
				break;
				/*return;*/
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		}
		else
		{
			k=(l+ir) >> 1;
			dble_temp=arr[k];
			arr[k]=arr[l+1];
			arr[l+1]=dble_temp;
			ul_temp=brr[k];
			brr[k]=brr[l+1];
			brr[l+1]=ul_temp;
			if (arr[l] > arr[ir])
			{
				dble_temp=arr[l];
				arr[l]=arr[ir];
				arr[ir]=dble_temp;
				ul_temp=brr[l];
				brr[l]=brr[ir];
				brr[ir]=ul_temp;
			}
			if (arr[l+1] > arr[ir])
			{
				dble_temp=arr[l+1];
				arr[l+1]=arr[ir];
				arr[ir]=dble_temp;
				ul_temp=brr[l+1];
				brr[l+1]=brr[ir];
				brr[ir]=ul_temp;
			}
			if (arr[l] > arr[l+1])
			{
				dble_temp=arr[l];
				arr[l]=arr[l+1];
				arr[l+1]=dble_temp;
				ul_temp=brr[l];
				brr[l]=brr[l+1];
				brr[l+1]=ul_temp;
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			b=brr[l+1];
			for (;;)
			{
				do
				{
					i++;
				}
				while(arr[i] < a);
				do
				{
					j--;
				}
				while (arr[j] > a);
				if (j < i)
				{
					break;
				}
				dble_temp=arr[i];
				arr[i]=arr[j];
				arr[j]=dble_temp;
				ul_temp=brr[i];
				brr[i]=brr[j];
				brr[j]=ul_temp;
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			brr[l+1]=brr[j];
			brr[j]=b;
			jstack += 2;
			if(jstack > NSTACK)
			{
				printf("NSTACK too small in quick_sort()\n");
				free(istack);
				exit(0);
			}
			if (ir-i+1 >= j-l)
			{
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			}
			else
			{
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
}

