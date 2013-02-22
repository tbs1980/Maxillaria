#ifndef MAXILLARIA_UTILS_H
#define MAXILLARIA_UTILS_H

#include<stdio.h>
#include"live_point.h"

void write_vector_to_txt_file(FILE* out_file,double* x,unsigned long num_dim);

void write_live_points_to_txt_file(char* file_name,live_point* lvs,unsigned num_live,unsigned num_dim);

void sort_live_points(live_point* lvs,unsigned num_live,unsigned num_dim);

void quick_sort(double *arr, unsigned long *brr, unsigned long n);

void write_live_points_to_txt_file_app(FILE* out_file,live_point* lvs,unsigned num_live,unsigned num_dim);

#endif/*MAXILLARIA_UTILS_H*/
