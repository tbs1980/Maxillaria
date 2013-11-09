#ifndef MAXILLARIA_NESTED_DIAGNOSTIC_DATA_H
#define MAXILLARIA_NESTED_DIAGNOSTIC_DATA_H

typedef struct
{
	double val;
	int keep_sampling;
	
} nested_diagnostic_data;

void init_diagnostic_data(nested_diagnostic_data* nested_diag_data);
void read_diagnostic_data_from_file(nested_diagnostic_data* nested_diag_data,char* diag_file_name);

#endif/*MAXILLARIA_NESTED_DIAGNOSTIC_DATA_H*/
