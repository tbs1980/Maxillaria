#include"nested_diagnostic.h"

void init_diagnostic_data(nested_diagnostic_data* nested_diag_data)
{
	nested_diag_data->val=0;
	nested_diag_data->keep_sampling=1;
}
void read_diagnostic_data_from_file(nested_diagnostic_data* nested_diag_data,char* diag_file_name)
{
	nested_diag_data->keep_sampling=1;
}
