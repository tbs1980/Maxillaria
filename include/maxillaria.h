#ifndef MAXILLARIA_MAXILLARIA_H
#define MAXILLARIA_MAXILLARIA_H

#define MAXILLARIA_VERSION_MAJOR 1
#define MAXILLARIA_VERSION_MINOR 0

void run_maxillaria_nested_sampler(unsigned num_live,unsigned num_dim,unsigned num_par,
	double z_tol,unsigned update_interval,unsigned long seed,
	char* prefix_to_files,int feedback,int resume,
	void (*log_likelihood)(double *cube, unsigned ndim, unsigned npar, double *lnew),
	void (*dumper)(double* log_z));

#endif /*MAXILLARIA_MAXILLARIA_H*/
