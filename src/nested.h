#ifndef CATTLEYA_NESTED_H
#define CATTLEYA_NESTED_H

#include<string.h>
#include<assert.h>

#ifdef __INTEL_COMPILER
	#define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__
	#if  __GNUC_MINOR__ > 4
		#define NESTRUN __nested_MOD_nestrun
	#else
		#define NESTRUN __nested__nestrun
	#endif
#else
       #error COMPIPLER_PREFIX_NOT_KNOWN
#endif

#define FSTRING_LENGTH 128

void NESTRUN(int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
	int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, char *root, int *seed,
	int *pWrap, int *fb, int *resume, void (*Loglike)(double *Cube, int *n_dim, int *n_par, double *lnew),
	void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *), int *root_len);
	
void run_multinest_sampler(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, int maxModes,
	int updInt, double Ztol, char* root_file_name, int seed, int *pWrap, int fb, int resume,
	void (*LogLike)(double *Cube, int *n_dim, int *n_par, double *lnew),
	void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *))
{
	int fstring_len;
	int cpylen;
	int i;
	int root_len;
	char f_chroot[FSTRING_LENGTH];
	fstring_len=FSTRING_LENGTH;
	
	cpylen=strlen(root_file_name);
	assert(cpylen <= fstring_len);
	memcpy(f_chroot,root_file_name,cpylen);
	for(i=cpylen;i <fstring_len;++i)
	{
		f_chroot[i]=' ';
	}
	root_len = strlen(f_chroot);
	
	NESTRUN(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
		&f_chroot[0], &seed, pWrap, &fb, &resume, LogLike, dumper, &root_len);
}


#endif/*CATTLEYA_NESTED_H*/
