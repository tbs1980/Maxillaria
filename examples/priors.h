#ifndef CATTLEYA_PRIORS_H
#define CATTLEYA_PRIORS_H

/* flat prior */
double flat_prior(double r, double x1, double x2);

/* log prior */
double log_prior(double r, double x1, double x2);

double apply_prior(int prior_type, double r, double prior_arg[]);

#endif/*CATTLEYA_PRIORS_H*/
