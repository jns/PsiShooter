#ifndef PSI_SHOOTER_H
#define PSI_SHOOTER_H

typedef struct {
	double energy;
	PS_DATA wavefunction;
	} PS_SOLUTION;

typedef struct {
	double energy_min; // The minimum energy to try
	double energy_max; // The maximum energy to try
	int n_iter_coarse; // The number of iterations to try on the coarse pass
	int n_iter_fine; // The number of iterations to try on the fine pass
	} PS_SOLVE_PARAMETERS;

#endif