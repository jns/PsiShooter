#ifndef PSI_SHOOTER_H
#define PSI_SHOOTER_H

typedef struct {
	double energy;
	PS_DATA wavefunction;
	} PS_SOLUTION;

typedef struct {
	double energy_min; // The minimum energy to try
	double energy_max; // The maximum energy to try
	int n_iter;
	} PS_SOLVE_PARAMETERS;

#endif