#ifndef PS_CONSTANTS_H
#define PS_CONSTANTS_H


/**
 * This header file defines necessary physical constants
 *
 * All constants are defined in terms of the following units:
 * Energy		meV, milli electron volts
 * Length		cm, centimeters
 * Time			s, seconds
 * Mass			?, is thus defined as [meV s^2 cm^-2]
 */ 


//[meV s^2 cm^-2], electron mass
#define MASS_ELECTRON = 5.6856269e-13 

//[meV s], Planck's Constant
#define H_PLANCK = 4.1356668e-12 

//[meV s], reduced Planck's constant, h_bar = h/(2*pi)
#define HBAR_PLANCK = 6.582118e-3 

// [dimensionless], the effective mass of electrons at the conduction band minimum in GaAs relative to the rest mass of an electrion
#define M_EFF_GAAS = 0.067




#endif 