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
//#define MASS_ELECTRON  5.6856269e-13 

// electron mass [g]
#define MASS_ELECTRON 9.10956e-28

//[meV s], Planck's Constant
//#define H_PLANCK  4.1356668e-12 

// Planck's Constant [erg s]
#define H_PLANCK 6.62620e-27

//[meV s], reduced Planck's constant, h_bar = h/(2*pi)
//#define HBAR_PLANCK  6.582118e-13 

//  Reduced Planck's Constant  [erg s]
#define HBAR_PLANCK 1.05459e-27

// Conversion factor ev -> ergs
#define EV_TO_ERGS 1.60217646e-12

// [dimensionless], the effective mass of electrons at the conduction band minimum in GaAs relative to the rest mass of an electrion
#define M_EFF_GAAS  0.067




#endif 