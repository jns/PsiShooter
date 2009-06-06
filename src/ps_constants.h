#ifndef PS_CONSTANTS_H
#define PS_CONSTANTS_H

//New units!
// electron mass []
#define MASS_ELECTRON 0.5d

//Planck's Constant, 2*pi*hbar
#define H_PLANCK 6.2831853071795864769232867665590037683936d

//reduced Planck's constant, h_bar = h/(2*pi)
#define HBAR_PLANCK  1d

// Reduced Planck's Constant Squared [eV s]^2
//#define HBAR_PLANCK_SQ 4.33242792353826e-31

// Mass of Electron / Hbar^2 [nm^-2 ev-1]
//#define G_COEFF 13.1235713792476

// Conversion factor ev -> ergs
#define EV_TO_ERGS 1.60217646e-12

// [dimensionless], the effective mass of electrons at the conduction band minimum in GaAs relative to the rest mass of an electrion
#define M_EFF_GAAS  0.067




/**
 * This header file defines necessary physical constants
 *
 * All constants are defined in terms of the following units:
 * Energy		eV, electron volts
 * Length		nm, nanometers
 * Time			s, seconds
 * Mass			?, is thus defined as [eV s^2 nm^-2] I guess.
 */ 


/* 
// electron mass [ev s^2 nm^-2]
#define MASS_ELECTRON 5.68569271e-12

//[eV s], Planck's Constant
#define H_PLANCK  4.1356668e-15

//[eV s], reduced Planck's constant, h_bar = h/(2*pi)
#define HBAR_PLANCK  6.582118e-16 

// Reduced Planck's Constant Squared [eV s]^2
#define HBAR_PLANCK_SQ 4.33242792353826e-31

// Mass of Electron / Hbar^2 [nm^-2 ev-1]
#define G_COEFF 13.1235713792476

// Conversion factor ev -> ergs
#define EV_TO_ERGS 1.60217646e-12

// [dimensionless], the effective mass of electrons at the conduction band minimum in GaAs relative to the rest mass of an electrion
#define M_EFF_GAAS  0.067
*/



// ... ah the good old erg days of PsiShooter. This is why physicists have a hard time making friends
// Planck's Constant [erg s] 
//#define H_PLANCK 6.62620e-27

//  Reduced Planck's Constant  [erg s]
// #define HBAR_PLANCK 1.05459e-27

// Reduced Planck's Constant Squared [erg s]^2
//#define HBAR_PLANCK_SQ 1.1121600681e-54

//[meV s^2 cm^-2], electron mass
//#define MASS_ELECTRON  5.6856269e-13 

// electron mass [g]
//#define MASS_ELECTRON 9.10956e-28

#endif 