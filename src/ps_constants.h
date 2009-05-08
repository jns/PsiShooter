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


const double m_e = 5.6856269e-13; //[meV s^2 cm^-2], electron mass
const double h = 4.1356668e-12; //[meV s], Planck's Constant
const double h_bar = 6.582118e-3; //[meV s], reduced Planck's constant, h_bar = h/(2*pi)

const double m_eff_GaAs = 0.067; // [dimensionless], the effective mass of electrons at the conduction band minimum in GaAs relative to the rest mass of an electrion




#endif 