#include <stdio.h>
#include "ps_constants.h"
#include "ps_errors.h"
#include "ps_data.h"

/**
 * PsiShooter program entry point
 */
int main(int argc, char **argv) {
	printf("PsiShooter -- a shooting method solver for the time independant Schrodinger equation under the effective mass approximation.\n");
	
	char LOG_FILENAME_V[] = "V.dat"; //Output file for the potential energy profile, V(x)
    char LOG_FILENAME_x[] = "x.dat"; //Output file for the indeptendant variable, x
    char LOG_FILENAME_F[] = "F.dat"; //Output file for the envelope function (wavefunction)
    char LOG_FILENAME_G[] = "G.dat"; //Output file for the auxillary function coupled to F
    char LOG_FILENAME_BS[] = "BS.dat"; //Output file for the envelope function (wavefunction) of bound states
    char LOG_FILENAME_E[] = "E.dat"; //Output file for the eigenvalues corresponding to bound states
    

	
	
	return 0;
}






/**
 * Generates a test potential. Not really a permenant feature, but rather its 
 * something I intend to use in order to do some useful solver engine work 
 * while the generation of potentials from files/user input is up in air.
 */
 PS_DATA test_potential_1D() {
 	
	/*    
	 * Lets define a single well as a test case for the BV/shooter method
	 *   1: V=Vb
	 *   2: V=0
	 *   3: V=Vb
	 *
	 * The thickness of well (region 2) is "well_width"
	 * The thickness of the barriers (regions 1 & 3) is "barrier_width"
	 *    
	 * V=Vb       _______     _______
	 *                   |   |    
	 *                   |   |        
	 * V=0               |___|    
	 * Region:      1      2     3  
	 */
	 
	 
	 
 	PS_DATA V;
 	
	// V.xsize = 1000;
 	//typedef struct {
	//unsigned int xsize;
	//unsigned int ysize;
	//double xstep;
	//double ystep;
	//double **data;
//} 

 	
 }