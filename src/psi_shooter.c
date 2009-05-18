#include <stdio.h>
#include "ps_constants.h"
#include "ps_errors.h"
#include "ps_data.h"
#include "math.h"

/*The Big To Do list:
 * Consider the supporting the simplifying case where m_eff(x,y) -> m_eff   (constant effective mass)
 *     Why? ... one effective mass for the entire system, this should cut out a bunch of derivatives that need to be calculated
 *     For many potentials electrons don't "see" much of the barriers compared to the wells, so m_eff(x,y) ~= m_eff_well_material 
 *     Maybe only support this case, it has a bunch less terms when you expand out the TISE with fintite differencing 
 *     note: A position dependant mass fucks up the hermicity of H (hamiltonian operator)
 */


//
// Solve for the bound energies given by the potential. 
// energies is a buffer of energies to test.
// bound_energies is where the actual bound energies will be stored.  It should be the same size as energies.
// both buffers should be of size buf_size.
// returns the number of bound energies found, 0 or negative for error
//
////////////////////////
//Solving the 1D TISE //
////////////////////////
//
//(0)    H * Psi = E * Psi
//
//E is the eignenvalue
//H is the hamiltonian operator: H = T + V = p^2/(2m) + V  
//    ... where p is the momentum operator, p = h/i * Div
//
//Moving to the effective mass regime we no longer solve for the wavefunction
//but instead we sove for an envelope function. Careful consideration shows that 
//if mass is allowed to vary with positon, m -> m(x), then simple subsitution of 
//m(x) for m in H causes H to become non-Hermitian. 
//Fortunately that can be worked around, the result is:
//
//(1)    T * F = -h_bar^2/2 * Div * ( m_eff(x,y,E)^-1 * Div * F )
//           ...note: F is the envelope function that is replacing the role of Psi
// 
//Put it all together to get the effective mass equation that we are going to
//use:
//    +--------------------------------------------------------------+                               
//(2) |    -h_bar^2/2 * Div * (1/m_eff * Div * F) + V * F = E * F    |  
//    +--------------------------------------------------------------+----> This is the main equation
//
//Note: Skip to equations (14) & (15) to see how we are going to proceed in the actual solver we are implenting.
//
	//Rearranging (2):
	//(3)    Div * (M * Div * F) = 2*(V-E)*F/h_bar^2
	//           ... note: Let 1/m_eff(x,y,E) = M  (For compactness of these notes only, it is still the effective mass that is position dependant)
	//take the coordinate system to be cartesisian and Div becomes: Div = d/dx + d/dy   -->   Div*F = dF/dx + dF/dy
	//(4)    Div * { M(dF/dx) + M(dF/dy) } = 2*(V-E)*F/h_bar^2
	//  
	//The LHS of (4) has a (2D) divergance of 2 terms that each have a product of two functions (M, dF), both of which depend on x and y. 
	//The chain rule will get invoked and there will be quite a few terms after we expand this
	//(5)    Div*(M(dF/dx)) + Div*(M(dF/dy)) = ...     ...note: I am leaving off the RHS for awhile
	//(6)    (d/dx + d/dy)*(M(dF/dx)) + (d/dx + d/dy)*(M(dF/dy)) = ...
	//(7)    (d/dx)*(M(dF/dx)) + (d/dy)*(M(dF/dx)) + (d/dx)*(M(dF/dy)) + (d/dy)*(M(dF/dy))= ...
	//(8)    (M)(d^2F/dx^2) + (dF/dx)(dM/dx)  + (M)(d^2F/(dxdy)) + (dF/dx)(dM/dy) + (M)(d^2F/(dydx)) + (dF/dy)(dM/dx) + (M)(d^2F/dy^2) + (dF/dy)(dM/dy) = ...
	//(9)    (M)(d^2F/dx^2) + (M)(d^2F/dy^2) + 2(M)(d^2F/(dxdy)) + (dF/dx)(dM/dx) + (dF/dx)(dM/dy) + (dF/dy)(dM/dx) + (dF/dy)(dM/dy) = ...
	//
	//The next step is to expand (9) out in terms of finite differences, here are some pieces that will be used:
	//(10)   dF/dx = (F[x+dx]-F[x-dx])/(2dx)
	//(11)   d^2F/dx^2 = ((F[x+dx+dx]-F[x+dx-dx])/(2dx) - (F[x-dx+dx]-F[x-dx-dx])/(2dx)) / (2dx)
	//                 = (F[x+2dx] - 2F[x] + F[x-2dx])*(2dx)^-2
	//(12)   d^2F/(dxdy) = (d/dy) * dF/dx
	//                   = (d/dy) * (F[x+dx]-F[x-dx])/(2dx)
	//                   = ((F[x+dx,y+dy]-F[x-dx,y+dy])/(2dx) - (F[x+dx,y-dy]-F[x-dx,y-dy])/(2dx)) / (2dy)
	//                   = ((F[x+dx,y+dy]-F[x-dx,y+dy]) - (F[x+dx,y-dy]-F[x-dx,y-dy])) / (2dx2dy)
	//
	//Using the finite differencing versions (10)-(12) of the derivatives in (9) and making M & F explcit functions of x & y: 
	//(9)    (M)(d^2F/dx^2) +         --> (13)    M[x,y] * (F[x+2dx,y] - 2F[x,y] + F[x-2dx,y])*(2dx)^-2 +  
	//       (M)(d^2F/dy^2) +         -->         M[x,y] * (F[x,y+2dy] - 2F[x,y] + F[x,y-2dy])*(2dy)^-2 + 
	//       2(M)(d^2F/(dxdy)) +      -->         2*M[x,y] * ((F[x+dx,y+dy]-F[x-dx,y+dy])-(F[x+dx,y-dy]-F[x-dx,y-dy]))/(2dx2dy) +
	//       (dF/dx)(dM/dx) +         -->         (F[x+dx,y]-F[x-dx,y])/(2dx) * (M[x+dx,y]-M[x-dx,y])/(2dx) + 
	//       (dF/dx)(dM/dy) +         -->         (F[x+dx,y]-F[x-dx,y])/(2dx) * (M[x,y+dy]-M[x,y-dy])/(2dy) +
	//       (dF/dy)(dM/dx) +         -->         (F[x,y+dy]-F[x,y-dy])/(2dy) * (M[x+dx,y]-M[x-dx,y])/(2dx) +
	//       (dF/dy)(dM/dy) = ...     -->         (F[x,y+dy]-F[x,y-dy])/(2dy) * (M[x,y+dy]-M[x,y-dy])/(2dy) = 2*(V[x,y]-E)*F[x,y]/h_bar^2
	//
	//The form in (13) can be rearranged so that you can use the shooting method. 
	//
	//Instead of proceeding directly with that we are going to first write the 
	//orginal 2nd order differential equation (3) as a system of coupled ODEs
	//(3)    Div * (M * Div * F) = 2*(V-E)*F/h_bar^2
//
//Let:  
//(14)    G = (M) * Div * F
//Then (3) becomes:
//(15)   Div * G = 2*(V-E)*F/h_bar^2
//
//Rearranging (14) and (15) slightly:
//      +--------------------------------------------------------------+
//(16)  |             Div*F = G / M                                    |
//(17)  |             Div*G = F * 2*(V-E)/h_bar^2                      |
//      +--------------------------------------------------------------+----> These are the useful equations the solver will implement
//
//Expanding (16) with finite differencing:
//       Div*F = G/M
//             = dF/dx + dF/dy = G/M
//(18)         = (F[x+dx,y]-F[x-dx,y])/(2dx) + (F[x,y+dy]-F[x,y-dy])/(2dy) = G[x,y]/M[x,y]
//
//Expanding (17) with finite differencing:
//       Div*G = F * 2*(V-E)/h_bar^2
//             = dG/dx + dG/dy = F * 2*(V-E)/h_bar^2
//(19)         = (G[x+dx,y]-G[x-dx,y])/(2dx) + (G[x,y+dy]-G[x,y-dy])/(2dy) = F[x,y]*2*(V[x,y]-E)/h_bar^2
//
//
/////////////////////
//Initial Conditons//
/////////////////////
//
//Take the structure to start in a barrier high enough enough to warrant using:
//(20)   F[0,0] = 0
//
//As for the initial value of G: G is porportional to the slope of F (the envelope function)
//solutions for the envelope function will be approx. exponentials inside barrier regions.
//Exponentials are fast functions so any reasonable value of G will get the value of F going 
//to where it is going to converge too. i.e. the overall tragectory will not be sensitive to
//the intial value of G if we are starting deep enough inside a tall barrier.
//
//(21)   G[0,0] = 1    
//Now, for the Feature Presentation:
// To find the lowest (quasi)bound state we can solve the infinite square well
// problem to get an upper bound on the energy.
// The solution for a partible in an infinite square well is known to be:
//(13) En = hbar^2*pi^2*n^2/(2*m*L^2)  ... where n=1 is the ground state
// The infinite square well solutions are used as a guess to find boundstates in the 
// finite wells
int ps_solve_1D(PS_DATA potential) {

	printf("PsiShooter -- a shooting method solver for the time independant Schrodinger equation under the effective mass approximation.\n");
	
	//TO DO: these will not stay hard coded. Someday we will have a slick way to handle file I/O.
	char LOG_FILENAME_V[] = "V.dat"; //Output file for the potential energy profile, V(x)
    char LOG_FILENAME_x[] = "x.dat"; //Output file for the indeptendant variable, x
    char LOG_FILENAME_F[] = "F.dat"; //Output file for the envelope function (wavefunction)
    char LOG_FILENAME_G[] = "G.dat"; //Output file for the auxillary function coupled to F
    char LOG_FILENAME_BS[] = "BS.dat"; //Output file for the envelope function (wavefunction) of bound states
    char LOG_FILENAME_E[] = "E.dat"; //Output file for the eigenvalues corresponding to bound states

    
    printf("Iterate through Energy Eigenvalues to find the lowest bound state.\n");      
	
	//TO DO: the energies being searched will not stay hard coded. We will probably define a structure to pass in that 
	//controls the energy range that will be searched for solutions
	double Emin = 0; //meV
	double Emax = 500; //meV
	double Estep = (Emax-Emin)/1000; //meV

	// wavefunction storage.  
	int N = potential->xsize;
	double dx = potential->xstep; //cm, size of differential length
	int N_threshold = 20; //point at which the threshold magnitude is taken. Its a kludgy way to do it, but for now its fine. To Do: Make this more general
	double F[N]; //the envelope function (wavefunction)
	double G[N]; //the aux function
	
	//F
	//intial eqn                 -->  Div*F = G/M	
	//convert to finite diff eqn -->  (F[x+dx,y]-F[x-dx,y])/(2dx) + (F[x,y+dy]-F[x,y-dy])/(2dy) = G[x,y]/M[x,y]	
	//convert to 1D              -->  (F[x+dx]-F[x-dx])/(2dx) = G[x]/M[x]	
	//rearrange                  -->  F[x+dx] = (2dx)*G[x]/M[x] + F[x-dx]
	//convert to c style         -->  F[i+1] = (2dx)*G[i]*m[i] + F[i-1]
	//                                F[i+1] = F_coeff*G[i]*m[i] + F[i-1]
	double F_coeff = 2*dx; //handy prefactor that would otherwise be for every point evaluated
	
	//G
	//intial eqn                 -->  Div*G = F * 2*(V-E)/h_bar^2
	//convert to finite diff eqn -->  (G[x+dx,y]-G[x-dx,y])/(2dx) + (G[x,y+dy]-G[x,y-dy])/(2dy) = F[x,y]*2*(V[x,y]-E[x,y])/h_bar^2
	//convert to 1D              -->  G[x+dx]-G[x-dx]/(2dx) = F[x]*2*(V[x]-E)/h_bar^2	
	//rearrange                  -->  G[x+dx] = (2dx)*F[x]*2*(V[x]-E)/h_bar^2 + G[x-dx]
	//convert to c style         -->  G[i+1] = 4*dx/hbar^2 * F[i]*(V[i]-E) + G[i-1]
	//                                G[i+1] = G_coeff * F[i]*(V[i]-E) + G[i-1]
	double G_coeff = 4*dx/(HBAR_PLANCK*HBAR_PLANCK); //handy prefactor that would otherwise be computed for every point evaluated 

    int bound_state_count = 0;	
    int threshold_set_flag = 0; 
    double F_threshold = 0;
	int i = 0;
	double E; //meV, Current energy
	double V; //meV, Current potential
	double m_eff = MASS_ELECTRON*M_EFF_GAAS;// To Do: make the electron mass be part of the structure that we are simulating. i.e. in general it can be a position dependant quantity just like the potential (think heterostructures with different band edge curvatures)
	int err;
		
	for(E = Emin; E < Emax; E+=Estep) {

        printf("\tE=%g meV\n", E);
        
        F[0] = 0;
        G[0] = 1;
		F[1] = 0;
        G[1] = 1;
        
		for(i=1; i<N-2; i++) {
			V = ps_data_value(potential, 0,i); //V[i]
			if (PS_OK != err) {
				printf("BADNESS.\n");
				goto END;
			}
			
			F[i+1] = F_coeff * G[i] * m_eff + F[i-1]; //subbed in m_eff for m[i], To Do: support a position dependant mass by storing different masses at different locations (add to the PS_DATA structure probably)			
			G[i+1] = G_coeff * F[i]*(V-E) + G[i-1];
        }
        
        printf("\tF[N-1]=%g\n", F[N-1]);
        
        //Is the boundry value near the threshold (and therefore a solution)?
        if(!threshold_set_flag) {
            F_threshold = F[N_threshold]; //Use the value that F gets to half way through the first barrier (the boundry barrier)
            printf("\tSetting boundry value threshold for F_threshold=%e, which is the value of F @ point N=%d (1/2 through the first barrier for the first \"shot\"), \n", F_threshold, N_threshold);
            threshold_set_flag = 1; 
        }   
        
        //printf("\tF[N]=%e\n", F[N/2]);
        //abs(F[N]) < F_threshold
		
        if(fabs(F[N-1]) < F_threshold) {
            //print the energy and the WFN envelope to a file 
			//bound_energies[bound_state_count++] = E;
            printf("\t\tBoundstate number %d with E=%e found, F[N]=%e < F_threshold=%e printing to log file, %s \n", bound_state_count, E, F[N], F_threshold, LOG_FILENAME_BS);
            
            FILE *pFile_E;
            char buffer_E[BUFSIZ];
            pFile_E = fopen(LOG_FILENAME_E, "a"); //append
            if (pFile_E != NULL) {
                setbuf(pFile_E, buffer_E);
                //fputs("#Output of F",pFile_F);
                fprintf(pFile_E, "%e\n", E);                            
                fclose(pFile_E); //fflush(pFile1); //closing flushes already
            }
            
            FILE *pFile_BS;
            char buffer_BS[BUFSIZ];
            pFile_BS = fopen(LOG_FILENAME_BS, "a"); //append
            if (pFile_BS != NULL) {
                setbuf(pFile_BS, buffer_BS);
                //fputs("#Output of F",pFile_F);
                for(i=0; i<N; i++) {
                    fprintf(pFile_BS, "%e ", F[i]);
                }            
                fprintf(pFile_BS, "\n");
                fclose(pFile_BS); //fflush(pFile1); //closing flushes already
            }
        }
        

		//Output printf files
        //print F[i] to log
        //print G[i] to log
		
        printf("\t\tPrinting F to log file, %s \n", LOG_FILENAME_F);
        FILE *pFile_F;
        char buffer_F[BUFSIZ];
        pFile_F = fopen(LOG_FILENAME_F, "a"); //append
        if (pFile_F != NULL) {
            setbuf(pFile_F, buffer_F);
            //fputs("#Output of F",pFile_F);
            for(i=0; i<N; i++) {
                fprintf(pFile_F, "%e\t%e\n", i*dx, F[i]);
            }            
            fclose(pFile_F); //fflush(pFile1); //closing flushes already
        }
		
        printf("\t\tPrinting G to log file, %s \n", LOG_FILENAME_G);
        FILE *pFile_G;
        char buffer_G[BUFSIZ];
        pFile_G = fopen(LOG_FILENAME_G, "a"); //append
        if (pFile_G != NULL) {
            setbuf(pFile_G, buffer_G);
            //fputs("#Output of G",pFile_G);
            for(i=0; i<N; i++) {
                fprintf(pFile_G, "%e\t%e\n", i*dx, G[i]);
            }            
            fclose(pFile_G); //fflush(pFile1); //closing flushes already
        }
    }
    

	END:
    	
	
	return bound_state_count;
}

//
// Generates a test potential. Not really a permenant feature, but rather its 
// something I intend to use in order to do some useful solver engine work 
// while the generation of potentials from files/user input is up in air.
//    
// Lets define a single well as a test case for the BV/shooter method
//      1: V=Vb     2: V=0     3: V=Vb
// 
// The thickness of well (region 2) is "well_width"
// The thickness of the barriers (regions 1 & 3) is "barrier_width"
// 
// V=Vb       _______     _______
//                   |   |    
//                   |   |        
// V=0               |___|    
// Region:      1      2     3  
PS_DATA test_potential_1D() {
	//For convience of trying different scenarios the test potential is defined here in terms of nanometers * cm/nm
	double well_width = 10 * 1e-7; //nm * cm/nm = cm, region 2
	double barrier_width = 10 * 1e-7; //nm * cm/nm , regions 1 and 3
	int number_of_points = 1000;
	 
	int xsize = number_of_points;
	double xstep = (well_width + barrier_width + barrier_width)/((double)number_of_points); //total width converted to cm divided by the number of points = width per point
	int ysize = 1; //1D for now so only 1 "row"
	double ystep = 0;
	double Vb = 500.0; // meV barrier
	
	PS_DATA potential = ps_create_data(xsize, ysize, xstep, ystep);
	
	//temporary local variables
	int i;
	int err;
	double V;
	double x = 0;
	
	//Generate the potential, 
 	for (i = 0; i < xsize; i++) {	
		if(xstep > barrier_width && xstep < (barrier_width+well_width)) 
			err = ps_data_set_value_at_row_column(potential, 0, 0, i); //in the well
        else 
			err = ps_data_set_value_at_row_column(potential, Vb, 0, i);	//in the barriers
		x += xstep;
	}
	
	
	return potential;	
 }


//
// PsiShooter program entry point
int main(int argc, char **argv) {
	//get a 1D test potential
	PS_DATA potential = test_potential_1D();
	
	// Solve
	int nfound = ps_solve_1D(potential);
	printf("Found %i Bound Energies. Have a nice day!\n", nfound);	
	
	//clean up
	ps_destroy_data(potential);
	
	return 0;
}

