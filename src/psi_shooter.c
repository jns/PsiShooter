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
    

	////////////////////////
    //Solving the 1D TISE://
    ////////////////////////
    //The hamiltonian is: H = T + V = p^2/(2m) + V  ... where p is the momentum operator
    //moving to the effective mass regime we no longer solve for the wavefunction
    //but instead we get an envelope function. Careful consideration also shows that 
    //if m -> m(x) then simple subsitution causes H to become non-Hermitian. 
    //Fortunately that can be worked around, the result is (apparently):
    //(1) T * F = -h_bar^2/2 * Del * ( m_eff(x,E)^-1 * Del * F )
	//     note: F is the envelope function. ... taking the role of Psi
    // 
    //Put it all together to get the effective mass equation that we are going to
    //use:
    //
    //(2) -h_bar^2/2 * Del * (m_eff * Del * F) + V * F = E * F   
    //   ... take the 1D direction to be x, then Del = d/dx
    //   ... finite diff. version of deriative: dF/dx = (F[x+dx]-F[x-dx]) / (2dx)    
    //(3) Del * (m_eff * (F[x+dx]-F[x-dx])/(2dx)) = 2*(V-E)*F / h_bar^2
    //(4) F[x+2dx]/m_eff[x+dx] = (2*(2dx/h_bar)^2*(V-E) + 1/m_eff[x+dx] + 1/m,_eff[x-dx])*F - F[x-2dx]/m_eff[x-dx]
    // 
    // The form in (4) is suitable for using in a shooting method because it has
    // F[x+1] as a function of of F[x] and F[x-1]. In order to simplify the math
    // we can write the second order equation (2) as 2 coupled first order 
    // equations instead. With the help of an intermediate function, G.
    // Subsitute (5) into (2) to get the other coupled equation, F:    
    // (5) Del*F = m_eff*G   
    // (6) Del*G = 2/h_bar^2*(V-E)*F
    // 
    // The two coupled equations (5) and (6) can then can have their differentials
    // turned into finite differences to yield (using the forward difference only):
    // (7) (F[x+dx]-F[x]) / dx = m_eff[x]*G[x]
    // (8) (G[x+dx]-G[x]) / dx = 2/h_bar^2*(V-E)*F[x]    
    // ...Rearranging:
    // (9) F[x+dx] = F[x] + dx*m_eff[x]*G[x]
    //(10) G[x+dx] = G[x] + 2dx/h_bar^2*(V-E)*F[x]
    //
    /////////////////////
    //Initial Conditons//
    /////////////////////
    // Take the structure to start in a barrier, deep enough to warrant using:
    //(11) F[0] = 0
    // As for the initial value of G, solutions for F will be exponentials inside 
    // the barrier regions and since exponentials are pretty fast functions any
    // value of G will get the value of F going to where it is going to converge too
    // reasonably well (well, any non zero G; G=0 gives us a trivial case).
    //(12) G[0] = 1
    
    //Now, for the Feature Presentation:
    // To find the lowest (quasi)bound state we can solve the infinite square well
    // problem to get an upper bound on the energy.
    // The solution for a partible in an infinite square well is known to be:
    //(13) En = hbar^2*pi^2*n^2/(2*m*L^2)  ... where n=1 is the ground state
    // The infinite square well solutions are used as a guess to find boundstates in the 
    // finite wells
    
    printf("Iterate through Energy Eigenvalues to find the lowest bound state.\n");      
    printf("\tGround state for the inf square well of the same width, E1=%g eV\n", E1_infsqwell);
    double E = E1_infsqwell;   
    int bound_state_count = 0;
	
    int threshold_set_flag = 0; 
    double F_threshold = 0;
	
    while(E > 0) {
        printf("\tE=%g eV\n", E);
        
        F[0] = 0;
        G[0] = 1;
        double F_coeff = dx*m_eff; //handy prefactor that would otherwise be computed N times in the following loop for F[i+1]
        double G_coeff = 2*dx/(hbar*hbar); //handy prefactor that would otherwise be computed N times in the following loop for G[i+1]       
        for(i=0; i<N-1; i++) {
            F[i+1] = F[i] + F_coeff*G[i]; // (9) F[x+dx] = F[x] + dx*m_eff[x]*G[x]
            G[i+1] = G[i] + G_coeff*(V[i]-E)*F[i]; //(10) G[x+dx] = G[x] + 2dx/h_bar^2*(V-E)*F[x]                
        }
        
        printf("\tF[N-1]=%e\n", F[N-1]);
        
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
            bound_state_count++;
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
        
        //progress to the next eigenstate (enegry) for the next time through the shooter loop
        E = E - dE;
        
		//Output log files
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
                fprintf(pFile_F, "%e\t%e\n", x[i], F[i]);
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
                fprintf(pFile_G, "%e\t%e\n", x[i], G[i]);
            }            
            fclose(pFile_G); //fflush(pFile1); //closing flushes already
        }
		/* */
    }
    
	
    
    
    //double Vb = 0.5; //eV, height of barriers
    //double well_width = 9e-9; //m, width of the well (region 3)
    //double barrier_width = 3e-9; //m, width of the barriers (regions 2 & 4)
    //const double xmin = -15e-9; //nm, lower bound of the computational domain in position space
    //const double xmax = 15e-9; //nm, upper bound of the computational domain in position space
    //const int N = 1e4; //dimensionless, a count of the number of discrete points to use. Larger number results in more computation but more accuracy 
    //double domain_x_length = xmax - xmin; //nm, length of the computational domain in position space
    //domain_x_length = fabs(domain_x_length);
    //double dx = domain_x_length/N; //nm, length of the discrete steps in x, "delta x"  
    //double x[N];
    
	/**/       
	
    
    //cleanup
    free(x);
    free(V);
    free(F);
    free(G);
	
	
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