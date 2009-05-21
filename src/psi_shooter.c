#include <stdio.h>
#include <stdlib.h>
#include "ps_constants.h"
#include "ps_errors.h"
#include "ps_data.h"
#include "ps_data_io.h"
#include "ps_list.h"
#include "math.h"

/*The Big To Do list:
 * Consider the supporting the simplifying case where m_eff(x,y) -> m_eff   (constant effective mass)
 *     Why? ... one effective mass for the entire system, this should cut out a bunch of derivatives that need to be calculated
 *     For many potentials electrons don't "see" much of the barriers compared to the wells, so m_eff(x,y) ~= m_eff_well_material 
 *     Maybe only support this case, it has a bunch less terms when you expand out the TISE with fintite differencing 
 *     note: A position dependant mass fucks up the hermicity of H (hamiltonian operator) 
 */

void ps_log(char *msg) {
	fprintf(stdout, msg);
}

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
PS_LIST ps_solve_1D(PS_DATA potential, PS_SOLVE_PARAMETERS *params) {

	char log_message[256];
	ps_log("PsiShooter -- a shooting method solver for the time independant Schrodinger equation under the effective mass approximation.\n");
	
    // This is a linked list for storing solutions
	PS_LIST solution_list = ps_list_create();

    ps_log("Iterate through Energy Eigenvalues to find the lowest bound state.\n");      
	

	// wavefunction storage.  
	int N = potential->xsize; // used for generating strings with sprintf for sending to ps_log
//	double dx = potential->xstep; //cm, size of differential length (TODO: Compute this inside loop)
	double dx = ps_data_dx_at(potential, 2); // 0 to 1 is wierd
	// int N_threshold = 100; //point at which the threshold magnitude is taken. Its a kludgy way to do it, but for now its fine. To Do: Make this more general
	double F[N]; //the envelope function (wavefunction)
	double f_cache[params->n_iter]; // We cache the last value of the envelope for each solution
	double G[N]; //the aux function
	
	//F
	//intial eqn                 -->  Div*F = G/M	
	//convert to finite diff eqn -->  (F[x+dx,y]-F[x-dx,y])/(2dx) + (F[x,y+dy]-F[x,y-dy])/(2dy) = G[x,y]/M[x,y]	
	//convert to 1D              -->  (F[x+dx]-F[x-dx])/(2dx) = G[x]/M[x]	
	//rearrange                  -->  F[x+dx] = (2dx)*G[x]/M[x] + F[x-dx]
	//convert to c style         -->  F[i+1] = (2dx)*G[i]*m[i] + F[i-1]
	//                                F[i+1] = F_coeff*G[i]*m[i] + F[i-1]
	double F_coeff = dx; //handy prefactor that would otherwise be for every point evaluated
	//from dx to 2*dx jh sometime in may
	
	
	//G
	//intial eqn                 -->  Div*G = F * 2*(V-E)/h_bar^2
	//convert to finite diff eqn -->  (G[x+dx,y]-G[x-dx,y])/(2dx) + (G[x,y+dy]-G[x,y-dy])/(2dy) = F[x,y]*2*(V[x,y]-E[x,y])/h_bar^2
	//convert to 1D              -->  G[x+dx]-G[x-dx]/(2dx) = F[x]*2*(V[x]-E)/h_bar^2	
	//rearrange                  -->  G[x+dx] = (2dx)*F[x]*2*(V[x]-E)/h_bar^2 + G[x-dx]
	//convert to c style         -->  G[i+1] = 4*dx/hbar^2 * F[i]*(V[i]-E) + G[i-1]
	//                                G[i+1] = G_coeff * F[i]*(V[i]-E) + G[i-1]
	double G_coeff = 2*dx/(HBAR_PLANCK_SQ); //handy prefactor that would otherwise be computed for every point evaluated 
	//from 4*dx to 2*dx jh sometime in may
    
	int bound_state_count = 0;	
    int threshold_set_flag = 0; 
    double F_threshold = 0;
	int i, iter;
	double E = params->energy_min;
	double Estep = (params->energy_max - params->energy_min)/(params->n_iter); 
	double V; //meV, Current potential
	double m_eff = MASS_ELECTRON*M_EFF_GAAS;// To Do: make the electron mass be part of the structure that we are simulating. i.e. in general it can be a position dependant quantity just like the potential (think heterostructures with different band edge curvatures)
		
	for (iter = 0; iter < params->n_iter; iter++) {
		
		//initial conditions
        F[0] = 0;
        G[0] = 1;
		F[1] = 0;
        G[1] = 1;
        
		for(i=1; i<N-1; i++) {
			V = ps_data_value(potential, 0,i); //V[i]
			F[i+1] = F_coeff * G[i] * m_eff + F[i]; //subbed in m_eff for m[i], To Do: support a position dependant mass by storing different masses at different locations (add to the PS_DATA structure probably)			
			G[i+1] = G_coeff * F[i] * (V-E) + G[i];
        }
		f_cache[iter] = F[N-1];
		
		// Consider removing this for speed
		//         sprintf(log_message, "\tE=%g eV\tF[N-1]=%g\n", E/EV_TO_ERGS, F[N-1]);
		// ps_log(log_message);
		
		//Try to detect which solutions are eigen states (bound states)
        if((iter > 0) && ((f_cache[iter] > 0 && f_cache[iter-1] < 0) || (f_cache[iter] < 0 && f_cache[iter-1] > 0))) { //if there was a change in sign between the last point of the prev and current envelope function then there was a zero crossing and there a solution
		// test for a sign crossing

			PS_DATA wavefunction = ps_data_copy(potential); // copy the potential
			ps_data_set_data(wavefunction, F); // Overwrite the potential with the wavefunction
			
			PS_SOLUTION *solution = (PS_SOLUTION*)malloc(sizeof(PS_SOLUTION));
			solution->energy = E;
			solution->wavefunction = wavefunction;
			
			// Add the solution to the list of bound states
			ps_list_add(solution_list, solution);

			// print a log message
            sprintf(log_message, "\tBoundstate number %d with E=%e found, F[N]=%e < F_threshold=%e\n", ++bound_state_count, solution->energy/EV_TO_ERGS, F[N], F_threshold);
			ps_log(log_message);
        }

		// Increment the energy
		E += Estep;
    }    	
	
	return solution_list;
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
	double Vb = 0.5*EV_TO_ERGS; // eV barrier
	
	PS_DATA potential = ps_data_create(xsize, ysize);
	potential->xstep = xstep; // This is temporary.  Jere and I have changed the file format to support non-uniform rectilinear grids
		  				  // which means that each x and y value is specified and the dx must be queried at every point.
	
	//temporary local variables
	int i;
	double x;
	
	//Generate the potential, 
 	for (i = 0; i < xsize; i++) {	
		if(x > barrier_width && x < (barrier_width+well_width)) {
			ps_data_set_value_at_row_column(potential, 0, 0, i); //in the well			
		} else {
			ps_data_set_value_at_row_column(potential, Vb, 0, i);	//in the barriers			
		} 
		ps_data_set_x_value_at(potential, i, x);
		x += xstep;
	}
	
	return potential;	
 }

//
// Generates a test potential. Not really a permenant feature, but rather its 
// something I intend to use in order to do some useful solver engine work 
// while the generation of potentials from files/user input is up in air.
// The size of the square well (region 0) is "well_width_x * well_width_y"
// The thickness of the barriers (regions 1 & 3) is "barrier_width"
//
//            well_width_x
//                +---+
//           
//           1 1 1 1 1 1 1 1 1
//           1 1 1 1 1 1 1 1 1  
//           1 1 1 1 1 1 1 1 1  
//           1 1 1 0 0 0 1 1 1  +
//           1 1 1 0 0 0 1 1 1  | well_width_y
//           1 1 1 0 0 0 1 1 1  +
//           1 1 1 1 1 1 1 1 1  
//           1 1 1 1 1 1 1 1 1  
//           1 1 1 1 1 1 1 1 1 
//
//           +---+       +---+
//             \          /
//             barrier_width  
//              
// Region 1: V=Vb
// Regoin 0: V=0
PS_DATA test_potential_2D() {
	//For convience of trying different scenarios the test potential is defined here in terms of nanometers * cm/nm
	double well_width_x = 10 * 1e-7; //nm * cm/nm = cm, region 2
	double barrier_width_x = 10 * 1e-7; //nm * cm/nm , regions 1 and 3
	double well_width_y = well_width_x; //square for now.
	double barrier_width_y = barrier_width_x; //square for now.

	int number_of_points_x = 50;
	int number_of_points_y = 50;
	
	int xsize = number_of_points_x;
	double xstep = (well_width_x + barrier_width_x + barrier_width_x)/((double)number_of_points_x); //total width converted to cm divided by the number of points = width per point
	int ysize = number_of_points_y;
	double ystep = (well_width_y + barrier_width_y + barrier_width_y)/((double)number_of_points_y); //total width converted to cm divided by the number of points = width per point

	double Vb = 0.5*EV_TO_ERGS; // eV barrier
	
	PS_DATA potential = ps_data_create(xsize, ysize);
	potential->xstep = xstep; // This is temporary.  Jere and I have changed the file format to support non-uniform rectilinear grids
	potential->ystep = ystep; // which means that each x and y value is specified and the dx must be queried at every point.
	
	//temporary local variables
	int i,j;
	int err;
	double V;
	double x = 0, y = 0;
	
	//Generate the potential, 
 	for (i = 0; i < xsize; i++) {	
		y = 0;
		for (j = 0; j < ysize; j++) {	

			if( (x > barrier_width_x && x < (barrier_width_x+well_width_x)) && 
			    (y > barrier_width_y && y < (barrier_width_y+well_width_y)) ) {				
				err = ps_data_set_value_at_row_column(potential, 0, j, i); //in the well			
			} else {
				err = ps_data_set_value_at_row_column(potential, Vb, j, i);	//in the barriers			
			} 
			
			if (0 == i) {
				ps_data_set_y_value_at(potential, j, y);				
			}
			
			y += ystep;
		}
		x += xstep;
		printf("\n");
		ps_data_set_x_value_at(potential, i, x);
	}	
	return potential;	
}


//
// PsiShooter program entry point
int main(int argc, char **argv) {

	char msg[256];
	PS_DATA potential;
	
	if (1 == argc) {

		sprintf(msg, "No file specified. Using builtin potential.\n");
		ps_log(msg);
		//get a 1D test potential
		potential = test_potential_2D();

		FILE *f = fopen("V_2d.dat", "w");
		ps_data_write_bin(potential, f);
		fclose(f);
		return 0;
		
	} else if (2 == argc) {
		// Interpret argument as file to process

		FILE *infile = fopen(argv[1], "r");
		if (infile == NULL) {
			sprintf(msg, "Error opening file '%s'\n", argv[1]);
			ps_log(msg);
			
			return PS_ERROR_FILE_NOT_FOUND;  // Exit abnormally
			
		} else {
			sprintf(msg, "Using potential from file '%s'\n", argv[1]);
			ps_log(msg);

			potential = ps_data_read_bin(infile);
			fclose(infile);
		}
	}

	// Setup solution parameters
	PS_SOLVE_PARAMETERS params;
	params.energy_min = ps_data_min_value(potential);
	params.energy_max = ps_data_max_value(potential);
	params.n_iter = 10000; // The number of energies to try
	double e_step = (params.energy_max - params.energy_min)/(params.n_iter);
		
	sprintf(msg, "Testing energies from %g to %g in %g increments\n", params.energy_min/EV_TO_ERGS, params.energy_max/EV_TO_ERGS, e_step/EV_TO_ERGS);
	ps_log(msg);
	
	// Create a list for solutions
	PS_LIST solutions = ps_list_create();
	
	// 1st pass coarse Solver
	PS_LIST coarse_solutions = ps_solve_1D(potential, &params);

	// For each solution, do a second pass with finer grained energy steps
	PS_SOLUTION *s = ps_list_front(coarse_solutions);
	while (s != NULL) {
		// I know that ps_solve_1d looks for a sign crossing as it increases
		// the test energy E and saves the solution for the larger energy
		params.energy_min = s->energy - e_step;
		params.energy_max = s->energy;
		params.n_iter = 10000;
		
		// Add the fine grain solutions to the final list
		PS_LIST e_solutions = ps_solve_1D(potential, &params);
		ps_list_add_all(solutions, e_solutions);
		ps_list_destroy(e_solutions); // Destroy the list without destroying the data.
		
		// Next loop iteration
		s = ps_list_next(coarse_solutions);
	}
	ps_list_destroy_all(coarse_solutions);
		
	// Write out solutions
	sprintf(msg, "Found %i solutions.  Writing to E.txt and BS.dat\n", ps_list_size(solutions));
	ps_log(msg);
	
	FILE *efile = fopen("E.txt", "w");
	fprintf(efile, "# Bound State Energies [eV]\n");
	FILE *bsfile = fopen("BS.dat", "w");
	s = ps_list_front(solutions);
	while (s != NULL) {
		// Write the energies to the energy file
		fprintf(efile, "%g\n", s->energy/EV_TO_ERGS);
		// Write the wavefunctions to the solutino file
		ps_data_write_bin(s->wavefunction, bsfile);
		s = ps_list_next(solutions);
	}
	fclose(bsfile);
	fclose(efile);

	//clean up
	ps_data_destroy(potential);
	ps_list_destroy_all(solutions);
		
	return PS_OK;
}

