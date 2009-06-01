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

PS_LIST g_solutions;

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
//(16)  |             Div*F = G / M                                    |     !!! DONT BE FOOLED,  M = 1/m_eff (JNS)
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
//	double G_coeff = 2*dx/(HBAR_PLANCK_SQ); //handy prefactor that would otherwise be computed for every point evaluated 
	//from 4*dx to 2*dx jh sometime in may

	double G_coeff = 2*dx*G_COEFF; // Redefined to use M_ELECTRON/HBAR^2 using assuming units of nm and eV
    
	int bound_state_count = 0;	
	int threshold_set_flag = 0; 
	double F_threshold = 0;
	int i, iter;
	double E = params->energy_min;
	double Estep = (params->energy_max - params->energy_min)/(params->n_iter); 
	double V; //meV, Current potential
//	double m_eff = MASS_ELECTRON*M_EFF_GAAS;// To Do: make the electron mass be part of the structure that we are simulating. i.e. in general it can be a position dependant quantity just like the potential (think heterostructures with different band edge curvatures)
	double m_eff = M_EFF_GAAS;
	
	for (iter = 0; iter < params->n_iter; iter++) {
		
	  //initial conditions
	  F[0] = 0;
	  G[0] = 1;
	  F[1] = 0;
	  G[1] = 1;
        
	  for(i=1; i<N-1; i++) {
	    V = ps_data_value(potential, 0,i); //V[i]
	    //Mix the bidirectional and forward derivatives to prevent two seperate solutions from forming. This
	    //Makes the upper solutions stable and oscillation free. (Probably only required on the i==1 index -jns) 
		// TODO compare results of this branching statement with IF (1==i) ... ELSE ...
	    if (i%2==1){
	      F[i+1] = F_coeff * G[i] * m_eff + F[i]; 
	      //subbed in m_eff for m[i], To Do: support a position dependant
	      // mass by storing different masses at different locations (add to the PS_DATA structure probably)
	      G[i+1] = G_coeff * F[i] * (V-E) + G[i];
	    }
	    if (i%2==0){
	      F[i+1] = 2*F_coeff * G[i] * m_eff + F[i-1]; 
	      //subbed in m_eff for m[i], To Do: support a position dependant
	      // mass by storing different masses at different locations (add to the PS_DATA structure probably)
	      G[i+1] = 2*G_coeff * F[i] * (V-E) + G[i-1];
	    }
	  }
	  f_cache[iter] = F[N-1];
		
	  // Consider removing this for speed
	  //         sprintf(log_message, "\tE=%g eV\tF[N-1]=%g\n", E/EV_TO_ERGS, F[N-1]);
	  // ps_log(log_message);
		
	  //Try to detect which solutions are eigen states (bound states)
	  if((iter > 0) && ((f_cache[iter] > 0 && f_cache[iter-1] < 0) || (f_cache[iter] < 0 && f_cache[iter-1] > 0))) { 
	    //if there was a change in sign between the last point of the prev
	    // and current envelope function then there was a zero crossing and there a solution
	    // test for a sign crossing

	    PS_DATA wavefunction = ps_data_copy(potential); // copy the potential
	    ps_data_set_data(wavefunction, F); // Overwrite the potential with the wavefunction
			
	    PS_SOLUTION *solution = (PS_SOLUTION*)malloc(sizeof(PS_SOLUTION));
	    solution->energy = E;
	    solution->wavefunction = wavefunction;
			
	    // Add the solution to the list of bound states
	    ps_list_add(solution_list, solution);

	    // print a log message
            sprintf(log_message, "\tBoundstate number %d with E=%e found, F[N]=%e < F_threshold=%e\n", ++bound_state_count, solution->energy, F[N], F_threshold);
	    ps_log(log_message);
	  }

	  // Increment the energy
	  E += Estep;
	}    	
	
	return solution_list;
}

//
//The 2D version of the solving engine.
//
//intial eqn                 
//-->  Div*F = G/M	
//
//convert to finite diff eqn, use a uniform mesh, Let, dx = dy = d  and Let P[x,y] = G[x,y]/M[x,y] = G*m
//-->  (F[x+d,y]-2F[x,y]+F[x-d,y])/d^2 + (F[x,y+d]-2F[x,y]+F[x,y-d])/d^2 = P[x,y]
//
//rearrange:
//-->  F[x+d,y] + F[x-d,y] + F[x,y+d] + F[x,y-d] - 4F[x,y] = d^2 * P[x,y]
//
//convert to c style: Let x = x0+j*d     and    y = y0+l*d 
//                  where j = 0,1,2...J  and    l = 0,1,2...L
//-->  F[j+1,l] + F[j-1,l] + F[j,l+1] + F[j,l-1] - 4F[j,l] = d^2 * P[j,l]y
//
//Goal: make a linear system of eqns in matrix form -- firstly we make a vector out of F. i.e. map j,l to a linear index 'i'
//     Let, i = j(L+1) + l   ...where j = 0,1,...J and l = 0,1,...L
//-->  F[i+L+1] + F[i-(L+1)] + F[i+1] + F[i-1] - 4F[i] = d^2 * P[i]
//    ...where this eqn holds for the interior points: j = 1,2,...,J and l = 1,2,...,L
//
//On the boundries F should be specified, i.e. where:
//    j = 0    or    i = 0, 1, 2, ..,L
//    j = J    or    i = J(L+1), J(L+1)+1, J(L+1)+2, ...,J(L+1)+L
//    l = 0    or    i = 0, (L+1), 2(L+1), ..., J(L+1)
//    l = L    or    i = L, (L+1)+L, 2(L+1)+L, ..., J(L+1)+L
PS_LIST ps_solve_2D(PS_DATA potential, PS_SOLVE_PARAMETERS *params) {
	
	char log_message[256];
	ps_log("PsiShooter -- a shooting method solver for the time independant Schrodinger equation under the effective mass approximation.\n");
	
    // This is a linked list for storing solutions
	PS_LIST solution_list = ps_list_create();

    ps_log("Iterate through Energy Eigenvalues to find the lowest bound state.\n");      
	
	// wavefunction storage.  
	int Nx = potential->xsize; // used for generating strings with sprintf for sending to ps_log
	int Ny = potential->ysize; // used for generating strings with sprintf for sending to ps_log
	
	//	double dx = potential->xstep; //cm, size of differential length 
	//(TODO: Compute this inside loop in order to support meshes with non-uniform differential steps)
	double dx = ps_data_dx_at(potential, 2); 
	double dy = ps_data_dy_at(potential, 2); // 0 to 1 is wierd, maybe. 
	double dx_dy = dx/dy; 
	double dy_dx = dy/dx;
	//F
	//intial eqn                 -->  Div*F = G * m_eff	
	//convert to finite diff eqn -->  (F[x+dx,y]-F[x,y])/dx + (F[x,y]-F[x,y-dy])/dy = G[x,y]M[x,y]	
	//rearrange                  -->  F[x+dx,y]-F[x,y] + dx*(F[x,y]-F[x,y-dy])/dy = dx*G[x,y]/M[x,y]
	//                           -->  F[x+dx,y] = dx*G[x,y]/M[x,y]+F[x,y]-dx/dy*(F[x,y]-F[x,y-dy])
	//                                What does this mean for our 2D simulations? The more closely spaced the y
	//                                Coordinate mesh points are, the more significant their effects will be. That
	//                                Is pretty terrible for us. 
	//
	//                                We need to decouple the coordinates so we don't need to know future F and G from
	//                                from future calculated y mesh points...
	//                                Is there some other way we can decouple them? We can't just look back instead of forward.
	//
	//                                How about we fill the array with 1D shots across the 1 direction, and then refill
	//                                with the 2D shots? Then the array will have some approximate knowledge of the future
	//                                and the algorithm might work.
	
	//convert to c style         -->  F[x+dx,y] = F[i][j+1]
	// *****README*******    	 --> (NOTE: x is iterated over with j index. x=columns, i=rows) -jns

	//in 2D                      -->  F[i][j+1] = dx*G[i][j]/M[i][j] + F[i][j] - dx/dy*(F[i][j]-F[i-1][j])

	
	//convert to 1D              -->  (F[x+dx]-F[x-dx])/(2dx) = G[x]/M[x]	
	//rearrange                  -->  F[x+dx] = (2dx)*G[x]/M[x] + F[x-dx]
	//convert to c style         -->  F[i+1] = (2dx)*G[i]*m[i] + F[i-1]
	//                                F[i+1] = F_coeff*G[i]*m[i] + F[i-1]
	double F_coeff = dx; //handy prefactor that would otherwise be for every point evaluated
	//from dx to 2*dx jh sometime in may
	
	//G
	//intial eqn                 -->  Div*G = F * 2*(V-E)/h_bar^2
	//convert to finite diff eqn -->  (G[x+dx,y]-G[x-dx,y])/(2dx) + (G[x,y+dy]-G[x,y-dy])/(2dy) = F[x,y]*2*(V[x,y]-E[x,y])/h_bar^2
	//in 2d:                     -->  G[x+dx,y] = (2dx)*F[x,y]*2*(V[x,y]-e[x,y])/h_bar^2+G[x-dx][y]-(2dx)/(2dy)(G[x,y+dy]-G[x,y-dy])
	//                           -->  G[i][j+1] = G_coeff * F[i][j]*(V[i][j]-E) + G[i][j-1] - dx/dy(G[i+1][j] + G[i-1][j])
	//
	// The equivalent FWD DIFFEQ -->  G[x+dx][y] = dx*F[x,y]*2*(V[x,y] - e[x,y])/h_bar^2 + G[x,y] - dx/dy(G[x,y] - G[x,y-dy])
	// 						     -->  G[i][j+1] = dx*F[i][j]*2*(V[i][j] - e[i][j])/h_bar^2 + G[i][j] - dx/dy(G[i][j] - G[i-1][j])


	//convert to 1D              -->  G[x+dx]-G[x-dx]/(2dx) = F[x]*2*(V[x]-E)/h_bar^2	
	//rearrange                  -->  G[x+dx] = (2dx)*F[x]*2*(V[x]-E)/h_bar^2 + G[x-dx]
	//convert to c style         -->  G[i+1] = 4*dx/hbar^2 * F[i]*(V[i]-E) + G[i-1]
	//                                G[i+1] = G_coeff * F[i]*(V[i]-E) + G[i-1]
//	double G_coeff = 2*dx/(HBAR_PLANCK_SQ); //handy prefactor that would otherwise be computed for every point evaluated 
	//from 4*dx to 2*dx jh sometime in may

	double G_coeff = 2*dx*G_COEFF;
	
	//dx = dy here.
	double F[Nx][Ny]; //the envelope function (wavefunction)
	double f_cache[params->n_iter]; // We cache the last value of the envelope for each solution
	double G[Nx][Ny]; //the aux function
	int bound_state_count = 0;	
	int threshold_set_flag = 0; 
	double F_threshold = 0;
	int i, j, iter;
	double E = params->energy_min;
	double Estep = (params->energy_max - params->energy_min)/(params->n_iter); 
	double V; //meV, Current potential
//	double m_eff = MASS_ELECTRON*M_EFF_GAAS;
	double m_eff = M_EFF_GAAS;

	// To Do: make the electron mass be part of the structure that we are simulating. i.e. 
	//in general it can be a position dependant quantity just like the potential 
	//(think heterostructures with different band edge curvatures)
		
	for (iter = 0; iter < params->n_iter; iter++) {

	  //Going to try filling the array with y-direction shots and then running the 2D shots based on those.
		
	  //initial conditions
	  for (j=0; j<Nx; j++){
	    F[0][j] = 0;
	    G[0][j] = 1;
		F[1][j] = 1;
		G[1][j] = 1;
	  }
	
	  for(i=1; i<Ny-1; i++) {

		    for(j=1; j<Nx-1; j++) {
		     	V = ps_data_value(potential, i,j); //V[i][j]
				// 5-point stencil in X and Y
				F[i+1][j] = dx*dx*m_eff*G[i][j]  - F[i][j+1] - F[i][j-1] - F[i-1][j] + 4*F[i][j];
				G[i+1][j] = 2*dx*dx*G_COEFF*(V-E)*F[i][j] - G[i][j+1] - G[i][j-1] - G[i-1][j] + 4*G[i][j];
		    }
			// For First and last point, Compute slope in x direction for F and G directly

			// F[x+3] = F[x] + 3h*F'[x] + (3h)^2/2!*F''[x]  (Taylor Expansion)
			// F'[x] = 1/(12*h)*[-F[x+2] + 8*F[x+1] - 8*F[x-1] + F[x-2]]  (5-point stencil)
			// F''[x] = 1/(12*h*h)*[-F[x+2] + 16*F[x+1] -30*F[x] +16*F[x+1] - F[x-2]] (5-point stencil)
			G[i+1][j] = G[i+1][j-3] + 0.25*(-G[i+1][j-1] + 8*G[i+1][j-2] - 8*G[i+1][j-4] + G[i+1][j-5]) + 9.0/24.0*(-G[i+1][j-5] + 16*G[i+1][j-4] - 30*G[i+1][j-3] + 16*G[i+1][j-2] - G[i+1][j-1]);
			F[i+1][j] = F[i+1][j-3] + 0.25*(-F[i+1][j-1] + 8*F[i+1][j-2] - 8*F[i+1][j-4] + F[i+1][j-5]) + 9.0/24.0*(-F[i+1][j-5] + 16*F[i+1][j-4] - 30*F[i+1][j-3] + 16*F[i+1][j-2] - F[i+1][j-1]);			

			// F[x-3] = F[x] - 3h*F'[x] + (3h)^2/2!*F''[x]  (Taylor Expansion)
			// F'[x] = 1/(12*h)*[-F[x+2] + 8*F[x+1] - 8*F[x-1] + F[x-2]]
			// F''[x] = 1/(12*h*h)*[-F[x+2] + 16*F[x+1] -30*F[x] +16*F[x+1] - F[x-2]]
			G[i+1][0] = G[i+1][3] - 0.25*(-G[i+1][5] + 8*G[i+1][4] - 8*G[i+1][2] + G[i+1][1]) + 9.0/24.0*(-G[i+1][5] + 16*G[i+1][4] - 30*G[i+1][3] + 16*G[i+1][2] - G[i+1][1]);
			F[i+1][0] = F[i+1][3] - 0.25*(-F[i+1][5] + 8*F[i+1][4] - 8*F[i+1][2] + F[i+1][1]) + 9.0/24.0*(-F[i+1][5] + 16*F[i+1][4] - 30*F[i+1][3] + 16*F[i+1][2] - F[i+1][1]);			
	  }
	  //Why not just look at the change on the very last point? This definitely won't be representative, but it
	  //might be useful during debugging.
	  // It needs to check the entire boundary. (jns)
	  f_cache[iter] = F[Nx-1][Ny-1];
		
	  // Consider removing this for speed
	  //         sprintf(log_message, "\tE=%g eV\tF[N-1]=%g\n", E/EV_TO_ERGS, F[N-1]);
	  // ps_log(log_message);
		
	  //Try to detect which solutions are eigen states (bound states)
	  //if((iter > 0) && ((f_cache[iter] > 0 && f_cache[iter-1] < 0) || (f_cache[iter] < 0 && f_cache[iter-1] > 0))) { 
	    //if there was a change in sign between the last point of the prev
	    // and current envelope function then there was a zero crossing and there a solution
	    // test for a sign crossing

	    PS_DATA wavefunction = ps_data_copy(potential); // copy the potential
	    ps_data_set_data(wavefunction, F); // Overwrite the potential with the wavefunction
			
	    PS_SOLUTION *solution = (PS_SOLUTION*)malloc(sizeof(PS_SOLUTION));
	    solution->energy = E;
	    solution->wavefunction = wavefunction;
			
	    // Add the solution to the list of bound states
	    ps_list_add(solution_list, solution);

		// Save the g solutions
		PS_DATA g = ps_data_copy(potential);
		ps_data_set_data(g, G);
		PS_SOLUTION *gsolution = (PS_SOLUTION*)malloc(sizeof(PS_SOLUTION));
		gsolution->energy = E;
		gsolution->wavefunction = g;
		ps_list_add(g_solutions, gsolution);
		
	    // print a log message
        sprintf(log_message, "\tBoundstate number %d with E=%e found, F[Nx][Ny]=%e < F_threshold=%e\n", ++bound_state_count, solution->energy/EV_TO_ERGS, F[Nx][Ny], F_threshold);
	    ps_log(log_message);
	    //}

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
	double well_width = 10; // * 1e-7; //nm * cm/nm = cm, region 2
	double barrier_width = 10; // * 1e-7; //nm * cm/nm , regions 1 and 3
	int number_of_points = 1000;
	 
	int xsize = number_of_points;
	double xstep = (well_width + barrier_width + barrier_width)/((double)number_of_points); //total width converted to cm divided by the number of points = width per point
	int ysize = 1; //1D for now so only 1 "row"
	double Vb = 0.5; //*EV_TO_ERGS; // eV barrier
	
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
	double well_width_x = 10; // * 1e-7; //nm * cm/nm = cm, region 2
	double barrier_width_x = 10; // * 1e-7; //nm * cm/nm , regions 1 and 3
	double well_width_y = well_width_x; //square for now.
	double barrier_width_y = barrier_width_x; //square for now.
	double total_width_x = well_width_x + barrier_width_x + barrier_width_x;
	double total_width_y = well_width_y + barrier_width_y + barrier_width_y;

	int number_of_points_x = 100;
	int number_of_points_y = 100;
	
	int xsize = number_of_points_x;
	double xstep = total_width_x/((double)number_of_points_x); //total width converted to cm divided by the number of points = width per point
	int ysize = number_of_points_y;
	double ystep = total_width_y/((double)number_of_points_y); //total width converted to cm divided by the number of points = width per point

	double Vb = 0.5; //*EV_TO_ERGS; // eV barrier
	
	PS_DATA potential = ps_data_create(xsize, ysize);
	potential->xstep = xstep; // This is temporary.  Jere and I have changed the file format to support non-uniform rectilinear grids
	potential->ystep = ystep; // which means that each x and y value is specified and the dx must be queried at every point.
	
	//temporary local variables
	int i,j;
	int err;
	double V;
	double x = -total_width_x/2.0, y = -total_width_y/2.0;
	// double x = 0,y=0;
	//Generate the potential, 
 	for (i = 0; i < ysize; i++) {	
		y = -total_width_y/2.0;
		// y = 0;
		for (j = 0; j < xsize; j++) {	
			
			// if( (x > barrier_width_x && x < (barrier_width_x+well_width_x)) && 
			//     (y > barrier_width_y && y < (barrier_width_y+well_width_y)) ) {				
			// 	err = ps_data_set_value_at_row_column(potential, 0, j, i); //in the well			
			// } else {
			// 	err = ps_data_set_value_at_row_column(potential, Vb, j, i);	//in the barriers			
			// } 

			V = Vb/(total_width_x*total_width_x)*(x*x + y*y);
			ps_data_set_value_at_row_column(potential, V, j, i);
			
			if (0 == i) {
				ps_data_set_y_value_at(potential, j, y);				
			}
			
			y += ystep;
		}
		x += xstep;
		ps_data_set_x_value_at(potential, i, x);
	}	
	return potential;	
}


//
// PsiShooter program entry point
int main(int argc, char **argv) {

	char msg[256];
	PS_DATA potential;
	int solver;
	PS_SOLUTION *s;
	
	if (1 == argc) {

		sprintf(msg, "No file specified. Using builtin potential.\n");
		ps_log(msg);
		//get a 1D test potential
		potential = test_potential_2D();
		
		FILE *f = fopen("V_2d.dat", "w");
		ps_data_write_bin(potential, f);
		fclose(f);
		solver = 2;
	} 
	else if (2 == argc) {
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
			solver = 1;
		}
	} 
	else if (3 == argc) {
		// lets say that if there is a third argument, it is a 2d request.
                //Interpret argument as file to process

		FILE *infile = fopen(argv[1], "r");
		if (infile == NULL) {
			sprintf(msg, "Error opening file '%s'\n", argv[1]);
			ps_log(msg);
			
			return PS_ERROR_FILE_NOT_FOUND;  // Exit abnormally
			
		} else {
			sprintf(msg, "Using *2D* potential from file '%s'\n", argv[1]);
			ps_log(msg);

			potential = ps_data_read_bin(infile);
			fclose(infile);
			solver = 2;
		}
        }


	// Setup solution parameters
	PS_SOLVE_PARAMETERS params;
	params.energy_min = ps_data_min_value(potential);
	params.energy_max = ps_data_max_value(potential);
	params.n_iter = 50; // The number of energies to try
	double e_step = (params.energy_max - params.energy_min)/(params.n_iter);
		
	sprintf(msg,"Testing energies from %g to %g in %g increments\n",params.energy_min/EV_TO_ERGS,params.energy_max/EV_TO_ERGS,e_step/EV_TO_ERGS);
	ps_log(msg);
	
	// Create a list for solutions
	PS_LIST solutions = ps_list_create();
	g_solutions = ps_list_create();
	
	  // 1st pass coarse Solver
	if (1 == solver) {
		// 1d solver with coarse/fine pass
		PS_LIST coarse_solutions = ps_solve_1D(potential, &params);	
		 // For each solution, do a second pass with finer grained energy steps
		 s = ps_list_front(coarse_solutions);
		 while (s != NULL) {
		   // I know that ps_solve_1d looks for a sign crossing as it increases
		   // the test energy E and saves the solution for the larger energy
		   params.energy_min = s->energy - e_step;
		   params.energy_max = s->energy;
		   params.n_iter = 100;

		   // Add the fine grain solutions to the final list
		   PS_LIST e_solutions = ps_solve_1D(potential, &params);
		   ps_list_add_all(solutions, e_solutions);
		   ps_list_destroy(e_solutions); // Destroy the list without destroying the data.

		   // Next loop iteration
		   s = ps_list_next(coarse_solutions);
		 }
		ps_list_destroy_all(coarse_solutions);

	} else {
		// 2d solver
		PS_LIST coarse_solutions = ps_solve_2D(potential, &params);
		ps_list_add_all(solutions, coarse_solutions);
		ps_list_destroy(coarse_solutions);
	}

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

	FILE *gfile = fopen("G.dat", "w");
	s = ps_list_front(g_solutions);
	while (s != NULL) {
		ps_data_write_bin(s->wavefunction, gfile);
		s = ps_list_next(g_solutions);
	}
	fclose(gfile);
	ps_list_destroy_all(g_solutions);
	
	//clean up
	ps_data_destroy(potential);
	ps_list_destroy_all(solutions);

	return PS_OK;
}
