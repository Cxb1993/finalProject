#include "boundary_val.h"
#include "helper.h"

/* ----------------------------------------------------------------------- */
/*                             Set Boundary Conditions                     */
/* ----------------------------------------------------------------------- */

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
                    int imax,
                    int jmax,
                    double **U,
                    double **V,
                    const int wl,
                    const int wr,
                    const int wt,
                    const int wb,
                    int **Flag
                    ) {
    
	int i,j;
	/*Initialize corners*/
	U[0][0]=0.0;
	U[0][jmax+1]=0.0;
	U[imax+1][0]=0.0;
	U[imax+1][jmax+1]=0.0;
	V[0][0]=0.0;
	V[0][jmax+1]=0.0;
	V[imax+1][0]=0.0;
	V[imax+1][jmax+1]=0.0;
    
    
    
	/* Set values for all the outside boundary depending on the value that
	 * the data file is set to. NO_SLIP = 1, FREE_SLIP=2 and OUTFLOW=3
	 * We start with the left boundary. */
	switch(wl){
        case NO_SLIP:
            for (j = 1; j < jmax + 1; j++){
                /*U velocity on left boundary */
                U[0][j] = 0;
                /*V velocity left boundary */
                V[0][j]=-1*V[1][j];
            }
            break;
        case FREE_SLIP:
            for (j = 1; j < jmax + 1; j++){
                /*U velocity on left boundary */
                U[0][j] = 0;
                /*V velocity left boundary */
                V[0][j] = V[1][j];
            }
            break;
        case OUTFLOW:
            for (j = 1; j < jmax + 1; j++){
                /*U velocity on left boundary */
                U[0][j] = U[1][j];
                /*V velocity left boundary */
                V[0][j]= V[1][j];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the right boundary*/
	switch(wr){
        case NO_SLIP:
            for (j = 1; j < jmax + 1; j++){
                U[imax][j] = 0;
                V[imax+1][j]=-1*V[imax][j];
            }
            break;
        case FREE_SLIP:
            for (j = 1; j < jmax + 1; j++){
                U[imax][j] = 0;
                V[imax+1][j] = V[imax][j];
            }
            break;
        case OUTFLOW:
            for (j = 1; j < jmax + 1; j++){
                U[imax][j] = U[imax-1][j];
                V[imax+1][j]= V[imax][j];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the top boundary*/
	switch(wt){
        case NO_SLIP:
            for (i = 1; i < imax + 1; i++){
                V[i][jmax] = 0;
                U[i][jmax+1]= -1*U[i][jmax];
            }
            break;
        case FREE_SLIP:
            for (i = 1; i < imax + 1; i++){
                V[i][jmax] = 0;
                U[i][jmax+1] = U[i][jmax];
            }
            break;
        case OUTFLOW:
            for (i = 1; i < imax + 1; i++){
                V[i][jmax] = V[i][jmax-1];
                U[i][jmax+1] = U[i][jmax];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the bottom boundary*/
	switch(wb){
        case NO_SLIP:
            for (i = 1; i < imax + 1; i++){
                V[i][0] = 0;
                U[i][0]=-1*U[i][1];
            }
            break;
        case FREE_SLIP:
            for (i = 1; i < imax + 1; i++){
                V[i][0] = 0;
                U[i][0] = U[i][1];
            }
            break;
        case OUTFLOW:
            for (i = 1; i < imax + 1; i++){
                V[i][0] = V[i][1];
                U[i][0]= U[i][1];
            }
            break;
        default:
            break;
	}
    
	/**
	 * Loop to check for boundary cells in the inner domain, and assign
	 * correct values of U and V.
	 */
	for(i = 1; i < imax+1; i++){
		for(j = 1; j < jmax+1; j++){
			if((Flag[i][j]&31)==B_N){
				V[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				U[i][j]=-1*U[i][j+1];
			}
			else if((Flag[i][j]&31)==B_S){
				V[i][j-1]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				U[i][j]=-1*U[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				U[i-1][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
				V[i][j]=-1*V[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				U[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
			}
			else if((Flag[i][j]&31)==B_NW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
			}
			else if((Flag[i][j]&31)==B_SO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				V[i][j-1]=0;
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_SW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j-1];
				V[i][j]=-1*V[i-1][j];
				V[i][j-1]=0;
			}
		}
	}
    
}

void spec_boundary_val(
                       char *problem,
                       int imax,
                       int jmax,
                       double **U,
                       double **V
                       ){
	int i,j;
	/**
	 * Special boundary condition for the cavity problem.
	 */
	if(strcmp(problem,"cavity")==0){
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1]= 2.0-1*U[i][jmax];
		}
	}
	/*
	 * Special boundary condition for the Karman Vortex street
	 * problem
	 */
	if(strcmp(problem,"KarmanVortexStreet")==0){
		for (j = 1; j < jmax + 1; j++){
			U[0][j] = 1.0;
			V[0][j]= 0.0;
		}
	}
	/*
	 * Special boundary condition for the Flow Over Step problem
	 */
	if(strcmp(problem,"FlowOverStep")==0){
		for (j = 1; j <= (jmax)/2; j++){
			U[0][j] = 0.0;
			V[0][j]= 0.0;
		}
		for (j = ((jmax)/2) + 1; j < jmax + 1; j++){
			U[0][j] = 1.0;
			V[0][j]= 0.0;
		}
	}
}

/* ----------------------------------------------------------------------- */
/*                     Set Boundary Conditions for Temperature             */
/* ----------------------------------------------------------------------- */


void TEMP_BoundaryCondition(
                            int imax,
                            int jmax,
                            double **TEMP,
                            const int wlt,
                            const int wrt,
                            const int wtt,
                            const int wbt,
                            const double TL,
                            const double TR,
                            const double TB,
                            const double TT,
                            int **Flag
                            ) {
    /* Set values for all the outside boundary depending on the value that
	 * the data file is set to. Temperature = 1 and Adiabatic = 2
     * We start with the left boundary. */
	int i,j;
    
    switch(wlt){
        case SET_TEMP:
            for (j = 0; j <= jmax + 1; j++){
                /*Temperature on left boundary */
                TEMP[0][j]=2.0*TL-1*TEMP[1][j];
            }
            break;
        case ADIABATIC:
            for (j = 0; j <= jmax + 1; j++){
                /*Temperature on left boundary */
                TEMP[0][j] = TEMP[1][j];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the right boundary*/
    
	switch(wrt){
        case SET_TEMP:
            for (j = 0; j <= jmax + 1; j++){
                TEMP[imax+1][j]=2.0*TR-1*TEMP[imax][j];
            }
            break;
        case ADIABATIC:
            for (j = 0; j <= jmax + 1; j++){
                TEMP[imax+1][j] = TEMP[imax][j];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the top boundary*/
	switch(wtt){
        case SET_TEMP:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][jmax+1]= 2.0*TT-1*TEMP[i][jmax];
            }
            break;
        case ADIABATIC:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][jmax+1] = TEMP[i][jmax];
            }
            break;
        default:
            break;
	}
    
	/*Set values for the bottom boundary*/
	switch(wbt){
        case SET_TEMP:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][0]=2.0*TB-1*TEMP[i][1];
            }
            break;
        case ADIABATIC:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][0] = TEMP[i][1];
            }
            break;
        default:
            break;
	}
    
    /*treatment in case we have boundary inside*/
}

void TEMP_SpecBoundaryCondition(
                                char *problem,
                                int imax,
                                int jmax,
                                double **TEMP
                                ) {
    
}

