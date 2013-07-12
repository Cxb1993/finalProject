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
	if(strcmp(problem,"2FinsT")==0){
		for (j = 1; j < jmax + 1; j++){
			U[0][j] = 1.0;
			V[0][j]= 0.0;
		}
	}
}

/* ----------------------------------------------------------------------- */
/*                     Set Boundary Conditions for Temperature for fluid   */
/* ----------------------------------------------------------------------- */


void TEMP_BoundaryCondition(
                            double dx,
                            double dy,
                            int imax,
                            int jmax,
                            double **TEMP,
                            double **TEMP_S,
                            const int wlt,
                            const int wrt,
                            const int wtt,
                            const int wbt,
                            const double TL,
                            const double TR,
                            const double TB,
                            const double TT,
                            const double QL,
                            const double QR,
                            const double QB,
                            const double QT,
                            double kratio,
                            int **Flag
                            ) {
    /* Set values for all the outside boundary depending on the value that
	 * the data file is set to. Dirichlet = 1 and Neumann = 2
     * We start with the left boundary. */
	int i,j;
    
	/*Set values for the left boundary*/

    switch(wlt){
        case Dirichlet:
            for (j = 0; j <= jmax + 1; j++){
                TEMP[0][j]=2.0*TL-1*TEMP[1][j];
            }
            break;
        case Neumann:
            for (j = 0; j <= jmax + 1; j++){
                if (( Flag[1][j]&B_C ) == B_C ){
                TEMP[0][j] = TEMP[1][j]+QL*dx;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the right boundary*/
    
	switch(wrt){
        case Dirichlet:
            for (j = 0; j <= jmax + 1; j++){
                TEMP[imax+1][j]=2.0*TR-1*TEMP[imax][j];
            }
            break;
        case Neumann:
            for (j = 0; j <= jmax + 1; j++){
                if (( Flag[imax][j]&B_C ) == B_C ){
                    TEMP[imax+1][j] = TEMP[imax][j]-QR*dx;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the top boundary*/
	switch(wtt){
        case Dirichlet:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][jmax+1]= 2.0*TT-1*TEMP[i][jmax];
            }
            break;
        case Neumann:
            for (i = 0; i <= imax + 1; i++){
                if (( Flag[i][jmax]&B_C ) == B_C ){
                    TEMP[i][jmax+1] = TEMP[i][jmax]-QT*dy;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the bottom boundary*/
	switch(wbt){
        case Dirichlet:
            for (i = 0; i <= imax + 1; i++){
                TEMP[i][0]=2.0*TB-1*TEMP[i][1];
            }
            break;
        case Neumann:
            for (i = 0; i <= imax + 1; i++){
                if (( Flag[i][1]&B_C ) == B_C ){
                    TEMP[i][0] = TEMP[i][1]+QB*dy;
                }
            }
            break;
        default:
            break;
	}
    
    /*treatment in case we have boundary inside*/
    /**
	 * Loop to check for boundary cells in the inner domain, and assign
	 * correct values of T.
	 */
	for(i = 1; i < imax+1; i++){
		for(j = 1; j < jmax+1; j++){
            /*
			 * If it's not a fluid cell but it has fluid cell neighbors, some values must still
			 * be calculated using the boundary flags.
			 */
			if((Flag[i][j]&31)==B_N){
				TEMP[i][j] = ( (kratio-1)*TEMP[i][j+1] + 2*TEMP_S[i][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==B_S){
				TEMP[i][j]=( (kratio-1)*TEMP[i][j-1] + 2*TEMP_S[i][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==B_W){
				TEMP[i][j]=( (kratio-1)*TEMP[i-1][j] + 2*TEMP_S[i][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==B_O){
				TEMP[i][j]=( (kratio-1)*TEMP[i+1][j] + 2*TEMP_S[i][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==B_NO){
				TEMP[i][j]= ( ( (kratio-1)*TEMP[i][j+1] + 2*TEMP_S[i][j] ) / (1+kratio) + ( (kratio-1)*TEMP[i+1][j] + 2*TEMP_S[i][j] ) / (1+kratio) ) / 2;
			}
			else if((Flag[i][j]&31)==B_NW){
				TEMP[i][j]= ( ( (kratio-1)*TEMP[i][j+1] + 2*TEMP_S[i][j] ) / (1+kratio) + ( (kratio-1)*TEMP[i-1][j] + 2*TEMP_S[i][j] ) / (1+kratio) ) / 2;
			}
			else if((Flag[i][j]&31)==B_SO){
				TEMP[i][j]= ( ( (kratio-1)*TEMP[i][j-1] + 2*TEMP_S[i][j] ) / (1+kratio) + ( (kratio-1)*TEMP[i+1][j] + 2*TEMP_S[i][j] ) / (1+kratio) ) / 2;
			}
			else if((Flag[i][j]&31)==B_SW){
				TEMP[i][j]= ( ( (kratio-1)*TEMP[i][j-1] + 2*TEMP_S[i][j] ) / (1+kratio) + ( (kratio-1)*TEMP[i-1][j] + 2*TEMP_S[i][j] ) / (1+kratio) ) / 2;;
			}
        }
	}
    
}

void TEMP_SpecBoundaryCondition(
                                char *problem,
                                int imax,
                                int jmax,
                                double **TEMP,
                                int **Flag
                                ) {
    int i, first, last;
    first = 0;
    last = 0;
    for (i = 1; i < imax + 1; i++){
        if((Flag[i][1]&B_C)!=B_C){
            if(first==0){
                first=i;
            }
            last = i;
        }
    }
    
    for(i = 1; i < first; i++){
        TEMP[i][0] = TEMP[i][1];
    }
    for(i = last+1; i < imax+1; i++){
        TEMP[i][0] = TEMP[i][1];
    }

    
}

/* ----------------------------------------------------------------------- */
/*                     Set Boundary Conditions for Temperature for solid   */
/* ----------------------------------------------------------------------- */


void TEMP_BoundaryCondition_Solid(
                            double dx,
                            double dy,
                            int imax,
                            int jmax,
                            double **TEMP,
                            double **TEMP_S,
                            const int wlt,
                            const int wrt,
                            const int wtt,
                            const int wbt,
                            const double TL,
                            const double TR,
                            const double TB,
                            const double TT,
                            const double QL,
                            const double QR,
                            const double QB,
                            const double QT,
                            double kratio,
                            int **Flag
                            ) {
    /* Set values for all the outside boundary depending on the value that
	 * the data file is set to. Dirichlet = 1 and Neumann = 2
     * We start with the left boundary. */
	int i,j;
    
    
    switch(wlt){
        case Dirichlet:
            for (j = 0; j <= jmax + 1; j++){
                TEMP_S[0][j]=2.0*TL-1*TEMP_S[1][j];
            }
            break;
        case Neumann:
            for (j = 0; j <= jmax + 1; j++){
                if (( Flag[1][j]&B_C ) != B_C ){
                    TEMP_S[0][j] = TEMP_S[1][j]+QL*kratio*dx;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the right boundary*/
    
	switch(wrt){
        case Dirichlet:
            for (j = 0; j <= jmax + 1; j++){
                TEMP_S[imax+1][j]=2.0*TR-1*TEMP_S[imax][j];
            }
            break;
        case Neumann:
            for (j = 0; j <= jmax + 1; j++){
                if (( Flag[imax][j]&B_C ) != B_C ){
                    TEMP_S[imax+1][j] = TEMP_S[imax][j]-QR*kratio*dx;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the top boundary*/
	switch(wtt){
        case Dirichlet:
            for (i = 0; i <= imax + 1; i++){
                TEMP_S[i][jmax+1]= 2.0*TT-1*TEMP_S[i][jmax];
            }
            break;
        case Neumann:
            for (i = 0; i <= imax + 1; i++){
                if (( Flag[i][jmax]&B_C ) != B_C ){
                    TEMP_S[i][jmax+1] = TEMP_S[i][jmax]-QT*kratio*dy;
                }
            }
            break;
        default:
            break;
	}
    
	/*Set values for the bottom boundary*/
	switch(wbt){
        case Dirichlet:
            for (i = 0; i <= imax + 1; i++){
                TEMP_S[i][0]=2.0*TB-1*TEMP_S[i][1];
            }
            break;
        case Neumann:
            for (i = 0; i <= imax + 1; i++){
                if (( Flag[i][1]&B_C ) != B_C ){
                    TEMP_S[i][0] = TEMP_S[i][1]+QB*kratio*dy;
                }
            }
            break;
        default:
            break;
	}
    
    /*treatment in case we have boundary inside*/
    /**
	 * Loop to check for boundary cells in the inner domain, and assign
	 * correct values of T.
	 */
	for(i = 1; i < imax+1; i++){
		for(j = 1; j < jmax+1; j++){
            /*
			 * If it's not a fluid cell but it has fluid cell neighbors, some values must still
			 * be calculated using the boundary flags.
			 */
			if((Flag[i][j]&31)==BS_N){
				TEMP_S[i][j] = ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j+1] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==BS_S){
				TEMP_S[i][j] = ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j-1] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==BS_W){
				TEMP_S[i][j] = ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i-1][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==BS_O){
				TEMP_S[i][j] = ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i+1][j] ) / (1+kratio);
			}
			else if((Flag[i][j]&31)==BS_NO){
				TEMP_S[i][j] = ( ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j+1] ) / (1+kratio) + ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i+1][j] ) / (1+kratio) )/2;
			}
			else if((Flag[i][j]&31)==BS_NW){
				TEMP_S[i][j] = ( ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j+1] ) / (1+kratio) + ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i-1][j] ) / (1+kratio) )/2;
			}
			else if((Flag[i][j]&31)==BS_SO){
				TEMP_S[i][j] = ( ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j-1] ) / (1+kratio) + ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i+1][j] ) / (1+kratio) )/2;
			}
			else if((Flag[i][j]&31)==BS_SW){
				TEMP_S[i][j] = ( ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i][j-1] ) / (1+kratio) + ( 2*kratio*TEMP[i][j]+(1-kratio)*TEMP_S[i-1][j] ) / (1+kratio) )/2;
			}
        }
	}
    
}

void TEMP_SpecBoundaryCondition_Solid(
                                char *problem,
                                int imax,
                                int jmax,
                                double **TEMP_S,
                                int **Flag                                      
                                ) {

}
