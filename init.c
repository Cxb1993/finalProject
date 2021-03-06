#include "helper.h"
#include "init.h"

/* ----------------------------------------------------------------------- */
/*                             Reading Parameters                          */
/* ----------------------------------------------------------------------- */


/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 * @param wl,wr,wt,wb boundary type
 * @param problem	 define problem to be solved
 * @param lp, rp, dp defines values of pressure at the left and right boundary, or the difference keeping right constant.
 * @param argv		 input argument for the problem
 * @param argc		 count there is only one input
 */

int read_parameters(
        double *Re,                /* reynolds number   */
		double *Pr,                /* prandtl number   */
		double *beta,              /* coefficient of thermal expansion */
        double *kratio,            /* Ratio of thermal conductivity of fluid to solid k_f/k_s*/
        double *Kappa,             /* coefficient of thermal diffusivity */
		double *UI,                /* velocity x-direction */
		double *VI,                /* velocity y-direction */
		double *PI,                /* pressure */
		double *TI,                /* initial temperature for fluid */
        double *TSI,                /* initial temperature for solid */
		double *GX,                /* gravitation x-direction */
		double *GY,                /* gravitation y-direction */
		double *t_end,             /* end time */
		double *xlength,           /* length of the domain x-dir.*/
		double *ylength,           /* length of the domain y-dir.*/
		double *dt,                /* time step */
		double *dx,                /* length of a cell x-dir. */
		double *dy,                /* length of a cell y-dir. */
		int  *imax,                /* number of cells x-direction*/
		int  *jmax,                /* number of cells y-direction*/
		double *alpha,             /* uppwind differencing factor*/
		double *omg,               /* relaxation factor */
		double *gamma,
		double *tau,               /* safety factor for time step*/
		int  *itermax,             /* max. number of iterations  */
                                    /* for pressure per time step */
		double *eps,               /* accuracy bound for pressure*/
		double *dt_value,			/* time for output */
		int *wl,						/* boundary type for left wall (1:no-slip 2: free-slip 3: outflow) */
		int *wr,						/* boundary type for right wall (1:no-slip 2: free-slip 3: outflow) */
		int *wt,						/* boundary type for top wall (1:no-slip 2: free-slip 3: outflow) */
		int *wb,						/* boundary type for bottom wall (1:no-slip 2: free-slip 3: outflow) */
		char *problem,               /* name of problem */
		double *lp,					/* pressure at left boundary */
		double *rp,					/* pressure at right boundary */
		double *dp,					/* pressure difference */
		int *wlt,
		int *wrt,
		int *wtt,
		int *wbt,
		double *TL,
		double *TR,
		double *TB,
		double *TT,
        double *QL,
        double *QR,
        double *QB,
        double *QT,
		int argc,
		char *argv,
		int *fins,
		int *finWidth,
		int *finHeight,
		int *finDistance,
		int *finCore
)
{
	char szFileName[200];
	strcpy(szFileName, argv);
	strcat(szFileName, ".dat");

	if(argc<4){
		READ_DOUBLE( szFileName, *xlength );
		READ_DOUBLE( szFileName, *ylength );

		READ_DOUBLE( szFileName, *Re      );
		READ_DOUBLE( szFileName, *Pr      );
		READ_DOUBLE( szFileName, *beta    );
        READ_DOUBLE( szFileName, *kratio  );
        READ_DOUBLE( szFileName, *Kappa   );



		READ_DOUBLE( szFileName, *t_end );
		READ_DOUBLE( szFileName, *dt    );

		READ_INT   ( szFileName, *imax );
		READ_INT   ( szFileName, *jmax );

		READ_DOUBLE( szFileName, *omg   );
		READ_DOUBLE( szFileName, *gamma );

		READ_DOUBLE( szFileName, *eps   );
		READ_DOUBLE( szFileName, *tau   );
		READ_DOUBLE( szFileName, *alpha );

		READ_INT   ( szFileName, *itermax );
		READ_DOUBLE( szFileName, *dt_value );

		READ_DOUBLE( szFileName, *UI );
		READ_DOUBLE( szFileName, *VI );
		READ_DOUBLE( szFileName, *GX );
		READ_DOUBLE( szFileName, *GY );
		READ_DOUBLE( szFileName, *PI );
		READ_DOUBLE( szFileName, *TI );
        READ_DOUBLE( szFileName, *TSI );


		READ_INT( szFileName, *wl );
		READ_INT( szFileName, *wr );
		READ_INT( szFileName, *wt );
		READ_INT( szFileName, *wb );

		READ_INT( szFileName, *wlt );
		READ_INT( szFileName, *wrt );
		READ_INT( szFileName, *wtt );
		READ_INT( szFileName, *wbt );

		READ_DOUBLE( szFileName, *TL );
		READ_DOUBLE( szFileName, *TR );
		READ_DOUBLE( szFileName, *TB );
		READ_DOUBLE( szFileName, *TT );

		READ_DOUBLE( szFileName, *QL );
		READ_DOUBLE( szFileName, *QR );
		READ_DOUBLE( szFileName, *QB );
		READ_DOUBLE( szFileName, *QT );
        
		strcpy(problem, argv);

		READ_DOUBLE ( szFileName, *lp );
		READ_DOUBLE ( szFileName, *rp );
		READ_DOUBLE ( szFileName, *dp );

		READ_INT( szFileName, *fins );
		READ_INT( szFileName, *finWidth);
		READ_INT( szFileName, *finHeight );
		READ_INT( szFileName, *finDistance );
		READ_INT( szFileName, *finCore );

		*dx = *xlength / (double)(*imax);
		*dy = *ylength / (double)(*jmax);

		return 1;
	}
	else{
		/* In case there was only one argument print an error and return 0*/
				ERROR("One and only one argument (data file) should be passed to the function" );
				return 0;
	}
}

/* ----------------------------------------------------------------------- */
/*                             Initializing U,V,P & T                       */
/* ----------------------------------------------------------------------- */

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */

void init_uvp(
		double UI,
		double VI,
		double PI,
		double TI,
        double TSI,
		int imax,
		int jmax,
		double **U,
		double **V,
		double **P,
		double **TEMP,
        double **TEMP_S,
		int **Flag
)

{
	int i;
	int j;

	for ( i = 0 ; i<= imax+1;  i++ )
	{
		for (j = 0; j<= jmax+1; j++ )
		{
			/*
			 * Initialize fluid cells with the values from the data file
			 */
			if((Flag[i][j]&B_C)==B_C){
				U[i][j] = UI ;
				V[i][j] = VI ;
				P[i][j] = PI ;
				TEMP[i][j] = TI ;
			}
            else if ((Flag[i][j]&B_C)!=B_C){
                TEMP_S[i][j] = TSI;
            }
			else{
				U[i][j] = 0.0 ;
				V[i][j] = 0.0 ;
				P[i][j] = 0.0 ;
				TEMP[i][j] = 0.0 ;
                TEMP_S[i][j] = 0.0 ;
			}
		}
	}
    
}

/* ----------------------------------------------------------------------- */
/*                             Initializing Flag array                     */
/* ----------------------------------------------------------------------- */

/*The array Flag is initialized with the flags C_F for fluid cells and C_B for obstacle cells as
 specified by the parameter problem. This must be followed by a loop over all cells where
 the boundary cells are marked with the appropriate flags B_xy depending on the direction, in
 which neighboring fluid cells lie.*/
void init_flag(
		const char *problem,
		int imax,
		int jmax,
		double lp,
		double rp,
		double dp,
		int **Flag
){

	char image[84];
	int i, j;
	int **tmp;
	strcpy(image, problem);
	strcat(image, ".pgm");

	tmp = imatrix(0, imax, 0, jmax);
	tmp = read_pgm(image);

	/* Outer boundaries will always be the same, so we assign those flags first*/
	/* Corners: */
	init_imatrix(Flag, 0, imax, 0, jmax, 0);


	/*Initialize corners*/
	Flag[0][0]=0;
	Flag[0][jmax+1]=0;
	Flag[imax+1][0]=0;
	Flag[imax+1][jmax+1]=0;

	/* Left boundary: */
	for(j = 1; j < jmax + 1; j++){
		if((tmp[1][j] & B_C)){
			Flag[0][j] = B_O;
		}
		/*
		 * In case there is a value for the pressure that should be assigned in the left boundary
		 * set the P_L flag to 1. In this case it could be directly because of lp being defined or
		 * dp not equal 0.
		 */
		if(lp>=0||dp!=0){
			Flag[0][j]|= P_L;
		}
	}

	/* Right boundary: */
	for(j = 1; j < jmax + 1; j++){
		if((tmp[imax][j] & B_C)){
			Flag[imax+1][j] = B_W;
		}
		/*
		 * Set flag P_R if pressure value at right boundary is defined
		 */
		if(rp>=0||dp!=0){
			Flag[imax+1][j]|= P_R;
		}
	}

	/* Top boundary: */
	for(i = 1; i < imax + 1; i++){
		if((tmp[i][jmax] & B_C)){
			Flag[i][jmax+1] = B_S;
		}	}

	/* Bottom boundary: */
	for(i = 1; i < imax + 1; i++){
		if((tmp[i][1] & B_C)){
			Flag[i][0] = B_N;
		}
	}

	/*
	 * Now loop over all inner cells checking the four neighbors (no corners)
	 * and check if it is a fluid cell and set the flags for the neighboring cells
	 * */
	for(i = 1; i < imax + 1; i++){
		for(j = 1; j < jmax + 1; j++){
			if(tmp[i][j]==C_F){
				Flag[i][j] = B_C;
			}
			else{
				Flag[i][j] = 0;
			}
			if((tmp[i-1][j]&B_C)==B_C||(tmp[i-1][j]&C_F)==C_F){
				Flag[i][j] |= B_W;
			}
			if((tmp[i+1][j]&B_C)==B_C||(tmp[i+1][j]&C_F)==C_F){
				Flag[i][j] |= B_O;
			}
			if((tmp[i][j+1]&B_C)==B_C||(tmp[i][j+1]&C_F)==C_F){
				Flag[i][j] |= B_N;
			}
			if((tmp[i][j-1]&B_C)==B_C||(tmp[i][j-1]&C_F)==C_F){
				Flag[i][j] |= B_S;
			}
			/*
			 * This is a check for possible forbidden boundaries.
			 */
			if((Flag[i][j]&31)==3||(Flag[i][j]&31)==7||((Flag[i][j]&31)>10&&(Flag[i][j]&31)<16)){
				printf("\nERROR! The flag field contains a forbidden boundary cell at i= %i j= %i\n",i,j);
			}
		}
	}
	free_imatrix(tmp, 0, imax, 0, jmax);
}


