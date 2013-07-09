#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


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
                    );
/**
 * Function to set special boundary values.
 */
void spec_boundary_val(
                       char *problem,
                       int imax,
                       int jmax,
                       double **U,
                       double **V);

/**
 * The temperature boundary values of the problem are set.
 */

void TEMP_BoundaryCondition(
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
                            double kratio,
                            int **Flag
                            );

/**
 * Function to set special boundary values for temperature.
 */

void TEMP_SpecBoundaryCondition(
                                char *problem,
                                int imax,
                                int jmax,
                                double **TEMP
                                );


#endif
