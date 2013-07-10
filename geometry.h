#ifndef GEOMETRY_H_
#define GEOMETRY_H_

/**
 * Create .pgm file for geometry
 */
void create_geometry(
		const char* szFileName,
		int imax,
		int jmax,
		double dx,
		double dy,
		int fins,
		int finWidth,
		int finHeight,
		int finDistance,
		int finCore,
		int option
);

#endif /* GEOMETRY_H_ */
