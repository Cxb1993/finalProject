#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include "geometry.h"
#include "helper.h"

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
){
	int i, j, k;
	double realWidth;
	double realHeight;
	double realDistance;
	int start;
	int end;
	int **Geometry;

	FILE * fh;
	fh = 0;
	start = 0;

	Geometry = imatrix(0, imax-1, 0, jmax-1);
	fh = fopen( szFileName, "w");	/* overwrite file/write new file */
	if( fh == NULL )			/* opening failed ? */
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
		ERROR( szBuff );
	}
	if(finDistance%2!=0){
		finDistance++;
		printf("Error! Fin distance must be an even number, value modified to %i\n", finDistance);
	}
	if(finWidth%2!=0){
		finWidth++;
		printf("Error! Fin distance must be an even number, value modified to %i\n", finWidth);
	}

	init_imatrix(Geometry, 0, imax-1, 0, jmax-1, 255);

	if(option==horizontal){
		realWidth = dx * (double)finWidth;
		realHeight = dy * (double)finHeight;
		realDistance = dx * (double)finDistance;
		while(imax < (fins*(finDistance+finWidth)-finDistance)&& fins > 0){
			fins--;
		}
		if(fins>0){
			start = (imax - (fins*(finDistance+finWidth)-finDistance))/2;
			end = start + fins*(finDistance+finWidth) - finDistance;
			printf("Geometry: %i fins of width:%f, height:%f, and distance in between of %f\n", fins, realWidth, realHeight, realDistance);
			fprintf(fh, "P3\n");
			fprintf(fh, "# %i fins of width:%f, height:%f, and distance in between of %f\n", fins, realWidth, realHeight, realDistance);
			fprintf(fh, "%i %i\n", imax, jmax);
			for(k = 0; k < fins; k++){
				for(i = 1; i <= finWidth; i++){
					for(j = finCore; j < finHeight+finCore; j++){
						Geometry[start+(k*(finWidth+finDistance)+i)][j]= 0;
					}
				}
			}
			for(i = start; i < end; i++){
				for(j = 0; j < finCore; j++){
					Geometry[i][j]= 0;
				}
			}
		}
		for(j = jmax-1; j >= 0; j--){
			for(i = 0; i < imax; i++){
				fprintf(fh, "%i\n", Geometry[i][j]);
			}
		}
	}

	if(option==vertical){
		realWidth = dx * (double)finHeight;
		realHeight = dy * (double)finWidth;
		realDistance = dy * (double)finDistance;
		while(jmax < (fins*(finDistance+finWidth)-finDistance)&& fins > 0){
			fins--;
		}
		if(fins>0){
			start = (imax/2)-finCore-finHeight;
			end = start + 2*finCore + 2*finHeight;
		}
		printf("Geometry: %i fins of width:%f, height:%f, and distance in between of %f\n", 2*fins, realWidth, realHeight, realDistance);
		fprintf(fh, "P3\n");
		fprintf(fh, "# %i fins of width:%f, height:%f, and distance in between of %f\n", 2*fins, realWidth, realHeight, realDistance);
		fprintf(fh, "%i %i\n", imax, jmax);
		for(k = 0; k < fins; k++){
			for(j = 0; j < finWidth; j++){
				for(i = start; i < end; i++){
					Geometry[i][j+k*(finWidth+finDistance)]= 0;
				}
			}
		}
		for(k = 0; k < fins-1; k++){
			for(j = 0; j < finDistance; j++){
				for(i = start+finHeight; i < start+finHeight+2*finCore; i++){
					Geometry[i][j+finWidth+k*(finWidth+finDistance)]= 0;
				}
			}
		}
		for(j = jmax-1; j >= 0; j--){
			for(i = 0; i < imax; i++){
				fprintf(fh, "%i\n", Geometry[i][j]);
			}
		}
	}

	if(fclose(fh))
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

	free_imatrix(Geometry, 0, imax-1, 0, jmax-1);

}
