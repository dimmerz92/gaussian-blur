#include "qdbmp.h"		/*Our Bitmap operations library */
#include <stdio.h>		/*For I/O */
#include <math.h>		/*For maths operations e.g pow,sqrt */
#include <stdlib.h>		/*For utils e.g. malloc */
#include <time.h>		/*For clock(...) operations */

#ifndef  M_PI
#define  M_PI  3.1415926535897932384626433
#endif


void generateGaussianKernel (double **kernel, int kernel_dim, double sd,
			     int origin, double *kernel_max,
			     double *colour_max);
void GroundColorMix (double *color, double x, double min, double max);
void bitmapFromSquareMatrix (float **mat, const char *filename, int mat_dim,
			     float mat_max, float colour_min,
			     float colour_max);
void applyConvolution (double **kernel, int kernel_dim, double kernel_origin,
		       double colour_max, BMP * old_bmp, BMP * new_bmp);
