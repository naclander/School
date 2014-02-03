/**
 *  cdf97.c - Fast discrete biorthogonal CDF 9/7 wavelet forward and inverse transform (lifting implementation)
 *  
 *  written by Di Wu - allwu@cs.ucla.edu
 *  based on the version from Gregoire Pau (2006) - gregoire.pau@ebi.ac.uk
 *
 *  removed global states for better parallelization
 *  only supports image sizes that are product of 2^k
 *
 *  2013.7
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cdf97.h"

#define a1 (-1.586134342)
#define a2 (-0.05298011854)
#define a3 (0.8829110762)
#define a4 (0.4435068522)

#define k1 1.1496043988602418
//#define k2 (1/1.1496043988602418)

#define x(i, j)			x[(i)*M+(j)]
#define tempbank(i, j)	tempbank[(i)*n+(j)]

/**
 * predict and update
 *
 * computation is an even-odd update:
 * first update all odd columns based on even columns, then update even columns based on odd columns
 */
void fwt97_pd(
		float* x, 
		int n, int m,
		int N, int M)
{
	int i, j, jj, ii, iii, iiii;
	#pragma omp parallel for private(i,j) 
	for (i=0; i<n; i+=1) {	
		ii = i+1;
		iii = ii+1;
		iiii = iii+1;


		j = 1;
		x(i, j) += a1*(x(i, j-1) + x(i, j+1));
		x(i, 0) += 2*a2*x(i, 1);
		j+=2;
		x(i, j) += a1*(x(i, j-1) + x(i, j+1));
		j=2;
		x(i, j) += a2*(x(i, j-1) + x(i, j+1));
		j+= 3;
		jj = 1;
		x(i, j) += a1*(x(i, j-1) + x(i, j+1));
		x(i, jj) += a3*(x(i, jj-1) + x(i, jj+1));
		x(i, 0) += 2*a4*x(i, 1);
		j = 4;
		x(i, j) += a2*(x(i, j-1) + x(i, j+1));
		j+=3;
		x(i, j) += a1*(x(i, j-1) + x(i, j+1));
		jj = 3;
		x(i, jj) += a3*(x(i, jj-1) + x(i, jj+1));
		j = 6;
		x(i, j) += a2*(x(i, j-1) + x(i, j+1));
		j+=3;
		x(i, j) += a1*(x(i, j-1) + x(i, j+1));


		//#pragma omp parallel for private(j)
		for(j = 8, jj = 2; j<m-4; j--, jj--){
			//printf("j : %d, jj: %d\n",j,jj);
			x(i, j) += a2*(x(i, j-1) + x(i, j+1));
			x(i, jj) += a4*(x(i, jj-1) + x(i, jj+1));
			
			j+=3;
			jj+=3;

			x(i, j) += a1*(x(i, j-1) + x(i, j+1));
			x(i, jj) += a3*(x(i,jj-1) + x(i, jj+1));
		}
		x(i, j) += a2*(x(i, j-1) + x(i, j+1));
		x(i, m-1) += 2*a1*x(i, m-2);
		x(i, m-2) += a2*(x(i, m-2-1) + x(i, m-2+1));	
		x(i, jj) += a4*(x(i, jj-1) + x(i, jj+1));
		jj+=3;
		x(i, jj) += a3*(x(i,jj-1) + x(i, jj+1));
		jj--;
		
		x(i, jj) += a4*(x(i, jj-1) + x(i, jj+1));
		jj+=3;
		x(i, jj) += a3*(x(i,jj-1) + x(i, jj+1));
		jj--;
		
		x(i, jj) += a4*(x(i, jj-1) + x(i, jj+1));
		jj+=3;
		x(i, jj) += a3*(x(i,jj-1) + x(i, jj+1));
		jj--;
		x(i, j) += a4*(x(i, j-1) + x(i, j+1));
		x(i, m-1) += 2*a3*x(i, m-2);
		x(i, m-2) += a4*(x(i, m-2-1) + x(i, m-2+1));

	}
}

/**
 * de-interleave with transpose
 *
 * two different function (foo and bar) are used to move data
 * between x and tempbank
 * the reason for this design is to avoid 'spilling' in x when data
 * is not square (n!=m)
 *
 **/
void fwt97_dl_foo (
		float* x, 
		float* tempbank, 
		int n, int m, 
		int N, int M) 
{
	int i, j;
	double k2 = 1/1.1496043988602418;
    // de-interleave and transpose
	#pragma omp parallel for private(i,j)
	for (i=0; i<n; i++) {
        // simultaneously transpose the matrix when deinterleaving
		for (j=0; j<m/2; j+= 2) {
			tempbank(j/2, i) = k1*x(i, j);
		}
		int jj;
		for (jj = 1; j<m; j+= 2, jj += 2) {
			tempbank(j/2, i) = k1*x(i, j);
			tempbank(jj/2+m/2, i) = k2*x(i, jj);
		}
		for(; jj<m; jj+= 2){
				tempbank(jj/2+m/2, i) = k2*x(i, jj);
		}
	}
}

void fwt97_dl_bar (
		float* x, 
		float* tempbank, 
		int n, int m, 
		int N, int M) 
{
	int i, j;
	double k2 = 1/1.1496043988602418;
    // de-interleave and transpose
	#pragma omp parallel for private(i,j)
	for (i=0; i<n; i++) {
        // simultaneously transpose the matrix when deinterleaving
		for (j=0; j<m/2; j+= 2) {
			x[(j/2)*N+i] = k1*tempbank[i*m+j];
		}
		int jj;
		for (j=m/2, jj = 1; j<m; j+= 2, jj+=2) {
			x[(j/2)*N+i] = k1*tempbank[i*m+j];
			x[(jj/2+m/2)*N+i] = k2*tempbank[i*m+jj];
		}
		for (; jj<m; jj+= 2) {
			x[(jj/2+m/2)*N+i] = k2*tempbank[i*m+jj];
		}
	}
}

/**
 *	CDF9/7 wavelet transform
 *
 *	this version assumes a column-based memory layout (Matlab default layout)
 *	so the wavelet along columns of the original matrix would be the transform
 *	along row in this code.
 *
 *	input matrix: nxm (n cols, m rows)
 *
 **/
void cdf97(
		float* x, 
		float* tempbank, 
		int n, int m, 
		int level
		) 
{

	int i;
	int nlevel = (level < 0 ? -level : level);

	int N = n;
	int M = m;

	int plevel = 1<<nlevel-1;

	// this means image is not big enough for transformation 
	// which we assume will not happen
	if (N/plevel < 2 || M/plevel < 2) 
		return;

	if (level > 0) {	
		// forward DWT (1. col; 2. row)
		//#pragma omp parallel for
		for (i=0; i<nlevel; i++) {

			// Stage 1
			fwt97_pd(x, m, n, M, N);

			// Stage 2
			fwt97_dl_foo(x, tempbank, m, n, M, N);

			// Stage 3
			fwt97_pd(tempbank, n, m, n, m);

			// Stage 4
			fwt97_dl_bar(x, tempbank, n, m, N, M);

			// work on the upper left image (low pass filtered) in the next level
			n /= 2;
			m /= 2;
		}
	}
	else {				
		// backward DWT, omitted for the project
		;	
	}
}


