/*
 * A file to reconstruct the original image from an edge file.
 *
 * Throughout the program, the following coordinate convention is assumed:
 
	  j	^
		|
		|
		|
		|
		|
		|
		|____________________>
							i

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MPI_Comms.h"
//#include "FAKE_Comms.h"
#include "pgmio.h"
#include "arralloc.h"
#include "precision.h"

#define MAXITER   1500
#define PRINTFREQ  20
#define MAX_DIMS 2

RealNumber boundaryval(int i, int m);

int main (int argc, char **argv) {


	RealNumber **old, **new, **edge, **masterbuf, **buf;

	int i, j, iter, maxiter;
	char *filename;

	int rank, comm;
	int dims[MAX_DIMS];
	int ***domains;

	initialise_MP(&comm, &rank, dims);
	
	// The filename should be passed in to the program
	filename = argv[1];

	int M, N;	
	pgmsize(filename, &M, &N);

	// Configures the maximum domaian size
	int MP = ceil((RealNumber)M / (RealNumber)dims[0]);
	int NP = ceil((RealNumber)N / (RealNumber)dims[1]);

	// Allocate space for the arrays
	masterbuf = (RealNumber **)arralloc(sizeof(RealNumber), 2, M, N);
	buf = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP, NP);

	new = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	old = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	edge = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);

	domains = (int ***)arraloc(sizeof(int), 3, dims[0], dims[1], 2);

	int i, j, iter, maxiter;
	RealNumber val;

	printf("Processing %d x %d image\n", M, N);
	printf("Number of iterations = %d\n", MAXITER);

	

	printf("\nReading <%s>\n", filename);
	pgmread(filename, masterbuf, M, N);
	printf("\n");

	// Configure the domain sizes, giving more work to the top and bottom processes
	// since they have one less communication direction (non-periodic)
	// but still ensuring that each row of processes has the same number of column pixels,
	// and each column of processes has the same number of row pixels

	int base_i = MP - 1;
	int base_j = NP - 1;
	int rem_i = M - (base_i*dims[0]);
	int rem_j = M - (base_j*dims[1]);

	for (i = 0; i < dims[0]; i++) {
		for (j = 0; j < dims[1]; j++) {
			domains[i][j][0] = base_i;
			domains[i][j][1] = base_j;
		}
	}

	// give the first extra row of pixels to the bottom set of processors
	if (rem_j > 0) {
		for (i = 0; i < dims[0]; i++) {
			domains[i][0][0]++;
		}
		rem_j--;
	}
	// Put the remaining extra rows from the first row down until we run out
	for (j = 0; j < rem_j; j++) { 
		for (i = 0; i < dims[0]; i++) {
			domains[i][j][0]++;
		}
	}

	// The remaining columns can be distributed to all but the first
	for (i = 1; i < rem_i; i++) {
		for (j = 1; j < dims[1]; j++) {
			domains[i][j][1]++;
		}
	}



	for (i=1;i<MP+1;i++) {
		for (j=1;j<NP+1;j++) {
			edge[i][j]=buf[i-1][j-1];
		}
	}

	for (i=0; i<M+2;i++) {
		for (j=0;j<N+2;j++) {
			old[i][j]=255.0;
		}
	}

	/* Set fixed boundary conditions on the bottom and top sides */

	for (i=1; i < M+1; i++) {
		/* compute sawtooth value */
     
		val = boundaryval(i, M);

		old[i][0]   = (int)(255.0*val);
		old[i][N+1] = (int)(255.0*(1.0-val));
	}

	for (iter=1;iter<=MAXITER; iter++) {
		if (iter%PRINTFREQ == 0) {
			printf("Iteration %d\n", iter);
		}


		for (i=1;i<M+1;i++) {
			for (j=1;j<N+1;j++) {
			new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
					- edge[i][j]);
			}
		}
	
		for (i=1;i<M+1;i++) {
			for (j=1;j<N+1;j++) {
				old[i][j]=new[i][j];
			}
		}
	}

	printf("\nFinished %d iterations\n", iter-1);

	for (i=1;i<M+1;i++) {
		for (j=1;j<N+1;j++) {
			buf[i-1][j-1]=old[i][j];
		}
	}

	filename="imagenew192x128.pgm";
	printf("\nWriting <%s>\n", filename); 
	pgmwrite(filename, buf, M, N);
} 

RealNumber boundaryval(int i, int m) {

	RealNumber val;

	val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
	if (i >= m/2+1) val = 2.0-val;
  
	return val;
}