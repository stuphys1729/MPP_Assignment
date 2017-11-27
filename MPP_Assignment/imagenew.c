/*
 * A file to reconstruct the original image from an edge file.
 *
 * Throughout the program, the following coordinate convention is assumed
 * For an array[i][j]:
 
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
#include <mpi.h>

#include "MPI_Comms.h"
//#include "FAKE_Comms.h"
#include "pgmio.h"
#include "arralloc.h"
#include "precision.h"

#define MAXITER   1500
#define PRINTFREQ  20
#define MAX_DIMS 2

#define TRUE 1
#define FALSE 0
/*
#define M 192
#define N 128

#define P 4

#define MP M/2
#define NP N/2
*/

RealNumber boundaryval(int i, int m);

int main(int argc, char **argv) {

	MPI_Init(NULL, NULL);

	RealNumber **old, **new, **edge, **masterbuf, **buf, **tempbuf;

	//RealNumber old[MP + 2][NP + 2], new[MP + 2][NP + 2], edge[MP + 2][NP + 2];
	int dims[MAX_DIMS];//, domains[2][2][2], counts[P], disps[P];

	//RealNumber masterbuf[M][N];
	//RealNumber buf[MP][NP];

	int i, j, iter;
	char *filename;

	int rank, size;
	
	int ***domains, *disps, *counts; //*hor_send_counts, *vert_send_counts, *hor_disps, *vert_disps;

	//initialise_MP(&comm, &rank, &size, dims);
	

	// Initialise a cartesian topology
	int periods[MAX_DIMS];

	for (i = 0; i < MAX_DIMS; i++) {
		dims[i] = 0;
		if (i == 0) periods[i] = TRUE;
		else periods[i] = FALSE;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm comm;
	MPI_Dims_create(size, MAX_DIMS, dims);
	MPI_Cart_create(MPI_COMM_WORLD, MAX_DIMS, dims, periods, TRUE, &comm);

	
	MPI_Comm_rank(comm, &rank);
	
	// The filename should be passed in to the program
	//filename = argv[1];
	filename = "edgenew192x128.pgm";

	/* Section for dynamic arrays */
	int M, N;	
	pgmsize(filename, &M, &N);

	// Configures the maximum domaian size
	int MP = ceil((RealNumber)M / (RealNumber)dims[0]);
	int NP = ceil((RealNumber)N / (RealNumber)dims[1]);

	//printf("MP: %d", MP);
	//printf("NP: %d", NP);

	// Allocate space for the arrays
	tempbuf = (RealNumber **)arralloc(sizeof(RealNumber), 2, M, NP);
	buf = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP, NP);

	new = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	old = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	edge = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);

	domains = (int ***)arralloc(sizeof(int), 3, dims[0], dims[1], 2);
	disps = (int *)arralloc(sizeof(int), 1, size);
	counts = (int *)arralloc(sizeof(int), 1, size);

	

	printf("rank %d is here\n",rank);

	if (rank == 0) {

		printf("Processing %d x %d image\n", M, N);
		printf("Number of iterations = %d\n", MAXITER);

		masterbuf = (RealNumber **)arralloc(sizeof(RealNumber), 2, M, N);

		printf("\nReading <%s>\n", filename);
		pgmread(filename, masterbuf, M, N);
		printf("Done.\n");
	}

	/*	Configure the domain sizes, giving more work to the top and bottom processes
		since they have one less communication direction (non-periodic)
		but still ensuring that each row of processes has the same number of column pixels,
		and each column of processes has the same number of row pixels
	*/
	printf("Rank %d Determining Domain Sizes\n", rank);

	int base_i = floor((RealNumber)M / (RealNumber)dims[0]);
	int base_j = floor((RealNumber)M / (RealNumber)dims[1]);
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
	for (i = dims[0]; i > 0; i--) {
		if (rem_i == 0) break;
		for (j = 1; j < dims[1]; j++) {
			domains[i][j][1]++;
		}
		rem_i--;
	}

	printf("Determined Domain Sizes\n");
	

	/* Now distribute the image across the processors */

	/* Stuff for sending as blocks
	MPI_Datatype Base_section;
	MPI_Datatype Extra_row;
	MPI_Datatype Extra_column;
	MPI_Datatype Max_section;

	MPI_Type_vector(MP-1, NP-1, N, MPI_REALNUMBER, &Base_section);
	MPI_Type_vector(MP-1, NP, N, MPI_REALNUMBER, &Extra_row);
	MPI_Type_vector(MP, NP-1, N, MPI_REALNUMBER, &Extra_column);
	MPI_Type_vector(MP, NP, N, MPI_REALNUMBER, &Max_section);
	*/

	/* Stuff for sending same-size blocks */
	/*	We can send each processor the same (maximum) sized blocks, but they can
		overwrite bits that are not part of their domain */

	

	MPI_Datatype Send_section;
	MPI_Type_vector(MP, NP, N, MPI_REALNUMBER, &Send_section);

	MPI_Datatype Recv_section;
	MPI_Type_vector(MP, NP, NP + 2, MPI_REALNUMBER, &Recv_section);

	for (i = 0; i < size; i++) {
		counts[i] = 1;
	}
	int offset = 0;
	for (i = 0; i < dims[0]; i++) {

		for (j = 0; j < dims[1]; j++) {
			disps[i + j] = offset;
			if (j < (dims[1]-1)) offset += domains[i][j][0];
		}
		if (i == (size - 1)) break;
		offset = 0;
		for (int x = 0; x < i + 1; x++) {
			offset += domains[x][0][1] * M;
		}
		
	}
	printf("About to scatter\n");
	MPI_Scatterv(&masterbuf[0][0], counts, disps, Send_section, buf, 1, Recv_section, 0, comm);

	

	/* Stuff for sending with cart_sub 

	hor_send_counts = (int *)arraloc(sizeof(int), 1, dims[0]);
	vert_send_counts = (int *)arraloc(sizeof(int), 1, dims[1]);
	vert_disps = (int *)arraloc(sizeof(int), 1, dims[1]);
	hor_disps = (int *)arraloc(sizeof(int), 1, dims[0]);

	int offset = 0;
	for (j = 0; j < dims[1]; j++) {
	vert_send_counts[j] = domains[0][j][1];
	vert_disps[j] = offset;
	offset += vert_send_counts[j];
	}
	offset = 0;
	for (i = 0; i < dims[0]; i++) {
	hor_send_counts[i] = domains[i][0][0];
	hor_disps = offset;
	offset += hor_send_counts[i];
	}

	// Sending down
	MPI_Datatype Base_row;

	MPI_Type_vector(M, 1, N, MPI_REALNUMBER, &Base_row);

	MPI_Comm vertical;
	MPI_Cart_sub(comm, (FALSE, TRUE), &vertical);

	MPI_Scatterv(&masterbuf[0][0], vert_send_counts, vert_disps, Base_row, &tempbuf, NP, Base_row, 0, vertical);


	// Sending across
	MPI_Datatype Base_column;
	MPI_Datatype Extra_column;

	MPI_Type_vector(1, NP - 1, N, MPI_REALNUMBER, &Base_column);
	MPI_Type_vector(1, NP, N, MPI_REALNUMBER, &Extra_column);
	
	MPI_Comm horizontal;
	MPI_Cart_sub(comm, (TRUE, FALSE), &horizontal);

	MPI_Scatterv(&tempbuf[0][0], hor_send_counts, hor_disps, )

	*/

	/*
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
	*/
	/* Set fixed boundary conditions on the bottom and top edges 
	for (i=1; i < M+1; i++) {
		// compute sawtooth value 
     
		val = boundaryval(i, M);

		old[i][0]   = (int)(255.0*val);
		old[i][N+1] = (int)(255.0*(1.0-val));
	}
	*/

	/*
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
	*/
	//printf("\nFinished %d iterations\n", iter-1);

	/* Gather the data back to process 0 */

	if (rank == 0) {

		filename = "imagenew192x128.pgm";
		printf("\nWriting <%s>\n", filename);
		//pgmwrite(filename, buf, M, N);

	}

	MPI_Finalize();
} 

RealNumber boundaryval(int i, int m) {

	RealNumber val;

	val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
	if (i >= m/2+1) val = 2.0-val;
  
	return val;
}