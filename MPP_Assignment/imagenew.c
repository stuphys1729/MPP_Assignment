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

#include "pgmio.h"
#include "arralloc.h"
#include "precision.h"
//#include "MPI_Comms.h"
 //#include "FAKE_Comms.h"

#define MAXITER 1500
#define PRINTFREQ  200
#define MAX_DIMS 2

#define TRUE 1
#define FALSE 0
#define LEFT_TAG 5
#define RIGHT_TAG 10
#define DEFAULT_TAG 0

RealNumber boundaryval(int i, int m);

int main(int argc, char **argv) {

	MPI_Init(NULL, NULL);

	RealNumber **old, **new, **edge, **masterbuf, **buf, **tempbuf;

	int i, j, iter;
	char *filename;

	RealNumber start, taken;

	int rank, size, dims[MAX_DIMS];
	
	int ***domains, *disps, *counts;
	
	/* Initialise a cartesian topology */
	int periods[MAX_DIMS];

	for (i = 0; i < MAX_DIMS; i++) {
		dims[i] = 0;
		// Only periodic boundary conditions horizontally
		if (i == 0) periods[i] = TRUE;
		else periods[i] = FALSE;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm comm;
	MPI_Dims_create(size, MAX_DIMS, dims);
	
	/* The filename should be passed in to the program */
	filename = argv[1];

	/* Section for dynamic arrays */
	int M, N;	
	pgmsize(filename, &M, &N);

	if (M % dims[0] || N % dims[1]) {
		// This configuration does not give equal domain sizes
		if (!(M % dims[1]) && !(N % dims[0])) {
			// Swapping the orientation works
			int temp = dims[0];
			dims[0] = dims[1];
			dims[1] = temp;
		}
		else {
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0) {
				printf("The passed image and processor combination is not compatible\n");
			}
			MPI_Finalize();
			exit(1);
		}
	}

	MPI_Cart_create(MPI_COMM_WORLD, MAX_DIMS, dims, periods, TRUE, &comm);
	MPI_Comm_rank(comm, &rank);

	/* Configures the maximum domaian size */
	int MP = M / dims[0];
	int NP = N / dims[1];

	/* Allocate space for the arrays */
	new = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	old = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	edge = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);

	domains = (int ***)arralloc(sizeof(int), 3, dims[0], dims[1], 2);
	disps = (int *)arralloc(sizeof(int), 1, size);
	counts = (int *)arralloc(sizeof(int), 1, size);


	if (rank == 0) {
		printf("Processing %d x %d image\n", M, N);
		printf("Number of iterations = %d\n", MAXITER);

		masterbuf = (RealNumber **)arralloc(sizeof(RealNumber), 2, M, N);

		printf("\nReading <%s>\n", filename);
		pgmread(filename, masterbuf, M, N);
	}

	int sizes[MAX_DIMS] = { M,N };
	int sub_sizes[MAX_DIMS] = { MP,NP };
	int starts[MAX_DIMS] = { 0,0 };

	MPI_Datatype Send_section;
	MPI_Datatype Small_send_section;
	MPI_Type_create_subarray(2, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_REALNUMBER, &Send_section);

	/*	We change MPI's understanding of where each section ends to be just one number after
		it starts to allow us to give precise locations of domain beginnings	*/
	MPI_Type_create_resized(Send_section, 0, NP*sizeof(RealNumber), &Small_send_section);
	MPI_Type_commit(&Small_send_section);

	MPI_Datatype Recv_section;
	MPI_Type_vector(MP, NP, NP + 2, MPI_REALNUMBER, &Recv_section);
	MPI_Type_commit(&Recv_section);

	for (i = 0; i < size; i++) {
		counts[i] = 1;
	}
	int offset = 0;
	for (i = 0; i < dims[0]; i++) {

		for (j = 0; j < dims[1]; j++) {
			disps[i*dims[1] + j] = offset;
			offset++;
		}
		offset = (i+1)*MP*dims[1];
		
	}
	if (rank == 0) {
		for (i = 0; i < size; i++) {
			printf("disp: %d\n", disps[i]);
		}
	}

	/* Begin timing the computation */
	start = MPI_Wtime();

	if (rank == 0) {
		printf("Scattering the original image.\n");
	}
	MPI_Scatterv(masterbuf, counts, disps, Small_send_section, &edge[1][1], 1, Recv_section, 0, comm);

	/* Initialise our 'guess' of the image to a bright square */
	for (i=0; i<MP+2;i++) {
		for (j=0;j<NP+2;j++) {
			old[i][j]=255.0;
		}
	}

	int up, down, left, right;
	RealNumber val;

	/* Find out where to send and receive halo swaps */
	MPI_Cart_shift(comm, 0, 1, &left, &right);
	MPI_Cart_shift(comm, 1, 1, &down, &up);

	/* For the sawtooth boundary conditions, we must find the overall
		pixel position, depending on the processor's domian */
	offset = ( rank / dims[1]) * MP;
	if (up == MPI_PROC_NULL) {
		for (i = 1; i < MP+1; i++) {

			val = boundaryval(offset + i, M);
			old[i][NP+1] = (int)(255.0*(1.0-val));
		}
	}

	if (down == MPI_PROC_NULL) {
		for (i = 1; i < MP+1; i++) {

			val = boundaryval(offset+i, M);
			old[i][0] = (int)(255.0*val);
		}
	}
	
	MPI_Datatype sides, top_bottom;
	MPI_Type_contiguous(NP, MPI_REALNUMBER, &sides);
	MPI_Type_vector(MP, 1, NP + 2, MPI_REALNUMBER, &top_bottom);
	
	MPI_Request requests[2*(2*MAX_DIMS)];
	MPI_Status statuses[2*(2*MAX_DIMS)];
	
	for (iter= 1;iter<=MAXITER; iter++) {
		if (iter%PRINTFREQ == 0 && rank==0) {
			printf("Iteration %d\n", iter);
		}
		
		/* Due to periodic boundaries left and right, we need a different tag for
		   the two different messages in case 'left' and 'right' are the same process */
		MPI_Isend(&old[1][1], 1, sides, left, RIGHT_TAG, comm, &requests[0]);// send left
		MPI_Irecv(&old[0][1], 1, sides, left, LEFT_TAG, comm, &requests[1]);// recv left

		MPI_Isend(&old[MP][1], 1, sides, right, LEFT_TAG, comm, &requests[2]);// send right
		MPI_Irecv(&old[MP+1][1], 1, sides, right, RIGHT_TAG, comm, &requests[3]);// recv right
		
		MPI_Isend(&old[1][1], 1, top_bottom, down, DEFAULT_TAG, comm, &requests[4]);// send down
		MPI_Irecv(&old[1][0], 1, top_bottom, down, DEFAULT_TAG, comm, &requests[5]);// recv down

		MPI_Isend(&old[1][NP], 1, top_bottom, up, DEFAULT_TAG, comm, &requests[6]);// send up
		MPI_Irecv(&old[1][NP+1], 1, top_bottom, up, DEFAULT_TAG, comm, &requests[7]);// recv up
		
		/* We can work on the internal pixels whilst we wait for the halo-swaps */
		for (i = 2; i < MP; i++) {
			for (j = 2; j < NP; j++) {
			new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
					- edge[i][j]);
			}
		}
		
		/* Make sure that all of the communications have completed */
		MPI_Waitall(2 * MAX_DIMS, requests, statuses);
		
		/* Calculate the sides */
		for (j = 1; j < NP + 1; j++) {
			new[1][j] = 0.25*(old[0][j] + old[2][j] + old[1][j - 1] + old[1][j + 1]
				- edge[1][j]);
			new[MP][j] = 0.25*(old[MP - 1][j] + old[MP + 1][j] + old[MP][j - 1] + old[MP][j + 1]
				- edge[MP][j]);
		}
		/* Calculate the top and bottom */
		for (i = 1; i < MP; i++) {
			new[i][1] = 0.25*(old[i - 1][1] + old[i + 1][1] + old[i][0] + old[i][2]
				- edge[i][1]);
			new[i][NP] = 0.25*(old[i - 1][NP] + old[i + 1][NP] + old[i][NP - 1] + old[i][NP + 1]
				- edge[i][NP]);
		}
		
		/* Copy the new array into the old for the next iteration */
		for (i = 1; i<MP+1; i++) {
			for (j = 1; j<NP+1; j++) {
				old[i][j] = new[i][j];
			}
		}

	}
	if (rank == 0) {
		printf("\nFinished %d iterations\n", iter - 1);
	}	
	
	/* Gather the data back to process 0 */
	MPI_Gatherv(&old[1][1], 1, Recv_section, masterbuf, counts, disps, Small_send_section, 0, comm);

	/* Stop timing the computation */
	taken = MPI_Wtime() - start;

	if (rank == 0) {

		printf("Total Computation Time: %f", taken);
		sprintf(filename, "imagenew%dx%d_%d.pgm", M, N, size);
		printf("\nWriting <%s>\n", filename);
		pgmwrite(filename, masterbuf, M, N);

	}

	MPI_Finalize();
} 

RealNumber boundaryval(int i, int m) {

	RealNumber val;

	val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
	if (i >= m/2+1) val = 2.0-val;
  
	return val;
}
