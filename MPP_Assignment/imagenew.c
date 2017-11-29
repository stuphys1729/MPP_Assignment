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

#define MAXITER 1500
#define PRINTFREQ  50
#define MAX_DIMS 2

#define TRUE 1
#define FALSE 0
#define DEFAULT_TAG 0
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
	

	/* Initialise a cartesian topology */
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
	//filename = "edgenew192x128.pgm";
	filename = "edgenew768x768.pgm";
	//filename = "edgenew512x384.pgm";

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
	//buf = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP, NP);

	new = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	old = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);
	edge = (RealNumber **)arralloc(sizeof(RealNumber), 2, MP + 2, NP + 2);

	domains = (int ***)arralloc(sizeof(int), 3, dims[0], dims[1], 2);
	disps = (int *)arralloc(sizeof(int), 1, size);
	counts = (int *)arralloc(sizeof(int), 1, size);

	

	//printf("rank %d is here\n",rank);

	if (rank == 0) {

		printf("Processing %d x %d image\n", M, N);
		printf("Number of iterations = %d\n", MAXITER);

		masterbuf = (RealNumber **)arralloc(sizeof(RealNumber), 2, M, N);

		printf("\nReading <%s>\n", filename);
		pgmread(filename, masterbuf, M, N);
		//printf("Done.\n");
	}

	/*	Configure the domain sizes, giving more work to the top and bottom processes
		since they have one less communication direction (non-periodic)
		but still ensuring that each row of processes has the same number of column pixels,
		and each column of processes has the same number of row pixels
	*/
	//printf("Rank %d Determining Domain Sizes\n", rank);

	int base_i = floor((RealNumber)M / (RealNumber)dims[0]);
	int base_j = floor((RealNumber)N / (RealNumber)dims[1]);
	int rem_i = M - (base_i*dims[0]);
	int rem_j = N - (base_j*dims[1]);

	for (i = 0; i < dims[0]; i++) {
		for (j = 0; j < dims[1]; j++) {
			domains[i][j][0] = base_i;
			domains[i][j][1] = base_j;
		}
	}
	
	// give the first extra row of pixels to the top set of processors
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

	if (rank == 0) {
		printf("Determined Domain Sizes\n");
		for (i = 0; i < dims[0]; i++) {
			for (j = 0; j < dims[1]; j++) {
				printf("(%d,%d) ", domains[i][j][0], domains[i][j][1]);
			}
			printf("\n");
		}
	}

	/* Stuff for sending same-size blocks */
	/*	We can send each processor the same (maximum) sized blocks, but they can
		overwrite bits that are not part of their domain */

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

	printf("About to scatter\n");
	//MPI_Scatterv(masterbuf, counts, disps, Small_send_section, buf, MP*NP, MPI_REALNUMBER, 0, comm);
	MPI_Scatterv(masterbuf, counts, disps, Small_send_section, &edge[1][1], 1, Recv_section, 0, comm);

	for (i=0; i<MP+2;i++) {
		for (j=0;j<NP+2;j++) {
			old[i][j]=255.0;
		}
	}

	int up, down, left, right;
	RealNumber val;

	MPI_Cart_shift(comm, 0, 1, &left, &right);
	MPI_Cart_shift(comm, 1, 1, &down, &up);

	if (rank == 0) {
		printf("Rank %d has neighbours:1n", rank);
		printf("up: %d\n", up);
		printf("down: %d\n", down);
		printf("left: %d\n", left);
		printf("right: %d\n", right);
	}

	if (up == MPI_PROC_NULL) {
		printf("Rank %d is doing top boundary conditions\n", rank);
		for (i = 1; i < MP; i++) {

			val = boundaryval(i, M);
			old[i][0] = (int)(255.0*val);
		}
	}

	if (down == MPI_PROC_NULL) {
		printf("Rank %d is doing bottom boundary conditions\n", rank);
		for (i = 1; i < MP; i++) {

			val = boundaryval(i, M);
			old[i][NP+1] = (int)(255.0*val);
		}
	}
	
	MPI_Datatype sides, top_bottom;
	MPI_Type_contiguous(MP, MPI_REALNUMBER, &sides);
	MPI_Type_vector(MP, 1, NP + 2, MPI_REALNUMBER, &top_bottom);
	
	//MPI_Request send_up, send_down, send_left, send_right;
	//MPI_Request recv_up, recv_down, recv_left, recv_right;
	MPI_Request requests[2*MAX_DIMS*MAX_DIMS];
	MPI_Status statuses[2*MAX_DIMS*MAX_DIMS];
	
	for (iter= 1;iter<=MAXITER; iter++) {
		if (iter%PRINTFREQ == 0 && rank==0) {
			printf("Iteration %d\n", iter);
		}
		
		MPI_Isend(&old[1][1], 1, sides, left, DEFAULT_TAG, comm, &requests[0]);// send_left);
		MPI_Isend(&old[MP][1], 1, sides, right, DEFAULT_TAG, comm, &requests[1]);//send_right);
		MPI_Isend(&old[1][1], 1, top_bottom, down, DEFAULT_TAG, comm, &requests[2]);//send_down);
		MPI_Isend(&old[1][NP], 1, top_bottom, up, DEFAULT_TAG, comm, &requests[3]);//send_up);

		MPI_Irecv(&old[0][1], 1, sides, left, DEFAULT_TAG, comm, &requests[4]);//recv_left);
		MPI_Irecv(&old[MP][1], 1, sides, right, DEFAULT_TAG, comm, &requests[5]);//recv_right);
		MPI_Irecv(&old[1][0], 1, top_bottom, down, DEFAULT_TAG, comm, &requests[6]);//recv_down);
		MPI_Irecv(&old[1][NP], 1, top_bottom, up, DEFAULT_TAG, comm, &requests[7]);//recv_up);
		

		for (i = 2; i < MP; i++) {
			for (j = 2; j < NP; j++) {
			new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
					- edge[i][j]);
			}
		}
		MPI_Waitall(2 * MAX_DIMS*MAX_DIMS, requests, statuses);
		
		i = 1;
		for (j = 1; j < NP + 1; j++) {
			new[i][j] = 0.25*(old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1]
				- edge[i][j]);
		}
		i = MP;
		for (j = 1; j < NP + 1; j++) {
			new[i][j] = 0.25*(old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1]
				- edge[i][j]);
		}
		j = 1;
		for (i = 1; i < MP + 1; i++) {
			new[i][j] = 0.25*(old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1]
				- edge[i][j]);
		}
		j = NP;
		for (i = 1; i < MP + 1; i++) {
			new[i][j] = 0.25*(old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1]
				- edge[i][j]);
		}
		/*
		// Calculate the sides
		for (j = 1; j < NP + 1; j++) {
			new[1][j] = 0.25*(old[0][j] + old[2][j] + old[1][j - 1] + old[1][j + 1]
				- edge[1][j]);
			new[MP][j] = 0.25*(old[MP - 1][j] + old[MP + 1][j] + old[MP][j - 1] + old[MP][j + 1]
				- edge[MP][j]);
		}
		// Calculate the top and bottom
		for (i = 1; i < MP; i++) {
			new[i][1] = 0.25*(old[i - 1][1] + old[i + 1][1] + old[i][0] + old[i][2]
				- edge[i][1]);
			new[i][NP] = 0.25*(old[i - 1][NP] + old[i + 1][NP] + old[i][NP - 1] + old[i][NP + 1]
				- edge[i][NP]);
		}
		*/

		for (i = 1; i<MP+1; i++) {
			for (j = 1; j<NP+1; j++) {
				old[i][j] = new[i][j];
			}
		}

	}
	
	//printf("\nFinished %d iterations\n", iter-1);

	/* Gather the data back to process 0 */
	//MPI_Gatherv(buf, MP*NP, MPI_REALNUMBER, masterbuf, counts, disps, Small_send_section, 0, comm);
	MPI_Gatherv(&old[1][1], 1, Recv_section, masterbuf, counts, disps, Small_send_section, 0, comm);

	if (rank == 0) {

		//filename = "imagenew192x128.pgm";
		filename = "imagenew768x768.pgm";
		//filename = "imagenew512x384.pgm";
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