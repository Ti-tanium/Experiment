#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define MAX_PROCESSOR_NUM 12
void MatrixMultiply(int n, double *a, double *b, double *c);
main (int argc, char *argv[])
{
    int i,j,k,m,p;
    int n, nlocal;
    double *a, *b, *c;
    int npes, dims[2], periods[2];
	int myrank, my2drank, mycoords[2];
	int shiftsource, shiftdest, rightrank;
    sint leftrank, downrank, uprank;
    MPI_Status status;
    MPI_Comm comm_2d;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if (argc !=2){
        if (myrank ==0)
            printf("Usage: %s <the dimension of  thematrix>\n", argv[0]);                
    	    MPI_Finalize();
         exit(0);
    }
    n=atoi(argv[1]);
    dims[0] = sqrt(npes);
    dims[1] = npes/dims[0];
    if (dims[0] != dims[1]){
        if (myrank ==0 )
            printf("The number of processes must be a perfect square.\n");
        MPI_Finalize();
        exit(0);
    }
    periods[0]=periods[1]=1; 
    MPI_Cart_create(MPI_COMM_WORLD,2, dims,periods, 0, &comm_2d);
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d,my2drank,2, mycoords);
    nlocal = n/dims[0];
    a=(double*)malloc(nlocal*nlocal* sizeof (double));
    b= (double *)malloc(nlocal*nlocal*sizeof (double));
    c= (double *)malloc(nlocal*nlocal*sizeof (double));
    srand48((long)myrank);
    for ( i=0; i<nlocal*nlocal; i++)  {
        a[i] = b[i] =(double)i; 
        c[i] = 0.0;
    }
    MPI_Cart_shift(comm_2d, 0, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a,nlocal*nlocal, MPI_DOUBLE, shiftdest,1,shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 1, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_DOUBLE,shiftdest, 1, shiftsource, 1, comm_2d, &status);    
    MPI_Cart_shift(comm_2d, 0,-1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 1, -1, &downrank, &uprank);
    for (i=0; i<dims[0]; i++) {
      	MatrixMultiply(nlocal, a, b, c);
 		MPI_Sendrecv_replace(a, nlocal*nlocal, 
        MPI_DOUBLE,leftrank,1,rightrank,  1, comm_2d, &status);
        MPI_Sendrecv_replace(b, nlocal*nlocal,MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);
    }
    MPI_Cart_shift(comm_2d, 0, +mycoords[0], &shiftsource,   &shiftdest);
    MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 1, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource,  1, comm_2d, &status);
	MPI_Cart_shift(comm_2d, 1, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal*nlocal, 	MPI_DOUBLE, shiftdest, 1, shiftsource, 1,comm_2d, &status);
    MPI_Comm_free(&comm_2d); 
    if (myrank ==0){
        puts("Random Matrix A");
        for(i = 0; i < nlocal; i ++){
            for(j = 0; j < nlocal; j ++)
 		        printf("%9.7f ", a[i*nlocal+j]);
        	printf("\n");
    	} 	
        puts("Random Matrix B");
        for(i = 0; i < nlocal; i ++){
            for(j = 0; j < nlocal; j ++)
                printf("%9.7f ", b[i*nlocal+j]);
            printf("\n");
        }
        puts("Matrix C = A*B");

        for(i = 0; i < nlocal; i ++){
            for(j = 0; j < nlocal; j ++) 
                printf("%9.7f ", c[i*nlocal+j]);
            printf("\n");
        }
        free(a);
	    free(b);
        free(c);
        MPI_Finalize();
    }
}
void MatrixMultiply(int n, double *a, double *b, double *c) {
   	int i, j, k;
   	for (i = 0; i < n; i++)
     	for (j = 0; j < n; j++)
       		for (k = 0; k < n; k++)
           		c[i*n+j]+=a[i*n+k]*b[k*n+j];
}

