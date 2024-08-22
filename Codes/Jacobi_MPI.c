//Parallel Jacobi Method using MPI 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

//Memory allocation for matrices
double** allocate2DArray(int rows, int cols) {
    double **arr = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        arr[i] = (double *)malloc(cols * sizeof(double));
    }
    return arr;
}

//Memory deallocation for matrices
void deallocate2DArray(double **arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}

int main(){
    MPI_Init(NULL,NULL);
    double t1,t2;
    t1 = MPI_Wtime();
    int i,j,k,my_id,size,tag1 = 0, tag2 = 1;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);

    //Mesh Parameters
    int xmin = 0, xmax = 1, ymin = 0, ymax = 1;
    int N = 33;

    int N_local = N/size;
    if(my_id == size-1){
        N_local += N%size;
    }

    double delta1 = (double)(xmax - xmin)/(N-1); 

    //Memory Allocation
    double **phi, **phi_exact, **phi_old;
    double **q;

    phi = allocate2DArray(N_local, N);
    phi_old = allocate2DArray(N_local, N);
    phi_exact = allocate2DArray(N_local, N);
    q = allocate2DArray(N_local, N);

    double ystart = ymin + delta1*my_id*(N/size);

    //Initial Guesses and Boundary Conditions
    double x,y;
    for(i = 0; i<N_local ;i++){
        for(j = 0; j<N; j++){
            x = xmin + delta1*j;
            y = ystart + delta1*i;
            phi[i][j] = 0;
            if (j == 0)
                phi[i][j] = exp(-2*y);
            if (j == N-1)
                phi[i][j] = exp(1-2*y);

        }
    }

    if(my_id == 0){
        for(j = 0; j<N; j++){
            x = xmin + delta1*j;
            phi[0][j] = exp(x);
        }          
    }

    if(my_id == size-1){
        for(j = 0; j<N; j++){
            x = xmin + delta1*j;
            phi[N_local-1][j] = exp(x-2);
        }          
    }

    //RHS of Poisson equation
    for(i = 0; i<N_local ;i++){
        for(j = 0; j<N; j++){
            x = xmin + j*delta1; 
            y = ystart + i*delta1;
            q[i][j] = 5*exp(x)*exp(-2*y);
       }
    }

    //Exact Solution
    for(i = 0; i<N_local; i++){
        for(j = 0; j<N; j++){
            x = xmin + delta1*j;
            y = ystart + delta1*i;  
            phi_exact[i][j] = exp(x)*exp(-2*y);
        }
    }

    //Norm of exact solution
    double norm_exact = 0;
    for(i = 0; i<N_local; i++){
        for(j = 0; j<N; j++){
            norm_exact += phi_exact[i][j]*phi_exact[i][j];
        }
    }
    double norm_exact_global;    
    MPI_Allreduce(&norm_exact,&norm_exact_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    norm_exact_global = sqrt(norm_exact_global);

    int iter = 0;
    double err = 1;

    //Buffers for send receive
    double phi_top[N];
    double phi_bot[N];
    
    //Main loop
    while(err > 0.0001){

        //Storing old values of phi
        for(i = 0; i<N_local; i++){
            for(j = 0; j<N; j++){
                phi_old[i][j] = phi[i][j];
            }
        }

        //Sending the ghost nodes required at the top from i-1th processor to ith processor
        if(my_id != size-1)
            MPI_Send(phi_old[N_local-1],N,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
        if(my_id != 0){
            MPI_Recv(&phi_top,N,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
            for(j = 1; j<N-1; j++){
                phi[0][j] = 0.25*(phi_top[j] + phi_old[1][j] + phi_old[0][j+1] + phi_old[0][j-1] - delta1*delta1*q[0][j]);
            }
        }

        //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
        if(my_id !=0)
            MPI_Send(phi_old[0],N,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
        if(my_id != size-1){
            MPI_Recv(&phi_bot,N,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
            for(j = 1; j<N-1; j++){
                    phi[N_local-1][j] = 0.25*(phi_old[N_local-2][j] + phi_bot[j] + phi_old[N_local-1][j+1] + phi_old[N_local-1][j-1] - delta1*delta1*q[N_local-1][j]);
            }
        }

        //Updating interior phi values
        for(i = 1; i<N_local-1; i++){
            for(j = 1; j<N-1; j++){
                    phi[i][j] = 0.25*(phi_old[i+1][j] + phi_old[i-1][j] + phi_old[i][j+1] + phi_old[i][j-1] - delta1*delta1*q[i][j]); 
            }
        }
    
        //Error Calculation
        double err_local = 0;
        err = 0;
        for (i = 0; i<N_local; i++){
            for(j = 0; j<N; j++){
                err_local += (phi_exact[i][j] - phi[i][j])*(phi_exact[i][j] - phi[i][j]);
            }
        }
             
        MPI_Allreduce(&err_local,&err,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        err = sqrt(err)/norm_exact_global;
        iter += 1;
    }

    t2 = MPI_Wtime();

    if(my_id ==0){
        printf("Number of iterations for Parallel Jacobi Method using MPI are : %d\n",iter);
        printf("The problem size is : %d x %d\n",N,N);
        printf("The error is : %lf\n",err);
        printf("The time taken with %d processors is : %lf\n",size,t2-t1);
    }

     //Gathering the results using MPI_Gatherv as each processor has different number of elements
    int *recvcounts = (int*)malloc(size*sizeof(int));
    for(i = 0; i<size; i++){
        recvcounts[i] = N*(N/size);}
    recvcounts[size-1] += N*(N%size);

    int *displs = (int*)malloc(size * sizeof(int));
    int disp = 0;
    for (int i = 0; i < size; i++) {
        displs[i] = disp;
        disp += recvcounts[i];
    }

    //Allocating memory for send receive buffers
    double *phi_send = (double*)malloc(N*N_local*sizeof(double));
    double *phi_global_recv = (double*)malloc(N*N*sizeof(double));
    double **phi_global = (double**)malloc(N*sizeof(double*));
    for(i = 0; i<N; i++){
        phi_global[i] = (double*)malloc(N*sizeof(double));
    }

    //Converting matrices to arrays before sending
    for(i = 0; i<N_local; i++){
        for(j = 0; j<N; j++){
            phi_send[i*N + j] = phi[i][j];
        }
    }

    //Collecting all arrays on processor 0
    MPI_Gatherv(phi_send,N*N_local,MPI_DOUBLE,phi_global_recv,recvcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(my_id == 0){
        //Converting the received array into a matrix
        for(i = 0; i<N; i++){
            for(j = 0; j<N; j++){
                phi_global[i][j] = phi_global_recv[i*N + j];
            }
        }
    }

    //Freeing the allocated memory
    deallocate2DArray(phi, N_local);
    deallocate2DArray(phi_old, N_local);
    deallocate2DArray(phi_exact, N_local);
    deallocate2DArray(q, N_local);

    MPI_Finalize();

    return 0;
}