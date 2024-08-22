//Parallel Red Black Method with Multigrid using MPI
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
    int N1 = 33, N2 = 17;
    int pre = 5;

    int N1_local = N1/size; int N2_local = N2/size;
    if(my_id == size-1){
        N1_local += N1%size;
        N2_local += N2%size;
    }

    double delta1 = (double)(xmax - xmin)/(N1-1); 
    double delta2 = (double)(xmax - xmin)/(N2-1);

    //Memory Allocation
    double **phi, **phi_exact;
    double **res1, **res2;
    double **delx1, **delx2;
    double **q;

    phi = allocate2DArray(N1_local, N1);
    phi_exact = allocate2DArray(N1_local, N1);
    res1 = allocate2DArray(N1_local, N1);
    res2 = allocate2DArray(N2_local, N2);
    delx1 = allocate2DArray(N1_local, N1);
    delx2 = allocate2DArray(N2_local, N2);
    q = allocate2DArray(N1_local, N1);

    //Row decomposition of fine mesh
    double ystart = ymin + delta1*my_id*(N1/size);

    //Initial Guesses and Boundary Conditions
    double x,y;
    for(i = 0; i<N1_local ;i++){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            y = ystart + delta1*i;
            phi[i][j] = 0;
            if (j == 0)
                phi[i][j] = exp(-2*y);
            if (j == N1-1)
                phi[i][j] = exp(1-2*y);

            res1[i][j] = 0;
            delx1[i][j] = 0;
        }
    }

    if(my_id == 0){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            phi[0][j] = exp(x);
        }          
    }

    if(my_id == size-1){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            phi[N1_local-1][j] = exp(x-2);
        }          
    }

    for(i = 0; i<N2_local ;i++){
        for(j = 0; j<N2; j++){
            res2[i][j] = 0;
            delx2[i][j] = 0;
        }
    }

    //RHS of Poisson equation
    for(i = 0; i<N1_local ;i++){
        for(j = 0; j<N1; j++){
            x = xmin + j*delta1; 
            y = ystart + i*delta1;
            q[i][j] = 5*exp(x)*exp(-2*y);
       }
    }
    //Exact Solution
    for(i = 0; i<N1_local; i++){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            y = ystart + delta1*i;  
            phi_exact[i][j] = exp(x)*exp(-2*y);
        }
    }

    //Norm of exact solution
    double norm_exact = 0;
    for(i = 0; i<N1_local; i++){
        for(j = 0; j<N1; j++){
            norm_exact += phi_exact[i][j]*phi_exact[i][j];
        }
    }
    double norm_exact_global;    
    MPI_Allreduce(&norm_exact,&norm_exact_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    norm_exact_global = sqrt(norm_exact_global);

    int iter = 0;
    double err = 1;

    //Buffers for send receive
    double phi_top[N1]; double phi_bot[N1];
    double res1_top[N1]; double res1_bot[N1];
    double delx2_top[N2]; double delx2_bot[N2];

    //Main loop
    while(err > 0.0001){

        //Pre Smoothing
        for(k = 0; k<pre; k++){

            //Odd points
            //Sending the ghost nodes required at the top from i-1th processor to ith processor
            if(my_id != size-1)
                MPI_Send(phi[N1_local-1],N1,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
            if(my_id != 0){
                MPI_Recv(&phi_top,N1,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + j) % 2 == 1){
                        phi[0][j] = 0.25*(phi_top[j] + phi[1][j] + phi[0][j+1] + phi[0][j-1] - delta1*delta1*q[0][j]);
                    }
                }
            }

            //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
            if(my_id !=0)
                MPI_Send(phi[0],N1,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
            if(my_id != size-1){
                MPI_Recv(&phi_bot,N1,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + N1_local - 1 + j) % 2 == 1){
                        phi[N1_local-1][j] = 0.25*(phi[N1_local-2][j] + phi_bot[j] + phi[N1_local-1][j+1] + phi[N1_local-1][j-1] - delta1*delta1*q[N1_local-1][j]);
                    }
                }
            }

            //Updating interior phi values
            for(i = 1; i<N1_local-1; i++){
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + i + j) % 2 == 1){
                        phi[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - delta1*delta1*q[i][j]);  
                    }
                }
            }

            //Even Points
            //Sending the ghost nodes required at the top from i-1th processor to ith processor
            if(my_id != size-1)
                MPI_Send(phi[N1_local-1],N1,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
            if(my_id !=0){
                MPI_Recv(&phi_top,N1,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + j) % 2 == 0 ){
                        phi[0][j] = 0.25*(phi_top[j] + phi[1][j] + phi[0][j+1] + phi[0][j-1] - delta1*delta1*q[0][j]);
                    }
                }
            }

            //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
            if(my_id !=0)
                MPI_Send(phi[0],N1,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
            if(my_id != size-1){
                MPI_Recv(&phi_bot,N1,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + N1_local - 1 +j) % 2 == 0){
                        phi[N1_local-1][j] = 0.25*(phi[N1_local-2][j] + phi_bot[j] + phi[N1_local-1][j+1] + phi[N1_local-1][j-1] - delta1*delta1*q[N1_local-1][j]);
                    }
                }
            }

            //Updating interior phi values
            for(i = 1; i<N1_local-1; i++){
                for(j = 1; j<N1-1; j++){
                    if( ((N1/size)*my_id + i + j) % 2 == 0){
                        phi[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - delta1*delta1*q[i][j]);  
                    }
                }
            }
        }


        //Residual 
        //Sending the ghost nodes required at the top from i-1th processor to ith processor
        if(my_id != size-1)
            MPI_Send(phi[N1_local-1],N1,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
        if(my_id !=0){
            MPI_Recv(&phi_top,N1,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
            for(j = 1; j<N1-1; j++){
                res1[0][j] =  q[0][j] - (phi[1][j] + phi_top[j] + phi[0][j+1] + phi[0][j-1] - 4*phi[0][j])/(delta1*delta1);
            }
        }

        //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
        if(my_id !=0)
            MPI_Send(phi[0],N1,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
        if(my_id != size-1){
            MPI_Recv(&phi_bot,N1,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
            for(j = 1; j<N1-1; j++){
                res1[N1_local-1][j] =  q[N1_local-1][j] - (phi_bot[j] + phi[N1_local-2][j] + phi[N1_local-1][j+1] + phi[N1_local-1][j-1] - 4*phi[N1_local-1][j])/(delta1*delta1);         
            }
        }

        for(i = 1; i<N1_local-1 ;i++){
            for(j = 1; j<N1-1; j++){
                res1[i][j] =  q[i][j] - (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - 4*phi[i][j])/(delta1*delta1);
            }
        }

        //Restriction
        //Sending the ghost nodes required at the top from i-1th processor to ith processor
        if(my_id != size-1)
            MPI_Send(res1[N1_local-1],N1,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
        if(my_id != 0){
            MPI_Recv(&res1_top,N1,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
            for(j = 1; j<N2-1; j++){
                res2[0][j] = 0.25*res1[0][2*j] + 0.125*(res1[1][2*j] + res1_top[2*j] + res1[0][2*j+1] + res1[0][2*j-1]) + 0.0625*(res1[1][2*j+1] + res1_top[2*j-1] + res1[1][2*j-1] + res1_top[2*j+1]);
            }
        }

        //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
        if(my_id != size-1){
            for(j = 1; j<N2-1; j++){
                k = N2_local-1;
                res2[k][j] = 0.25*res1[N1_local-2][2*j] + 0.125*(res1[N1_local-1][2*j] + res1[N1_local-3][2*j] + res1[N1_local-2][2*j+1] + res1[N1_local-2][2*j-1]) + 0.0625*(res1[N1_local-1][2*j+1] + res1[N1_local-3][2*j-1] + res1[N1_local-1][2*j-1] + res1[N1_local-3][2*j+1]);
            }
        }

        for(i = 1; i<N2_local-1 ;i++){
            for(j = 1; j<N2-1; j++){
                res2[i][j] = 0.25*res1[2*i][2*j] + 0.125*(res1[2*i+1][2*j] + res1[2*i-1][2*j] + res1[2*i][2*j+1] + res1[2*i][2*j-1]) + 0.0625*(res1[2*i+1][2*j+1] + res1[2*i-1][2*j-1] + res1[2*i+1][2*j-1] + res1[2*i-1][2*j+1]);
            }
        }   

        //Smoothing
        for(k = 0; k<pre; k++){

            //Odd Points
            //Sending the ghost nodes required at the top from i-1th processor to ith processor
            if(my_id != size-1)
                MPI_Send(delx2[N2_local-1],N2,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
            if(my_id != 0){
                MPI_Recv(&delx2_top,N2,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + j) % 2 == 1){
                        delx2[0][j] = 0.25*(delx2_top[j] + delx2[1][j] + delx2[0][j+1] + delx2[0][j-1] - delta2*delta2*res2[0][j]);
                    }
                }
            }

            //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
            if(my_id !=0)
                MPI_Send(delx2[0],N2,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
            if(my_id != size-1){
                MPI_Recv(&delx2_bot,N2,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + N2_local - 1 + j) % 2 == 1){
                        delx2[N2_local-1][j] = 0.25*(delx2[N2_local-2][j] + delx2_bot[j] + delx2[N2_local-1][j+1] + delx2[N2_local-1][j-1] - delta2*delta2*res2[N2_local-1][j]);
                    }
                }
            }

            //Updating interior phi values
            for(i = 1; i<N2_local-1; i++){
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + i + j) % 2 == 1){
                        delx2[i][j] = 0.25*(delx2[i+1][j] + delx2[i-1][j] + delx2[i][j+1] + delx2[i][j-1] - delta2*delta2*res2[i][j]);  
                    }
                }
            }

            //Even Points
            //Sending the ghost nodes required at the top from i-1th processor to ith processor
            if(my_id != size-1)
                MPI_Send(delx2[N2_local-1],N2,MPI_DOUBLE,my_id+1,tag1,MPI_COMM_WORLD);
            if(my_id != 0){
                MPI_Recv(&delx2_top,N2,MPI_DOUBLE,my_id-1,tag1,MPI_COMM_WORLD,&status);
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + j) % 2 == 0){
                        delx2[0][j] = 0.25*(delx2_top[j] + delx2[1][j] + delx2[0][j+1] + delx2[0][j-1] - delta2*delta2*res2[0][j]);
                    }
                }
            }

            //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
            if(my_id !=0)
                MPI_Send(delx2[0],N2,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
            if(my_id != size-1){
                MPI_Recv(&delx2_bot,N2,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + N2_local - 1 + j) % 2 == 0){
                        delx2[N2_local-1][j] = 0.25*(delx2[N2_local-2][j] + delx2_bot[j] + delx2[N2_local-1][j+1] + delx2[N2_local-1][j-1] - delta2*delta2*res2[N2_local-1][j]);
                    }
                }
            }

            //Updating interior phi values
            for(i = 1; i<N2_local-1; i++){
                for(j = 1; j<N2-1; j++){
                    if( ((N2/size)*my_id + i + j) % 2 == 0){
                        delx2[i][j] = 0.25*(delx2[i+1][j] + delx2[i-1][j] + delx2[i][j+1] + delx2[i][j-1] - delta2*delta2*res2[i][j]);  
                    }
                }
            }
        }


        //Interpolation
        //Sending the ghost nodes required at the bottom from i+1th processor to ith processor
        if(my_id !=0)
            MPI_Send(delx2[0],N2,MPI_DOUBLE,my_id-1,tag2,MPI_COMM_WORLD);
        if(my_id != size-1){
            MPI_Recv(&delx2_bot,N2,MPI_DOUBLE,my_id+1,tag2,MPI_COMM_WORLD,&status);
            for(j = 1; j<N2-1; j++){
                    delx1[N1_local-1][2*j] = 0.5*(delx2[N2_local-1][j] + delx2_bot[j]);
                    delx1[N1_local-1][2*j+1] = 0.25*(delx2[N2_local-1][j] + delx2_bot[j+1] + delx2_bot[j] + delx2[N2_local-1][j+1]);
                    delx1[N1_local-1][2*j-1] = 0.25*(delx2[N2_local-1][j] + delx2_bot[j-1] + delx2_bot[j] + delx2[N2_local-1][j-1]);
            }
        }

        for(i = 0; i<N2_local ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i][2*j] = delx2[i][j]; 
            }
        }

        //Vertical Sweep
        for(i = 0; i<N2_local-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i+1][2*j] = 0.5*(delx2[i][j] + delx2[i+1][j]);
            }
        }

        //Horizontal Sweep
        for(i = 0; i<N2_local ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i][2*j+1] = 0.5*(delx2[i][j] + delx2[i][j+1]);
                if(j == 1)
                    delx1[2*i][2*j-1] = 0.5*(delx2[i][j] + delx2[i][j-1]);
            }
        }

        //4 point interpolation
        for(i = 0; i<N2_local-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i+1][2*j+1] = 0.25*(delx2[i][j] + delx2[i+1][j+1] + delx2[i+1][j] + delx2[i][j+1]);
                delx1[2*i+1][2*j-1] = 0.25*(delx2[i][j] + delx2[i+1][j-1] + delx2[i+1][j] + delx2[i][j-1]);
            }
        }

        //Adding the correction
        for(i = 0; i<N1_local ;i++){
            for(j = 0; j<N1; j++){
                phi[i][j] += (delx1[i][j]);
            }
        }

        //Error Calculation
        double err_local = 0;
        err = 0;
        for (i = 0; i<N1_local; i++){
            for(j = 0; j<N1; j++){
                err_local += (phi_exact[i][j] - phi[i][j])*(phi_exact[i][j] - phi[i][j]);
            }
        }
             
        MPI_Allreduce(&err_local,&err,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        err = sqrt(err)/norm_exact_global;

        iter += 1;
    }

    t2 = MPI_Wtime();

    if(my_id ==0){
        printf("Number of iterations for Parallel Red Black Method with Multigrid using MPI are : %d\n",2*pre*iter);
        printf("The problem size is : %d x %d\n",N1,N1);
        printf("The error is : %lf\n",err);
        printf("The time taken with %d processors is : %lf\n",size,t2-t1);
    }

     //Gathering the results using MPI_Gatherv as each processor has different number of elements
    int *recvcounts = (int*)malloc(size*sizeof(int));
    for(i = 0; i<size; i++){
        recvcounts[i] = N1*(N1/size);}
    recvcounts[size-1] += N1*(N1%size);

    int *displs = (int*)malloc(size * sizeof(int));
    int disp = 0;
    for (int i = 0; i < size; i++) {
        displs[i] = disp;
        disp += recvcounts[i];
    }

    //Allocating memory for send receive buffers
    double *phi_send = (double*)malloc(N1*N1_local*sizeof(double));
    double *phi_global_recv = (double*)malloc(N1*N1*sizeof(double));
    double **phi_global = (double**)malloc(N1*sizeof(double*));
    for(i = 0; i<N1; i++){
        phi_global[i] = (double*)malloc(N1*sizeof(double));
    }

    //Converting matrices to arrays before sending
    for(i = 0; i<N1_local; i++){
        for(j = 0; j<N1; j++){
            phi_send[i*N1 + j] = phi[i][j];
        }
    }

    //Collecting all arrays on processor 0
    MPI_Gatherv(phi_send,N1*N1_local,MPI_DOUBLE,phi_global_recv,recvcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(my_id == 0){

        //Converting the received array into a matrix
        for(i = 0; i<N1; i++){
            for(j = 0; j<N1; j++){
                phi_global[i][j] = phi_global_recv[i*N1 + j];
            }
        }
    }

    //Freeing the allocated memory
    deallocate2DArray(phi, N1_local);
    deallocate2DArray(phi_exact, N1_local);
    deallocate2DArray(res1, N1_local);
    deallocate2DArray(res2, N2_local);
    deallocate2DArray(delx1, N1_local);
    deallocate2DArray(delx2, N2_local);
    deallocate2DArray(q, N1_local);

    MPI_Finalize();

    return 0;
}