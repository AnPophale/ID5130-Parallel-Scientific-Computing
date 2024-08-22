//Parallel Jacobi Method using OpenMP
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

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
    double t1,t2;
    t1 = omp_get_wtime();

    //Mesh Parameters
    int xmin = 0, xmax = 1, ymin = 0, ymax = 1, thread_count;
    int N = 33;

    //No. of threads
    thread_count = 2;

    double delta = (double)(xmax - xmin)/(N-1); 

    //Memory Allocation
    double **phi, **phi_exact, **phi_old;
    double **q;

    phi = allocate2DArray(N, N);
    phi_old = allocate2DArray(N, N);
    phi_exact = allocate2DArray(N, N);
    q = allocate2DArray(N, N);

    //Norm of exact solution
    double norm_exact = 0;

    //Initial Guesses and Boundary Conditions
    int i,j,k; double x,y;

    #pragma omp parallel num_threads(thread_count) default(none) private(i,j,x,y) shared(N,xmin,delta,ymin,phi,q,phi_exact,norm_exact)
    {   
        #pragma omp for collapse(2)
        for(i = 0; i<N ;i++){
        for(j = 0; j<N; j++){
            x = xmin + delta*j;
            y = ymin + delta*i;
            phi[i][j] = 0;
            if (i == 0)
                phi[i][j] = exp(x);
            if (i == N-1)
                phi[i][j] = exp(x-2);
            if (j == 0)
                phi[i][j] = exp(-2*y);
            if (j == N-1)
                phi[i][j] = exp(1-2*y);
        }
    }

    //RHS of Poisson equation
    #pragma omp for collapse(2)
    for(i = 1; i<N-1 ;i++){
        for(j = 1; j<N-1; j++){
            x = xmin + j*delta; 
            y = ymin + i*delta;
            q[i][j] = 5*exp(x)*exp(-2*y);
       }
    }

    //Exact Solution
    #pragma omp for collapse(2)
    for(i = 0; i<N; i++){
        for(j = 0; j<N; j++){
            x = xmin + delta*j;
            y = ymin + delta*i;  
            phi_exact[i][j] = exp(x)*exp(-2*y);
        }
    }

    //Norm of exact solution
    #pragma omp for collapse(2) reduction(+:norm_exact)
    for(i = 0; i<N; i++){
        for(j = 0; j<N; j++){
            norm_exact += phi_exact[i][j]*phi_exact[i][j];
            }
        }

    }
    
    norm_exact = sqrt(norm_exact);

    int iter = 0;
    double err = 1;

    //Main loop
    while(err>0.0001){
            
            err = 0;

            #pragma omp parallel num_threads(thread_count) default(none) private(i,j) shared(N,phi_old,phi,phi_exact,delta,err,q)
            {
            //Storing old values of phi
            #pragma omp for collapse(2)
            for(i = 0; i<N; i++){
                for(j = 0; j<N; j++){
                    phi_old[i][j] = phi[i][j];
                }
            }
            
            #pragma omp for collapse(2)
            for(i = 1; i<N-1; i++){
                for(j = 1; j<N-1; j++){
                    phi[i][j] = 0.25*(phi_old[i][j+1] + phi_old[i][j-1] + phi_old[i+1][j] + phi_old[i-1][j] - delta*delta*q[i][j]);
                }
            }
            

        //Error Calculation
        #pragma omp for collapse(2) reduction(+:err)
        for (i = 0; i<N; i++){
            for(j = 0; j<N; j++){
                err += (phi_exact[i][j] - phi[i][j])*(phi_exact[i][j] - phi[i][j]);
            }
        }

        }
            
        err = sqrt(err)/norm_exact;
        iter += 1;
    }

    t2 = omp_get_wtime();

    printf("Number of iterations for Parallel Jacobi Method using OpenMP are : %d\n",iter);
    printf("The problem size is : %d x %d\n",N,N);
    printf("The error is : %lf\n",err);
    printf("The time taken with %d threads is : %lf\n",thread_count,t2-t1);

    //Freeing the allocated memory
    deallocate2DArray(phi, N);
    deallocate2DArray(phi_exact, N);
    deallocate2DArray(phi_old,N);
    deallocate2DArray(q, N);


    return 0;
}