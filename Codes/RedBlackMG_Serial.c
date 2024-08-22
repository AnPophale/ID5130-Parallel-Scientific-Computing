//Serial Red Black Method with Multigrid
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
    int xmin = 0, xmax = 1, ymin = 0, ymax = 1;
    int N1 = 33, N2 = 17;
    int pre = 5;

    double delta1 = (double)(xmax - xmin)/(N1-1); 
    double delta2 = (double)(xmax - xmin)/(N2-1);

    //Memory Allocation
    double **phi, **phi_exact;
    double **res1, **res2;
    double **delx1, **delx2;
    double **q;

    phi = allocate2DArray(N1, N1);
    phi_exact = allocate2DArray(N1, N1);
    res1 = allocate2DArray(N1, N1);
    res2 = allocate2DArray(N2, N2);
    delx1 = allocate2DArray(N1, N1);
    delx2 = allocate2DArray(N2, N2);
    q = allocate2DArray(N1, N1);

    //Initial Guesses and Boundary Conditions
    int i,j,k; double x,y;
    for(i = 0; i<N1 ;i++){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            y = ymin + delta1*i;
            phi[i][j] = 0;
            if (i == 0)
                phi[i][j] = exp(x);
            if (i == N1-1)
                phi[i][j] = exp(x-2);
            if (j == 0)
                phi[i][j] = exp(-2*y);
            if (j == N1-1)
                phi[i][j] = exp(1-2*y);

            res1[i][j] = 0;
            delx1[i][j] = 0;
        }
    }
    
    for(i = 0; i<N2 ;i++){
        for(j = 0; j<N2; j++){
            res2[i][j] = 0;
            delx2[i][j] = 0;
        }
    }

    //RHS of Poisson equation
    for(i = 1; i<N1-1 ;i++){
        for(j = 1; j<N1-1; j++){
            x = xmin + j*delta1; 
            y = ymin + i*delta1;
            q[i][j] = 5*exp(x)*exp(-2*y);
       }
    }

    //Exact Solution
    for(i = 0; i<N1; i++){
        for(j = 0; j<N1; j++){
            x = xmin + delta1*j;
            y = ymin + delta1*i;  
            phi_exact[i][j] = exp(x)*exp(-2*y);
        }
    }

    //Norm of exact solution
    double norm_exact = 0;
    for(i = 0; i<N1; i++){
        for(j = 0; j<N1; j++){
            norm_exact += phi_exact[i][j]*phi_exact[i][j];
            }
        }
    norm_exact = sqrt(norm_exact);

    int iter = 0;
    double err = 1;

    //Main loop
    while(err>0.0001){

        //Pre Smoothing
        for(k = 0; k<pre; k++){

            //Updating odd elements
            for(i = 1; i<N1-1; i++){
                for(j = 1; j<N1-1; j++){
                    if((i+j)%2 == 1)
                        phi[i][j] = 0.25*(phi[i][j+1] + phi[i][j-1] + phi[i+1][j] + phi[i-1][j] - delta1*delta1*q[i][j]);
                }
            }
            
            //Updating even elements
            for(i=1; i<N1-1; i++){
                for(j=1; j<N1-1; j++){
                    if((i+j)%2 == 0)
                        phi[i][j] = 0.25*(phi[i][j+1] + phi[i][j-1] + phi[i+1][j] + phi[i-1][j] - delta1*delta1*q[i][j]);
                }
            }
        }

        //Residual 
        for(i = 1; i<N1-1 ;i++){
            for(j = 1; j<N1-1; j++){
                res1[i][j] =  q[i][j] - (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - 4*phi[i][j])/(delta1*delta1);
            }
        }

        //Restricton
        for(i = 1; i<N2-1 ;i++){
            for(j = 1; j<N2-1; j++){
                res2[i][j] = 0.25*res1[2*i][2*j] + 0.125*(res1[2*i+1][2*j] + res1[2*i-1][2*j] + res1[2*i][2*j+1] + res1[2*i][2*j-1]) + 0.0625*(res1[2*i+1][2*j+1] + res1[2*i-1][2*j-1] + res1[2*i+1][2*j-1] + res1[2*i-1][2*j+1]);
            }
        }   
    
        //Smoothing
        for(k = 0; k<pre; k++){

            //Updating odd elements
            for(i = 1; i<N2-1; i++){
                for(j = 1; j<N2-1; j++){
                    if((i+j)%2 == 1)
                        delx2[i][j] = 0.25*(delx2[i][j+1] + delx2[i][j-1] + delx2[i+1][j] + delx2[i-1][j] - delta2*delta2*res2[i][j]);
                }
            }
            
            //Updating even elements
            for(i = 1; i<N2-1; i++){
                for(j = 1; j<N2-1; j++){
                    if((i+j)%2 ==0)
                        delx2[i][j] = 0.25*(delx2[i][j+1] + delx2[i][j-1] + delx2[i+1][j] + delx2[i-1][j] - delta2*delta2*res2[i][j]);
                }
            }
        }

        //Interpolation
        for(i = 1; i<N2-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i][2*j] = delx2[i][j]; 
            }
        }

        //Vertical Sweep
        for(i = 1; i<N2-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i+1][2*j] = 0.5*(delx2[i][j] + delx2[i+1][j]);
                if(i == 1)
                    delx1[2*i-1][2*j] = 0.5*(delx2[i][j] + delx2[i-1][j]);
            }
        }

        //Horizontal Sweep
        for(i = 1; i<N2-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i][2*j+1] = 0.5*(delx2[i][j] + delx2[i][j+1]);
                if(j == 1)
                    delx1[2*i][2*j-1] = 0.5*(delx2[i][j] + delx2[i][j-1]);
            }
        }

        //4 point interpolation
        for(i = 1; i<N2-1 ;i++){
            for(j = 1; j<N2-1; j++){
                delx1[2*i+1][2*j+1] = 0.25*(delx2[i][j] + delx2[i+1][j+1] + delx2[i+1][j] + delx2[i][j+1]);
                delx1[2*i+1][2*j-1] = 0.25*(delx2[i][j] + delx2[i+1][j-1] + delx2[i+1][j] + delx2[i][j-1]);
                delx1[2*i-1][2*j-1] = 0.25*(delx2[i][j] + delx2[i-1][j-1] + delx2[i-1][j] + delx2[i][j-1]);
                delx1[2*i-1][2*j+1] = 0.25*(delx2[i][j] + delx2[i-1][j+1] + delx2[i-1][j] + delx2[i][j+1]);
            }
        }

        //Adding the correction
        for(i = 1; i<N1-1 ;i++){
            for(j = 1; j<N1-1; j++){
                phi[i][j] += (delx1[i][j]);
            }
        }

        //Error Calculation
        err = 0;
        for (i = 0; i<N1; i++){
            for(j = 0; j<N1; j++){
                err += (phi_exact[i][j] - phi[i][j])*(phi_exact[i][j] - phi[i][j]);
            }
        }
        err = sqrt(err)/norm_exact;
        iter += 1;
    }

    t2 = omp_get_wtime();

    printf("Number of iterations for Serial Red Black with Multigrid are : %d\n",2*pre*iter);
    printf("The problem size is : %d x %d\n",N1,N1);
    printf("The error is : %lf\n",err);
    printf("The time taken is : %lf\n",t2-t1);

    //Freeing the allocated memory
    deallocate2DArray(phi, N1);
    deallocate2DArray(phi_exact, N1);
    deallocate2DArray(res1, N1);
    deallocate2DArray(res2, N2);
    deallocate2DArray(delx1, N1);
    deallocate2DArray(delx2, N2);
    deallocate2DArray(q, N1);


    return 0;
}