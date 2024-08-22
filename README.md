# ID5130-Parallel-Scientific-Computing
### Parallel implementation of the Multigrid method

This project was part of the course ID5130 Parallel Scientific Computing taken during the Jan-May 2024 semester. The main objective of this work is to parallelize the multigrid method using OpenMP for shared memory systems and MPI for distributed memory systems. Further, the performance of standard iterative solvers such as the Jacobi and Gauss Seidel Methods is compared with the Multigrid method along with parallel implementations of both.

A detailed explanation for work can be found in the [project report](https://github.com/AnPophale/ID5130-Parallel-Scientific-Computing/blob/f3e6d66c87d151d2d5ce7e71c112c3df2c1f00d9/ID5130%20Project%20Report.pdf) and the C codes to implement all the methods are given in the [Codes](Codes) folder. The data collected for different methods, such as the number of iterations, runtime and speedup are tabulated in the [supplementary information](https://github.com/AnPophale/ID5130-Parallel-Scientific-Computing/blob/f3e6d66c87d151d2d5ce7e71c112c3df2c1f00d9/Supplementary%20Data.xlsx) 

**Motivation:**  

**Multigrid Method:**  


**Codes:**  
We consider an example of a 2D Poisson equation discretized using the Finite Difference method which leads to a system of linear equations to be solved. To solve these equations we implement the Jacobi and Gauss Seidel methods with and without applying the Multigrid principle as well as serial and parallel versions of each. An analytical solution exists for the chosen Poisson equation which is used to determine the convergence criteria.

Jacobi Method:
* [Jacobi_Serial.c](Codes/Jacobi_Serial.c) - Serial Jacobi Method
* [Jacobi_OMP.c](Codes/Jacobi_OMP.c) - Parallel Jacobi Method using OpenMP
* [Jacobi_MPI.c](Codes/Jacobi_MPI.c) - Parallel Jacobi Method using MPI
* [JacobiMG_Serial.c](Codes/JacobiMG_Serial.c) - Serial Jacobi Method with Multigrid
* [JacobiMG_OMP.c](Codes/JacobiMG_OMP.c) - Parallel Jacobi Method with Multigrid using OpenMP
* [JacobiMG_MPI.c](Codes/JacobiMG_MPI.c) - Parallel Jacobi Method with Multigrid using MPI

Red Black Gauss Seidel Method:  
The Red Black Gauss Seidel Method is a modification of the standard Gauss Seidel algorithm allowing it to run in parallel which is not possible otherwise due to the data depedency inherent to the method. More details about this method can be found [here](https://ocw.mit.edu/courses/16-920j-numerical-methods-for-partial-differential-equations-sma-5212-spring-2003/2351cb5ce7f15d89fa4cb3fa17eb6f64_lec6_notes.pdf)  
* [RedBlack_Serial.c](Codes/RedBlack_Serial.c) - Serial Red Black Gauss Seidel Method
* [RedBlack_OMP.c](Codes/RedBlack_OMP.c) - Parallel Red Black Method using OpenMP
* [RedBlack_MPI.c](Codes/RedBlack_MPI.c) - Parallel Red Black Method using MPI
* [RedBlackMG_Serial.c](Codes/RedBlackMG_Serial.c) - Serial Red Black Method with Multigrid
* [RedBlackMG_OMP.c](Codes/RedBlackMG_OMP.c) - Parallel Red Black Method with Multigrid using OpenMP
* [RedBlackMG_MPI.c](Codes/RedBlackMG_MPI.c) - Parallel Red Black Method with Multigrid using MPI

To run the parallel codes it is necessary to install the OpenMP and MPI libraries for C without which the codes may not compile and run properly.

**Results:**  
Some of the results for the Jacobi and Gauss Seidel methods with Multigrid are shown below in Fig. 3 and 4. We analyze the variation of the speed up with the problem size for a fixed number of processors as well as the variation of speedup with number of processors used while keeping the problem size fixed.

<p align="center">
  <img src="https://github.com/user-attachments/assets/000fd2fa-0ac3-4a78-a99a-c254e56ede4a" alt="Fig 3: Speed up vs Problem Size for parallel multigrid methods using 4 processors" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 3: Speed up vs Problem Size for parallel multigrid methods using 4 processors</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/7370662f-b3de-4e6d-b10b-cc3a5be7ed9b" alt="Fig 4: Speed up vs Number of processors for parallel multigrid methods using a problem size of 513" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 4: Speed up vs Number of processors for parallel multigrid methods using a problem size of 513</em>
</p>

**References:**  
[1] Mazumder, S. (2015). Numerical methods for partial differential equations: Finite Difference and Finite Volume Methods. Academic Press.
