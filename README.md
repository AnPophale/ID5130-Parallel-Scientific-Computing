# ID5130-Parallel-Scientific-Computing
### Parallel implementation of the Multigrid method

This project was part of the course ID5130 Parallel Scientific Computing taken during the Jan-May 2024 semester. The main objective of this work is to parallelize the multigrid method using OpenMP for shared memory systems and MPI for distributed memory systems. Further, the performance of standard iterative solvers such as the Jacobi and Gauss Seidel Methods is compared with the Multigrid method along with parallel implementations of both.

A detailed explanation for work can be found in the [project report](ID5130 Project Report.pdf) and the C codes to implement all the methods are given in the [Codes](Codes) folder. All the data collected for different methods, such as the number of iterations, run time and speedup are tabulated in the [supplementary information](Supplementary Data.xlsx) 

**Motivation**  

**Multigrid Method**  

**Codes:**  
Jacobi Method:
* [Jacobi_Serial.c](Codes/Jacobi_Serial.c) - Serial Jacobi Method
* [Jacobi_OMP.c](Codes/Jacobi_OMP.c) - Parallel Jacobi Method using OpenMP
* [Jacobi_MPI.c](Codes/Jacobi_MPI.c) - Parallel Jacobi Method using MPI
* [JacobiMG_Serial.c](Codes/JacobiMG_Serial.c) - Serial Jacobi Method with Multigrid
* [JacobiMG_OMP.c](Codes/JacobiMG_OMP.c) - Parallel Jacobi Method with Multigrid using OpenMP
* [JacobiMG_MPI.c](Codes/JacobiMG_MPI.c) - Parallel Jacobi Method with Multigrid using MPI

Red Black Gauss Seidel Method:  
* [RedBlack_Serial.c](Codes/RedBlack_Serial.c) - Serial Red Black Method
* [RedBlack_OMP.c](Codes/RedBlack_OMP.c) - Parallel Red Black Method using OpenMP
* [RedBlack_MPI.c](Codes/RedBlack_MPI.c) - Parallel Red Black Method using MPI
* [RedBlackMG_Serial.c](Codes/RedBlackMG_Serial.c) - Serial Red Black Method with Multigrid
* [RedBlackMG_OMP.c](Codes/RedBlackMG_OMP.c) - Parallel Red Black Method with Multigrid using OpenMP
* [RedBlackMG_MPI.c](Codes/RedBlackMG_MPI.c) - Parallel Red Black Method with Multigrid using MPI

To run the parallel codes it is necessary to install the OpenMP and MPI libraries for C without which the codes may not compile and run properly.

**Results**  

**References:**  
[1] Mazumder, S. (2015). Numerical methods for partial differential equations: Finite Difference and Finite Volume Methods. Academic Press.
