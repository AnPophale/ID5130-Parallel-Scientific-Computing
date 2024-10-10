# ID5130-Parallel-Scientific-Computing
### Parallel implementation of the Multigrid method

This project was part of the course ID5130 Parallel Scientific Computing taken during the Jan-May 2024 semester. The main objective of this work is to parallelize the multigrid method using OpenMP for shared memory systems and MPI for distributed memory systems. Further, the performance of standard iterative solvers such as the Jacobi and Gauss Seidel methods is compared with the Multigrid method, along with parallel implementations of both.

Presented below is a short summary of the project, while a detailed explanation of the work can be found in the [project report](https://github.com/AnPophale/ID5130-Parallel-Scientific-Computing/blob/f3e6d66c87d151d2d5ce7e71c112c3df2c1f00d9/ID5130%20Project%20Report.pdf) along with the C codes to implement all the methods, which have been uploaded in the [codes](Codes) folder. The data collected for different methods, such as the number of iterations, runtime and speedup is tabulated in the [supplementary information](https://github.com/AnPophale/ID5130-Parallel-Scientific-Computing/blob/f3e6d66c87d151d2d5ce7e71c112c3df2c1f00d9/Supplementary%20Data.xlsx).

**Motivation:**  
Tackling engineering problems with a computational approach often involves solving partial differential equations. These equations when solved numerically using methods such as the finite volume or finite difference methods generate a system of linear equations that need to be solved. The solution to these equations is generally obtained using iterative methods such as the Jacobi and Gauss Seidel methods. It has been shown that as the error decreases, the rate of convergence of these methods decreases making them slow.  

The Multigrid method is a principle which can be applied to such iterative methods to accelerate their convergence and obtain solutions with a reduced computational effort which is especially important for larger problems. In this work, we parallelize the multigrid method to investigate the potential speed up benefits when the method is distributed across several processors.  

<p align="center">
  <img src="https://github.com/user-attachments/assets/baae457e-326a-4da0-81b2-cb54e8bb9401" alt="Fig 1: The Multigrid Algorithm" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 1: The Multigrid Algorithm</em>
</p>

**Multigrid Method:**    
Decomposing the error field over the discretized domain for iterative solvers using a Fourier series reveals some interesting and useful properties of the error. It is observed that high frequency error components are eliminated quickly but the low frequency components slow down the convergence, requiring a larger number of equations to smoothen out. The multigrid method makes use of this property of the error field. 

The key idea in the multigrid method is that low frequency error components on a fine mesh, have a high frequency on a coarse mesh which can be eliminated quickly. Hence, in multigrid method, the errors are transferred from a fine mesh to a coarser mesh where they can be smoothed out faster. This methodology can be applied to any of the iterative solvers which are known as smoothers in this context and we use the Jacobi and Gauss Seidel methods as smoothers in this study. This algorithm for the Multigrid method is depicted in Fig. 1 which shows an example of the V cycle multigrid with 3 mesh levels.  

<p align="center">
  <img src="https://github.com/user-attachments/assets/c85762aa-2678-4a7a-8877-985d9b91a2ff" alt="Fig 2: Types of Multigrid cycles" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 2: Types of Multigrid cycles</em>
</p>

Based on the number of levels of the coarse and fine meshes as well as the patterns in which the solution data is transferred across various meshes, there are 3 types of multigrid methods which are the V cycle, W cycle and the F cycle as shown in Fig. 2. Further the method can be classified as Geometric Multigrid (GMG) or Algebraic Multigrid (AMG), where GMG requires the geometric information of all mesh levels while AMG does not need such data and hence is useful for more complicated cases such as with unstructured meshes.  

In this project, we use 2 level V cycle geometric multigrid for simplicity as the main aim is to parallelize the code and compare its performance with other methods.

**Codes:**  
We consider an example of a 2D Poisson equation discretized using the Finite Difference method which leads to a system of linear equations to be solved. To solve these equations we implement the Jacobi and Gauss Seidel methods with and without applying the Multigrid principle as well as serial and parallel versions of each. An analytical solution exists for the chosen Poisson equation which is used to determine the convergence criteria.  

_Jacobi Method_
* [Jacobi_Serial.c](Codes/Jacobi_Serial.c) - Serial Jacobi Method
* [Jacobi_OMP.c](Codes/Jacobi_OMP.c) - Parallel Jacobi Method using OpenMP
* [Jacobi_MPI.c](Codes/Jacobi_MPI.c) - Parallel Jacobi Method using MPI
* [JacobiMG_Serial.c](Codes/JacobiMG_Serial.c) - Serial Jacobi Method with Multigrid
* [JacobiMG_OMP.c](Codes/JacobiMG_OMP.c) - Parallel Jacobi Method with Multigrid using OpenMP
* [JacobiMG_MPI.c](Codes/JacobiMG_MPI.c) - Parallel Jacobi Method with Multigrid using MPI

_Red Black Gauss Seidel Method:_  
The Red Black Gauss Seidel Method is a modification of the standard Gauss Seidel algorithm allowing it to run in parallel which is not possible otherwise due to the data dependency inherent to the method. More details about this method can be found [here](https://ocw.mit.edu/courses/16-920j-numerical-methods-for-partial-differential-equations-sma-5212-spring-2003/2351cb5ce7f15d89fa4cb3fa17eb6f64_lec6_notes.pdf).  
* [RedBlack_Serial.c](Codes/RedBlack_Serial.c) - Serial Red Black Gauss Seidel Method
* [RedBlack_OMP.c](Codes/RedBlack_OMP.c) - Parallel Red Black Method using OpenMP
* [RedBlack_MPI.c](Codes/RedBlack_MPI.c) - Parallel Red Black Method using MPI
* [RedBlackMG_Serial.c](Codes/RedBlackMG_Serial.c) - Serial Red Black Method with Multigrid
* [RedBlackMG_OMP.c](Codes/RedBlackMG_OMP.c) - Parallel Red Black Method with Multigrid using OpenMP
* [RedBlackMG_MPI.c](Codes/RedBlackMG_MPI.c) - Parallel Red Black Method with Multigrid using MPI

To run the parallel codes it is necessary to install the OpenMP and MPI libraries for C without which the codes may not compile and run properly.

**Results:**  
Some of the results for the Jacobi and Gauss Seidel methods with Multigrid are shown below in Fig. 3 and 4. We analyze the variation of the speed up with the problem size for a fixed number of processors as well as the variation of speed up with number of processors used while keeping the problem size fixed.

<p align="center">
  <img src="https://github.com/user-attachments/assets/000fd2fa-0ac3-4a78-a99a-c254e56ede4a" alt="Fig 3: Speed up vs Problem Size for parallel multigrid methods using 4 processors" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 3: Speedup vs Problem Size for Parallel Multigrid Methods using 4 processors</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/7370662f-b3de-4e6d-b10b-cc3a5be7ed9b" alt="Fig 4: Speed up vs Number of processors for parallel multigrid methods using a problem size of 513" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 4: Speed up vs Number of processors for Parallel Multigrid Methods using a problem size of 513</em>
</p>


**Conclusions:**  
Some of the key conclusions from this study are as follows:
* The Multigrid method requires a significantly lower number of iterations compared to Jacobi or Gauss Seidel methods as expected.
* The Multigrid method with the Red Black smoother was found to be the fastest among all the methods.
* Along with this the OpenMP parallelization of this method provides the best performance when compared to MPI. For OpenMP, 16 threads give the best performance, while for MPI 8 processors provide the best speedup.
* This work can be further extended using various levels and cycles for the multigrid method as well as using GPU acceleration to parallelize the multigrid algorithm. 


**References:**  
[1] Mazumder, S. (2015). Numerical methods for partial differential equations: Finite Difference and Finite Volume Methods. Academic Press.
