HOW TO USE MY CODE 

Type 'make' in the cmd to compile both 'newton' and 'quasi_newton' implementations 
     
     - Alterantively Type 'make newton' or 'make quasi_newton' to compile individually

To run 'newton'/'quasi_newton', type './newton' or './quasi_newton' in cmd. 

To remove any temporary files type 'make clean' in cmd 

when specifying variables, 

     - m>1 and integer
     - tolerance must be real 
     - kmax must be integer

'Memory avaliable to store computed r_k? (yes-1,no-0):
	- input 1 to save all computed r_k and produce table required in q2 and 4c
	- input 0 if no memory avaliable, or want to optimise solution computation speed 

Print solution U? (yes-1,no-0)
      - input 1 to print the computed solution for U
      - input 0 to ignore the computed solution for U

Note: 
      - U(1/2,1/2) will only display if it exists
      - avergae time for iterative loop will only display if at least 1 iteration of the loop is completed 
      
FUNCTIONS 

quasi_newton.f90: computes the numerical solution to the nonlinear equation using efficient practical quasi_newton method (4a)

newton.f90: computes the numerical solution to the non-linear equation using newtons method

laplace.f90: produces upper banded matrix for A as described in assignment1

func.f90: computes F(U)
	  
jacobian.f90: computes F'(U)

bzero.f90: computes F'(U) and stores the cholesky decomposition in same m x (m-1)^2 matrix

OTHER

Makefile: compiles both quasi_newton and newton, with linking to BLAS and LAPACK 	
