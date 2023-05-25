F95 = gfortran
OPT = -O3 -Wall -fimplicit-none

# Link to BLAS and LAPACK

LIB_LIST = -L${BLASDIR} ${BLASLIB} -lpthread

all: newton quasi_newton

func.o: func.f90
	$(F95) $(OPT) -c func.f90

laplace.o: laplace.f90 
	$(F95) $(OPT) -c laplace.f90

jacobian.o: jacobian.f90
	$(F95) $(OPT) -c jacobian.f90

newton.o: newton.f90
	$(F95) $(OPT) -c newton.f90

quasi_newton.o: quasi_newton.f90
	$(F95) $(OPT) -c quasi_newton.f90

bzero.o: bzero.f90
	$(F95) $(OPT) -c bzero.f90

quasi_newton: quasi_newton.o laplace.o bzero.o func.o   
	$(F95) -o quasi_newton quasi_newton.o laplace.o bzero.o func.o $(LIB_LIST)

newton: newton.o laplace.o jacobian.o func.o 
	$(F95) -o  newton newton.o laplace.o jacobian.o func.o $(LIB_LIST) 

clean: 
	rm -f *.o core.* newton quasi_newton
# creates a comand to remove any temp files
