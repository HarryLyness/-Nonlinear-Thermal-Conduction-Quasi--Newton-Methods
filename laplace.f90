subroutine laplace(AB,m) 

  implicit none

  ! fortran95  function  to assemble the banded matrix (upper triangle only)
  ! arising from the finite difference discretisation of the Laplacian
  ! partial differential operator.
  !
  ! Recall the structure for m=4 is
  ! *  *  *   -1 -1 -1    -1 -1 -1   } -I 
  ! *  *  0    0  0  0     0  0  0
  ! * -1 -1    0 -1 -1     0 -1 -1   }  
  ! 4  4  4    4  4  4     4  4  4   }  B
  

  integer, intent(in) :: m
  real(kind=8),  intent(out),  dimension(m,(m-1)**2) :: AB
  integer :: i

  ! initialise with zeros
  AB(:,:) = 0.0_8

  ! Main diagonal
  AB(m,:) = 4.0_8
  ! Outer off-diagonal. Since ghost elements are not used, just set all to -1
  AB(1,:) = -1.0_8
  ! Inner off-diagonal. Start with filling the majority of elements (-1)
  AB(m-1,:) = -1.0_8
  ! Set zeros separating the blocks B in the off-diagonal
  do i=1,m-1
     AB(m-1, 1+(i-1)*(m-1)) = 0.0_8
  end do

  AB(:,:) = (m**2)*AB(:,:)   ! 1/h^2 factor

end subroutine laplace
