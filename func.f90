subroutine func(A,u,m,r)
   implicit none
   ! Evaluates the nonlinear function F(U)=AU-G(U)
   integer, intent(in) :: m
   real(kind=8), intent(in), dimension((m-1)**2) :: u
   real(kind=8),intent(in), dimension(m,(m-1)**2) :: A
   real(kind=8),intent(out), dimension((m-1)**2) :: r
   integer :: i, j, n
   real (kind=8) :: x1, x2, uij, h
   real(kind=8) :: beta, lambda
   real(kind=8), parameter :: pi = 3.1415926535897931_8
   ! Total number of unknowns
   n = (m-1)**2
   beta = 0.12_8
   lambda = 0.19_8
   ! Cell size
   h = 1.0_8/real(m,8)
   ! Compute the nonlinear part G(U) first
   do i = 1,m-1
      do j = 1,m-1
         ! Coordinates, and
         x1  = i*h
         x2  = j*h
         ! Solution at (i,j)th point
         uij = u((j-1)*(m-1)+i)
         ! We map (i,j) indices lexicographically into the global index
         ! G_{i,j}
         r((j-1)*(m-1)+i) = lambda*exp(uij/(1.0_8+beta*uij)) + &
             1.0d2*sin(pi*x1)*sin(pi*x2)
      end do
   end do
   ! Compute AU-G(U) and stores the answer in r
   call dsbmv('U', n, m-1, 1.0_8,A,m,u,1,-1.0_8,r,1)
end subroutine func


