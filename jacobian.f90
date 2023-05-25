subroutine jacobian(A,u,m,J)
   implicit none
   ! Evaluate the Jacobian F'(U)
   integer, intent(in) :: m
   real(kind=8),intent(in),dimension(m,(m-1)**2) :: A
   real(kind=8),intent(in),dimension((m-1)**2) :: u
   real(kind=8),intent(out),dimension(m,(m-1)**2) :: J
   integer :: i
   real(kind=8) :: y, beta, lambda
   beta = 0.12_8
   lambda = 0.19_8
   J = A
   ! Calculates F'(U) = A - diag(G'(U))
   do i=1,(m-1)**2
      y = 1.0_8 + beta*u(i)
      ! since J in banded format need mth row since take elements of diag
      J(m,i) = J(m,i) - lambda*exp(u(i)/y)/(y*y)
   end do
end subroutine jacobian
