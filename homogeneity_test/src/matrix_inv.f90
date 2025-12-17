! ! test the matrix inverse function inv
! program test_inv
!   implicit none
!   real(8) :: aa(3,3),bb(3,3)
!   integer :: i
!   interface
!      function inv(A) result(Ainv)
!        real(8), dimension(:,:), intent(in) :: A
!        real(8), dimension(size(A,1),size(A,2)) :: Ainv
!      end function inv
!   end interface
!   aa(1,1)=1.d0
!   aa(1,2)=0.d0
!   aa(1,3)=2.d0
!   aa(2,1)=3.d0
!   aa(2,2)=5.d0
!   aa(2,3)=1.d0
!   aa(3,1)=0.d0
!   aa(3,2)=2.d0
!   aa(3,3)=6.d0
!   bb=inv(aa)
!   write(*,*) 'inv of aa is :', (bb(i,:),i=1,3)
!   return
! end program test_inv





! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  include 'info.f'
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n!, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
!     stop 'Matrix is numerically singular!'
     write(*,*) 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
!     stop 'Matrix inversion failed!'
     write(*,*) 'Matrix inversion failed!'
  end if
end function inv
