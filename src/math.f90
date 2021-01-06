module math
  implicit none

  public :: pi, skewSymmetric, matinv3
  private :: matmatmul_operator, matvecmul_operator, &
             vecmatmul_operator, cross_operator

  real*8, parameter :: pi = 4*atan(1d0)

  interface operator (.matmul.)
    module procedure matmatmul_operator
    module procedure matvecmul_operator
    module procedure vecmatmul_operator
  end interface

  interface operator (.cross.)
    module procedure cross_operator
  end interface 

  contains

    pure function matmatmul_operator(A, B) result(C)
      real*8, intent(in), dimension(:,:) :: A, B
      real*8, dimension(size(A,1),size(B,2)) :: C

      C = matmul(A,B)
    end function matmatmul_operator

    pure function matvecmul_operator(A, B) result(C)
      real*8, intent(in)  :: A(:,:), B(:)
      real*8, dimension(size(A,1)) :: C

      C = matmul(A,B)
    end function matvecmul_operator

    pure function vecmatmul_operator(A, B) result(C)
      real*8, intent(in)  :: A(:), B(:,:)
      real*8, dimension(size(A)) :: C
      !^^^^ Wrong shape

      C = matmul(A,B)
    end function vecmatmul_operator

    pure function cross_operator(x, y) result(z)
      real*8, intent(in) :: x(3), y(3)
      real*8             :: z(3)

      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
    end function cross_operator

    pure subroutine skewSymmetric(vector, S)
      real*8, intent(in) :: vector(3)
      real*8, intent(out) :: S(3,3)

      S(1,:) = (/        0d0, -vector(3),  vector(2) /)
      S(2,:) = (/  vector(3),        0d0, -vector(1) /)
      S(3,:) = (/ -vector(2),  vector(1),        0d0 /)

    end subroutine skewSymmetric

    pure function matinv3(A) result(B)
      ! From http://fortranwiki.org/fortran/show/Matrix+inversion
      !! Performs a direct calculation of the inverse of a 3×3 matrix.
      real*8, intent(in) :: A(3,3)   !! Matrix
      real*8             :: B(3,3)   !! Inverse matrix
      real*8             :: detinv

      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

      ! Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function


end module math
