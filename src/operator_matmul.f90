module operator_matmul
  implicit none

  private :: matmatmul_operator, matvecmul_operator, &
             vecmatmul_operator, cross_operator

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
      real*8, dimension(size(A,2),size(B,1)) :: C

      C = matmul(A,B)
    end function matmatmul_operator

    pure function matvecmul_operator(A, B) result(C)
      real*8, intent(in)  :: A(:,:), B(:)
      real*8, dimension(size(A,2)) :: C

      C = matmul(A,B)
    end function matvecmul_operator

    pure function vecmatmul_operator(A, B) result(C)
      real*8, intent(in)  :: A(:), B(:,:)
      real*8, dimension(size(A)) :: C

      C = matmul(A,B)
    end function vecmatmul_operator

    pure function cross_operator(x, y) result(z)
      real*8, intent(in) :: x(3), y(3)
      real*8             :: z(3)

      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
    end function cross_operator

end module operator_matmul
