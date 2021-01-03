module aero
  implicit none

  private :: calcA_solution, calcA_YU

  interface calcA
    module procedure calcA_YU
    module procedure calcA_Solution
  end interface

  contains

    pure subroutine getAeroCoeffs(alpha3, Cf, Cm)
      use transformation
      real*8, intent(out), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      real*8, intent(in) :: alpha3 !< Angle of Attack

      Cf(1) = -0.1*alpha3**2 ! drag
      Cf(2) = 0d0  ! side
      Cf(3) = -3.1*alpha3**2 ! lift

      ! Cm(1) = -0.1*alpha3 ! rolling moment
      ! Cm(2) = -0.1*alpha3 ! pitching moment
      ! Cm(3) = -0.1*alpha3 ! yawing moment (ie. resisting spin)

      ! Cf = 0d0
      Cm = 0d0

    end subroutine getAeroCoeffs

    pure subroutine calcAeroMomentForce(Cf, Cm, F4, M4, problem, U)
      use types
      real*8, intent(in), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      type(problemData), intent(in) :: problem
      real*8, intent(in), dimension(6) :: U
      real*8, intent(out), dimension(3) :: F4, M4 !< Force and moment vectors
      real*8 :: q

      q = 0.5*problem%env%density*norm2(U(1:3) - problem%env%wind)**2! dynamic pressure

      F4(:) = q*problem%disc%A*Cf(:)
      M4(:) = q*problem%disc%A*problem%disc%D*Cm(:)

    end subroutine calcAeroMomentForce

    pure subroutine calcA_YU(Y, U, problem, A, dAdU)
      use types
      use transformation
      use math

      real*8, dimension(6), intent(in) :: Y, U
      type(problemData), intent(in):: problem
      real*8, intent(out) :: A(6), dAdU(6,6)

      real*8, dimension(3) :: Cf, Cm, Iinv
      real*8, dimension(3) :: F4, M4, F2, M2
      real*8, dimension(3,3) :: Ta12, Ta23, Ta34, Ta42, Tr12
      real*8, dimension(3,3) :: F_factor, M_factor
      real*8 :: beta2, alpha3
      real*8 :: wind2(3), wind3(3)
      real*8 :: vel2(3), vel3(3)
      real*8 :: omega2(3)
      integer :: j

      Ta12 = T_a12(Y(4:6))
      wind2 = transpose(Ta12).matmul.problem%env%wind
      vel2 = transpose(Ta12).matmul.U(1:3)
      beta2 = ATAN( (vel2(2) - wind2(2)) / (vel2(1) - wind2(1)) )

      Ta23 = T_a23(beta2)
      wind3 = transpose(Ta23).matmul.wind2
      vel3 = transpose(Ta23).matmul.vel2
      alpha3 = ATAN( (vel3(2) - wind3(2)) / (vel3(1) - wind3(1)) )

      Ta34 = T_a34(alpha3)

      call getAeroCoeffs(alpha3, Cf, Cm)
      call calcAeroMomentForce(Cf, Cm, F4, M4, problem, U)

      Ta42 = transpose(Ta23).matmul.transpose(Ta34)
      F2 = Ta42.matmul.F4 + (transpose(Ta12).matmul.problem%env%gravity)*problem%disc%m
      M2 = Ta42.matmul.M4

      Iinv = (problem%disc%I**(-1))
      A(:3) = F2/problem%disc%m
      A(4:) = M2*Iinv

      ! Calculating A contribution to Tangent Matrix (ie. dAdU)
      F_factor = 0.5d0*problem%env%density*problem%disc%A*Ta42
      M_factor = F_factor

      F_factor = F_factor/problem%disc%m
      M_factor = M_factor*problem%disc%D

      dAdU = 0
      forall(j=1:3) dAdU(j,j) = 2*U(j)*Cf(j)
      forall(j=1:3) dAdU(j+3,j+3) = 2*U(j)*Cm(j)*Iinv(j)

      dAdU(1:3,1:3) = F_factor.matmul.dAdU(1:3,1:3)
      dAdU(4:6,4:6) = F_factor.matmul.dAdU(4:6,4:6)

    end subroutine

    pure subroutine calcA_solution(solution, problem, A, dAdU)
      use types

      type(solution_state), intent(in) :: solution
      type(problemData), intent(in):: problem
      real*8, intent(out) :: A(6), dAdU(6,6)
      real*8, dimension(6):: Y, U

      Y = solution%Y
      U = solution%U

      call calcA_YU(Y, U, problem, A, dAdU)

    end subroutine

end module aero
