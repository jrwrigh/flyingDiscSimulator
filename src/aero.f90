module aero
  implicit none

  private :: calcA_solution, calcA_YU

  interface calcA
    module procedure calcA_YU
    module procedure calcA_Solution
  end interface

  contains

    ! pure subroutine getAeroCoeffs(alpha3, U, D, Cf, Cm, Cm_damping)
    pure subroutine getAeroCoeffs(alpha3, omega3, apparent_velocity, Cf, Cm, problem)
      use transformation
      use types
      use math
      real*8, intent(in), dimension(3) :: omega3 !< Angular velocity solution
      real*8, intent(in) :: alpha3 !< Angle of Attack
      real*8, intent(in) :: apparent_velocity !< Linear velocity Magnitude
      type(problemData), intent(in) :: problem
      real*8, intent(out), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      real*8 :: damping_factor

      Cf = 0d0
      Cm = 0d0

      Cf(1) = -(3.3d0*(alpha3 + 0.052)**2 + 0.085d0) ! drag
      Cf(2) = 0d0  ! side
      Cf(3) = -(3.09d0*alpha3 + 0.13d0) ! lift

      damping_factor = problem%disc%D/(2*apparent_velocity)
      Cm(1) = problem%disc%Cm_damping(1)*omega3(1)*damping_factor ! rolling moment
      Cm(2) = 0.057*alpha3 - 0.01 + problem%disc%Cm_damping(2)*omega3(2)*damping_factor ! pitching moment
      ! Cm(2) = 1d-3 ! pitching moment
      Cm(3) = problem%disc%Cm_damping(3)*omega3(3)*damping_factor ! yawing moment (ie. resisting spin)
    end subroutine getAeroCoeffs

    pure subroutine calcA_YU(Y, U, problem, A, dAdU)
      use types
      use transformation
      use math

      real*8, dimension(6), intent(in), target :: Y, U
      type(problemData), intent(in):: problem
      real*8, intent(out) :: A(6), dAdU(6,6)

      real*8, dimension(3) :: Cf, Cm, Iinv
      real*8, dimension(3) :: F4, M4, F2, M2
      real*8, dimension(3,3) :: Ta12, Ta34, Ta42, Tr12
      real*8, dimension(3,3) :: Ta32
      real*8, dimension(3,3) :: F_factor, M_factor
      real*8, dimension(3,6) :: J_Cm, J_Cf
      real*8 :: beta2, alpha3
      real*8 :: wind2(3), wind3(3)
      real*8, pointer :: vel2(:)
      real*8 :: vel3(3)
      integer :: i, j
      
      Ta12 = T_a12(Y(4:6))
      wind2 = transpose(Ta12).matmul.problem%env%wind
      beta2 = ATAN2( (U(2) - wind2(2)), (U(1) - wind2(1)) )

      Ta32 = transpose(T_a23(beta2))
      wind3 = Ta32.matmul.wind2
      vel3 = Ta32.matmul.U(1:3)
      alpha3 = ATAN2( (vel3(3) - wind3(3)), (vel3(1) - wind3(1)) )

      Ta34 = T_a34(alpha3)
      Ta42 = Ta32.matmul.transpose(Ta34)

      associate(omega3 => wind2, apparent_velocity => beta2)
        omega3 = Ta32.matmul.U(4:6)
        apparent_velocity = norm2(U(1:3) - problem%env%wind)

        call getAeroCoeffs(alpha3, omega3, apparent_velocity, Cf, Cm, problem)

        dynpress: associate(q=> alpha3)
          q = 0.5*problem%env%density*apparent_velocity**2! dynamic pressure

          F4(:) = q*problem%disc%A*Cf(:)
          M4(:) = q*problem%disc%A*problem%disc%D*Cm(:)
        end associate dynpress


        F2 = (Ta42.matmul.F4) + (transpose(Ta12).matmul.problem%env%gravity)*problem%disc%m
        M2 = Ta42.matmul.M4

        Iinv = (problem%disc%I**(-1))
        A(:3) = F2/problem%disc%m
        A(4:) = M2*Iinv

        ! Calculating A contribution to Tangent Matrix (ie. dAdU)
        F_factor = 0.5d0*problem%env%density*problem%disc%A*Ta42
        M_factor = F_factor

        F_factor = F_factor/problem%disc%m
        M_factor = M_factor*problem%disc%D

        call calcCfJacobian(Ta32, problem, U, J_Cf)
        call calcCmJacobian(Ta32, apparent_velocity, omega3, problem, U, problem%disc%Cm_damping, J_Cm)
      end associate

      associate(Fkj => J_Cf, Mkj => J_Cm, velocity => beta2)
        velocity = norm2(U(1:3))
        Fkj = Fkj*velocity**2
        Mkj = Mkj*velocity**2

        forall(i=1:3, j=1:3) Fkj(i,j) = Fkj(i,j) + 2*Cf(i)*U(j)
        forall(i=1:3, j=1:3) Mkj(i,j) = Mkj(i,j) + 2*Cm(i)*U(j)

        dAdU(1:3,:) = F_factor.matmul.Fkj
        dAdU(4:6,:) = M_factor.matmul.Mkj
      end associate

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

    pure subroutine calcCmJacobian(Ta32, apparent_velocity, omega3, problem, U, Cm_damping, J_Cm)
      use math
      use types, only: problemData
      real*8, intent(in) :: U(6), Cm_damping(3)
      type(problemData), intent(in) :: problem
      real*8, intent(in):: omega3(3), apparent_velocity
      real*8, dimension(3,3), intent(in) :: Ta32
      real*8, intent (out) :: J_Cm(3,6)
      integer :: i

      forall(i=1:3) J_Cm(i,1:3) = U(1:3)*Cm_damping(i)*omega3(i)*problem%disc%D/apparent_velocity**3
      forall(i=1:3) J_Cm(i,4:6) = 0.5d0*Ta32(i,:)*Cm_damping(i)*problem%disc%D*apparent_velocity
    end subroutine

    pure subroutine calcCfJacobian(Ta32, problem, U, J_Cf)
      use math
      use types, only: problemData
      real*8, intent(in) :: Ta32(3,3), U(6)
      type(problemData), intent(in) :: problem
      real*8, intent (out) :: J_Cf(3,6)

      J_Cf = 0
    end subroutine

end module aero
