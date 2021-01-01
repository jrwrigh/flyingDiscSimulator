module aero
  implicit none

  contains

    pure subroutine getAeroCoeffs(alpha3, Cf, Cm)
      use transformation
      real*8, intent(out), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      real*8, intent(in) :: alpha3 !< Angle of Attack

      Cf(1) = -0.1*alpha3 ! drag
      Cf(2) = 0.1*alpha3  ! side
      Cf(3) = -0.1*alpha3**2 ! drag

      Cm(1) = -0.1*alpha3 ! rolling moment
      Cm(2) = -0.1*alpha3 ! pitching moment
      Cm(3) = -0.1*alpha3 ! yawing moment (ie. resisting spin)

    end subroutine getAeroCoeffs

    pure subroutine calcAeroMomentForce(Cf, Cm, F4, M4, problem, solution)
      use types
      real*8, intent(in), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      type(problemData), intent(in) :: problem
      type(solution_state), intent(in) :: solution
      real*8, intent(out), dimension(3) :: F4, M4 !< Force and moment vectors
      real*8 :: q

      q = 0.5*problem%env%density*norm2(solution%U(1:3) - problem%env%wind)**2! dynamic pressure

      F4(:) = q*problem%disc%A*Cf(:)
      M4(:) = q*problem%disc%A*problem%disc%D*Cm(:)

    end subroutine calcAeroMomentForce

    pure subroutine calcA(solution, problem, A)
      use types
      use transformation
      use operator_matmul

      type(solution_state), intent(in) :: solution
      type(problemData), intent(in):: problem
      real*8, intent(out) :: A(6)

      real*8, dimension(3) :: Cf, Cm
      real*8, dimension(3) :: F4, M4, F2, M2
      real*8, dimension(3,3) :: Ta12, Ta23, Ta34, Ta42, Tr12
      real*8, dimension(3,3) :: omega_tilde
      real*8 :: beta2, alpha3
      real*8 :: wind2(3), wind3(3)
      real*8 :: vel2(3), vel3(3)
      real*8 :: omega2(3)

      Ta12 = T_a12(solution%Y(4:6))
      wind2 = transpose(Ta12).matmul.problem%env%wind
      vel2 = transpose(Ta12).matmul.solution%U(1:3)
      beta2 = ATAN( (vel2(2) - wind2(2)) / (vel2(1) - wind2(1)) )

      Ta23 = T_a23(beta2)
      wind3 = transpose(Ta23).matmul.wind2
      vel3 = transpose(Ta23).matmul.vel2
      alpha3 = ATAN( (vel3(2) - wind3(2)) / (vel3(1) - wind3(1)) )

      Ta34 = T_a34(alpha3)

      call getAeroCoeffs(alpha3, Cf, Cm)
      call calcAeroMomentForce(Cf, Cm, F4, M4, problem, solution)

      Ta42 = transpose(Ta23).matmul.transpose(Ta34)
      F2 = Ta42.matmul.F4 + (transpose(Ta12).matmul.problem%env%gravity)*problem%disc%m
      M2 = Ta42.matmul.M4

      A(:3) = F2
      A(4:) = M2

    end subroutine

end module aero
