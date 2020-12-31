module aero
  implicit none

  contains

    pure subroutine getAeroCoeffs(alpha3, Cf, Cm)
      use core
      real*8, intent(out), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      real*8, intent(in) :: alpha3 !< Angle of Attack

      Cf(1) = -0.1*alpha3 ! drag
      Cf(2) = 0.1*alpha3  ! side
      Cf(3) = -0.1*alpha3**2 ! drag

      Cm(1) = -0.1*alpha3 ! rolling moment
      Cm(2) = -0.1*alpha3 ! pitching moment
      Cm(3) = -0.1*alpha3 ! yawing moment (ie. resisting spin)

    end subroutine getAeroCoeffs

    pure subroutine calcAeroMomentForce(Cf, Cm, F4, M4, problem, state)
      use types
      real*8, intent(in), dimension(3) :: Cf, Cm !< Vector of force and moment coefficients
      type(problemData), intent(in) :: problem
      type(disc_status), intent(in) :: state
      real*8, intent(out), dimension(3) :: F4, M4 !< Force and moment vectors
      real*8 :: q

      q = 0.5*problem%env%density*norm2(state%vel - problem%env%wind)**2! dynamic pressure

      F4(:) = q*problem%disc%A*Cf(:)
      M4(:) = q*problem%disc%A*problem%disc%D*Cm(:)

    end subroutine calcAeroMomentForce




end module aero
