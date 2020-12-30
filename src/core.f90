module core
  implicit none

  private :: T_a12_IndvAngles, T_a12_ArrayAngles

  interface T_a12
    module procedure T_a12_IndvAngles
    module procedure T_a12_ArrayAngles
  end interface T_a12

  contains
    subroutine tester(input)
      use types
      type(disc_status), intent(in) :: input
      print *, input%phi
    end subroutine tester

    function test() result(j)
      real*8 :: j
      j = 2
    end function test

    pure function T_a12_IndvAngles(phi, theta, psi) result(T_a)
      real*8 :: T_a(3,3)
      real*8, intent(in) :: phi, theta, psi
      real*8 :: cosphi, costheta, cospsi, &
                sinphi, sintheta, sinpsi
      cosphi = cos(phi)
      costheta = cos(theta)
      cospsi = cos(psi)
      sinphi = sin(phi)
      sintheta = sin(theta)
      sinpsi = sin(psi)

      T_a(1,1) =  cosphi*cospsi
      T_a(1,2) =  sinphi*sintheta*cospsi + cosphi*sinpsi
      T_a(1,3) = -cosphi*sintheta*cospsi + sinphi*sinpsi

      T_a(2,1) =  costheta*sinpsi
      T_a(2,2) =  sinphi*sintheta*sinpsi + cosphi*cospsi
      T_a(2,3) =  cosphi*sintheta*sinpsi + sinphi*cospsi

      T_a(3,1) =  sintheta
      T_a(3,2) = -sinphi*costheta
      T_a(3,3) =  cosphi*costheta
    end function T_a12_IndvAngles

    pure function T_a12_ArrayAngles(angles) result(T_a)
      real*8 :: T_a(3,3)
      real*8, intent(in) :: angles(3)
      real*8 :: phi, theta, psi
      real*8 :: cosphi, costheta, cospsi, &
                sinphi, sintheta, sinpsi
      cosphi = cos(angles(1))
      costheta = cos(angles(2))
      cospsi = cos(angles(3))
      sinphi = sin(angles(1))
      sintheta = sin(angles(2))
      sinpsi = sin(angles(3))

      T_a(1,1) =  cosphi*cospsi
      T_a(1,2) =  sinphi*sintheta*cospsi + cosphi*sinpsi
      T_a(1,3) = -cosphi*sintheta*cospsi + sinphi*sinpsi

      T_a(2,1) =  costheta*sinpsi
      T_a(2,2) =  sinphi*sintheta*sinpsi + cosphi*cospsi
      T_a(2,3) =  cosphi*sintheta*sinpsi + sinphi*cospsi

      T_a(3,1) =  sintheta
      T_a(3,2) = -sinphi*costheta
      T_a(3,3) =  cosphi*costheta
    end function T_a12_ArrayAngles

    pure function T_a23(beta) result(T_a)
      real*8 :: T_a(3,3)
      real*8, intent(in) :: beta
      real*8 :: cosbeta, sinbeta

      cosbeta = cos(beta)
      sinbeta = sin(beta)

      T_a(1,:) = (/  cosbeta, -sinbeta, 0D0 /)
      T_a(2,:) = (/ -sinbeta,  cosbeta, 0D0 /)
      T_a(3,:) = (/      0D0,      0D0, 1D0 /)

    end function T_a23

    pure function T_a34(alpha) result(T_a)
      real*8 :: T_a(3,3)
      real*8, intent(in) :: alpha
      real*8 :: cosalpha, sinalpha

      cosalpha = cos(alpha)
      sinalpha = sin(alpha)

      T_a(1,:) = (/  cosalpha, 0D0, -sinalpha /)
      T_a(2,:) = (/       0D0, 1D0,       0D0 /)
      T_a(3,:) = (/  sinalpha, 0D0,  cosalpha /)

    end function T_a34





end module core
