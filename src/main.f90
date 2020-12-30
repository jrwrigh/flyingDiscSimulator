program flyingDiscSimulator

use core
use types
implicit none

type(disc_status) :: disc

disc%phi(:) = (/ 1.1D0, 1.2D0, 1.3D0 /)

call tester(disc)

print *, 'it works'

print *, T_a12(0D0, 0D0, 0D0)
print *, T_a12(disc%phi)
! print *, T_a23(0D0)
! print *, T_a34(0D0)

call calcuation()

contains
  subroutine calcuation()
    use types
    use core
    use aero

    real*8, dimension(3) :: Cf, Cm
    real*8, dimension(3) :: F4, M4
    real*8, dimension(3,3) :: Ta12, Ta23, Ta34
    real*8 :: beta2, alpha3
    real*8 :: wind2(3), wind3(3)
    real*8 :: vel2(3), vel3(3)
    type(problemData) :: problem
    type(disc_status) :: state

    ! Initializing problem data
    problem%disc%I = (/ 0.1D0, 0.1D0, 0.2D0 /)
    problem%disc%m = 200
    problem%disc%A = 1D-5
    problem%disc%D = 0.2D0

    problem%env%wind = (/ 0d0, 1d0, 0d0 /)
    problem%env%density = 1.223d0
    problem%env%gravity = (/ 0d0, 0d0, 9.81d0 /)
    problem%env%viscosity = 1.5d-5

    ! Initialize state data
    state%phi(:) = (/ 0d0, 0d0, 0d0 /)
    state%loc(:) = (/ 0d0, 0d0, 0d0 /)
    state%vel(:) = (/ 0d0, 0d0, 0d0 /)
    state%omega(:) = (/ 0d0, 0d0, 0d0 /)

    Ta12 = T_a12(state%phi(:))
    wind2 = matmul(transpose(Ta12), problem%env%wind)
    vel2 = matmul(transpose(Ta12), state%vel)
    beta2 = ATAN( (vel2(2) - wind2(2)) / (vel2(1) - wind2(1)) )

    Ta23 = T_a23(beta2)
    wind3 = matmul(transpose(Ta23), wind2)
    vel3 = matmul(transpose(Ta23), vel2)
    alpha3 = ATAN( (vel3(2) - wind2(2)) / (vel3(1) - wind2(1)) )

    Ta34 = T_a34(alpha3)

    call getAeroCoeffs(alpha3, Cf, Cm)
    call calcAeroMomentForce(Cf, Cm, F4, M4, problem)
  end subroutine calcuation

end program flyingDiscSimulator
