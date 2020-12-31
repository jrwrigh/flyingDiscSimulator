program flyingDiscSimulator

use core
use types
implicit none

type(disc_status) :: disc

disc%theta(:) = (/ 1.1D0, 1.2D0, 1.3D0 /)

call tester(disc)

print *, 'it works'

print *, T_a12(0D0, 0D0, 0D0)
print *, T_a12(disc%theta)
! print *, T_a23(0D0)
! print *, T_a34(0D0)

call calcuation()

contains
  subroutine calcuation()
    use types
    use core
    use aero
    use operator_matmul

    real*8, dimension(3) :: Cf, Cm
    real*8, dimension(3) :: F4, M4, F2, M2
    real*8, dimension(3,3) :: Ta12, Ta23, Ta34, Ta42, Tr12
    real*8, dimension(3,3) :: omega_tilde
    real*8 :: beta2, alpha3
    real*8 :: wind2(3), wind3(3)
    real*8 :: vel2(3), vel3(3)
    real*8 :: omega2(3)
    type(problemData) :: problem
    type(disc_status) :: state

    real*8 :: results(2,3)

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
    state%theta(:) = (/ 0d0, 0d0, 0d0 /)
    state%loc(:) = (/ 0d0, 0d0, 0d0 /)
    state%vel(:) = (/ 0d0, 0d0, 0d0 /)
    state%omega(:) = (/ 0d0, 0d0, 0d0 /)

    Ta12 = T_a12(state%theta(:))
    wind2 = transpose(Ta12).matmul.problem%env%wind
    vel2 = transpose(Ta12).matmul.state%vel
    beta2 = ATAN( (vel2(2) - wind2(2)) / (vel2(1) - wind2(1)) )

    Ta23 = T_a23(beta2)
    wind3 = transpose(Ta23).matmul.wind2
    vel3 = transpose(Ta23).matmul.vel2
    alpha3 = ATAN( (vel3(2) - wind2(2)) / (vel3(1) - wind2(1)) )

    Ta34 = T_a34(alpha3)

    call getAeroCoeffs(alpha3, Cf, Cm)
    call calcAeroMomentForce(Cf, Cm, F4, M4, problem, state)

    Ta42 = transpose(Ta23).matmul.transpose(Ta34)
    F2 = Ta42.matmul.F4 + (transpose(Ta12).matmul.problem%env%gravity)*problem%disc%m
    M2 = Ta42.matmul.M4

    Tr12 = T_r12(state%theta)
    omega2 = transpose(Tr12).matmul.state%omega

    results(1,:) = F2/problem%disc%m - omega2.cross.state%vel
    results(2,:) = problem%disc%I**(-1)  * (M2 - omega2.cross.(problem%disc%I*omega2))


    ! print *, vel2**(-1)
    print *, 'it ran!'
  end subroutine calcuation

end program flyingDiscSimulator
