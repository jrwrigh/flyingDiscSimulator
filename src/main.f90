program flyingDiscSimulator

use transformation
use types
use aero
use time_integrator
use math, only: pi
implicit none

! type(disc_status) :: disc
type(problemData) :: problem
type(solution_state) :: soln, solnp1
integer, parameter :: nsteps=9
type(solution_state) :: solution(nsteps)
type(solver_settings) :: solver
real*8, dimension(6) :: A
real*8 :: Imat(3,3)
integer :: k, j

! Initializing problem data
problem%disc%I = (/ 1.2d-3, 1.2d-3, 2.4d-3 /)
problem%disc%m = 0.175d0
problem%disc%D = 0.275d0
problem%disc%A = pi*(problem%disc%D*0.5)**2

problem%env%wind = (/ 0d0, 0d0, 0d0 /)
problem%env%density = 1.223d0
problem%env%gravity = (/ 0d0, 0d0, 9.81d0 /)
problem%env%viscosity = 1.5d-5

soln%Y(:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln%U(:) = (/ 20d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln%Udot(:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)

solution(1) = soln

solver%rho_infty = 0.5
solver%nsteps = 3
solver%delta_t = 1d-4
solver%tolerance = 1d-16

call calc_alphas_gamma(solver)

! call predictor_UdotConstant(soln, solnp1, solver)
! call iterate(soln, solnp1, solver, problem)

print *, 'STEP: ', 1
call predictor_UdotConstant(solution(1), solution(2), solver)
call iterate(solution(1), solution(2), solver, problem)
do k=2,nsteps-1
  print *, 'STEP: ', k
  call predictor_deltaUConstant(solution(k), solution(k-1), solution(k+1), solver)
  call iterate(solution(k), solution(k+1), solver, problem)
end do

! Imat = 0
! forall(j=1:3) Imat(j,j) = problem%disc%I(j)
! print *, Imat

contains
  ! subroutine calcuation()
  !   use types
  !   use transformation
  !   use aero
  !   use operator_matmul

  !   real*8, dimension(3) :: Cf, Cm
  !   real*8, dimension(3) :: F4, M4, F2, M2
  !   real*8, dimension(3,3) :: Ta12, Ta23, Ta34, Ta42, Tr12
  !   real*8, dimension(3,3) :: omega_tilde
  !   real*8 :: beta2, alpha3
  !   real*8 :: wind2(3), wind3(3)
  !   real*8 :: vel2(3), vel3(3)
  !   real*8 :: omega2(3)
  !   type(problemData) :: problem
  !   type(disc_status) :: state

  !   real*8 :: results(2,3)

  !   ! Initializing problem data
  !   problem%disc%I = (/ 0.1D0, 0.1D0, 0.2D0 /)
  !   problem%disc%m = 200
  !   problem%disc%A = 1D-5
  !   problem%disc%D = 0.2D0

  !   problem%env%wind = (/ 0d0, 1d0, 0d0 /)
  !   problem%env%density = 1.223d0
  !   problem%env%gravity = (/ 0d0, 0d0, 9.81d0 /)
  !   problem%env%viscosity = 1.5d-5

  !   ! Initialize state data
  !   state%theta(:) = (/ 0d0, 0d0, 0d0 /)
  !   state%loc(:) = (/ 0d0, 0d0, 0d0 /)
  !   state%vel(:) = (/ 0d0, 0d0, 0d0 /)
  !   state%omega(:) = (/ 0d0, 0d0, 0d0 /)

  !   Ta12 = T_a12(state%theta(:))
  !   wind2 = transpose(Ta12).matmul.problem%env%wind
  !   vel2 = transpose(Ta12).matmul.state%vel
  !   beta2 = ATAN( (vel2(2) - wind2(2)) / (vel2(1) - wind2(1)) )

  !   Ta23 = T_a23(beta2)
  !   wind3 = transpose(Ta23).matmul.wind2
  !   vel3 = transpose(Ta23).matmul.vel2
  !   alpha3 = ATAN( (vel3(2) - wind2(2)) / (vel3(1) - wind2(1)) )

  !   Ta34 = T_a34(alpha3)

  !   call getAeroCoeffs(alpha3, Cf, Cm)
  !   call calcAeroMomentForce(Cf, Cm, F4, M4, problem, state)

  !   Ta42 = transpose(Ta23).matmul.transpose(Ta34)
  !   F2 = Ta42.matmul.F4 + (transpose(Ta12).matmul.problem%env%gravity)*problem%disc%m
  !   M2 = Ta42.matmul.M4.matmul.Ta42

  !   Tr12 = T_r12(state%theta)
  !   omega2 = transpose(Tr12).matmul.state%omega

  !   results(1,:) = F2/problem%disc%m - omega2.cross.state%vel
  !   results(2,:) = problem%disc%I**(-1)  * (M2 - omega2.cross.(problem%disc%I*omega2))


    ! print *, vel2**(-1)
    ! print *, shape(F4(2:3))
    ! print *, 'it ran!'
  ! end subroutine calcuation

end program flyingDiscSimulator
