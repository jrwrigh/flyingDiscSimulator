program flyingDiscSimulator

use transformation
use types
use aero
use time_integrator
use math, only: pi
implicit none

! type(disc_status) :: disc
type(problemData) :: problem
! type(solution_state) :: soln, solnp1
type(solution_state_ptr) :: soln_ptr, solnp1_ptr, solnm1_ptr
integer, parameter :: nsteps=1e4
! type(solution_state) :: solution(nsteps)
real*8, allocatable, target :: solution(:,:)
type(solver_settings) :: solver
real*8, dimension(6) :: A
real*8 :: Imat(2,2), vec(2)
integer :: k, j, outu

vec = [ 1, 2 ]
Imat(1,:) = [2, 3]
Imat(2,:) = [4, 5]

! print*, Imat*vec

allocate(solution(18,nsteps))
solution = 0

! Initializing problem data
problem%disc%I = (/ 1.2d-3, 1.2d-3, 2.4d-3 /)
problem%disc%m = 0.175d0
problem%disc%D = 0.275d0
problem%disc%A = pi*(problem%disc%D*0.5)**2

problem%env%wind = (/ 0d0, 0d0, 0d0 /)
problem%env%density = 1.223d0
problem%env%gravity = (/ 0d0, 0d0, 9.81d0 /)
problem%env%viscosity = 1.5d-5

! forall(j=1:18) solution(j,1) = j

call soln_ptr%setptr(solution, 1)
soln_ptr%Y(:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln_ptr%U(:) = (/ 20d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln_ptr%Udot(:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln_ptr%Udot(1:3) = soln_ptr%Udot(1:3) + problem%env%gravity

solver%rho_infty = 1d0
solver%niters = 5
solver%delta_t = 1d-4
solver%tolerance = 1d-16

call iterations()
contains

  ! subroutine iterations(solver, soln_ptr, solnp1_ptr, outu, solution)
  subroutine iterations()
    use types, only:solver_settings, solution_state_ptr
    ! type(solver_settings), intent(inout) :: solver
    ! type(solution_state_ptr), intent(in) :: soln_ptr, solnp1_ptr
    ! integer, intent(inout) :: outu
    ! real*8, intent(inout) :: solution(:,:)

    call solver%calc_alphas_gamma()

    ! print *, 'STEP: ', 1
    call soln_ptr%setptr(solution, 1); call solnp1_ptr%setptr(solution, 2)

    call predictor_UdotConstant(soln_ptr, solnp1_ptr, solver)
    call iterate(soln_ptr, solnp1_ptr, solver, problem)
    do k=2,nsteps-1
      ! print *, 'STEP: ', k
      call soln_ptr%setptr(solution, k)
      call solnm1_ptr%setptr(solution, k-1); call solnp1_ptr%setptr(solution, k+1)

      call predictor_deltaUConstant(soln_ptr, solnm1_ptr, solnp1_ptr, solver)
      call iterate(soln_ptr, solnp1_ptr, solver, problem)
    end do

    outu = 0
    open (unit=outu,form="unformatted",file="test.fbin",action="write")
    write(outu) solution
    close(outu)
  end subroutine

end program flyingDiscSimulator
