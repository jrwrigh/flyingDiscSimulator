program flyingDiscSimulator

use transformation
use types
use aero
use time_integrator
! use math, only: pi
use math
implicit none

! type(disc_status) :: disc
type(problemData) :: problem
! type(solution_state) :: soln, solnp1
type(solution_state_ptr) :: soln_ptr, solnp1_ptr, solnm1_ptr
integer, parameter :: nsteps=1e4
! type(solution_state) :: solution(nsteps)
real*8, allocatable, target :: solution(:,:)
real*8, allocatable, target :: residual(:)
type(solver_settings) :: solver
integer :: k, j, outu
! print*, Imat*vec

allocate(solution(18,nsteps))
solution = 0

! Initializing problem data
problem%disc%I = (/ 1.2d-3, 1.2d-3, 2.4d-3 /)
problem%disc%Cm_damping = (/ -1.3d0, -1.4d0, -1.2d-1/)
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
soln_ptr%U(:) = (/ 20d0, 0d0, 0d0, 0d0, 0d0, 52.85d0  /)
! soln_ptr%U(:) = (/ 20d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln_ptr%Udot(:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0  /)
soln_ptr%Udot(1:3) = soln_ptr%Udot(1:3) + problem%env%gravity

solver%rho_infty = 1d0
solver%niters = 5
solver%delta_t = 1d-4
solver%tolerance = 1d-16

allocate(residual(nsteps))
residual = 0

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
    call iterate(soln_ptr, solnp1_ptr, solver, problem, residual(1))
    do k=2,nsteps-1
      ! print *, 'STEP: ', k
      call soln_ptr%setptr(solution, k)
      call solnm1_ptr%setptr(solution, k-1); call solnp1_ptr%setptr(solution, k+1)

      call predictor_deltaUConstant(soln_ptr, solnm1_ptr, solnp1_ptr, solver)
      call iterate(soln_ptr, solnp1_ptr, solver, problem, residual(k))
    end do

    outu = 0
    open (unit=outu,form="unformatted",file="test.fbin",action="write")
    write(outu) solution
    close(outu)

    outu = 1
    open (unit=outu,form="unformatted",file="residual.fbin",action="write")
    write(outu) residual
    close(outu)
  end subroutine

end program flyingDiscSimulator
