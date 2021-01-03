module time_integrator
  ! use types, only: solver_settings, solution_state
  implicit none

  ! type :: solution_state
  !   real*8, dimension(6) :: Y !< axis1 positions and attitude
  !   real*8, dimension(6) :: U !< axis2 linear and angular velocity
  !   real*8, dimension(6) :: Udot !< axis2 linear and angular accelerations
  ! end type solution_state

  ! type(solver_settings) :: solver

  contains

    pure subroutine calc_alphas_gamma(solver)
      use types, only: solver_settings
      type(solver_settings), intent(inout) :: solver

      solver%alpha_m = (3 - solver%rho_infty)/(1 + solver%rho_infty)
      solver%alpha_f = 1/(1 + solver%rho_infty)
      solver%gamma = 0.5d0 + solver%alpha_m - solver%alpha_f
    end subroutine

    pure subroutine predictor_deltaUConstant(soln, solnm1, solnp1, solver)
      use types, only: solver_settings, solution_state
      use transformation
      use math
      type(solution_state), intent(in) :: soln, solnm1
      type(solver_settings), intent(in) :: solver
      type(solution_state), intent(out) :: solnp1
      real*8 :: Ty12(6,6)

      solnp1%U = soln%U + (soln%U - solnm1%U)
      solnp1%Udot = (soln%U - solnm1%U)/(solver%delta_t * solver%gamma) &
                    + ((solver%gamma - 1)*soln%Udot)/solver%gamma

      Ty12 = 0d0
      Ty12(1:3,1:3) = T_a12(soln%Y(1:3))
      Ty12(4:6,4:6) = transpose(T_r12(soln%Y(4:6)))
      solnp1%Y = soln%Y + (Ty12.matmul.soln%U)*solver%delta_t
      ! Old implementation without correct axis transformations
      ! solnp1%Y = soln%Y + soln%U*solver%delta_t + 0.5d0*soln%Udot*solver%delta_t**2 &
      !            + ( ((soln%Udot - soln%Udot)/solver%delta_t) * solver%delta_t**3 )/6d0
    end subroutine

    pure subroutine predictor_UdotConstant(soln, solnp1, solver)
      use types, only: solver_settings, solution_state
      use transformation
      use math
      type(solution_state), intent(in) :: soln
      type(solver_settings), intent(in) :: solver
      type(solution_state), intent(out) :: solnp1
      real*8 :: Ty12(6,6)

      solnp1%Udot = soln%Udot
      solnp1%U = soln%U + solver%delta_t*soln%Udot

      Ty12 = 0d0
      Ty12(1:3,1:3) = T_a12(soln%Y(1:3))
      Ty12(4:6,4:6) = transpose(T_r12(soln%Y(4:6)))
      solnp1%Y = soln%Y + (Ty12.matmul.soln%U)*solver%delta_t
      ! Old implementation without correct axis transformation
      ! solnp1%Y = soln%Y + soln%U*solver%delta_t + 0.5d0*soln%Udot*solver%delta_t**2
    end subroutine

    subroutine iterate(soln, solnp1, solver, problem)
      use aero
      use types, only: solver_settings, solution_state, problemData
      use transformation
      use math
      type(solution_state), intent(inout) :: soln
      type(solver_settings), intent(in) :: solver
      type(problemData), intent(in) :: problem
      type(solution_state), intent(out) :: solnp1
      ! TODO create subroutine to calculate A and B from U and Y
      !  [ ] Calculate them based on soln and solnp1
      !  [ ] Then calculate A_naf and B_naf from An, Anp1, Bn, and Bnp1
      !  [ ] Make calcA take Y and U as direct inputs rather than solution
      !  [ ] Use Tr12 and (subsequent) omega2 to generate omega_tilde
      !  [x] Make Y predictor use Tr12 since Y is in axis1 and U in axis2

      real*8, dimension(6) :: Udoti_nam, Ui_naf, Yi_naf, Gi, Ai_naf
      real*8, dimension(6,6) :: Bi_naf, Ki, Ty12_naf, dAi_nafdUi_naf
      real*8, dimension(3,3) :: omega_tilde, Imat
      real*8 :: time_factor
      real*8, dimension(3,3) :: F_factor, M_factor
      integer :: i, j, k

      time_factor = solver%alpha_m/(solver%gamma*solver%delta_t*solver%alpha_f)

      Udoti_nam = soln%Udot + solver%alpha_m*(solnp1%Udot - soln%Udot)
      Ui_naf = soln%U + solver%alpha_f*(solnp1%U - soln%U)
      Yi_naf = soln%Y + solver%alpha_f*(solnp1%Y - soln%Y)

      do i=1,solver%nsteps
        call calcA(Yi_naf, Ui_naf, problem, Ai_naf, dAi_nafdUi_naf)
        Bi_naf = 0
        Imat = 0
        forall(j=1:3) Imat(j,j) = problem%disc%I(j)
        call skewSymmetric(Ui_naf(4:6), omega_tilde)

        Bi_naf(1:3, 1:3) = omega_tilde
        Bi_naf(4:6, 4:6) = transpose(Imat).matmul.omega_tilde.matmul.Imat

        Gi = Udoti_nam - Ai_naf + (Bi_naf.matmul.Ui_naf)
        print *, 'Max Residual', maxval(abs(Gi))

        if (maxval(abs(Gi)) < solver%tolerance) then
          solnp1%U = soln%U + (Ui_naf - soln%U)/solver%alpha_f
          solnp1%Udot = soln%Udot + (Udoti_nam - soln%Udot)/solver%alpha_m
          solnp1%Y = soln%Y + (Yi_naf - soln%Y)/solver%alpha_f
          exit
        end if

        ! Create tangent matrix
        Ki = Bi_naf
        forall(j=1:6) Ki(j,j) = Ki(j,j) + time_factor


        ! invert block diagonal tangent matrix
        Ki(1:3,1:3) = matinv3(Ki(1:3,1:3))
        Ki(4:6,4:6) = matinv3(Ki(4:6,4:6))

        Ui_naf = Ui_naf - (Ki.matmul.Gi)
        Udoti_nam = (1 - solver%alpha_m/solver%gamma)*soln%Udot + time_factor*(Ui_naf - soln%U)

        Ty12_naf = 0d0
        Ty12_naf(1:3,1:3) = T_a12(soln%Y(1:3))
        Ty12_naf(4:6,4:6) = transpose(T_r12(soln%Y(4:6)))
        Yi_naf = soln%Y - (Ty12_naf.matmul.Ui_naf) * (solver%delta_t*(solver%alpha_f - 1))
        ! print *, solnp1%Y(1:3)

      end do

      print *, solnp1%Y(1:3)
    end subroutine



end module time_integrator
