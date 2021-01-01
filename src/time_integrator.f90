module time_integrator
  use types, only: solver_settings
  implicit none

  type :: solution_state
    real*8, dimension(6) :: Y !< axis1 positions and attitude
    real*8, dimension(6) :: U !< axis2 linear and angular velocity
    real*8, dimension(6) :: Udot !< axis2 linear and angular accelerations
  end type solution_state

  ! type(solver_settings) :: solver

  contains

    pure subroutine calc_alphas_gamma(solver)
      type(solver_settings), intent(inout) :: solver

      solver%alpha_m = (3 - solver%rho_infty)/(1 + solver%rho_infty)
      solver%alpha_f = 1/(1 + solver%rho_infty)
      solver%gamma = 0.5d0 + solver%alpha_m - solver%alpha_f
    end subroutine

    pure subroutine predictor_deltaUConstant(soln, solnm1, solnp1, solver)
      type(solution_state), intent(in) :: soln, solnm1
      type(solver_settings), intent(in) :: solver
      type(solution_state), intent(out) :: solnp1

      solnp1%U = soln%U + (soln%U - solnm1%U)
      solnp1%Udot = (soln%U - solnm1%U)/(solver%delta_t * solver%gamma) &
                    + ((solver%gamma - 1)*soln%Udot)/solver%gamma
      solnp1%Y = soln%Y + soln%U*solver%delta_t + 0.5d0*soln%Udot*solver%delta_t**2 &
                 + ( ((soln%Udot - soln%Udot)/solver%delta_t) * solver%delta_t**3 )/6d0
    end subroutine

    pure subroutine predictor_UdotConstant(soln, solnp1, solver)
      type(solution_state), intent(in) :: soln
      type(solver_settings), intent(in) :: solver
      type(solution_state), intent(out) :: solnp1

      solnp1%Udot = soln%Udot
      solnp1%U = soln%U + solver%delta_t*soln%Udot
      solnp1%Y = soln%Y + soln%U*solver%delta_t + 0.5d0*soln%Udot*solver%delta_t**2
    end subroutine

    pure subroutine iterate(soln, solnp1, solver)
      type(solution_state), intent(inout) :: soln
      type(solver_settings), intent(in) :: solver
      type(solution_state), intent(out) :: solnp1
      ! TODO create subroutine to calculate A and B from U and Y
      !  - Calculate them based on soln and solnp1
      !  - Then calculate A_naf and B_naf from An, Anp1, Bn, and Bnp1
      !  - Make calcA take Y and U as direct inputs rather than solution
      !  - Use Tr12 and (subsequent) omega2 to generate omega_tilde

      real*8, dimension(6) :: Udoti_nam, Ui_naf, Yi_naf, Gi
      real*8 :: omega_tilde(3,3)
      integer :: i

      do i=1,solver%nsteps
        Udoti_nam = soln%Udot + solver%alpha_m*(solnp1%Udot - soln%Udot)
        Ui_naf = soln%U + solver%alpha_f*(solnp1%U - soln%U)
        Yi_naf = soln%Y + solver%alpha_f*(solnp1%Y - soln%Y)



        Tr12 = T_r12(solution%Y(4:6))
        omega2 = transpose(Tr12).matmul.solution%U(4:6)

        ! Gi = Udoti_nam - A
      end do
    end subroutine



end module time_integrator
