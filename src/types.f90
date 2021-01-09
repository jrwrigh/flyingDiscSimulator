module types
  implicit none
  private

  type, public :: solution_state
    real*8, dimension(6) :: Y !< axis1 positions and attitude
    real*8, dimension(6) :: U !< axis2 linear and angular velocity
    real*8, dimension(6) :: Udot !< axis2 linear and angular accelerations
  end type solution_state

  type, public :: solution_state_ptr
    real*8, pointer :: Y(:) !< axis1 positions and attitude
    real*8, pointer :: U(:) !< axis2 linear and angular velocity
    real*8, pointer :: Udot(:) !< axis2 linear and angular accelerations
    contains
      procedure, private :: associatePointers_solution_state_ptr
      procedure, private :: associatePointersArrays_solution_state_ptr
      generic :: setptr => associatePointersArrays_solution_state_ptr, &
                           associatePointers_solution_state_ptr
      ! procedure :: setptr => associatePointers
  end type solution_state_ptr

  type, public :: disc_status
    real*8 :: loc(3)   !<Disc physical location
    real*8 :: theta(3)   !<Disc attitude
    real*8 :: vel(3)   !<Disc velocity
    real*8 :: omega(3) !<Disc angular velocity
  end type disc_status

  type, public :: disc_props
    real*8 :: I(3) !<Disc moment of inertia
    real*8 :: m   !<Disc mass
    real*8 :: A !< Disc plan-view area
    real*8 :: D !< Disc diameter
    real*8 :: Cm_damping(3) !<Moment damping coefficients
  end type disc_props

  type, public :: environment
    real*8 :: density !< Air density
    real*8 :: viscosity !< Air dynamic viscosity
    real*8 :: wind(3) !< External Wind
    real*8 :: gravity(3) !< gravitational vector
  end type environment

  type, public :: problemData
    type(environment) :: env
    type(disc_props) :: disc
  end type problemData

  type, public :: solver_settings
    real*8 :: rho_infty
    real*8 :: delta_t
    integer :: niters
    real*8 :: tolerance
    real*8 :: alpha_m = 0
    real*8 :: alpha_f = 0
    real*8 :: gamma = 0
    contains
      procedure :: calc_alphas_gamma
  end type solver_settings

  contains
    subroutine associatePointers_solution_state_ptr(self, Y, U, Udot)
      class(solution_state_ptr), intent(inout) :: self
      real*8, intent(in), dimension(:), target :: Y, U, Udot

      self%Y => Y
      self%U => U
      self%Udot => Udot
    end subroutine

    subroutine associatePointersArrays_solution_state_ptr(self, array, j)
      class(solution_state_ptr), intent(inout) :: self
      real*8, intent(in), dimension(:,:), target :: array
      integer, intent(in) :: j

      self%Y    => array(1:6,j)
      self%U    => array(7:12,j)
      self%Udot => array(13:18,j)
    end subroutine

    pure subroutine calc_alphas_gamma(self)
      class(solver_settings), intent(inout) :: self

      self%alpha_m = (3 - self%rho_infty)/(1 + self%rho_infty)
      self%alpha_f = 1/(1 + self%rho_infty)
      self%gamma = 0.5d0 + self%alpha_m - self%alpha_f
    end subroutine
end module types
