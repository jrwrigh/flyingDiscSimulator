module types
  implicit none

  private
  ! logical, parameter :: na = .false.

  type, public :: solution_state
    real*8, dimension(6) :: Y !< axis1 positions and attitude
    real*8, dimension(6) :: U !< axis2 linear and angular velocity
    real*8, dimension(6) :: Udot !< axis2 linear and angular accelerations
  end type solution_state

  type, public :: disc_status
    sequence

    real*8 :: loc(3)   !<Disc physical location
    real*8 :: theta(3)   !<Disc attitude
    real*8 :: vel(3)   !<Disc velocity
    real*8 :: omega(3) !<Disc angular velocity

  end type disc_status

  type, public :: disc_props
    sequence

    real*8 :: I(3) !<Disc moment of inertia
    real*8 :: m   !<Disc mass
    real*8 :: A !< Disc plan-view area
    real*8 :: D !< Disc diameter

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
    integer :: nsteps
    real*8 :: tolerance
    real*8 :: alpha_m = 0
    real*8 :: alpha_f = 0
    real*8 :: gamma = 0
  end type solver_settings
end module types
