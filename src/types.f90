module types
  implicit none

  logical, parameter :: na = .false.

  type :: disc_status
    sequence

    real*8 :: loc(3)   !<Disc physical location
    real*8 :: phi(3)   !<Disc attitude
    real*8 :: vel(3)   !<Disc velocity
    real*8 :: omega(3) !<Disc angular velocity

  end type disc_status

  type :: disc_props
    sequence

    real*8 :: I(3) !<Disc moment of inertia
    real*8 :: m   !<Disc mass
    real*8 :: A !< Disc plan-view area
    real*8 :: D !< Disc diameter

  end type disc_props

  type :: environment
    real*8 :: density !< Air density
    real*8 :: viscosity !< Air dynamic viscosity
    real*8 :: wind(3) !< External Wind
    real*8 :: gravity(3) !< gravitational vector
  end type environment

  type :: problemData
    type(environment) :: env
    type(disc_props) :: disc
  end type problemData

end module types
