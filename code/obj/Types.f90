! This module defines the types for atoms, bonds, angles, dihedrals, and nonbonded interactions
module Types
    use Constants
    implicit none

    type Atom
        character(len=1)    :: type                         ! Atom type ('C' for carbon, 'H' for hydrogen, etc.)
        real(kind=8)        :: x, y, z                      ! Cartesian coordinates (Å)
        real(kind=8)        :: R_nb                         ! Lennard-Jones radius parameter (Å)
        real(kind=8)        :: epsilon_nb                   ! Epsilon value (kcal / mol)
        real(kind=8)        :: q_nb                         ! Partial charge (e; dimensionless)
    end type Atom

    type Bond
        integer             :: atom1, atom2                 ! Indices of the two atoms forming the bond
        real(kind=8)        :: k                            ! Force constant for the bond (kcal / (mol * Å^2))
        real(kind=8)        :: r0                           ! Equilibrium bond length (Å)
    end type Bond

    type Angle
        integer             :: atom1, atom2, atom3          ! Indices of the three atoms forming the angle
        real(kind=8)        :: k                            ! Force constant for the bond (kcal / (mol * rad^2))
    end type Angle

    type Dihedral
        integer             :: atom1, atom2, atom3, atom4   ! Indices of the four atoms forming the dihedral angle
    end type Dihedral

    type NonBonded
        integer             :: atom1, atom2                 ! Indices of the two nonbonded atoms
    end type NonBonded
end module Types 