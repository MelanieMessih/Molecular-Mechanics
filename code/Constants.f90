! This module defines physical and molecular constants for molecular simulation calculations
module Constants
    implicit none

    ! Physical constants for calculations
    real(kind=8), parameter         :: pi = 2.0 * acos(0.0)
    real(kind=8), parameter         :: kB = 1.985875e-3             ! Boltzmann constant in kcal/(mol*K)

    ! Constants for bond identification and stretch energy calculation
    real(kind=8), parameter         :: r0_CC = 1.526                ! Equilibrium C-C bond length (Å)
    real(kind=8), parameter         :: r0_CH = 1.090                ! Equilibrium C-H bond length (Å)
    real(kind=8), parameter         :: kCC = 310.0                  ! Force constant for C-C bonds (kcal / (mol * Å^2))
    real(kind=8), parameter         :: kCH = 340.0                  ! Force constant for C-H bonds (kcal / (mol * Å^2))
    real(kind=8), parameter         :: dist_tolerance = 0.1         ! Tolerance for bond length variations (Å)

    ! Constants for angles identification and bending energy calculation
    real(kind=8), parameter         :: kCCC = 40.0                  ! Force constant for CCC angles (kcal / (mol * rad^2))
    real(kind=8), parameter         :: kCCH = 50.0                  ! Force constant for CCH angles (kcal / (mol * rad^2))
    real(kind=8), parameter         :: kHCH = 35.0                  ! Force constant for HCH angles (kcal / (mol * rad^2))
    real(kind=8), parameter         :: theta0 = 109.5 * pi / 180.0  ! Equilibrium bond angle (rad)

    ! Constants for dihedrals identification and torsional energy calculation
    real(kind=8), parameter         :: halfVn = 1.40                ! Half of the barrier height (kcal / mol)
    real(kind=8), parameter         :: gamma = 0                    ! Phase shift in dihedral angle calculations (rad)
    real(kind=8), parameter         :: n = 3.0                      ! Periodicity

    ! Constants for nonbonded identification and nonbonded energy calculation
    real(kind=8), parameter         :: R_C = 1.9080                 ! Lennard-Jones C radius parameter (Å)
    real(kind=8), parameter         :: R_H = 1.4870                 ! Lennard-Jones H radius parameter (Å)
    real(kind=8), parameter         :: epsilonC = 0.1094            ! Epsilon value for C (kcal / mol)
    real(kind=8), parameter         :: epsilonH = 0.0157            ! Epsilon value for H (kcal / mol)
    real(kind=8), parameter         :: qC = -0.344                  ! Partial charge for C (e; dimensionless)
    real(kind=8), parameter         :: qH = 0.078                   ! Partial charge for H (e; dimensionless)
    real(kind=8), parameter         :: scalingFactor = 332.0636     ! To obtain desired energy unit (kcal / mol * Å^-1 * e^-2)
end module Constants