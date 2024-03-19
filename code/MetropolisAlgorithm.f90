! This module contains a subroutine to minimize an energy value using a Metropolis algorithm
module MetropolisAlgorithm
    use EnergyCalculations
    implicit none

contains
 
    ! Performs energy minimization using the Metropolis algorithm
    subroutine MetropolisMinimization(atoms, bonds, angles, dihedrals, separationCount, T, r, maxIter, tolerance, &
                                      hasConverged, iterTaken)
        type(Atom), allocatable, intent(inout)  :: atoms(:)
        type(Atom), allocatable                 :: oldAtoms(:)
        type(Bond), allocatable, intent(in)     :: bonds(:)
        type(Angle), allocatable, intent(in)    :: angles(:)
        type(Dihedral), allocatable, intent(in) :: dihedrals(:)
        integer, allocatable, intent(inout)     :: separationCount(:,:)
        integer, intent(in)                     :: maxIter
        real(kind=8), intent(in)                :: T, r, tolerance
        real(kind=8)                            :: i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy
        real(kind=8)                            :: f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy
        real(kind=8)                            :: i_totalEnergy, f_totalEnergy, deltaE, pA
        real(kind=8)                            :: randVec(3)
        logical, intent(out)                    :: hasConverged
        integer, intent(out)                    :: iterTaken
        integer                                 :: i, j

        ! Allocate storage for old configuration to revert if needed
        allocate(oldAtoms(size(atoms)))

        ! Calculate the total energy of the initial configuration
        i_stretchEnergy = CalculateStretchEnergy(atoms, bonds)
        i_bendEnergy = CalculateBendEnergy(atoms, angles)
        i_torsionEnergy = CalculateTorsionEnergy(atoms, dihedrals)
        i_nbEnergy = CalculateNonBondedEnergy(atoms, separationCount)
        i_totalEnergy = CalculateTotalEnergy(i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy)

        ! Iteration counter and convergence flag
        iterTaken = 0
        hasConverged = .false.

        ! Metropolis algorithm
        do i = 1, maxIter
            ! Store current configuration
            oldAtoms = atoms
            iterTaken = i
        
            ! Apply random displacements to each atom to generate a new configuration
            do j = 1, size(atoms)
                call random_number(randVec)         ! Generates random numbers in [0,1)
                randVec = 2.0d0 * randVec - 1.0d0   ! Scale to [-1,1)
                atoms(j)%x = atoms(j)%x + r * randVec(1)
                atoms(j)%y = atoms(j)%y + r * randVec(2)
                atoms(j)%z = atoms(j)%z + r * randVec(3)
            end do
        
            ! Calculate the total energy for the new configuration
            f_stretchEnergy = CalculateStretchEnergy(atoms, bonds)
            f_bendEnergy = CalculateBendEnergy(atoms, angles)
            f_torsionEnergy = CalculateTorsionEnergy(atoms, dihedrals)
            f_nbEnergy = CalculateNonBondedEnergy(atoms, separationCount)
            f_totalEnergy = CalculateTotalEnergy(f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy)
        
            ! Calculate energy difference and acceptance probability
            deltaE = f_totalEnergy - i_totalEnergy
            pA = min(1.0d0, exp(-deltaE / (kB * T)))
        
            ! Generate a random number for acceptance check
            call random_number(randVec(1))

            ! Decide whether to accept the new configuration
            if (deltaE < 0.0d0 .or. randVec(1) < pA) then
                i_totalEnergy = f_totalEnergy       ! Accept new configuration
            else
                atoms = oldAtoms                    ! Revert atoms to previous configuration
            endif

            ! Check if the algorithm has converged based on the energy difference
            if (abs(deltaE) < tolerance) then
                print *
                print *, "The Metropolis algorithm has converged after ", i, " iterations."
                hasConverged = .true.
                exit
            endif
        end do

        ! Check afterwards if it did not converge
        if (.not. hasConverged) then
            print *, "The maximum number of iterations of ", maxIter, " has been reached without convergence."
            print *, "Consider adjusting the r value or increasing the temperature."
        endif

        ! Deallocate temporary storage
        deallocate(oldAtoms)
    end subroutine MetropolisMinimization
end module MetropolisAlgorithm