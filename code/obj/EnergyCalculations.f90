! This module contains procedures to calculate energy components from a molecule's structure
module EnergyCalculations
    use HelperCalculations
    implicit none

contains

    ! Calculates the stretch (bond) energy of a molecule
    function CalculateStretchEnergy(atoms, bonds) result(stretchEnergy)
        type(Atom), intent(in)                  :: atoms(:)
        type(Bond), allocatable, intent(in)     :: bonds(:)
        real(kind=8)                            :: stretchEnergy, bondLength, energy
        integer                                 :: i

        stretchEnergy = 0.0

        ! Check if there are any bonds to process
        if (allocated(bonds) .and. size(bonds) > 0) then
            do i = 1, size(bonds)
                ! Calculate the length of the identified bond
                bondLength = CalculateDistance(atoms(bonds(i)%atom1), atoms(bonds(i)%atom2))

                ! Calculate stretch energy contribution
                energy = bonds(i)%k * (bondLength - bonds(i)%r0)**2
                stretchEnergy = stretchEnergy + energy
            end do
        endif
    end function CalculateStretchEnergy

    ! Calculates the bending energy of a molecule
    function CalculateBendEnergy(atoms, angles) result(bendEnergy)
        type(Atom), intent(in)                  :: atoms(:)
        type(Angle), allocatable, intent(in)    :: angles(:)
        real(kind=8)                            :: bendEnergy, angleABC, energy
        integer                                 :: i
    
        bendEnergy = 0.0

        ! Check if there are any angles to process
        if (allocated(angles) .and. size(angles) > 0) then
            do i = 1, size(angles)
                ! Calculate the actual angle for the identified angle
                angleABC = CalculateAngle(atoms(angles(i)%atom1), atoms(angles(i)%atom2), atoms(angles(i)%atom3))
                
                ! Calculate bending energy contribution
                energy = angles(i)%k * (angleABC - theta0)**2
                bendEnergy = bendEnergy + energy
            end do
        endif
    end function CalculateBendEnergy

    ! Calculates the torsional (dihedral) energy of a molecule
    function CalculateTorsionEnergy(atoms, dihedrals) result(torsionEnergy)
        type(Atom), intent(in)                  :: atoms(:)
        type(Dihedral), allocatable, intent(in) :: dihedrals(:)
        type(Atom)                              :: atomA, atomB, atomC, atomD
        real(kind=8)                            :: torsionEnergy, dihedralAngle, energy
        integer                                 :: i
    
        torsionEnergy = 0.0

        ! Check if there are any dihedral angles to process
        if (allocated(dihedrals) .and. size(dihedrals) > 0) then
            do i = 1, size(dihedrals)
                ! Extract atom data for each atom involved
                atomA = atoms(dihedrals(i)%atom1)
                atomB = atoms(dihedrals(i)%atom2)
                atomC = atoms(dihedrals(i)%atom3)
                atomD = atoms(dihedrals(i)%atom4)

                ! Calculate the actual angle for the identified dihedral
                dihedralAngle = CalculateDihedralAngle(atomA, atomB, atomC, atomD)

                ! Calculate torsional energy contribution
                energy = 0.5 * halfVn * (1 + cos(n * dihedralAngle - gamma))
                torsionEnergy = torsionEnergy + energy
            end do
        endif
    end function CalculateTorsionEnergy

    ! Calculates the nonbonded (Van der Waals and electrostatic) energy of a molecule
    function CalculateNonBondedEnergy(atoms, separationCount) result(nonbondedEnergy)
        use Constants
        type(Atom), intent(in)                  :: atoms(:)
        integer, intent(in), allocatable        :: separationCount(:,:) ! Bonded states: 1 for nonbonded, 0 otherwise
        real(kind=8)                            :: nonbondedEnergy, vdwEnergy, elecEnergy
        real(kind=8)                            :: distance, epsilon_ij, R_ij, A_ij, B_ij
        integer                                 :: i, j
    
        nonbondedEnergy = 0.0
    
        do i = 1, size(atoms) - 1
            do j = i + 1, size(atoms)
                if (separationCount(i, j) == 1) then
                    ! Calculate the distance between the nonbonded atoms
                    distance = CalculateDistance(atoms(i), atoms(j))
    
                    ! Calculate Van der Waals energy contribution
                    R_ij = atoms(i)%R_nb + atoms(j)%R_nb
                    epsilon_ij = sqrt(atoms(i)%epsilon_nb * atoms(j)%epsilon_nb)
                    A_ij = epsilon_ij * R_ij**12
                    B_ij = 2.0 * epsilon_ij * R_ij**6
                    vdwEnergy = ((A_ij / distance**12) - (B_ij / distance**6))
    
                    ! Calculate Coulomb electrostatic energy contribution
                    elecEnergy = scalingFactor * (atoms(i)%q_nb * atoms(j)%q_nb) / distance
                    
                    ! Calculate total nonbonded energy contribution
                    nonbondedEnergy = nonbondedEnergy + vdWEnergy + elecEnergy
                endif
            end do
        end do
    end function CalculateNonBondedEnergy

    ! Calculates the total energy of a molecule by summing all energy components
    function CalculateTotalEnergy(stretchEnergy, bendEnergy, torsionEnergy, nonbondedEnergy) result(totalEnergy)
        implicit none
        real(kind=8), intent(in)            :: stretchEnergy, bendEnergy, torsionEnergy, nonbondedEnergy
        real(kind=8)                        :: totalEnergy
    
        ! Sum up all the energy components
        totalEnergy = stretchEnergy + bendEnergy + torsionEnergy + nonbondedEnergy
    end function CalculateTotalEnergy
end module EnergyCalculations