! This module coordinates a comprehensive molecular minimization experiment by optimizing energy configurations
module MinimizationExperiment
    use InputOutput
    use IdentifyInteractions
    use MetropolisAlgorithm
    implicit none

contains

    ! Runs a molecular energy minimization experiment by applying the Metropolis algorithm
    subroutine RunMinimizationExperiment(inputFilename, T, r, tolerance, maxIter)
        character(len=20), intent(in)   :: inputFilename
        real(kind=8), intent(in)        :: T, r, tolerance
        integer, intent(in)             :: maxIter
        type(Atom), allocatable         :: atoms(:)
        type(Bond), allocatable         :: bonds(:)
        type(Angle), allocatable        :: angles(:)
        type(Dihedral), allocatable     :: dihedrals(:)
        integer, allocatable            :: separationCount(:,:)
        integer                         :: nAtoms, bondCount, angleCount, dihedralCount
        integer                         :: CC_BondCount, CH_BondCount
        integer                         :: CCC_Count, CCH_Count, HCH_Count
        integer                         :: CCCC_Count, CCCH_Count, HCCH_Count
        integer                         :: CC_nbCount, HC_nbCount, HH_nbCount
        integer                         :: iterTaken
        real(kind=8)                    :: i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy, i_totalEnergy
        real(kind=8)                    :: f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy, f_totalEnergy
        logical                         :: hasConverged

        ! Read molecule data from XYZ file
        call ReadMolecule(inputFilename, atoms)
        nAtoms = size(atoms)

        ! Initial identification and energy calculation
        call IdentifyBonds(atoms, bonds, nAtoms, bondCount, CC_BondCount, CH_BondCount)
        call IdentifyAngles(atoms, bonds, angles, angleCount, CCC_Count, CCH_Count, HCH_Count)
        call IdentifyDihedrals(atoms, bonds, dihedrals, dihedralCount, CCCC_Count, CCCH_Count, HCCH_Count)
        call IdentifyNonBondedSeparations(atoms, bonds, angles, separationCount, CC_nbCount, HC_nbCount, HH_nbCount)

        ! Calculate initial energy values
        i_stretchEnergy = CalculateStretchEnergy(atoms, bonds)
        i_bendEnergy = CalculateBendEnergy(atoms, angles)
        i_torsionEnergy = CalculateTorsionEnergy(atoms, dihedrals)
        i_nbEnergy = CalculateNonBondedEnergy(atoms, separationCount)
        i_totalEnergy = CalculateTotalEnergy(i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy)

        ! Minimize energy using the Metropolis algorithm
        call MetropolisMinimization(atoms, bonds, angles, dihedrals, separationCount, T, r, maxIter, &
                                    tolerance, hasConverged, iterTaken)

        ! Recalculate energies after minimization
        f_stretchEnergy = CalculateStretchEnergy(atoms, bonds)
        f_bendEnergy = CalculateBendEnergy(atoms, angles)
        f_torsionEnergy = CalculateTorsionEnergy(atoms, dihedrals)
        f_nbEnergy = CalculateNonBondedEnergy(atoms, separationCount)
        f_totalEnergy = CalculateTotalEnergy(f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy)

        ! Write updated geometry to a file
        call WriteOptimizedGeometry(inputFilename, atoms)

        ! Write output energy values to a file
        call WriteEnergyOutput(inputFilename, CC_BondCount, CH_BondCount, &
                               CCC_Count, CCH_Count, HCH_Count, &
                               CCCC_Count, CCCH_Count, HCCH_Count, &
                               CC_nbCount, HC_nbCount, HH_nbCount, &
                               i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy, i_totalEnergy, &
                               f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy, f_totalEnergy, &
                               iterTaken, maxIter, hasConverged)

        ! Clean up
        if (allocated(atoms)) then
            deallocate(atoms)
        endif
        if (allocated(bonds)) then
            deallocate(bonds)
        endif
        if (allocated(angles)) then
            deallocate(angles)
        endif
        if (allocated(dihedrals)) then
            deallocate(dihedrals)
        endif
        if (allocated(separationCount)) then
            deallocate(separationCount)
        endif
                    
    end subroutine RunMinimizationExperiment
end module MinimizationExperiment