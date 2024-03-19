! This module contains subroutines for reading molecular data from files and writing output
module InputOutput
    use Types, only: Atom
    implicit none

contains

    ! Reads molecular data from a file and initializes atom properties
    subroutine ReadMolecule(filename, atoms)
        use Constants
        use Types
        character(len=*), intent(in)            :: filename
        type(Atom), allocatable, intent(out)    :: atoms(:)
        integer                                 :: i, unit, nAtoms
        character(len=2)                        :: atomType

        ! Open the file for reading
        open(newunit=unit, file=filename, status='old', action='read')

        ! Read the number of atoms
        read(unit, *) nAtoms
        allocate(atoms(nAtoms))

        ! Skip the blank line
        read(unit, *)

        ! Read each atom's data
        do i = 1, nAtoms
            read(unit, *) atomType, atoms(i)%x, atoms(i)%y, atoms(i)%z
            atoms(i)%type = atomType(1:1)

            ! Assign properties based on atom type
            select case(atoms(i)%type)
            case('C')
                atoms(i)%q_nb = qC
                atoms(i)%R_nb = R_C
                atoms(i)%epsilon_nb = epsilonC
            case('H')
                atoms(i)%q_nb = qH
                atoms(i)%R_nb = R_H
                atoms(i)%epsilon_nb = epsilonH
            case default
                print *, "Encountered an unknown atom type."
            end select
        end do

        close(unit)
    end subroutine ReadMolecule

    ! Writes simulation results to a file, including interaction type counts and energy components
    subroutine WriteEnergyOutput(inputFilename, CC_bondCount, CH_bondCount, &
                           CCC_count, CCH_count, HCH_count, &
                           CCCC_count, CCCH_count, HCCH_count, &
                           CC_nbCount, HC_nbCount, HH_nbCount, &
                           i_stretchEnergy, i_bendEnergy, i_torsionEnergy, i_nbEnergy, i_totalEnergy, &
                           f_stretchEnergy, f_bendEnergy, f_torsionEnergy, f_nbEnergy, f_totalEnergy, &
                           iterationsTaken, maxIterations, hasConverged)

        character(len=*), intent(in)            :: inputFilename
        integer, intent(in)                     :: CC_bondCount, CH_bondCount
        integer, intent(in)                     :: CCC_count, CCH_count, HCH_count
        integer, intent(in)                     :: CCCC_count, CCCH_count, HCCH_count
        integer, intent(in)                     :: CC_nbCount, HC_nbCount, HH_nbCount
        real(kind=8), intent(in)                :: i_stretchEnergy, f_stretchEnergy
        real(kind=8), intent(in)                :: i_bendEnergy, f_bendEnergy
        real(kind=8), intent(in)                :: i_torsionEnergy, f_torsionEnergy
        real(kind=8), intent(in)                :: i_nbEnergy, f_nbEnergy
        real(kind=8), intent(in)                :: i_totalEnergy, f_totalEnergy
        integer, intent(in)                     :: iterationsTaken, maxIterations
        logical, intent(in)                     :: hasConverged
        integer                                 :: unit
        character(len=100)                      :: moleculeName, outputFilename

        ! Extract the molecule name
        call ExtractMoleculeName(inputFilename, moleculeName)

        ! Construct the energy output filename based on the molecule name
        outputFilename = "energy_minimization_output_" // trim(moleculeName) // ".txt"

        ! Open the file for writing
        open(newunit=unit, file=outputFilename, status='replace', action='write')

        ! Write the molecule name
        call ConvertToUpperCase(moleculeName)
        write(unit, '(A, A)') "Molecule Name: ", trim(moleculeName)

        ! Interaction counts with formatted output
        write(unit, '(A, T30, I3)') "CC bonded interactions:", CC_bondCount
        write(unit, '(A, T30, I3)') "HC bonded interactions:", CH_bondCount
        write(unit, *)
        write(unit, '(A, T30, I3)') "CCC bending interactions:", CCC_count
        write(unit, '(A, T30, I3)') "CCH bending interactions:", CCH_count
        write(unit, '(A, T30, I3)') "HCH bending interactions:", HCH_count
        write(unit, *)
        write(unit, '(A, T30, I3)') "CCCC torsional interactions:", CCCC_count
        write(unit, '(A, T30, I3)') "CCCH torsional interactions:", CCCH_count
        write(unit, '(A, T30, I3)') "HCCH torsional interactions:", HCCH_count
        write(unit, *)
        write(unit, '(A, T30, I3)') "CC nonbonded interactions:", CC_nbCount
        write(unit, '(A, T30, I3)') "HC nonbonded interactions:", HC_nbCount
        write(unit, '(A, T30, I3)') "HH nonbonded interactions:", HH_nbCount
        write(unit, *)

        ! Write the algorithm convergence information with formatted output for alignment
        write(unit, '(A, T20, I5, A, I5)') "Iterations taken:", iterationsTaken, " out of ", maxIterations
        if (hasConverged) then
            write(unit, '(A, T30, A)') "Algorithm has converged: ", " Yes"
        else
            write(unit, '(A, T30, A)') "Algorithm has converged: ", " No"
        endif

        ! Initial energy values with formatted output
        write(unit, *)
        write(unit, '(A)') "Initial Energy Values:"
        write(unit, '(A, T30, F12.4, A)') "Stretch Energy: ", i_stretchEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Bending Energy: ", i_bendEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Torsional Energy: ", i_torsionEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Nonbonded Energy: ", i_nbEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Total Energy: ", i_totalEnergy, " kcal/mol"
        write(unit, *)

        ! Final energy values with formatted output
        write(unit, '(A)') "Final Energy Values After Minimization:"
        write(unit, '(A, T30, F12.4, A)') "Stretch Energy: ", f_stretchEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Bending Energy: ", f_bendEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Torsional Energy: ", f_torsionEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Nonbonded Energy: ", f_nbEnergy, " kcal/mol"
        write(unit, '(A, T30, F12.4, A)') "Total Energy: ", f_totalEnergy, " kcal/mol"

        ! Close the file
        close(unit)
    end subroutine WriteEnergyOutput

    ! Writes the optimized geometry to an XYZ file
    subroutine WriteOptimizedGeometry(inputFilename, atoms)
        character(len=*), intent(in)            :: inputFilename
        type(Atom), allocatable, intent(in)     :: atoms(:)
        character(len=100)                      :: moleculeName, optimizedFilename
        integer                                 :: i, unit
    
        ! Extract the molecule name
        call ExtractMoleculeName(inputFilename, moleculeName)
    
        ! Construct the optimized geometry filename
        optimizedFilename = 'optimized_geometry_' // trim(moleculeName) // '.xyz'
    
        ! Open the file for writing
        open(newunit=unit, file=optimizedFilename, status='replace', action='write')
    
        ! Write the number of atoms and a blank line
        write(unit, *) size(atoms)
        write(unit, *)
    
        ! Write each atom's type and coordinates
        do i = 1, size(atoms)
            write(unit, '(A,3F12.6)') atoms(i)%type, atoms(i)%x, atoms(i)%y, atoms(i)%z
        end do
    
        ! Close the file
        close(unit)
    end subroutine WriteOptimizedGeometry

    ! Extracts the molecule's name from the input filename
    subroutine ExtractMoleculeName(inputFilename, moleculeName)
        character(len=*), intent(in)            :: inputFilename
        character(len=100), intent(out)         :: moleculeName
        integer                                 :: i
    
        ! Extract molecule name by removing the extension
        i = index(inputFilename, '.')           ! Find the dot
        if (i > 0) then
            moleculeName = inputFilename(1:i-1) ! Extract name part before the dot

        else
            moleculeName = inputFilename
        endif
    end subroutine ExtractMoleculeName

    ! Converts a string to uppercase
    subroutine ConvertToUpperCase(str)
        character(len=*), intent(inout)         :: str
        integer                                 :: i

        ! Go over each character of the input string
        do i = 1, len(trim(str))
            ! Check if dealing with lowercase character
            if (ichar(str(i:i)) >= ichar('a') .and. ichar(str(i:i)) <= ichar('z')) then
                str(i:i) = char(ichar(str(i:i)) - 32) ! Convert to uppercase using ASCII code
            endif
        end do
    end subroutine ConvertToUpperCase

end module InputOutput