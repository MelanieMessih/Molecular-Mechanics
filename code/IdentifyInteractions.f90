! This module contains procedures to identify existent bond patterns in a molecule
module IdentifyInteractions
    use HelperCalculations
    implicit none

contains

    ! Identifies bonds and counts the different types of bonds (CC, CH)
    subroutine IdentifyBonds(atoms, bonds, nAtoms, bondCount, CC_BondCount, CH_BondCount)
        type(Atom), allocatable, intent(in)         :: atoms(:)
        type(Bond), allocatable, intent(inout)      :: bonds(:)
        integer, intent(in)                         :: nAtoms
        integer, intent(out)                        :: bondCount, CC_BondCount, CH_BondCount
        integer                                     :: i, j
        real(kind=8)                                :: distance
        character(len=1)                            :: typeA, typeB

        ! Initialize counters for bond types
        bondCount = 0
        CC_BondCount = 0
        CH_BondCount = 0

        ! Iterate through atom pairs to find covalent bonds based on distance and atom types
        do i = 1, nAtoms-1
            do j = i+1, nAtoms
                distance = CalculateDistance(atoms(i), atoms(j))

                ! Determine atom types involved
                typeA = atoms(i)%type
                typeB = atoms(j)%type

                ! Update bond information and dynamic bonds allocation based on atom types
                if (typeA == 'C' .and. typeB == 'C') then
                    if (abs(distance - r0_CC) <= dist_tolerance) then
                        call UpdateBondsAllocation(bonds, bondCount)
                        call UpdateBondInfo(bonds, bondCount, i, j, kCC, r0_CC)
                        CC_BondCount = CC_BondCount + 1
                    endif
                elseif ((typeA == 'C' .and. typeB == 'H') .or. (typeA == 'H' .and. typeB == 'C')) then
                    if (abs(distance - r0_CH) <= dist_tolerance) then
                        call UpdateBondsAllocation(bonds, bondCount)
                        call UpdateBondInfo(bonds, bondCount, i, j, kCH, r0_CH)
                        CH_BondCount = CH_BondCount + 1
                    endif
                endif
            end do
        end do
    end subroutine IdentifyBonds

    ! Identifies angles and counts the different types of angles (CCC, CCH, HCH)
    subroutine IdentifyAngles(atoms, bonds, angles, angleCount, CCC_count, CCH_count, HCH_count)
        use Constants
        type(Atom), intent(in)                      :: atoms(:)
        type(Bond), intent(in)                      :: bonds(:)
        type(Angle), allocatable, intent(out)       :: angles(:)
        integer, intent(out)                        :: angleCount, CCC_count, CCH_count, HCH_count
        integer                                     :: i, bondIndex1, bondIndex2
        logical                                     :: alreadyExists
        integer                                     :: atomA, atomB, atomC
        character(len=1)                            :: typeA, typeB, typeC
    
        ! Initialize counters for angle types
        angleCount = 0
        CCC_count = 0
        CCH_count = 0
        HCH_count = 0
    
        ! Iterate through atom triplets to find angles
        do i = 1, size(atoms)
    
            ! Examine each bond pair connected to potential vertex atom i to identify angles
            do bondIndex1 = 1, size(bonds) - 1
                if (bonds(bondIndex1)%atom1 == i .or. bonds(bondIndex1)%atom2 == i) then
                    do bondIndex2 = bondIndex1 + 1, size(bonds)
                        if (bonds(bondIndex2)%atom1 == i .or. bonds(bondIndex2)%atom2 == i) then
    
                            ! Determine the atoms forming the angle
                            atomB = i               ! Vertex of the angle
                            atomA = SelectAngleAtom(bonds(bondIndex1), i)
                            atomC = SelectAngleAtom(bonds(bondIndex2), i)
    
                            ! Check if this angle is already listed to avoid adding duplicates
                            alreadyExists = CheckForDuplicateAngle(angles, angleCount, atomA, atomB, atomC) 
                            
                            if (.not. alreadyExists) then
                                ! Update dynamic angles allocation
                                call UpdateAnglesAllocation(angles, angleCount)
    
                                ! Assign the atom indices for the new angle
                                angles(angleCount)%atom1 = atomA
                                angles(angleCount)%atom2 = atomB
                                angles(angleCount)%atom3 = atomC
    
                                ! Classify angle type based on atom types
                                typeA = atoms(atomA)%type
                                typeB = atoms(atomB)%type
                                typeC = atoms(atomC)%type
    
                                ! Assign the correct k value and increment counter based on angle type
                                if (typeA == 'C' .and. typeB == 'C' .and. typeC == 'C') then
                                    call UpdateAngleInfo(angles, angleCount, atomA, atomB, atomC, kCCC)
                                    CCC_count = CCC_count + 1
                                elseif (typeA == 'H' .and. typeB == 'C' .and. typeC == 'H') then
                                    call UpdateAngleInfo(angles, angleCount, atomA, atomB, atomC, kHCH)
                                    HCH_count = HCH_count + 1
                                else
                                    call UpdateAngleInfo(angles, angleCount, atomA, atomB, atomC, kCCH)
                                    CCH_count = CCH_count + 1
                                endif
                            endif
                        endif
                    end do
                endif
            end do
        end do
    end subroutine IdentifyAngles

    ! Identifies dihedral angles and counts the different types of angles (CCCC, CCCH, HCCH)
    subroutine IdentifyDihedrals(atoms, bonds, dihedrals, dihedralCount, CCCC_count, CCCH_count, HCCH_count)
        type(Atom), intent(in)                      :: atoms(:)
        type(Bond), intent(in)                      :: bonds(:)
        type(Dihedral), allocatable, intent(out)    :: dihedrals(:)
        integer, intent(out)                        :: dihedralCount
        integer, intent(out)                        :: CCCC_count, CCCH_count, HCCH_count
        logical                                     :: bondExistsAB, bondExistsBC, bondExistsCD
        integer                                     :: i, j, k, l
        logical                                     :: alreadyExists
        character(len=1)                            :: typeA, typeB, typeC, typeD
    
        ! Initialize counters for dihedral types
        dihedralCount = 0
        CCCC_count = 0
        CCCH_count = 0
        HCCH_count = 0
    
        ! Iterate through atom quadruplets to find dihedral angles
        do i = 1, size(atoms)
            do j = 1, size(atoms)
                if (j == i) cycle
    
                do k = 1, size(atoms)
                    if (k == i .or. k == j) cycle
    
                    do l = 1, size(atoms)
                        if (l == i .or. l == j .or. l == k) cycle
    
                        ! Check for continuous i-j-k-l chain forming a dihedral angle
                        bondExistsAB = CheckIfBondExists(i, j, bonds)
                        bondExistsBC = CheckIfBondExists(j, k, bonds)
                        bondExistsCD = CheckIfBondExists(k, l, bonds)
    
                        ! If there are bonds forming a continuous chain, we have a potential dihedral
                        if (bondExistsAB .and. bondExistsBC .and. bondExistsCD) then
                            ! Check if this dihedral is already listed to avoid adding duplicates
                            alreadyExists = CheckForDuplicateDihedral(dihedrals, dihedralCount, i, j, k, l)
    
                            if (.not. alreadyExists) then
                                call UpdateDihedralsAllocation(dihedrals, dihedralCount)
    
                                ! Assign the atom indices for the new dihedral
                                dihedrals(dihedralCount)%atom1 = i
                                dihedrals(dihedralCount)%atom2 = j
                                dihedrals(dihedralCount)%atom3 = k
                                dihedrals(dihedralCount)%atom4 = l

                                ! Classify dihedral type based on atom types
                                typeA = atoms(i)%type
                                typeB = atoms(j)%type
                                typeC = atoms(k)%type
                                typeD = atoms(l)%type
    
                                ! Increment the appropriate counter based on dihedral type
                                if (typeA == 'C' .and. typeB == 'C' .and. typeC == 'C' .and. atoms(l)%type == 'C') then
                                    CCCC_count = CCCC_count + 1
                                elseif (typeA == 'H' .and. typeB == 'C' .and. typeC == 'C' .and.  typeD == 'H') then
                                    HCCH_count = HCCH_count + 1
                                else
                                    CCCH_count = CCCH_count + 1
                                endif
                            endif
                        endif
                    end do
                end do
            end do
        end do
    end subroutine IdentifyDihedrals

    ! Identifies nonbonded interactions between atom pairs and counts the different types (CC, HC, HH)
    subroutine IdentifyNonBondedSeparations(atoms, bonds, angles, separationMatrix, CC_nbCount, HC_nbCount, HH_nbCount)
        type(Atom), intent(in)                      :: atoms(:)
        type(Bond), intent(in)                      :: bonds(:)
        type(Angle), intent(in)                     :: angles(:)
        integer, allocatable, intent(out)           :: separationMatrix(:,:)
        integer, intent(out)                        :: CC_nbCount, HC_nbCount, HH_nbCount
        integer                                     :: i, j
        logical                                     :: areNonBonded
        character(len=1)                            :: typeA, typeB
    
        ! Initialize separation matrix and counters for nonbonded interaction types
        allocate(separationMatrix(size(atoms), size(atoms)))
        separationMatrix = 0
        HH_nbCount = 0
        HC_nbCount = 0
        CC_nbCount = 0
    
        ! Iterate through pairs of atoms to identify nonbonded interactions
        do i = 1, size(atoms) - 1
            do j = i + 1, size(atoms)
                ! Determine if atoms i and j are considered nonbonded
                areNonBonded = CheckIfNonBonded(i, j, bonds, angles)
    
                ! Atoms that are three or more bonds apart are considered nonbonded
                if (areNonBonded) then
                    separationMatrix(i, j) = 1
                    separationMatrix(j, i) = 1

                    ! Classify nonbonded interaction type based on atom types
                    typeA = atoms(i)%type
                    typeB = atoms(j)%type
    
                    ! Increment the appropriate counter based on the nonbonded interaction type
                    if (typeA == 'C' .and. typeB == 'C') then
                        CC_nbCount = CC_nbCount + 1
                    elseif (typeA == 'H' .and. typeB == 'H') then
                        HH_nbCount = HH_nbCount + 1
                    else
                        HC_nbCount = HC_nbCount + 1
                    endif
                endif
            end do
        end do
    end subroutine IdentifyNonBondedSeparations

    ! Updates bond information with new atom indices, force constant, and equilibrium distance
    subroutine UpdateBondInfo(bonds, bondCount, atomA, atomB, kValue, r0Value)
        type(Bond), allocatable, intent(inout)  :: bonds(:)
        integer, intent(inout)                  :: bondCount
        integer, intent(in)                     :: atomA, atomB
        real(kind=8), intent(in)                :: kValue, r0Value
    
        ! Set the atom indices for the new bond and the corresponding k and r0 values
        bonds(bondCount)%atom1 = atomA
        bonds(bondCount)%atom2 = atomB
        bonds(bondCount)%k = kValue
        bonds(bondCount)%r0 = r0Value
    end subroutine UpdateBondInfo

    ! Updates angle information with new atom indices and force constant
    subroutine UpdateAngleInfo(angles, angleCount, atomA, atomB, atomC, kValue)
        type(Angle), allocatable, intent(inout) :: angles(:)
        integer, intent(inout)                  :: angleCount
        integer, intent(in)                     :: atomA, atomB, atomC
        real(kind=8), intent(in)                :: kValue
    
        ! Set the atom indices for the new angle and the corresponding k value
        angles(angleCount)%atom1 = atomA
        angles(angleCount)%atom2 = atomB
        angles(angleCount)%atom3 = atomC
        angles(angleCount)%k = kValue
    end subroutine UpdateAngleInfo

    ! Allocates or reallocates the bonds array based on the current bond count
    subroutine UpdateBondsAllocation(bonds, bondCount)
        type(Bond), allocatable, intent(inout)      :: bonds(:)
        integer, intent(inout)                      :: bondCount
        type(Bond), allocatable                     :: newBonds(:)
    
        bondCount = bondCount + 1

        ! Check if reallocation is needed
        if (.not. allocated(bonds)) then
            allocate(bonds(1))
        else
            allocate(newBonds(bondCount))
            newBonds(:bondCount-1) = bonds          ! Copy existing bonds
            call move_alloc(from=newBonds, to=bonds)
        endif
    end subroutine UpdateBondsAllocation

    ! Allocates or reallocates the bonds array based on the current angle count
    subroutine UpdateAnglesAllocation(angles, angleCount)
        type(Angle), allocatable, intent(inout)     :: angles(:)
        integer, intent(inout)                      :: angleCount
        type(Angle), allocatable                    :: newAngles(:)
    
        angleCount = angleCount + 1

        ! Check if reallocation is needed
        if (.not. allocated(angles)) then
            allocate(angles(1))
        else
            allocate(newAngles(angleCount))
            newAngles(:angleCount-1) = angles       ! Copy existing angles
            call move_alloc(from=newAngles, to=angles)
        endif
    end subroutine UpdateAnglesAllocation

    ! Allocates or reallocates the dihedrals array based on the current dihedral count
    subroutine UpdateDihedralsAllocation(dihedrals, dihedralCount)
        type(Dihedral), allocatable, intent(inout)  :: dihedrals(:)
        integer, intent(inout)                      :: dihedralCount
        type(Dihedral), allocatable                 :: newDihedrals(:)
    
        dihedralCount = dihedralCount + 1

        ! Check if reallocation is needed
        if (.not. allocated(dihedrals)) then
            allocate(dihedrals(1))
        else
            allocate(newDihedrals(dihedralCount))
            newDihedrals(:dihedralCount-1) = dihedrals  ! Copy existing dihedrals
            call move_alloc(from=newDihedrals, to=dihedrals)
        endif
    end subroutine UpdateDihedralsAllocation

    ! Returns the index of the non-vertex atom in a bond relative to a specified vertex atom of the angle
    function SelectAngleAtom(bondParameter, vertexAtom) result(nonVertexAtom)
        type(Bond), intent(in)                      :: bondParameter
        integer, intent(in)                         :: vertexAtom
        integer                                     :: nonVertexAtom
    
        ! Check which atom in the bond is not the vertex atom and return its index
        if (bondParameter%atom1 == vertexAtom) then
            nonVertexAtom = bondParameter%atom2
        else
            nonVertexAtom = bondParameter%atom1
        endif
    end function SelectAngleAtom

    ! Checks if a given angle is already listed in the angles array
    function CheckForDuplicateAngle(angles, angleCount, atomA, atomB, atomC) result(alreadyExists)
        type(Angle), intent(in)                     :: angles(:)
        integer, intent(in)                         :: angleCount
        integer, intent(in)                         :: atomA, atomB, atomC
        integer                                     :: i
        logical                                     :: alreadyExists

        alreadyExists = .false.

        do i = 1, angleCount
            ! Check if the current angle matches the one formed by atomA, atomB, atomC in any order
            if ((angles(i)%atom1 == atomA .and. angles(i)%atom2 == atomB .and. angles(i)%atom3 == atomC) .or. &
                (angles(i)%atom1 == atomC .and. angles(i)%atom2 == atomB .and. angles(i)%atom3 == atomA)) then
                alreadyExists = .true.
                exit
            endif
        end do
    end function CheckForDuplicateAngle

    ! Checks for the existence of a dihedral angle to prevent duplicates
    function CheckForDuplicateDihedral(dihedrals, dihedralCount, atomA, atomB, atomC, atomD) result(alreadyExists)
        type(Dihedral), allocatable, intent(in)     :: dihedrals(:)
        integer, intent(in)                         :: dihedralCount
        integer, intent(in)                         :: atomA, atomB, atomC, atomD
        integer                                     :: i
        logical                                     :: alreadyExists

        alreadyExists = .false.

        do i = 1, dihedralCount
            ! Check if the current dihedral angle matches the one formed by atomA, atomB, atomC and atomD in any order
            if ((dihedrals(i)%atom1 == atomA .and. &
                dihedrals(i)%atom2 == atomB .and. &
                dihedrals(i)%atom3 == atomC .and. &
                dihedrals(i)%atom4 == atomD) .or. &
                (dihedrals(i)%atom1 == atomD .and. &
                dihedrals(i)%atom2 == atomC .and. &
                dihedrals(i)%atom3 == atomB .and. &
                dihedrals(i)%atom4 == atomA)) then
                alreadyExists = .true.
                exit
            endif
        end do
    end function CheckForDuplicateDihedral

    ! Determines if two atoms are nonbonded considering direct bonds and two-bond separations
    function CheckIfNonbonded(atom1Index, atom2Index, bonds, angles) result(areNonbonded)
        integer, intent(in)                         :: atom1Index, atom2Index
        type(Bond), intent(in)                      :: bonds(:)
        type(Angle), intent(in)                     :: angles(:)
        integer                                     :: i
        logical                                     :: areNonBonded
    
        areNonBonded = .true.
    
        ! Check if a direct bond exists between the atoms
        do i = 1, size(bonds)
            if ((bonds(i)%atom1 == atom1Index .and. bonds(i)%atom2 == atom2Index) .or. &
                (bonds(i)%atom1 == atom2Index .and. bonds(i)%atom2 == atom1Index)) then
                areNonBonded = .false.
                return
            endif
        end do
    
        ! Check if the atoms form an angle, implying a two-bond separation
        do i = 1, size(angles)
            if ((angles(i)%atom1 == atom1Index .and. angles(i)%atom3 == atom2Index) .or. &
                (angles(i)%atom3 == atom1Index .and. angles(i)%atom1 == atom2Index)) then
                areNonBonded = .false.
                return
            endif
        end do
    end function CheckIfNonbonded
end module IdentifyInteractions