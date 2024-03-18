! This module contains helper functions for geometric and bonding calculations
module HelperCalculations
    use VectorMath
    implicit none

contains 

    ! Calculates the Eucledian distance between two atoms
    function CalculateDistance(atomA, atomB) result(distanceAB)
        type(Atom), intent(in)  :: atomA, atomB
        real(kind=8)            :: distanceAB

        distanceAB = sqrt((atomA%x - atomB%x)**2 + (atomA%y - atomB%y)**2 + (atomA%z - atomB%z)**2)
    end function CalculateDistance

    ! Calculates the angle formed by three atoms
    function CalculateAngle(atomA, atomB, atomC) result(angleABC)
        type(Atom), intent(in)  :: atomA, atomB, atomC
        real(kind=8)            :: angleABC, normBA, normBC
        real(kind=8)            :: vecBA(3), vecBC(3)
    
        ! Calculate bond vectors
        vecBA = Create3DVector(atomB, atomA)
        vecBC = Create3DVector(atomB, atomC)
    
        ! Calculate magnitudes (norms) of vectors BA and BC
        normBA = sqrt(sum(vecBA**2))
        normBC = sqrt(sum(vecBC**2))
    
        ! Calculate the angle using the dot product formula
        angleABC = acos(DotProduct(vecBA, vecBC) / (normBA * normBC))    
    end function CalculateAngle

    ! Calculates the dihedral angle formed by four atoms
    function CalculateDihedralAngle(atomA, atomB, atomC, atomD) result(dihedralAngle)
        type(Atom), intent(in)  :: atomA, atomB, atomC, atomD
        real(kind=8)            :: dihedralAngle, norm1, norm2
        real(kind=8)            :: vecAB(3), vecBC(3), vecCD(3), crossVec1(3), crossVec2(3)
    
        ! Calculate bond vectors
        vecAB = Create3DVector(atomA, atomB)
        vecBC = Create3DVector(atomB, atomC)
        vecCD = Create3DVector(atomC, atomD)
    
        ! Calculate normal vectors to the planes defined by ABC and BCD
        crossVec1 = CrossProduct(vecAB, vecBC)
        crossVec2 = CrossProduct(vecBC, vecCD)
    
        ! Normalize the normal vectors
        norm1 = sqrt(sum(crossVec1**2))
        norm2 = sqrt(sum(crossVec2**2))
        crossVec1 = crossVec1 / norm1
        crossVec2 = crossVec2 / norm2
    
        ! Calculate the dihedral angle using the dot product formula
        dihedralAngle = acos(DotProduct(crossVec1, crossVec2))
    
        ! Adjust the angle based on the sign to get the correct orientation
        if (DotProduct(CrossProduct(crossVec1, crossVec2), vecBC) < 0.0) then
            dihedralAngle = -dihedralAngle
        endif
    end function CalculateDihedralAngle

    ! Checks if a covalent bond is known to exist between two atoms
    logical function CheckIfBondExists(atomIndex1, atomIndex2, bonds)
        integer, intent(in)     :: atomIndex1, atomIndex2
        type(Bond), intent(in)  :: bonds(:)
        integer                 :: i

        CheckIfBondExists = .false.
        do i = 1, size(bonds)
            if ((bonds(i)%atom1 == atomIndex1 .and. bonds(i)%atom2 == atomIndex2) .or. &
                (bonds(i)%atom1 == atomIndex2 .and. bonds(i)%atom2 == atomIndex1)) then
                CheckIfBondExists = .true.
                exit
            endif
        end do
    end function CheckIfBondExists
end module HelperCalculations