! This module contains procedures to perform calculations between two vectors
module VectorMath
    use Types
    implicit none

contains

    ! Creates a 3D vector pointing from atom1 to atom2
    function Create3DVector(atom1, atom2) result(vector)
        type(Atom), intent(in)      :: atom1, atom2
        real(kind=8)                :: vector(3)

        vector(1) = atom2%x - atom1%x
        vector(2) = atom2%y - atom1%y
        vector(3) = atom2%z - atom1%z
    end function Create3DVector

    ! Calculates and returns the dot product of two vectors
    function DotProduct(vecA, vecB) result(dotProd)
        real(kind=8), intent(in)    :: vecA(3), vecB(3)
        real(kind=8)                :: dotProd
    
        dotProd = sum(vecA * vecB)
    end function DotProduct

    ! Calculates and returns the cross product of two vectors
    function CrossProduct(vecA, vecB) result(crossVec)
        real(kind=8), intent(in)    :: vecA(3), vecB(3)
        real(kind=8)                :: crossVec(3)
    
        crossVec(1) = vecA(2)*vecB(3) - vecA(3)*vecB(2)
        crossVec(2) = vecA(3)*vecB(1) - vecA(1)*vecB(3)
        crossVec(3) = vecA(1)*vecB(2) - vecA(2)*vecB(1)
    end function CrossProduct
end module VectorMath