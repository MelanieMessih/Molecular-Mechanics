! ========================================================================================================
! Program: OptimizeMoleculeEnergy
! Author: Melanie Messih (13362933)
! Date: 19-03-2024
! Course: Scientific Computing & Programming - week 6, 7 & 8
! Project 3 - Molecular Mechanics
! 
! Description:
!   This program models and minimizes the energy of a molecular structure using a parametrized energy 
!   function, focusing on molecules composed of both carbon and hydrogen atoms exclusively linked by
!   single bonds. The program employs the Metropolis algorithm, which adjusts atomic positions to find the 
!   minimum energy geometry. 
!   
!   This approach simplifies the complex task of electronic energy calculation, providing a practical 
!   technique for studying molecules through force field methods. This simulation offers insights into 
!   molecular stability and structure, both essential for understanding chemical properties and reactions.
! ========================================================================================================

program OptimizeMoleculeEnergy
    use MinimizationExperiment
    implicit none

    character(len=20)       :: inputFilename = 'ch4.xyz'
    real(kind=8), parameter :: T = 20, r = 1.0e-3, tolerance = 1.0e-5 ! Metropolis algorithm parameters
    integer                 :: maxIter = 10000

    call RunMinimizationExperiment(inputFilename, T, r, tolerance, maxIter)
end program OptimizeMoleculeEnergy