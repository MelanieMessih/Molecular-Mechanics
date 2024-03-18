! ========================================================================================================
! Program: OptimizeMoleculeEnergies
! Author: Melanie Messih (13362933)
! Date: 05-03-2024
! Course: Scientific Computing & Programming - week 6, 7 & 8
! Project 3 - Molecular Mechanics
! 
! Description: !!!!!!!!!!
!   This program demonstrates the application of finite differentiation to model energy functions useful 
!   for chemical systems. It utilizes three modules, EnergyFunctions, Differentiation and PrintModule, 
!   to compute and neatly print the first (gradient) and second (Hessian) derivatives for two types of 
!   potential energy models:
!   - Polyene - to model polyenes with alternating single and double bonds;
!   - Lennard-Jones - to model interactions between noble gas atoms.
! ========================================================================================================

program OptimizeMoleculeEnergy
    use MinimizationExperiment
    implicit none

    character(len=20)   :: inputFilename = 'c4h10.xyz'
    character(len=20)   :: outputFilename = 'c4h10_output.txt'
    real(kind=8)        :: T = 20, r = 1.0e-3, tolerance = 1.0e-5 ! Metropolis algorithm parameters
    integer             :: maxIter = 10000

    call RunMinimizationExperiment(inputFilename, outputFilename, T, r, tolerance, maxIter)
end program OptimizeMoleculeEnergy