# Molecular-Mechanics

### Scientific Computing & Programming - final project 
- University of Amsterdam & Vrije Universiteit Amsterdam
- March 2024

## Introduction to the project
...

### Requirements

This codebase is written in GNU Fortran (Rev4, Built by MSYS2 project) version 13.2.0.


### Structure

The following list describes the directories and files related to this project, and where those can be found:

- **/code**: Contains all code to perform this project. Each module is stored in a different file and the filename corresponds with the module name. The `main` program is stored in the `main.f90` file.
- **/documentation**: Contains all internal and external documentation needed to perform this project. 
  - **/overview.pdf**: Contains information about the modules, procedures, user-derived data types, dependencies, and the input and output relevant to this project.
  - **/molecular mechanics paper.pdf**: Contains the force field parameters needed for the energy calculations in this project.
  - **/project description.pdf**: Contains the project description and requirements.
- **/input & output**: Contains examples of input and output files for two different molecules: c4h10 and ch4.
  - **/{...}.xyz**: Contains the initial geometry of a molecule {...} as Cartesian coordinates, e.g. c4h10.
  - **/optimized_geometry_{...}.xyz**: Contains the final geometry of the molecule as Cartesian coordinates.
  - **/energy_minimization_output_{...}.txt**: Contains the simulation results, including the counts for the different interaction types (e.g. CC bonded interactions) and the initial and final values of the energy components. 


## Author
- Melanie Messih (13362933)

## Supervisor
- Chima Chibueze
