# Parallelizing iGEMDOCK - a fundamental virtual screening software for recognizing pharmacological interactions in drug design

## Description

Brief description of your project goes here. Explain the purpose of the program and any other relevant information.

## Getting Started

### Prerequisites

- GCC compiler for the serial version.
- CUDA Toolkit for the CUDA-accelerated version.

### Data preparation
Please prepare the .pdb file of protein cavity (Using PyMOL) and .mol files of drug molecules before docking.

Put your protein cavity in `WCavPDB` directory. (We have already prepared a pre-processed protein cavity profile in this directory : 1a4w_8A_cavity.pdb)

Put your drug molecules in `Drug` directory. (We have already prepared a few drug files in this directory : DB06699.mol, DB00080.mol, DB09099.mol)

### Compiling and Running

This project contains two versions: a serial version and a CUDA-accelerated version. Follow the instructions below to compile and run the version of your choice.

#### Serial Version
1. Navigate to the `src_serial` directory:
   ```sh
   cd src_serial
   ```
2. Prepare a list of drugs to dock, list the drugs you just prepared in the file. (ex : list.txt)

3. Compile the program using the `make` command:
   ```sh
   make
   ```
4. Execute the serial version of iGEMDOCK and you can set the population number yourself ( In this example population is 2. The higher the population, the longer the execution time. ) :
   ```sh
   ./bestdock_gemdock -lsm 2 1a4w_8A_cavity.pdb list.txt
   ```
#### CUDA Version
1. Navigate to the `src_CUDA` directory:
   ```sh
   cd src_CUDA
   ```
2. Prepare a list of drugs to dock, list the drugs you just prepared in the file. (ex : list.txt)

3. Compile the program using the `make` command:
   ```sh
   make
   ```
4. Execute the serial version of iGEMDOCK and you can set the population number yourself ( In this example population is 2. The higher the population, the longer the execution time. ) :
   ```sh
   ./bestdock_gemdock_CUDA -lsm 2 1a4w_8A_cavity.pdb list.txt
   ```

### Result
The final docking results will be stored in PrePDB/docked_Pose in the serial and CUDA version folders. (ex : 1a4w_8A_cavity-DB00014-0.pdb, 1a4w_8A_cavity-DB00014-1.pdb)

These docking results can be viewed using PyMOL to view the 3D structure results after molecular simulation.