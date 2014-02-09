
A soil-plant-atmosphere model based on MAESTRA and SPA
===================================================
  
Compiles with Intel Visual Fortran Compiler (version >10). 

To compile with gfortran, comment the following line: 
```
USE IFPORT
```
in `getmet.f90` and `inout.f90`. 
A Makefile is provided to compile Maes* on a Mac (thanks to Martin de Kauwe and Alejandro Morales).