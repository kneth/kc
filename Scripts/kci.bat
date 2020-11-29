@echo off
rem ***********************************************************
rem * Kinetic Compiler and Integrator                         *
rem * CopyWrong by Kenneth Geisshirt (kneth@osc.kiku.dk)      *
rem * Last updated: 4 July 1994                               *
rem ***********************************************************

go32 d:\science\chemcomp\kc\src\kc -q -m3 < %1
copy d:\science\chemcomp\kc\solvers\kksolver.c 
gcc -O2 -o _simul kksolver.c -lm
go32 _simul
del _simul
del model.c
del kksolver.c
