#!/usr/bin/bash


mpirun -np 8 ./lmp_mpi -in in.production

egrep '[0-9]{1}[.][0-9]{6}' log.lammps | grep -v [a-z] >> data

