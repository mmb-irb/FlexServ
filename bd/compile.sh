#!/bin/tcsh
gfortran -c -O -ffixed-line-length-none bdrandom.f
gfortran -O -o bd -ffixed-line-length-none bd2.f bdrandom.o
