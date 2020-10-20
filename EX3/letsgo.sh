#!/bin/bash
max_size=1500
gfortran MatTest.f03 -o MatTest.out -O0
./MatTest.out time_resultsO0 $max_size
gfortran MatTest.f03 -o MatTest.out -O1
./MatTest.out time_resultsO1 $max_size
gfortran MatTest.f03 -o MatTest.out -O2
./MatTest.out time_resultsO2 $max_size
gfortran MatTest.f03 -o MatTest.out -O3
./MatTest.out time_resultsO3 $max_size

gnuplot sin.gp