#!/bin/bash
max_size=1500
ifort MatTest1.f90 -o MatTest1_intel.out -parallel -heap-arrays -O0
./MatTest1_intel.out time_results_parallel_O0 $max_size
ifort MatTest1.f90 -o MatTest1_intel.out -parallel -heap-arrays -O1
./MatTest1_intel.out time_results_parallel_O1 $max_size
ifort MatTest1.f90 -o MatTest1_intel.out -parallel -heap-arrays -O2
./MatTest1_intel.out time_results_parallel_O2 $max_size
ifort MatTest1.f90 -o MatTest1_intel.out -parallel -heap-arrays -O3
./MatTest1_intel.out time_results_parallel_O3 $max_size


gnuplot sin.gp