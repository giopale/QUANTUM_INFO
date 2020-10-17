#!/bin/bash
max_size=1500
gfortran DMatrixCODE.f03 -o DMatrixCODE.out -framework Accelerate
./DMatrixCODE.out