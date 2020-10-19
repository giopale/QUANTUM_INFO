#!/bin/bash
gfortran DMatrixCODE.f03 -o DMatrixCODE.out -framework Accelerate
./DMatrixCODE.out