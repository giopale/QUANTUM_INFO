! **********************************************************************
! Project           : Quantum information, Ex1
! 
! Program name      : sum.f90
! 
! Author            : Giorgio Palermo
! 
! Date created      : 20201007
! 
! Purpose           : Test the capability of various data types
! 
! Revision History  :
!
! Date        Author      Ref    Revision (Date in YYYYMMDD format) 
!
! **********************************************************************

program summation
implicit none
integer*2 :: sum2, a2=2000000, b2=1
integer*4 :: sum4, a4=2000000, b4=1

real*4 :: pis=3.14159265359, rs=1.41421356237e21
real*4 :: sumsingle=0.
real*8 :: pid=3.14159265359, rd=1.41421356237e21
real*8 :: sumdouble=0.

sum2=a2+b2
sum4=a4+b4

print*, "		Sum of integers:		"

print*,"This is the sum of 2e6 and 1 as integer2 numbers:"
print*, sum2
print*, "This is the sum of 2e6 and 1 as integer4 numbers:"
print*, sum4

sumsingle= pis+rs
sumdouble = pid+rd

print*, "		Sum of reals:		"
print*,"This is the sum of Pi and sqrt(2)e21 in single precision:"
print*,sumsingle
print*,"This is the sum of Pi and sqrt(2)e21 in double precision:"
print*,sumdouble



print*, "End of the program"
 
 end


