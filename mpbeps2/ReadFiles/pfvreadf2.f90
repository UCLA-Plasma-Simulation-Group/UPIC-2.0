!-----------------------------------------------------------------------
! This program reads real 2d velocity data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppfvreadf2
      use pread2, only: pfvreadf2
      call pfvreadf2(-1)
      end program
