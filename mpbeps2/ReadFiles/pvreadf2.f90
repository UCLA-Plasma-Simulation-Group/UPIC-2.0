!-----------------------------------------------------------------------
! This program reads real periodic 2d vector data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppvreadf2
      use pread2, only: pvreadf2
      call pvreadf2(-1)
      end program
