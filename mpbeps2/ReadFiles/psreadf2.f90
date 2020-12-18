!-----------------------------------------------------------------------
! This program reads real 2d phase space data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppsreadf2
      use pread2, only: psreadf2
      call psreadf2(-1)
      end program
