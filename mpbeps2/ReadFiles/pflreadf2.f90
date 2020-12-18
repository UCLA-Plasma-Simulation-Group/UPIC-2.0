!-----------------------------------------------------------------------
! This program reads real periodic 2d fluid data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppflreadf2
      use pread2, only: pflreadf2
      call pflreadf2(-1)
      end program
