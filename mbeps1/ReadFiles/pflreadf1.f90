!-----------------------------------------------------------------------
! This program reads real periodic 1d fluid data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppflreadf1
      use pread1, only: pflreadf1
      call pflreadf1(-1)
      end program
