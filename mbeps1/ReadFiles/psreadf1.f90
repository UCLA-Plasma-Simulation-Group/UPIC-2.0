!-----------------------------------------------------------------------
! This program reads real 1d phase space data
! written for 1D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppsreadf1
      use pread1, only: psreadf1
      call psreadf1(-1)
      end program
