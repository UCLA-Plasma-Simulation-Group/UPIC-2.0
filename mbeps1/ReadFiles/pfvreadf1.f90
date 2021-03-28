!-----------------------------------------------------------------------
! This program reads real 1d velocity data
! written for 1D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppfvreadf1
      use pread1, only: pfvreadf1
      call pfvreadf1(-1)
      end program
