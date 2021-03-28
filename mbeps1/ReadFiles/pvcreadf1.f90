!-----------------------------------------------------------------------
! This program reads compressed complex periodic 1d vector data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppvcreadf1
      use pread1, only: pvcreadf1
      call pvcreadf1(-1)
      end program

