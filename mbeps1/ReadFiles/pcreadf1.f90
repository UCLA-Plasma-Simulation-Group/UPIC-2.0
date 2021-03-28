!-----------------------------------------------------------------------
! This program reads compressed complex periodic 1d scalar data
! written for 1D OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppcreadf1
      use pread1, only: pcreadf1
      call pcreadf1(-1)
      end program
