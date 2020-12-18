!-----------------------------------------------------------------------
! This program reads compressed complex periodic 2d scalar data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppcreadf2
      use pread2, only: pcreadf2
      call pcreadf2(-1)
      end program
