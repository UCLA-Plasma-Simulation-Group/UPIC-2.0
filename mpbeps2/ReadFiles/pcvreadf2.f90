!-----------------------------------------------------------------------
! This program reads compressed complex periodic 2d vector data
! written for 2D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppcvreadf2
      use pread2, only: pcvreadf2
      call pcvreadf2(-1)
      end program
