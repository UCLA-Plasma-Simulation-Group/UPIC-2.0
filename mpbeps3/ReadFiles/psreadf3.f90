!-----------------------------------------------------------------------
! This program reads real 3d phase space data
! written for 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppsreadf3
      use pread3, only: psreadf3
      call psreadf3(-1)
      end program
