!-----------------------------------------------------------------------
! This program reads real periodic 3d velocity data
! written for 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppfvreadf3
      use pread3, only: pfvreadf3
      call pfvreadf3(-1)
      end program
