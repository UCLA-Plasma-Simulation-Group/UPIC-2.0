!-----------------------------------------------------------------------
! This program reads real periodic 3d vector data
! written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppvreadf3
      use pread3, only: pvreadf3
      call pvreadf3(-1)
      end program
