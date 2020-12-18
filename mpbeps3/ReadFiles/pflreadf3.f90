!-----------------------------------------------------------------------
! This program reads real periodic 3d fluid data
! written by 3D MPI/OpenMP PIC codes
! written by Viktor K. Decyk, UCLA
      program ppflreadf3
      use pread3, only: pflreadf3
      call pflreadf3(-1)
      end program
