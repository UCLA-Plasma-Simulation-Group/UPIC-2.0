!-----------------------------------------------------------------------
! Null MPI library constants
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: january 31, 2017
      module mpi
      implicit none
!
      integer, parameter :: MPI_STATUS_SIZE = 8
      integer, parameter :: MPI_COMM_WORLD = 0
      integer, parameter :: MPI_INTEGER=18
      integer, parameter :: MPI_REAL=19
      integer, parameter :: MPI_DOUBLE_PRECISION=20
      integer, parameter :: MPI_COMPLEX=22
      integer, parameter :: MPI_DOUBLE_COMPLEX=23
      integer, parameter :: MPI_MAX=0
      integer, parameter :: MPI_SUM=2
      double precision, external :: MPI_WTIME
      save
!
      end module
