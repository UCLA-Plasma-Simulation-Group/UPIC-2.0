!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications with OpenMP
! npmpplib2.f90 contains basic communications procedures no MPI
! PPMAX performs parallel maximum of a real vector.
! PPBICAST broadcasts integer data from node 0
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: august 13, 2018
      module mpplib2
      use mpi
      implicit none
!
! common data for parallel processing
! lstat = length of status array
      integer, parameter :: lstat = MPI_STATUS_SIZE
! nproc = number of real or virtual processors obtained
      integer :: nproc = 1
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
      integer :: lgrp, mreal, mint, mcplx, mdouble, lworld
! msum = MPI_SUM
! mmax = MPI_MAX
      integer :: msum, mmax
      save
!
      private
      public :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      public :: PPMAX, PPBICAST
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PPMAX(f,g,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, f(j,k) = maximum as a function of k of f(j,k)
! at the end, all processors contain the same maximum.
! f = input and output real data
! g = scratch real array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! mmax = MPI_MAX
! local data
      integer j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mreal,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBICAST(if,nxp)
! this subroutine broadcasts integer data from node 0
! if = input and output integer data
! nxp = number of data values
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mint = default datatype for integers
! local data
      integer :: ierr
! return if only one processor
      if (nproc==1) return
! broadcast integer
      call MPI_BCAST(if,nxp,mint,0,lgrp,ierr)
      end subroutine
!
      end module
