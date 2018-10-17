!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications with OpenMP
! mpplib2.f90 contains basic communications procedures for 1d partitions
! PPINIT2 initializes parallel processing for Fortran90, returns
!         number of processors and processor id.
! PPEXIT terminates parallel processing.
! PPABORT aborts parallel processing.
! PWTIMERA performs parallel local wall clock timing.
! PPSUM performs parallel sum of a real vector.
! PPDSUM performs parallel sum of a double precision vector.
! PPMAX performs parallel maximum of a real vector.
! PPIMAX performs parallel maximum of an integer vector.
! PPDMAX performs parallel maximum of a double precision vector.
! PPDSCAN performs parallel prefix reduction of a double precision
!         vector
! PPBICAST broadcasts integer data from node 0
! PPBDCAST broadcasts double precision data from node 0
! PPISHFTR moves an integer array to the right.
! PPNCGUARD2L copies data to guard cells in y for scalar data, linear
!             interpolation, and distributed data with non-uniform
!             partition.
! PPNAGUARD2L adds guard cells in y for scalar array, linear
!             interpolation, and distributed data with non-uniform
!             partition.
! PPNACGUARD2L adds guard cells in y for vector array, linear
!              interpolation, and distributed data with non-uniform
!              partition.
! PPFMOVE2 moves fields into appropriate spatial regions, between
!          non-uniform and uniform partitions.
! PPTPOSE performs a transpose of a complex scalar array, distributed
!         in y, to a complex scalar array, distributed in x.
! PPNTPOSE performs a transpose of an n component complex vector array,
!          distributed in y, to an n component complex vector array,
!          distributed in x.
! PPMOVE2 moves particles into appropriate spatial regions with periodic
!         boundary conditions.  Assumes ihole list has been found.
! PPPMOVE2 moves particles into appropriate spatial regions for tiled
!          distributed data.
! PPWRITE2 collects distributed real 2d scalar data f and writes to a
!          direct access binary file with patial decomposition
! PPREAD2 reads real 2d scalar data f from a direct access binary file
!         and distributes it with spatial decomposition
! PPVWRITE2 collects distributed real 2d vector data f and writes to a
!           direct access binary file with spatial decomposition
! PPVREAD2 reads real 2d vector data f from a direct access binary file
!          and distributes it with spatial decomposition
! PPCWRITE2 collects distributed complex 2d scalar data f and writes to
!           a direct access binary file with patial decomposition
! PPCREAD2 reads complex 2d scalar data f from a direct access binary
!          file and distributes it with spatial decomposition
! PPVCWRITE2 collects distributed complex 2d vector data f and writes to
!            a direct access binary file with spatial decomposition
! PPVCREAD2 reads complex 2d vector data f from a direct access binary
!           file and distributes it with spatial decomposition
! PPWRPART2 collects distributed particle data part and writes to a
!           fortran unformatted file with spatial decomposition
! PPRDPART2 reads particle data part from a fortran unformatted file and
!           distributes it with spatial decomposition
! PPWRDATA2 collects distributed periodic real 2d scalar data f and
!           writes to a fortran unformatted file
! PPRDDATA2 reads periodic real 2d scalar data f from a fortran
!           unformatted file and distributes it
! PPWRVDATA2 collects distributed periodic real 2d vector data f and
!            writes to a fortran unformatted file with spatial
!            decomposition
! PPRDVDATA2 reads periodic real 2d vector data f from a fortran
!            unformatted file and distributes it with spatial
!            decomposition
! PPWRVCDATA2 collects distributed periodic complex 2d vector data f and
!             writes to a fortran unformatted file with spatial
!             decomposition
! PPRDVCDATA2 reads periodic complex 2d vector data f from a fortran
!             unformatted file and distributes it with spatial
!             decomposition
! PPARTT2 collects distributed test particle data
! PPADJFVS2 adjusts 3d velocity distribution in different regions of
!           space, so that partial regions have equal grid points
! PPWRVNDATA2 collects distributed real 2d vector non-uniform data f and
!             writes to a fortran unformatted file
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
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! msum = MPI_SUM
! mmax = MPI_MAX
      integer :: msum, mmax
      save
!
      private
      public :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      public :: PPINIT2, PPEXIT, PPABORT, PWTIMERA
      public :: PPSUM, PPDSUM, PPMAX, PPIMAX, PPDMAX, PPDSCAN
      public :: PPBICAST, PPBDCAST
      public :: PPISHFTR, PPNCGUARD2L, PPNAGUARD2L, PPNACGUARD2L
      public :: PPFMOVE2, PPTPOSE, PPNTPOSE, PPMOVE2, PPPMOVE2
      public :: PPWRITE2, PPREAD2, PPVWRITE2, PPVREAD2
      public :: PPCWRITE2, PPCREAD2, PPVCWRITE2, PPVCREAD2
      public :: PPWRPART2, PPRDPART2, PPWRDATA2, PPRDDATA2
      public :: PPWRVDATA2, PPRDVDATA2, PPWRVCDATA2, PPRDVCDATA2
      public :: PPARTT2, PPADJFVS2, PPWRVNDATA2
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
! this subroutine initializes parallel processing
! lgrp communicator = MPI_COMM_WORLD
! output: idproc, nvp
! idproc = processor id in lgrp communicator
! nvp = number of real or virtual processors obtained
      implicit none
      integer, intent(inout) :: idproc, nvp
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
! msum = MPI_SUM
! mmax = MPI_MAX
! local data
      integer :: ierror, ndprec, idprec
      integer :: iprec
      logical :: flag
      real :: prec
! ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (digits(prec) > 24) then
         ndprec = 1
      else
         ndprec = 0
      endif
! idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
      if (digits(iprec) > 31) then
         idprec = 1
      else
         idprec = 0
      endif
! this segment is used for mpi computers
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
! initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror /= 0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
! set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
! single precision real
      if (ndprec==0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
! double precision real
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
! single precision integer
!     if (idprec==0) then
!        mint = MPI_INTEGER
! double precision integer
!     else
!        mint = MPI_INTEGER8
!     endif
! operators
      msum = MPI_SUM
      mmax = MPI_MAX
      nvp = nproc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT()
! this subroutine terminates parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
! synchronize processes
         call MPI_BARRIER(lworld,ierror)
! terminate MPI execution environment
         call MPI_FINALIZE(ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT()
! this subroutine aborts parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: errorcode, ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
         errorcode = 1
! terminate MPI execution environment
         call MPI_ABORT(lworld,errorcode,ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
! this subroutine performs local wall clock timing
! input: icntrl, dtime
! icntrl = (-1,0,1) = (initialize,ignore,read) clock
! clock should be initialized before it is read!
! time = elapsed time in seconds
! dtime = current time
! written for mpi
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
! local data
      double precision :: jclock
! initialize clock
      if (icntrl==(-1)) then
         dtime = MPI_WTIME()
! read clock and write time difference from last clock initialization
      else if (icntrl==1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output real data
! g = scratch real array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! perform sum
      call MPI_ALLREDUCE(f,g,nxp,mreal,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! perform sum
       call MPI_ALLREDUCE(f,g,nxp,mdouble,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
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
      subroutine PPIMAX(if,ig,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, if(j,k) = maximum as a function of k of if(j,k)
! at the end, all processors contain the same maximum.
! if = input and output integer data
! ig = scratch integer array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mint = default datatype for integers
! mmax = MPI_MAX
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(if,ig,nxp,mint,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         if(j) = ig(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, f(j,k) = maximum as a function of k of f(j,k)
! at the end, all processors contain the same maximum.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! mmax = MPI_MAX
! local data
      integer j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSCAN(f,g,nxp)
! this subroutine performs a parallel prefix reduction of a vector,
! that is: f(j,k) = sum over k of f(j,k), where the sum is over k values
! less than idproc.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! performs a parallel prefixm sum
       call MPI_SCAN(f,g,nxp,mdouble,msum,lgrp,ierr)
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
!-----------------------------------------------------------------------
      subroutine PPBDCAST(f,nxp)
! this subroutine broadcasts double precision data from node 0
! f = input and output double precision data
! nxp = number of data values
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mdouble = default double precision type
! local data
      integer :: ierr
! return if only one processor
      if (nproc==1) return
! broadcast integer
      call MPI_BCAST(f,nxp,mdouble,0,lgrp,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPISHFTR(if,ig,nxp)
! this subroutine moves an integer array to the right
! with periodic boundary conditions
! if = input and output integer data
! ig = scratch integer array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
! lgrp = current communicator
! mint = default datatype for integers
! local data
      integer :: j, idproc, nvp, kr, kl, ltag
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! find processor id
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! find number of processors
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! return if only one processor
      if (nvp==1) return
! determine neighbors
      kr = idproc + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = idproc - 1
      if (kl.lt.0)  kl = kl + nvp
      ltag = nxp + 1
! perform right shift
      call MPI_IRECV(ig,nxp,mint,kl,ltag,lgrp,msid,ierr)
      call MPI_SEND(if,nxp,mint,kr,ltag,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! copy output from scratch array
      do 10 j = 1, nxp
      if(j) = ig(j)
   10 continue
      return
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! this subroutine copies data to guard cells in non-uniform partitions
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f
! nyp = number of primary gridpoints in field partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cell.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! special case for one processor
      if (nvp==1) then
         do j = 1, nxv
            f(j,nyp+1) = f(j,1)
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 2
! copy to guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0)  kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, nx1, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) then
         do j = 1, nx1
            f(j,1) = f(j,1) + f(j,nyp+1)
            f(j,nyp+1) = 0.
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
! add guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0) kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,ks),nxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
      do j = 1, nx1
         f(j,1) = f(j,1) + scr(j)
         f(j,nyp+1) = 0.0
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! f(ndim,j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(ndim,j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! ndim = leading dimension of array f
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) ::  nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, n, nx1, ks, moff, kl, kr
      integer :: nnxv
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) then
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,1) = f(n,j,1) + f(n,j,nyp+1)
               f(n,j,nyp+1) = 0.0
            enddo
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
      nnxv = ndim*nxv
! add guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0) kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(scr,nnxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,ks),nnxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
      do j = 1, nx1
         do n = 1, ndim
            f(n,j,1) = f(n,j,1) + scr(n,j)
            f(n,j,nyp+1) = 0.0
         enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPFMOVE2(f,g,noff,nyp,noffs,nyps,noffd,nypd,isign,kyp, &
     &ny,kstrt,nvp,nxv,nypmx,mter,ierr)
! this subroutine moves fields into appropriate spatial regions,
! between non-uniform and uniform partitions
! f(j,k) = real data for grid j,k in field partition.
! the grid is non-uniform and includes extra guard cells.
! output: f, g, ierr, and possibly mter
! g(j,k) = scratch data for grid j,k in field partition.
! noff = lowermost global gridpoint in field partition
! nyp = number of primary gridpoints in field partition
! noffs/nyps = source or scratch arrays for field partition
! noffd/nypd = destination or scratch arrays for field partition
! isign = -1, move from non-uniform (noff/nyp) to uniform (kyp) fields
! isign = 1, move from uniform (kyp) to non-uniform (noff/nyp) fields
! if isign = 0, the noffs/nyps contains the source partition, noffd/nypd
!    contains the destination partition, and noff/nyp, kyp are not used.
!    the partitions noffs/nyps and noffd/nypd are modified.
! kyp = number of complex grids in each uniform field partition.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, must be >= kyp+1
! mter = number of shifts required
! if mter = 0, then number of shifts is determined and returned
! ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer, intent(in) :: noff, nyp, isign, kyp, ny, kstrt, nvp
      integer, intent(in) :: nxv, nypmx
      integer, intent(inout) :: noffs, nyps, noffd, nypd
      integer, intent(inout) :: mter, ierr
      real, dimension(nxv,nypmx), intent(inout) :: f, g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: j, k, nbsize, ks, kyps, iter, npr, nps, nter, koff
      integer :: kl, kr, kk, nypmn, itermax
      integer, dimension(2) :: jsl, jsr, ibflg, iwork
      integer :: msid
      integer, dimension(lstat) :: istatus
! exit if certain flags are set
      if ((mter < 0).or.(nvp==1)) return
      ks = kstrt - 1
      kyps = min(kyp,max(0,ny-kyp*ks))
      nbsize = nxv*nypmx
      nypmn = nypmx
      iter = 2
      itermax = 200
      ierr = 0
! move from non-uniform to uniform fields
      koff = min(kyp*ks,ny)
! kl = processor location for grid point ny
      kl = (ny - 1)/kyp
      if (isign < 0) then
! copy non-uniform partition parameters
         noffs = noff
         nyps = nyp
         noffd = koff
         nypd = kyps
! extend partition to include ny+1 grid
! for non-uniform partition, append on last processor in y
         if (ks==(nvp-1)) nyps = nyps + 1
! for uniform partition, append just after ny grid point
         if (ks==kl) nypd = nypd + 1
         if (ks > kl) noffd = noffd + 1
! move from uniform to non-uniform fields
      else if (isign > 0) then
! set uniform partition parameters
         noffs = koff
         nyps = kyps
         noffd = noff
         nypd = nyp
! extend partition to include ny+1 grid
! for non-uniform partition, append on last processor in y
         if (ks==(nvp-1)) nypd = nypd + 1
! for uniform partition, append just after ny grid point
         if (ks==kl) nyps = nyps + 1
         if (ks > kl) noffs = noffs + 1
! move from non-uniform to non-uniform fields
      else
! extend partitions to include (ny+1) grid
         if (ks==(nvp-1)) then
            nyps = nyps + 1
            nypd = nypd + 1
         endif
      endif
! main iteration loop
! determine number of outgoing grids
   10 kl = noffd
      kr = kl + nypd
      jsl(1) = 0
      jsr(1) = 0
      do k = 1, nyps
         kk = k + noffs
! fields going right
         if (kk > kr) then
            jsr(1) = jsr(1) + 1
! fields going left
         else if (kk <= kl) then
            jsl(1) = jsl(1) + 1
         endif
      enddo
! copy fields
      iter = iter + 1
      npr = 0
      nter = 0
!
! get fields from left
      kr = ks + 1
      kl = ks - 1
      jsl(2) = 0
      jsr(2) = 0
      nps = nxv*jsr(1)
      koff = min(nyps-jsr(1)+1,nypmx)
! this segment is used for mpi computers
! post receive from left
      if (kl >= 0) then
         call MPI_IRECV(g,nbsize,mreal,kl,iter,lgrp,msid,ierr)
      endif
! send fields to right
      if (kr < nvp) then
         call MPI_SEND(f(1,koff),nps,mreal,kr,iter,lgrp,ierr)
      endif
! wait for fields to arrive
      if (kl >= 0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2) = nps/nxv
      endif
! adjust field size
      nyps = nyps - jsr(1)
! do not allow move to overflow field array
      jsr(1) = max0((nyps+jsl(2)-nypmn),0)
      if (jsr(1) > 0) then
         nyps = nyps - jsr(1)
         npr = max0(npr,jsr(1))
! save whatever is possible into end of g
         kk = min0(jsr(1),nypmn-jsl(2))
         do k = 1, kk
            do j = 1, nxv
               g(j,nypmn-kk+k) = f(j,nyps+k)
            enddo
         enddo
      endif
! shift data which is staying, if necessary
      if ((nyps > 0).and.(jsl(2) > 0)) then
         do k = 1, nyps
            kk = nyps - k + 1
            do j = 1, nxv
               f(j,kk+jsl(2)) = f(j,kk)
            enddo
         enddo
      endif
! insert data coming from left
      do k = 1, jsl(2)
         do j = 1, nxv
            f(j,k) = g(j,k)
         enddo
      enddo
! adjust field size and offset
      nyps = nyps + jsl(2)
      noffs = noffs - jsl(2)
!
! get fields from right
      kr = ks + 1
      kl = ks - 1
      nps = nxv*jsl(1)
      iter = iter + 1
! this segment is used for mpi computers
! post receive from right
      if (kr < nvp) then
         call MPI_IRECV(g,nbsize,mreal,kr,iter,lgrp,msid,ierr)
      endif
! send fields to left
      if (kl >= 0) then
         call MPI_SEND(f,nps,mreal,kl,iter,lgrp,ierr)
      endif
! wait for fields to arrive
      if (kr < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2) = nps/nxv
      endif
! adjust field size
      nyps = nyps - jsl(1)
      noffs = noffs + jsl(1)
! shift data which is staying, if necessary
      if ((nyps > 0).and.(jsl(1) > 0)) then
         do k = 1, nyps
            do j = 1, nxv
               f(j,k) = f(j,k+jsl(1))
            enddo
        enddo
      endif
! do not allow move to overflow field array
      jsl(1) = max0((nyps+jsr(2)-nypmn),0)
      if (jsl(1) > 0) then
         npr = max0(npr,jsl(1))
         jsr(2) = jsr(2) - jsl(1)
      endif
! process if no prior error
      if ((jsl(1) > 0).or.(jsr(1) <= 0)) then
! insert data coming from right
         do k = 1, jsr(2)
            do j = 1, nxv
               f(j,k+nyps) = g(j,k)
            enddo
         enddo
! adjust field size and offset
         nyps = nyps + jsr(2)
      endif
! check if new partition is uniform
      nter = nter + abs(nyps-nypd) + abs(noffs-noffd)
! calculate number of iterations
      nps = iter/2 - 1
      if (nps <= mter) then
! process errors
! should not happen, other processors do not know about this error
         if (npr /= 0) then
            ierr = npr
            write (2,*) kstrt, 'local field overflow error, ierr=', ierr
            return
         endif
! continue iteration
         if (nps < mter) go to 10
         return
      endif
! check errors, executed only first time, when mter = 0
      ibflg(1) = npr
      ibflg(2) = nter
      call PPIMAX(ibflg,iwork,2)
! field overflow error
      if (ibflg(1) /= 0) then
         ierr = ibflg(1)
         if (kstrt==1) then
            write (2,*) 'global field overflow error, ierr = ', ierr
         endif
         return
      endif
! check if any fields have to be passed further
      if (ibflg(2) > 0) then
         if (kstrt==1) then
            write (2,*) 'Info: fields being passed further = ', ibflg(2)
         endif
! continue iteration
         if (iter < itermax) go to 10
         ierr = -((iter-2)/2)
         if (kstrt==1) then
            write (2,*) 'Iteration overflow, iter = ', ierr
         endif
      endif
      mter = nps
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,  &
     &kypd)
! this subroutine performs a transpose of a matrix f, distributed in y,
! to a matrix g, distributed in x, that is,
! g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv/nyv = first dimension of f/g
! kypd/kxpd = second dimension of f/g
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      complex, dimension(nxv,kypd), intent(in) :: f
      complex, dimension(nyv,kxpd), intent(inout) :: g
      complex, dimension(kxp*kyp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyp
            do j = 1, kxp
               g(k,j) = f(j,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        do k = 1, min(kyp,max(0,ny-koff))
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              do j = 1, min(kxp,max(0,nx-joff))
!                 g(k+koff,j+joff) = f(j+joff,k+koff)
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyps
            do j = 1, ld
               s(j+ld*(k-1)) = f(j+joff,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         ld = ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ld
            do j = 1, kxps
               g(k+koff,j) = t(j+kxps*(k-1))
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv, &
     &kxpd,kypd)
! this subroutine performs a transpose of a matrix f, distributed in y,
! to a matrix g, distributed in x, that is,
! g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! ndim = leading dimension of arrays f and g
! nxv/nyv = first dimension of f/g
! kypd/kxpd = second dimension of f/g
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      complex, dimension(ndim,nxv,kypd), intent(in) :: f
      complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
      complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = ndim*kxp*kyp
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kyp
            do j = 1, kxp
               do i = 1, ndim
                  g(i,k,j) = f(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        do k = 1, min(kyp,max(0,ny-koff))
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              do j = 1, min(kxp,max(0,nx-joff))
!                 do i = 1, ndim
!                    g(i,k+koff,j+joff) = f(i,j+joff,k+koff)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id.lt.0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kyps
            do j = 1, ld
               do i = 1, ndim
                  s(i,j+ld*(k-1)) = f(i,j+joff,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         ld = ndim*ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ld
            do j = 1, kxps
               do i = 1, ndim
                  g(i,k+koff,j) = t(i,j+kxps*(k-1))
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny&
     &,kstrt,nvp,ndim,nc,idimp,npmax,idps,nbmax,ntmax,info)
! this subroutine moves particles into appropriate spatial regions
! periodic boundary conditions
! output: part, npp, sbufr, sbufl, rbufr, rbufl, info
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! edges(1:2) = lower:lower boundary of particle partition
! npp = number of particles in partition
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ihole = location of holes left in particle arrays
! ny = system length in y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! ndim = number of velocity dimensions = 2 or 3
! nc = (1,2) = (normal,tagged) partitioned co-ordinate to sort by
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition.
! idps = number of partition boundaries
! nbmax =  size of buffers for passing particles between processors
! ntmax =  size of hole array for particles leaving processors
! info = status information
! info(1) = ierr = (0,N) = (no,yes) error condition exists
! info(2) = maximum number of particles per processor
! info(3) = minimum number of particles per processor
! info(4) = maximum number of buffer overflows
! info(5) = maximum number of particle passes required
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, ndim, nc, idimp, npmax
      integer, intent(in) :: idps, nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(ntmax+1), intent(in) :: ihole
      integer, dimension(5), intent(inout) :: info
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
! iy = partitioned co-ordinate
!     integer, parameter :: iy = 2
      integer :: ierr, ks, ih, iter, nps, itg, kl, kr, j, j1, j2, i
      integer :: joff, jin, nbsize, nter, mter, itermax, iy
      real :: any, yt
      integer, dimension(4) :: msid
      integer, dimension(lstat) :: istatus
      integer, dimension(2) :: jsl, jsr, jss
      integer, dimension(4) :: ibflg, iwork
      if ((nc.lt.1).or.(nc.gt.2)) return
! determine co-ordinate to sort by
      if (nc.eq.1) then
         iy = 2
      else if (idimp.gt.(ndim+2)) then
         iy = ndim + 3
      else
         return
      endif
      any = real(ny)
      ks = kstrt - 1
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
      info(1) = 0
      info(5) = 0
      ih = ihole(1)
      joff = 1
      jin = 1
      itermax = 2000
! ih = number of particles extracted from holes
! joff = next hole location for extraction
! jss(1) = number of holes available to be filled
! jin = next hole location to be filled
! start loop
   10 mter = 0
      nps = 0
! buffer outgoing particles
      jsl(1) = 0
      jsr(1) = 0
! load particle buffers
      do j = 1, ih
         j1 = ihole(j+joff)
         yt = part(iy,j1)
! particles going down
         if (yt < edges(1)) then
            if (ks==0) yt = yt + any
            if (jsl(1) < nbmax) then
               jsl(1) = jsl(1) + 1
               do i = 1, idimp
                  sbufl(i,jsl(1)) = part(i,j1)
               enddo
               sbufl(iy,jsl(1)) = yt
            else
               nps = 1
               exit
            endif
! particles going up
         else
            if (ks==(nvp-1)) yt = yt - any
            if (jsr(1) < nbmax) then
               jsr(1) = jsr(1) + 1
               do i = 1, idimp
                  sbufr(i,jsr(1)) = part(i,j1)
               enddo
               sbufr(iy,jsr(1)) = yt
            else
               nps = 1
               exit
            endif
         endif
      enddo
      jss(1) = jsl(1) + jsr(1)
      joff = joff + jss(1)
      ih = ih - jss(1)
! check for full buffer condition
      ibflg(3) = nps
! copy particle buffers
   60 iter = iter + 2
      mter = mter + 1
! special case for one processor
      if (nvp==1) then
         jsl(2) = jsr(1)
         do j = 1, jsl(2)
            do i = 1, idimp
               rbufl(i,j) = sbufr(i,j)
            enddo
         enddo
         jsr(2) = jsl(1)
         do j = 1, jsr(2)
            do i = 1, idimp
               rbufr(i,j) = sbufl(i,j)
            enddo
         enddo
! this segment is used for mpi computers
      else
! get particles from below and above
         kr = ks + 1
         if (kr >= nvp) kr = kr - nvp
         kl = ks - 1
         if (kl < 0) kl = kl + nvp
! post receive
         itg = iter - 1
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg,lgrp,msid(1),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,iter,lgrp,msid(2),ierr)
! send particles
         jsr(1) = idimp*jsr(1)
         call MPI_ISEND(sbufr,jsr(1),mreal,kr,itg,lgrp,msid(3),ierr)
         jsl(1) = idimp*jsl(1)
         call MPI_ISEND(sbufl,jsl(1),mreal,kl,iter,lgrp,msid(4),ierr)
! wait for particles to arrive
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2) = nps/idimp
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2) = nps/idimp
      endif
! check if particles must be passed further
! check if any particles coming from above belong here
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, jsr(2)
         if (rbufr(iy,j) < edges(1)) jsl(1) = jsl(1) + 1
         if (rbufr(iy,j) >= edges(2)) jsr(1) = jsr(1) + 1
      enddo
!     if (jsr(1) /= 0) write (2,*) ks+1, 'Info: particles returning up'
! check if any particles coming from below belong here
      do j = 1, jsl(2)
         if (rbufl(iy,j) >= edges(2)) jsr(1) = jsr(1) + 1
         if (rbufl(iy,j) < edges(1)) jss(2) = jss(2) + 1
      enddo
!     if (jss(2) /= 0) write (2,*) ks+1,'Info: particles returning down'
      nps = jsl(1) + jsr(1) + jss(2)
      ibflg(2) = nps
! make sure sbufr and sbufl have been sent
      if (nvp /= 1) then
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
      if (nps==0) go to 180
! remove particles which do not belong here
! first check particles coming from above
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, jsr(2)
         yt = rbufr(iy,j)
! particles going down
         if (yt < edges(1)) then
            jsl(1) = jsl(1) + 1
            if (ks==0) yt = yt + any
            rbufr(iy,j) = yt
            do i = 1, idimp
               sbufl(i,jsl(1)) = rbufr(i,j)
            enddo
! particles going up, should not happen
         else if (yt >= edges(2)) then
            jsr(1) = jsr(1) + 1
            if (ks==(nvp-1)) yt = yt - any
            rbufr(iy,j) = yt
            do i = 1, idimp
               sbufr(i,jsr(1)) = rbufr(i,j)
            enddo
! particles staying here
         else
            jss(2) = jss(2) + 1
            do i = 1, idimp
               rbufr(i,jss(2)) = rbufr(i,j)
            enddo
         endif
      enddo
      jsr(2) = jss(2)
! next check particles coming from below
      jss(2) = 0
      do j = 1, jsl(2)
         yt = rbufl(iy,j)
! particles going up
         if (yt >= edges(2)) then
            if (jsr(1) < nbmax) then
               jsr(1) = jsr(1) + 1
               if (ks==(nvp-1)) yt = yt - any
               rbufl(iy,j) = yt
               do i = 1, idimp
                  sbufr(i,jsr(1)) = rbufl(i,j)
               enddo
            else
               jss(2) = 2*npmax
               exit
            endif
! particles going down, should not happen
         else if (yt < edges(1)) then
            if (jsl(1) < nbmax) then
               jsl(1) = jsl(1) + 1
               if (ks==0) yt = yt + any
               rbufl(iy,j) = yt
               do i = 1, idimp
                  sbufl(i,jsl(1)) = rbufl(i,j)
               enddo
            else
               jss(2) = 2*npmax
               exit
            endif
! particles staying here
         else
            jss(2) = jss(2) + 1
            do i = 1, idimp
               rbufl(i,jss(2)) = rbufl(i,j)
            enddo
         endif
      enddo
      jsl(2) = jss(2)
! check if move would overflow particle array
  180 nps = npp + jsl(2) + jsr(2) - jss(1)
      ibflg(1) = nps
      ibflg(4) = -min0(npmax,nps)
      call PPIMAX(ibflg,iwork,4)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1)
      if (ierr > npmax) then
!        write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
! distribute incoming particles from buffers
! distribute particles coming from below into holes
      jss(2) = min0(jss(1),jsl(2))
      do j = 1, jss(2)
         j1 = ihole(j+jin)
         do i = 1, idimp
            part(i,j1) = rbufl(i,j)
         enddo
      enddo
      jin = jin + jss(2)
      if (jss(1) > jsl(2)) then
         jss(2) = min0(jss(1)-jsl(2),jsr(2))
      else
         jss(2) = jsl(2) - jss(1)
      endif
      do j = 1, jss(2)
! no more particles coming from below
! distribute particles coming from above into holes
         if (jss(1) > jsl(2)) then
            j1 = ihole(j+jin)
            do i = 1, idimp
               part(i,j1) = rbufr(i,j)
            enddo
! no more holes
! distribute remaining particles from below into bottom
         else
            do i = 1, idimp
               part(i,j+npp) = rbufl(i,j+jss(1))
            enddo
         endif
      enddo
      if (jss(1) > jsl(2)) jin = jin + jss(2)
      nps = jsl(2) + jsr(2)
      if (jss(1) <= jsl(2)) then
         npp = npp + (jsl(2) - jss(1))
         jss(1) = jsl(2)
      endif
! no more holes
! distribute remaining particles from above into bottom
      jsr(2) = max0(0,nps-jss(1))
      jss(1) = jss(1) - jsl(2)
      do j = 1, jsr(2)
         do i = 1, idimp
            part(i,j+npp) = rbufr(i,j+jss(1))
         enddo
      enddo
      npp = npp + jsr(2)
! holes left over
! fill up remaining holes in particle array with particles from bottom
      if (ih.eq.0) then
         jsr(2) = max0(0,ihole(1)-jin+1)
         do j = 1, jsr(2)
            j1 = npp - j + 1
            j2 = ihole(jsr(2)-j+jin+1)
            if (j1 > j2) then
! move particle only if it is below current hole
               do i = 1, idimp
                  part(i,j2) = part(i,j1)
               enddo
            endif
         enddo
         jin = jin + jsr(2)
         npp = npp - jsr(2)
      endif
      jss(1) = 0
! check if any particles have to be passed further
      if (ibflg(3) > 0) ibflg(3) = 1
      info(5) = max0(info(5),mter)
      if (ibflg(2) > 0) then
!        write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (iter < itermax) go to 60
         ierr = -((iter-2)/2)
!        write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
! check if buffer overflowed and more particles remain to be checked
      if (ibflg(3) > 0) then
         nter = nter + 1
         go to 10
      endif
      info(4) = nter
!     if (nter > 0) then
!        write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
!     endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
! this subroutine moves particles into appropriate spatial regions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! output: rbufr, rbufl, mcll, mclr
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ncll = particle number being sent to lower processor
! nclr = particle number being sent to upper processor
! mcll = particle number being received from lower processor
! mclr = particle number being received from upper processor
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4 or 5
! nbmax =  size of buffers for passing particles between processors
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax), intent(in) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(3,mx1), intent(in) :: ncll, nclr
      integer, dimension(3,mx1), intent(inout) :: mcll, mclr
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: ierr, ks, kl, kr, i, j, jsl, jsr
      integer :: nbsize, ncsize
      integer, dimension(8) :: msid
      integer, dimension(4) :: itg
      integer, dimension(lstat) :: istatus
      data itg /3,4,5,6/
      ks = kstrt - 1
      nbsize = idimp*nbmax
      ncsize = 3*mx1
! copy particle buffers: update rbufl, rbufr, mcll, mclr
! special case for one processor
      if (nvp==1) then
         do j = 1, mx1
            do i = 1, 3
               mcll(i,j) = nclr(i,j)
            enddo
         continue
         enddo
         do j = 1, mx1
            do i = 1, 3
               mclr(i,j) = ncll(i,j)
            enddo
         enddo
         do j = 1, nclr(3,mx1)
            do i = 1, idimp
               rbufl(i,j) = sbufr(i,j)
            enddo
         enddo
         do j = 1, ncll(3,mx1)
            do i = 1, idimp
               rbufr(i,j) = sbufl(i,j)
            enddo
         enddo
! this segment is used for mpi computers
      else
! get particles from below and above
         kr = ks + 1
         if (kr >= nvp) kr = kr - nvp
         kl = ks - 1
         if (kl < 0) kl = kl + nvp
! post receives
         call MPI_IRECV(mcll,ncsize,mint,kl,itg(1),lgrp,msid(1),ierr)
         call MPI_IRECV(mclr,ncsize,mint,kr,itg(2),lgrp,msid(2),ierr)
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg(3),lgrp,msid(3),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,itg(4),lgrp,msid(4),ierr)
! send particle number offsets
         call MPI_ISEND(nclr,ncsize,mint,kr,itg(1),lgrp,msid(5),ierr)
         call MPI_ISEND(ncll,ncsize,mint,kl,itg(2),lgrp,msid(6),ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
! send particles
         jsr = idimp*nclr(3,mx1)
         call MPI_ISEND(sbufr,jsr,mreal,kr,itg(3),lgrp,msid(7),ierr)
         jsl = idimp*ncll(3,mx1)
         call MPI_ISEND(sbufl,jsl,mreal,kl,itg(4),lgrp,msid(8),ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
! make sure sbufr, sbufl, ncll, and nclr have been sent
      if (nvp /= 1) then
         do i = 1, 4
            call MPI_WAIT(msid(i+4),istatus,ierr)
         enddo
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRITE2(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec)
! this subroutine collects distributed periodic real 2d scalar data f
! and writes to a direct access binary file with spatial decomposition
! data must have a uniform partition
! f = input data to be written
! g = scratch data
! nx/ny = system length in x/y direction
! kyp = number of data values per block
! nxv = first dimension of data array f, must be >= nx
! nypmx = second dimension of data array f, must be >= kyp
! iunit = fortran unit number
! nrec = current record number for write, if nrec > 0
! input: all, output: nrec
      implicit none
      integer, intent(in) :: nx, ny, kyp, nxv, nypmx, iunit
      integer, intent(inout) :: nrec
      real, intent(in), dimension(nxv,nypmx) :: f
      real, intent(inout), dimension(nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, igo, nyp, kypp, id, i, j, k, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      nyp = nxv*kyp
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! kypp = actual size to send in y direction
      kypp = min(kyp,max(0,ny-kyp*idproc))
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
         nrec = nrec + 1
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nyp,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypp,ierr)
         kypp = kypp/nxv
         if (kypp > 0) then
            write (unit=iunit,rec=nrec) ((g(j,k),j=1,nx),k=1,kyp)
            nrec = nrec + 1
         endif
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         if (kypp==0) nyp = 0
         call MPI_SEND(f,nyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPREAD2(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec,irc)
! this subroutine reads periodic real 2d scalar data f from a direct
! access binary file and distributes it with spatial decomposition
! data must have a uniform partition
! f = output data to be read
! g = scratch data
! nx/ny = system length in x/y direction
! kyp = number of data values per block
! nxv = first dimension of data array f, must be >= nx
! nypmx = second dimension of data array f, must be >= kyp
! iunit = fortran unit number
! nrec = current record number for read, if nrec > 0
! irc = error indicator
! input: all, output: f, nrec, irc
      implicit none
      integer, intent(in) :: nx, ny, kyp, nxv, nypmx, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout), dimension(nxv,nypmx) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nyp, kypp, id, i, j, k, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nyp = nxv*kyp
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,rec=nrec,iostat=ios) ((f(j,k),j=1,nx),k=1,kyp)
!     if (ios /= 0) nrc(1) = 1
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,rec=nrec,iostat=ios) ((f(j,k),j=1,nx),k=1,kyp)
         if (ios /= 0) nrc(1) = 1
         nrec = nrec + 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
! kypp = actual size to send in y direction
            kypp = min(kyp,max(0,ny-kyp*id))
            if (kypp > 0) then
               read (unit=iunit,rec=nrec,iostat=ios) ((g(j,k),j=1,nx),  &
     &k=1,kyp)
               if (ios /= 0) then
                  if (nrc(1) /= 0) nrc(1) = i
               endif
               nrec = nrec + 1
            endif
! send data from node 0
            if (kypp==0) nyp = 0
            call MPI_SEND(g,nyp,mreal,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nyp,mreal,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVWRITE2(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec)
! this subroutine collects distributed periodic real 2d vector data f
! and writes to a direct access binary file with spatial decomposition
! data must have a uniform partition
! f = input data to be written
! g = scratch data
! nx/ny = system length in x/y direction
! kyp = number of data values per block
! ndim = first dimension of data array f
! nxv = second dimension of data array f, must be >= nx
! nypmx = third dimension of data array f, must be >= kyp
! iunit = fortran unit number
! nrec = current record number for write, if nrec > 0
! input: all, output: nrec
      implicit none
      integer, intent(in) :: nx, ny, kyp, ndim, nxv, nypmx, iunit
      integer, intent(inout) :: nrec
      real, intent(in), dimension(ndim,nxv,nypmx) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, igo, nnyp, kypp, id, i, j, k, n, nnxv
      integer :: ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      nnxv = ndim*nxv
      nnyp = nnxv*kyp
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit,rec=nrec) (((f(n,j,k),n=1,ndim),j=1,nx),k=1,kyp)
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! kypp = actual size to send in y direction
      kypp = min(kyp,max(0,ny-kyp*idproc))
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit,rec=nrec) (((f(n,j,k),n=1,ndim),j=1,nx),k=1, &
     &kyp)
         nrec = nrec + 1
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nnyp,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypp,ierr)
         kypp = kypp/nnxv
         if (kypp > 0) then
            write (unit=iunit,rec=nrec) (((g(n,j,k),n=1,ndim),j=1,nx),  &
     &k=1,kyp)
            nrec = nrec + 1
         endif
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         if (kypp==0) nnyp = 0
         call MPI_SEND(f,nnyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVREAD2(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec,irc)
! this subroutine reads periodic real 2d vector data f from a direct
! access binary file and distributes it with spatial decomposition
! data must have a uniform partition
! f = output data to be read
! g = scratch data
! nx/ny = system length in x/y direction
! kyp = number of data values per block
! ndim = first dimension of data array f
! nxv = second dimension of data array f, must be >= nx
! nypmx = third dimension of data array f, must be >= kyp
! iunit = fortran unit number
! nrec = current record number for read, if nrec > 0
! irc = error indicator
! input: all, output: f, nrec, irc
      implicit none
      integer, intent(in) :: nx, ny, kyp, ndim, nxv, nypmx, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout), dimension(ndim,nxv,nypmx) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nnyp, kypp, id, i, j, k, n, nnxv, ios
      integer :: ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nnxv = ndim*nxv
      nnyp = nnxv*kyp
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,rec=nrec,iostat=ios) (((f(n,j,k),n=1,ndim),j=1,  &
!    &nx),k=1,kyp)
!     if (ios /= 0) nrc(1) = 1
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,rec=nrec,iostat=ios) (((f(n,j,k),n=1,ndim),   &
     &j=1,nx),k=1,kyp)
         if (ios /= 0) nrc(1) = 1
         nrec = nrec + 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
! kypp = actual size to send in y direction
            kypp = min(kyp,max(0,ny-kyp*id))
            if (kypp > 0) then
               read (unit=iunit,rec=nrec,iostat=ios) (((g(n,j,k),n=1,   &
     &ndim),j=1,nx),k=1,kyp)
               if (ios /= 0) then
                  if (nrc(1) /= 0) nrc(1) = i
               endif
               nrec = nrec + 1
            endif
! send data from node 0
            if (kypp==0) nnyp = 0
            call MPI_SEND(g,nnyp,mreal,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nnyp,mreal,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPCWRITE2(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec)
! this subroutine collects distributed complex 2d scalar data f and
! writes to a direct access binary file with spatial decomposition
! data must have a uniform partition
! f = input data to be written
! g = scratch data
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! nyv = first dimension of data array f, must be >= ny
! kxpd = second dimension of data array f, must be >= kxp
! iunit = fortran unit number
! nrec = current record number for write, if nrec > 0
! input: all, output: nrec
      implicit none
      integer, intent(in) :: nx, ny, kxp, nyv, kxpd, iunit
      integer, intent(inout) :: nrec
      complex, intent(in), dimension(nyv,kxpd) :: f
      complex, intent(inout), dimension(nyv,kxpd) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, igo, mxp, nxp, kxpp, id, i, j, k, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      mxp = min(nx,kxp)
      nxp = nyv*mxp
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit,rec=nrec) ((f(k,j),k=1,ny),j=1,mxp)
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! kxpp = actual size to send in x direction
      kxpp = min(mxp,max(0,nx-kxp*idproc))
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         if (kxpd <= kxp) then
            write (unit=iunit,rec=nrec) ((f(k,j),k=1,ny),j=1,mxp)
! mode kx = NX/2 is stored at location kxp+1 when idproc=0
         else
            write (unit=iunit,rec=nrec) ((f(k,j),k=1,ny),j=1,kxp+1)
         endif
         nrec = nrec + 1
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nxp,mcplx,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
         kxpp = kxpp/nyv
         if (kxpp > 0) then
            write (unit=iunit,rec=nrec) ((g(k,j),k=1,ny),j=1,mxp)
            nrec = nrec + 1
         endif
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         if (kxpp==0) nxp = 0
         call MPI_SEND(f,nxp,mcplx,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPCREAD2(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec,irc)
! this subroutine reads complex 2d scalar data f from a direct access
! binary file and distributes it with spatial decomposition
! data must have a uniform partition
! f = output data to be read
! g = scratch data
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! nyv = first dimension of data array f, must be >= ny
! kxpd = second dimension of data array f, must be >= kxp
! iunit = fortran unit number
! nrec = current record number for read, if nrec > 0
! irc = error indicator
! input: all, output: f, nrec, irc
      implicit none
      integer, intent(in) :: nx, ny, kxp, nyv, kxpd, iunit
      integer, intent(inout) :: nrec, irc
      complex, intent(inout), dimension(nyv,kxpd) :: f, g
! lgrp = current communicator
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, mxp, nxp, kxpp, id, i, j, k, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      mxp = min(nx,kxp)
      nxp = nyv*kxp
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,rec=nrec,iostat=ios) ((f(k,j),k=1,ny),j=1,mxp)
!     if (ios /= 0) nrc(1) = 1
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         if (kxpd <= kxp) then
            read (unit=iunit,rec=nrec,iostat=ios) ((f(k,j),k=1,ny),j=1, &
     &mxp)
! mode kx = NX/2 is stored at location kxp+1 when idproc=0
         else
            read (unit=iunit,rec=nrec,iostat=ios) ((f(k,j),k=1,ny),j=1, &
     &kxp+1)
         endif
         if (ios /= 0) nrc(1) = 1
         nrec = nrec + 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
! kxpp = actual size to send in x direction
            kxpp = min(mxp,max(0,nx-kxp*id))
            if (kxpp > 0) then
               read (unit=iunit,rec=nrec,iostat=ios) ((g(k,j),k=1,ny),  &
     &j=1,mxp)
               if (ios /= 0) then
                  if (nrc(1) /= 0) nrc(1) = i
               endif
               nrec = nrec + 1
            endif
! send data from node 0
            if (kxpp==0) nxp = 0
            call MPI_SEND(g,nxp,mcplx,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nxp,mcplx,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVCWRITE2(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec)
! this subroutine collects distributed complex 2d vector data f and
! writes to a direct access binary file with spatial decomposition
! data must have a uniform partition
! f = input data to be written
! g = scratch data
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! ndim = first dimension of data array f
! nyv = second dimension of data array f, must be >= ny
! kxpd = third dimension of data array f, must be >= kxp
! iunit = fortran unit number
! nrec = current record number for write, if nrec > 0
! input: all, output: nrec
      implicit none
      integer, intent(in) :: nx, ny, kxp, ndim, nyv, kxpd, iunit
      integer, intent(inout) :: nrec
      complex, intent(in), dimension(ndim,nyv,kxpd) :: f
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, igo, mxp, nnxp, kxpp, id, i, j, k, n, nnyv
      integer :: ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      mxp = min(nx,kxp)
      nnyv = ndim*nyv
      nnxp = nnyv*mxp
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit,rec=nrec) (((f(n,k,j),n=1,ndim),k=1,ny),j=1,mxp)
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! kxpp = actual size to send in x direction
      kxpp = min(mxp,max(0,nx-kxp*idproc))
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         if (kxpd <= kxp) then
            write (unit=iunit,rec=nrec) (((f(n,k,j),n=1,ndim),k=1,ny),  &
     &j=1,mxp)
! mode kx = NX/2 is stored at location kxp+1 when idproc=0
         else
            write (unit=iunit,rec=nrec) (((f(n,k,j),n=1,ndim),k=1,ny),  &
     &j=1,kxp+1)
         endif
         nrec = nrec + 1
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nnxp,mcplx,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
         kxpp = kxpp/nnyv
         if (kxpp > 0) then
            write (unit=iunit,rec=nrec) (((g(n,k,j),n=1,ndim),k=1,ny),  &
     &j=1,mxp)
            nrec = nrec + 1
         endif
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         if (kxpp==0) nnxp = 0
         call MPI_SEND(f,nnxp,mcplx,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVCREAD2(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec,irc)
! this subroutine reads complex 2d vector data f from a direct access
! binary file and distributes it with spatial decomposition
! data must have a uniform partition
! f = output data to be read
! g = scratch data
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! ndim = first dimension of data array f
! nyv = second dimension of data array f, must be >= ny
! kxpd = third dimension of data array f, must be >= kxp
! iunit = fortran unit number
! nrec = current record number for read, if nrec > 0
! irc = error indicator
! input: all, output: f, nrec, irc
      implicit none
      integer, intent(in) :: nx, ny, kxp, ndim, nyv, kxpd, iunit
      integer, intent(inout) :: nrec, irc
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: f, g
! lgrp = current communicator
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, mxp, nnxp, kxpp, id, i, j, k, n, nnyv, ios
      integer :: ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      mxp = min(nx,kxp)
      nnyv = ndim*nyv
      nnxp = nnyv*mxp
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,rec=nrec,iostat=ios) (((f(n,k,j),n=1,ndim),k=1,  &
!    &ny),j=1,mxp)
!     if (ios /= 0) nrc(1) = 1
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         if (kxpd <= kxp) then
            read (unit=iunit,rec=nrec,iostat=ios) (((f(n,k,j),n=1,ndim),&
     &k=1,ny),j=1,mxp)
! mode kx = NX/2 is stored at location kxp+1 when idproc=0
         else
            read (unit=iunit,rec=nrec,iostat=ios) (((f(n,k,j),n=1,ndim),&
     &k=1,ny),j=1,kxp+1)
         endif
         if (ios /= 0) nrc(1) = 1
         nrec = nrec + 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
! kxpp = actual size to send in x direction
            kxpp = min(mxp,max(0,nx-kxp*id))
            if (kxpp > 0) then
               read (unit=iunit,rec=nrec,iostat=ios) (((g(n,k,j),n=1,   &
     &ndim),k=1,ny),j=1,mxp)
               if (ios /= 0) then
                  if (nrc(1) /= 0) nrc(1) = i
               endif
               nrec = nrec + 1
            endif
! send data from node 0
            if (kxpp==0) nnxp = 0
            call MPI_SEND(g,nnxp,mcplx,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nnxp,mcplx,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRPART2(part,npp,idimp,npmax,iunit,iscr)
! this subroutine collects distributed particle data part and writes to
! a fortran unformatted file with spatial decomposition
! part = input data to be written
! npp = number of particles in partition
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! iunit = fortran unit number
! iscr = unused unit number available for scratch file
      implicit none
      integer, intent(in) :: npp, idimp, npmax, iunit, iscr
      real, intent(inout), dimension(idimp,npmax) :: part
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, igo, ndp, id, i, j, k, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit) ((part(j,k),j=1,idimp),k=1,npp)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write particle data for node 0 to scratch array
         if (nvp > 1) then
            open(unit=iscr,form='unformatted',status='scratch')
            if (npp > 0) then
               write (iscr) ((part(j,k),j=1,idimp),k=1,npp)
            endif
         endif
! write particle data for node 0 to restart file
         write (unit=iunit) npp
         if (npp > 0) then
            write (unit=iunit)  ((part(j,k),j=1,idimp),k=1,npp)
         endif
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(part,idimp*npmax,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,ndp,ierr)
         ndp = ndp/idimp
         write (unit=iunit) ndp
         if (ndp > 0) then
            write (unit=iunit) ((part(j,k),j=1,idimp),k=1,ndp)
         endif
         enddo
! read back particle data for node 0 from scratch array
         if (nvp > 1) then
            rewind iscr
            if (npp > 0) then
               read (iscr) ((part(j,k),j=1,idimp),k=1,npp)
            endif
            close (iscr)
         endif
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         ndp = idimp*npp
         call MPI_SEND(part,ndp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDPART2(part,npp,idimp,npmax,iunit,iscr,irc)
! this subroutine reads particle data part from a fortran unformatted
! file and distributes it with spatial decomposition
! part = output data to be read
! npp = number of particles in partition
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! iunit = fortran unit number
! iscr = unused unit number available for scratch file
! irc = error indicator
! input: all except part, npp, irc, output: part, npp, irc
      implicit none
      integer, intent(in) :: idimp, npmax, iunit, iscr
      integer, intent(inout) :: npp, irc
      real, intent(inout), dimension(idimp,npmax) :: part
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, ndp, id, i, j, k, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,iostat=ios) (((part(j,k),j=1,idimp),k=1,npp)
!     if (ios /= 0) nrc(1) = 1
!     nrec = nrec + 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,iostat=ios) npp
         if (ios /= 0) nrc(1) = 1
         if (npp > 0) then
            read (unit=iunit,iostat=ios) ((part(j,k),j=1,idimp),k=1,npp)
         endif
         if (ios /= 0) nrc(1) = 1
! then write particle data for node 0 to scratch array
         if (nvp > 1) then
            open(unit=iscr,form='unformatted',status='scratch')
            if (npp > 0) then
               write (iscr) ((part(j,k),j=1,idimp),k=1,npp)
            endif
         endif
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
            read (unit=iunit,iostat=ios) ndp
            if (ios /= 0) nrc(1) = 1
            if (ndp > 0) then
               read (unit=iunit,iostat=ios) ((part(j,k),j=1,idimp),     &
     &k=1,ndp)
               if (ios /= 0) then
                  if (nrc(1) /= 0) nrc(1) = i
               endif
            endif
! send data from node 0
            ndp = idimp*ndp
            call MPI_SEND(part,ndp,mreal,id,98,lgrp,ierr)
         enddo
! read back particle data for node 0 from scratch array
         if (nvp > 1) then
            rewind iscr
            if (npp > 0) then
               read (iscr) ((part(j,k),j=1,idimp),k=1,npp)
            endif
            close (iscr)
         endif
! other nodes receive data from node 0
      else 
         call MPI_RECV(part,idimp*npmax,mreal,0,98,lgrp,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,ndp,ierr)
         npp = ndp/idimp
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRDATA2(f,g,nxv,nypmx,iunit)
! this subroutine collects distributed periodic real 2d scalar data f
! and writes to a fortran unformatted file with spatial decomposition
! f = input data to be written
! g = scratch data
! nxv = first dimension of data array f
! nypmx = second dimension of data array f
! iunit = fortran unit number
! input: all
      implicit none
      integer, intent(in) :: nxv, nypmx, iunit
      real, intent(in), dimension(nxv,nypmx) :: f
      real, intent(inout), dimension(nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, igo, nyp, id, i, j, k, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      nyp = nxv*nypmx
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit) ((f(j,k),j=1,nxv),k=1,nypmx)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit) ((f(j,k),j=1,nxv),k=1,nypmx)
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nyp,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         write (unit=iunit) ((g(j,k),j=1,nxv),k=1,nypmx)
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_SEND(f,nyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDDATA2(f,g,nxv,nypmx,iunit,irc)
! this subroutine reads periodic real 2d scalar data f from a fortran
! unformatted file and distributes it with spatial decomposition
! f = output data to be read
! g = scratch data
! nxv = first dimension of data array f
! nypmx = second dimension of data array f
! iunit = fortran unit number
! irc = error indicator
! input: all, output: f, irc
      implicit none
      integer, intent(in) :: nxv, nypmx, iunit
      integer, intent(inout) :: irc
      real, intent(inout), dimension(nxv,nypmx) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nyp, id, i, j, k, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nyp = nxv*nypmx
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,iostat=ios) ((f(j,k),j=1,nxv),k=1,nypmx)
!     if (ios /= 0) nrc(1) = 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,iostat=ios) ((f(j,k),j=1,nxv),k=1,nypmx)
         if (ios /= 0) nrc(1) = 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
            read (unit=iunit,iostat=ios) ((g(j,k),j=1,nxv),k=1,nypmx)
            if (ios /= 0) then
               if (nrc(1) /= 0) nrc(1) = i
            endif
! send data from node 0
            call MPI_SEND(g,nyp,mreal,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nyp,mreal,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVDATA2(f,g,ndim,nxv,nypmx,iunit)
! this subroutine collects distributed periodic real 2d vector data f
! and writes to a fortran unformatted file with spatial decomposition
! f = input data to be written
! g = scratch data
! ndim = first dimension of data array f
! nxv = second dimension of data array f
! nypmx = third dimension of data array f
! iunit = fortran unit number
! input: all
      implicit none
      integer, intent(in) :: ndim, nxv, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nypmx) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, igo, nnyp, id, i, j, k, n, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      nnyp = ndim*nxv*nypmx
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nypmx)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nypmx)
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nnyp,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         write (unit=iunit) (((g(n,j,k),n=1,ndim),j=1,nxv),k=1,nypmx)
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_SEND(f,nnyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDVDATA2(f,g,ndim,nxv,nypmx,iunit,irc)
! this subroutine reads periodic real 2d vector data f from a fortran
! unformatted file and distributes it with spatial decomposition
! f = output data to be read
! g = scratch data
! ndim = first dimension of data array f
! nxv = second dimension of data array f
! nypmx = third dimension of data array f
! iunit = fortran unit number
! irc = error indicator
! input: all, output: f, irc
      implicit none
      integer, intent(in) :: ndim, nxv, nypmx, iunit
      integer, intent(inout) :: irc
      real, intent(inout), dimension(ndim,nxv,nypmx) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nnyp, id, i, j, k, n, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nnyp = ndim*nxv*nypmx
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,iostat=ios) (((f(n,j,k),n=1,ndim),j=1,nxv),      &
!    &k=1,nypmx)
!     if (ios /= 0) nrc(1) = 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,iostat=ios) (((f(n,j,k),n=1,ndim),j=1,nxv),   &
     &k=1,nypmx)
         if (ios /= 0) nrc(1) = 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
            read (unit=iunit,iostat=ios) (((g(n,j,k),n=1,ndim),j=1,nxv),&
     &k=1,nypmx)
           if (ios /= 0) then
               if (nrc(1) /= 0) nrc(1) = i
            endif
! send data from node 0
            call MPI_SEND(g,nnyp,mreal,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nnyp,mreal,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVCDATA2(f,g,ndim,nyv,kxpd,iunit)
! this subroutine collects distributed periodic complex 2d vector data f
! and writes to a fortran unformatted file with spatial decomposition
! f = input data to be written
! g = scratch data
! ndim = first dimension of data array f
! nyv = second dimension of data array f
! kxpd = third dimension of data array f
! iunit = fortran unit number
! input: all
      implicit none
      integer, intent(in) :: ndim, nyv, kxpd, iunit
      complex, intent(in), dimension(ndim,nyv,kxpd) :: f
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, igo, nnxp, id, i, j, k, n, ierr
      integer :: msid
      integer, dimension(lstat) :: istatus
      nnxp = ndim*nyv*kxpd
      igo = 1
! this segment is used for shared memory computers
!     write (unit=iunit) (((f(n,k,j),n=1,ndim),k=1,nyv),j=1,kxpd)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit) (((f(n,k,j),n=1,ndim),k=1,nyv),j=1,kxpd)
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
! send go signal to sending node
         call MPI_SEND(igo,1,mint,id,98,lgrp,ierr)
         call MPI_IRECV(g,nnxp,mcplx,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         write (unit=iunit) (((g(n,k,j),n=1,ndim),k=1,nyv),j=1,kxpd)
         enddo
! other nodes send data to node 0 after receiving go signal
      else
         call MPI_IRECV(igo,1,mint,0,98,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_SEND(f,nnxp,mcplx,0,99,lgrp,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDVCDATA2(f,g,ndim,nyv,kxpd,iunit,irc)
! this subroutine reads periodic complex 2d vector data f from a fortran
! unformatted file and distributes it with spatial decomposition
! f = output data to be read
! g = scratch data
! ndim = first dimension of data array f
! nyv = second dimension of data array f
! kxpd = third dimension of data array f
! iunit = fortran unit number
! irc = error indicator
! input: all, output: f, irc
      implicit none
      integer, intent(in) :: ndim, nyv, kxpd, iunit
      integer, intent(inout) :: irc
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: f, g
! lgrp = current communicator
! mcplx = default datatype for complex type
! local data
      integer :: nvp, idproc, nnxp, id, i, j, k, n, ios, ierr
      integer, dimension(1) :: nrc, iwrk1
      integer, dimension(lstat) :: istatus
      nnxp = ndim*nyv*kxpd
      nrc(1) = 0
! this segment is used for shared memory computers
!     read (unit=iunit,iostat=ios) (((f(n,k,j),n=1,ndim),k=1,nyv),      &
!    &j=1,kxpd)
!     if (ios /= 0) nrc(1) = 1
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first read data for node 0
         read (unit=iunit,iostat=ios) (((f(n,k,j),n=1,ndim),k=1,nyv),   &
     &j=1,kxpd)
         if (ios /= 0) nrc(1) = 1
! then read data on node 0 to send to remaining nodes
         do i = 2, nvp
            id = i - 1
            read (unit=iunit,iostat=ios) (((g(n,k,j),n=1,ndim),k=1,nyv),&
     &j=1,kxpd)
            if (ios /= 0) then
               if (nrc(1) /= 0) nrc(1) = i
            endif
! send data from node 0
            call MPI_SEND(g,nnxp,mcplx,id,98,lgrp,ierr)
         enddo
! other nodes receive data from node 0
      else 
         call MPI_RECV(f,nnxp,mcplx,0,98,lgrp,istatus,ierr)
      endif
! check for error condition
      call PPIMAX(nrc,iwrk1,1)
      irc = nrc(1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPARTT2(partt,numtp,idimp,nprobt,irc)
! this subroutine collects distributed test particle data
! partt = tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 5 or 6
! nprobt = number of test charges whose trajectories will be stored.
! irc = error indicator
! input: all, output: partt, irc
      implicit none
      integer, intent(in) :: numtp, idimp, nprobt
      integer, intent(inout) :: irc
      real, intent(inout), dimension(idimp,nprobt) :: partt
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, noff, npbt, nntp, i, id, ltag, ierr
      integer :: msid
      integer, dimension(1) :: iwork
      integer, dimension(lstat) :: istatus
      nntp = idimp*numtp
      ltag = nprobt
      irc = 0
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
         noff = numtp + 1
         npbt = idimp*(nprobt - numtp)
         do i = 2, nvp
         id = i - 1
         call MPI_IRECV(partt(1,noff),npbt,mreal,id,ltag,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nntp,ierr)
         nntp = nntp/idimp
         noff = noff + nntp
         npbt = npbt - nntp
         enddo
! incorrect number of trajectories received
         noff = noff - 1
         if (noff.ne.nprobt) irc = noff - 1
         iwork(1) = irc
         call PPBICAST(iwork,1)
         irc = iwork(1)
! other nodes send data to node 0
      else
         call MPI_SEND(partt,nntp,mreal,0,ltag,lgrp,ierr)
         call PPBICAST(iwork,1)
         irc = iwork(1)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPADJFVS2(fvs,gvs,noff,nyp,ndim,nmv,mvy,nxb,nyb,nmvf)
! for 2-1/2d code, this subroutine adjusts 3d velocity distribution, in
! different regions of space, so that partial regions have equal amounts
! of spatial grid points
! input: all except gvs, output: fvs, gvs
! fvs = spatially resolved distribution function
! gvs = scratch array
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! ndim = number of velocity dimensions = 2 or 3
! mvy = number of grids in y for phase space aggregation
! nxb/nyb = number of segments in x/y for velocity distribution
! nmvf = first dimension of fvs
      implicit none
      integer, intent(in) :: noff, nyp
      integer, intent(in) :: ndim, nmv, mvy, nxb, nyb, nmvf
      real, dimension(nmvf,ndim,nxb,nyb+1), intent(inout) :: fvs
      real, dimension(nmvf,ndim,nxb), intent(inout) :: gvs
! local data
      integer :: i, j, k, l, ns, ne, nmv21, ndata, idproc, nvp, id
      integer :: msid, ltag, ierr
      integer, dimension(lstat) :: istatus
      real, dimension(ndim) :: scale
      nmv21 = 2*nmv + 1
      ns = noff - mvy*(noff/mvy)
      ne = noff + nyp
      ne = ne - mvy*(ne/mvy)
      ndata = nmvf*ndim*nxb
      ltag = nmvf*ndim
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! receive data from right if last region on this node has partial grid
      if (ne > 0) then
         id = idproc + 1
         if (id < nvp) then
            call MPI_IRECV(gvs,ndata,mreal,id,ltag,lgrp,msid,ierr)
         endif
      endif
! send data to right if last region on previous node has partial grid
      if (ns > 0) then
         id = idproc - 1
         if (id >= 0) then
            call MPI_SEND(fvs,ndata,mreal,id,ltag,lgrp,ierr)
         endif
! save scale
         do j = 1, ndim
         scale(j) = fvs(nmv21+1,j,1,1)
         enddo
! left shift data in y
         do l = 1, nyb
         do k = 1, nxb
         do j = 1, ndim
         do i = 1, nmvf
         fvs(i,j,k,l) = fvs(i,j,k,l+1)
         enddo
         enddo
         enddo
         enddo
! restore scale
         do j = 1, ndim
         fvs(nmv21+1,j,1,1) = scale(j)
         enddo
! zero out extra element
         do k = 1, nxb
         do j = 1, ndim
         do i = 1, nmvf
         fvs(i,j,k,nyb+1) = 0.0
         enddo
         enddo
         enddo
      endif
! receive data from right if last region on this node has partial grid
      if (ne > 0) then
         id = idproc + 1
         if (id < nvp) then
            call MPI_WAIT(msid,istatus,ierr)
            do k = 1, nxb
            do j = 1, ndim
            do i = 1, nmvf
            fvs(i,j,k,nyb) = fvs(i,j,k,nyb) + gvs(i,j,k)
            enddo
            enddo
            enddo
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVNDATA2(f,g,ndim,nxv,nyp,nypmx,iunit)
! this subroutine collects distributed periodic real 2d vector
! non-uniform data f and writes to a fortran unformatted file
! f = input data to be written
! g = scratch data
! ndim = first dimension of data array f
! nyp = actual data written for third dimension
! nxv = second dimension of data array f
! nypmx = third dimension of data array f
! iunit = fortran unit number
! input: all
      implicit none
      integer, intent(in) :: ndim, nxv, nyp, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nyp) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nps, nnxv, nnyp, nnypx, id, i, j, k, n
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nnxv = ndim*nxv
      nnyp = nnxv*nyp
      nnypx = nnxv*nypmx
! this segment is used for shared memory computers
!     write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nyp)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nyp)
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
         call MPI_IRECV(g,nnypx,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         nps = nps/nnxv
         write (unit=iunit) (((g(n,j,k),n=1,ndim),j=1,nxv),k=1,nps)
         enddo
! other nodes send data to node 0
      else
         call MPI_SEND(f,nnyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
      use mpplib2, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use mpplib2, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use mpplib2, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use mpplib2, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use mpplib2, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use mpplib2, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMAX(f,g,nxp)
      use mpplib2, only: SUB => PPMAX
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use mpplib2, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use mpplib2, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSCAN(f,g,nxp)
      use mpplib2, only: SUB => PPDSCAN
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBICAST(if,nxp)
      use mpplib2, only: SUB => PPBICAST
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if
      call SUB(if,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBDCAST(f,nxp)
      use mpplib2, only: SUB => PPBDCAST
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      call SUB(f,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPISHFTR(if,ig,nxp)
      use mpplib2, only: SUB => PPISHFTR
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNCGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      call SUB(f,nyp,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNAGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNACGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPFMOVE2(f,g,noff,nyp,noffs,nyps,noffd,nypd,isign,kyp, &
     &ny,kstrt,nvp,nxv,nypmx,mter,ierr)
      use mpplib2, only: SUB => PPFMOVE2
      implicit none
      integer, intent(in) :: noff, nyp, isign, kyp, ny, kstrt, nvp
      integer, intent(in) :: nxv, nypmx
      integer, intent(inout) :: noffs, nyps, noffd, nypd
      integer, intent(inout) :: mter, ierr
      real, dimension(nxv,nypmx), intent(inout) :: f, g
      call SUB(f,g,noff,nyp,noffs,nyps,noffd,nypd,isign,kyp,ny,kstrt,nvp&
     &,nxv,nypmx,mter,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,  &
     &kypd)
      use mpplib2, only: SUB => PPTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      complex, dimension(nxv,kypd), intent(in) :: f
      complex, dimension(nyv,kxpd), intent(inout) :: g
      complex, dimension(kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv, &
     &kxpd,kypd)
      use mpplib2, only: SUB => PPNTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      complex, dimension(ndim,nxv,kypd), intent(in) :: f
      complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
      complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny&
     &,kstrt,nvp,ndim,nc,idimp,npmax,idps,nbmax,ntmax,info)
      use mpplib2, only: SUB => PPMOVE2
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, ndim, nc, idimp, npmax
      integer, intent(in) :: idps, nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(ntmax+1), intent(in) :: ihole
      integer, dimension(5), intent(inout) :: info
      call SUB(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp&
     &,ndim,nc,idimp,npmax,idps,nbmax,ntmax,info)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
      use mpplib2, only: SUB => PPPMOVE2
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax), intent(in) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(3,mx1), intent(in) :: ncll, nclr
      integer, dimension(3,mx1), intent(inout) :: mcll, mclr
      call SUB(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,   &
     &idimp,nbmax,mx1)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRITE2(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec)
      use mpplib2, only: SUB => PPWRITE2
      implicit none
      integer, intent(in) :: nx, ny, kyp, nxv, nypmx, iunit
      integer, intent(inout) :: nrec
      real, intent(in), dimension(nxv,nypmx) :: f
      real, intent(inout), dimension(nxv,nypmx) :: g
      call SUB(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPREAD2(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec,irc)
      use mpplib2, only: SUB => PPREAD2
      implicit none
      integer, intent(in) :: nx, ny, kyp, nxv, nypmx, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout), dimension(nxv,nypmx) :: f, g
      call SUB(f,g,nx,ny,kyp,nxv,nypmx,iunit,nrec,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVWRITE2(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec)
      use mpplib2, only: SUB => PPVWRITE2
      implicit none
      integer, intent(in) :: nx, ny, kyp, ndim, nxv, nypmx, iunit
      integer, intent(inout) :: nrec
      real, intent(in), dimension(ndim,nxv,nypmx) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
      call SUB(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVREAD2(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec,irc)
      use mpplib2, only: SUB => PPVREAD2
      implicit none
      integer, intent(in) :: nx, ny, kyp, ndim, nxv, nypmx, iunit
      integer, intent(inout) :: nrec, irc
      real, intent(inout), dimension(ndim,nxv,nypmx) :: f, g
      call SUB(f,g,nx,ny,kyp,ndim,nxv,nypmx,iunit,nrec,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPCWRITE2(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec)
      use mpplib2, only: SUB => PPCWRITE2
      implicit none
      integer, intent(in) :: nx, ny, kxp, nyv, kxpd, iunit
      integer, intent(inout) :: nrec
      complex, intent(in), dimension(nyv,kxpd) :: f
      complex, intent(inout), dimension(nyv,kxpd) :: g
      call SUB(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPCREAD2(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec,irc)
      use mpplib2, only: SUB => PPCREAD2
      implicit none
      integer, intent(in) :: nx, ny, kxp, nyv, kxpd, iunit
      integer, intent(inout) :: nrec, irc
      complex, intent(inout), dimension(nyv,kxpd) :: f, g
      call SUB(f,g,nx,ny,kxp,nyv,kxpd,iunit,nrec,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVCWRITE2(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec)
      use mpplib2, only: SUB => PPVCWRITE2
      implicit none
      integer, intent(in) :: nx, ny, kxp, ndim, nyv, kxpd, iunit
      integer, intent(inout) :: nrec
      complex, intent(in), dimension(ndim,nyv,kxpd) :: f
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: g
      call SUB(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPVCREAD2(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec,irc)
      use mpplib2, only: SUB => PPVCREAD2
      implicit none
      integer, intent(in) :: nx, ny, kxp, ndim, nyv, kxpd, iunit
      integer, intent(inout) :: nrec, irc
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: f, g
      call SUB(f,g,nx,ny,kxp,ndim,nyv,kxpd,iunit,nrec,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRPART2(part,npp,idimp,npmax,iunit,iscr)
      use mpplib2, only: SUB => PPWRPART2
      implicit none
      integer, intent(in) :: npp, idimp, npmax, iunit, iscr
      real, intent(inout), dimension(idimp,npmax) :: part
      call SUB(part,npp,idimp,npmax,iunit,iscr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDPART2(part,npp,idimp,npmax,iunit,iscr,irc)
      use mpplib2, only: SUB => PPRDPART2
      implicit none
      integer, intent(in) :: idimp, npmax, iunit, iscr
      integer, intent(inout) :: npp, irc
      real, intent(inout), dimension(idimp,npmax) :: part
      call SUB(part,npp,idimp,npmax,iunit,iscr,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRDATA2(f,g,nxv,nypmx,iunit)
      use mpplib2, only: SUB => PPWRDATA2
      implicit none
      integer, intent(in) :: nxv, nypmx, iunit
      real, intent(in), dimension(nxv,nypmx) :: f
      real, intent(inout), dimension(nxv,nypmx) :: g
      call SUB(f,g,nxv,nypmx,iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDDATA2(f,g,nxv,nypmx,iunit,irc)
      use mpplib2, only: SUB => PPRDDATA2
      implicit none
      integer, intent(in) :: nxv, nypmx, iunit
      integer, intent(inout) :: irc
      real, intent(inout), dimension(nxv,nypmx) :: f, g
      call SUB(f,g,nxv,nypmx,iunit,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVDATA2(f,g,ndim,nxv,nypmx,iunit)
      use mpplib2, only: SUB => PPWRVDATA2
      implicit none
      integer, intent(in) :: ndim, nxv, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nypmx) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
      call SUB(f,g,ndim,nxv,nypmx,iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDVDATA2(f,g,ndim,nxv,nypmx,iunit,irc)
      use mpplib2, only: SUB => PPRDVDATA2
      implicit none
      integer, intent(in) :: ndim, nxv, nypmx, iunit
      integer, intent(inout) :: irc
      real, intent(inout), dimension(ndim,nxv,nypmx) :: f, g
      call SUB(f,g,ndim,nxv,nypmx,iunit,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVCDATA2(f,g,ndim,nyv,kxpd,iunit)
      use mpplib2, only: SUB => PPWRVCDATA2
      implicit none
      integer, intent(in) :: ndim, nyv, kxpd, iunit
      complex, intent(in), dimension(ndim,nyv,kxpd) :: f
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: g
      call SUB(f,g,ndim,nyv,kxpd,iunit)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRDVCDATA2(f,g,ndim,nyv,kxpd,iunit,irc)
      use mpplib2, only: SUB => PPRDVCDATA2
      implicit none
      integer, intent(in) :: ndim, nyv, kxpd, iunit
      integer, intent(inout) :: irc
      complex, intent(inout), dimension(ndim,nyv,kxpd) :: f, g
      call SUB(f,g,ndim,nyv,kxpd,iunit,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPARTT2(partt,numtp,idimp,nprobt,irc)
      use mpplib2, only: SUB => PPARTT2
      implicit none
      integer, intent(in) :: numtp, idimp, nprobt
      integer, intent(inout) :: irc
      real, intent(inout), dimension(idimp,nprobt) :: partt
      call SUB(partt,numtp,idimp,nprobt,irc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPADJFVS2(fvs,gvs,noff,nyp,ndim,nmv,mvy,nxb,nyb,nmvf)
      use mpplib2, only: SUB => PPADJFVS2
      implicit none
      integer, intent(in) :: noff, nyp
      integer, intent(in) :: ndim, nmv, mvy, nxb, nyb, nmvf
      real, dimension(nmvf,ndim,nxb,nyb+1), intent(inout) :: fvs
      real, dimension(nmvf,ndim,nxb), intent(inout) :: gvs
      call SUB(fvs,gvs,noff,nyp,ndim,nmv,mvy,nxb,nyb,nmvf)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVNDATA2(f,g,ndim,nxv,nyp,nypmx,iunit)
      use mpplib2, only: SUB => PPWRVNDATA2
      implicit none
      integer, intent(in) :: ndim, nxv, nyp, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nyp) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
      call SUB(f,g,ndim,nxv,nyp,nypmx,iunit)
      end subroutine
