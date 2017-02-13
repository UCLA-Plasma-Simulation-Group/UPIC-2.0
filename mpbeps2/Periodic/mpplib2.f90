!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications with OpenMP
! pplib2.f90 contains basic communications procedures for 1d partitions:
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
! PPPMOVE2 moves particles into appropriate spatial regions for tiled
!          distributed data.
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: february 4, 2017
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
      public :: PPSUM, PPDSUM, PPMAX, PPIMAX, PPDMAX, PPBICAST, PPBDCAST
      public :: PPISHFTR, PPNCGUARD2L, PPNAGUARD2L, PPNACGUARD2L
      public :: PPFMOVE2, PPTPOSE, PPNTPOSE, PPPMOVE2
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
      if (nproc.eq.1) return
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
      if (nproc.eq.1) return
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
      if (nproc.eq.1) return
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
      if (nproc.eq.1) return
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
      if (nproc.eq.1) return
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,mmax,lgrp,ierr)
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
      if (nproc.eq.1) return
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
      if (nproc.eq.1) return
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
      if (nvp.eq.1) return
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


