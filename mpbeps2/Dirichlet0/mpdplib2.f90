!-----------------------------------------------------------------------
! Parallel PIC library for MPI communications with OpenMP
! for non-periodic field boundary condtions
! PPNLCGUARD2L copies data to interior guard cells in y for scalar data,
!              linear interpolation, and distributed data with
!              non-uniform partition.
! PPNLAGUARD2L adds interior guard cells in y for scalar array, linear
!              interpolation, and distributed data with non-uniform
!              partition.
! PPNLACGUARD2L adds interior guard cells in y for vector array, linear
!               interpolation, and distributed data with non-uniform
!               partition.
! PPDBLSIN2C creates a doubled 2 component vector array to perform 2d
!            sine/cosine transforms with real to complex ffts, for
!            distributed data.
! PPDBLSIN2D creates a doubled scalar array to perform 2d sine transform
!            with real to complex fft, for distributed data.
! PPDBLSIN2B creates a doubled 3 component vector array to perform 2d
!            sine/cosine transforms with real to complex ffts, for
!            distributed data.
! PPDBLSIN2M creates a doubled 4 component tensor array to perform 2d
!            sine/cosine transforms with real to complex ffts, for
!            distributed data.
! PPDBLSIN22M creates a doubled 2 component tensor array to perform 2d
!             sine/cosine transforms with real to complex ffts, for
!             distributed data.
! PPHAFDBL2D copies data from a doubled scalar array q2 to a regular
!            scalar array q, for distributed data.
! PPHAFDBL2C copies data from a double vector array fxy2 to regular
!            vector array fxy, for distributed data.
! PPRTPOSE performs a transpose of a real scalar array, distributed
!          in y, to a real scalar array, distributed in x.
! PPRNTPOSE performs a transpose of a real vector array, distributed
!          in y, to a real vector array, distributed in x.
! mplcguard2 copies scalar data to interior guard cells
!            calls PPNLCGUARD2L
! mpnlcguard2 copies vector data to interior guard cells
!             calls PPNLCGUARD2L
! mpnlaguard2 adds scalar data from interior guard cells
!             calls PPNLAGUARD2L
! mpnlacguard2 adds vector data from interior guard cells
!              calls PPNLACGUARD2L
! mpdblsin2d creates an odd array q2 from an array q
!            calls PPDBLSIN2D
! mpdblsin2c creates a mixed even/odd array cu2 from an array cu
!            calls PPDBLSIN2C or PPDBLSIN2B
! mpdblsin2m creates mixed even/odd tensor array amu2 from an array amu
!            calls PPDBLSIN2M
! mphafdbl2d copies data from a double scaler array q2 to regular scalar
!            array q
!            calls PPHAFDBL2D
! mphafdbl2c copies data from a double vector array fxy2 to regular
!            vector array fxy
!            calls PPHAFDBL2C
! written by viktor k. decyk, ucla
! copyright 1999-2017, regents of the university of california
! update: march 29, 2017
      module mpdplib2
      use mpplib2
      implicit none
!
! scr = guard cell buffer received from nearby processors
      real, dimension(:), allocatable  :: scr
      integer :: szscr = -1
! qs = scratch arrays for q shift
      real, dimension(:), allocatable  :: qs
      integer :: szqs = -1
! q2s = scratch arrays for q2s shift
      real, dimension(:,:), allocatable  :: q2s
      integer :: szq2s = -1
! cus = scratch arrays for cu shift
      real, dimension(:,:), allocatable  :: cus
      integer :: szcus = -1
! cu2s = scratch arrays for cu2s shift
      real, dimension(:,:,:), allocatable  :: cu2s
      integer :: szcu2s = -1
! amus = scratch arrays for amu shift
      real, dimension(:,:), allocatable  :: amus
      integer :: szamus = -1
! amu2s = scratch arrays for amu2s shift
      real, dimension(:,:,:), allocatable  :: amu2s
      integer :: szamu2s = -1
      save
!
      private :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! this subroutine copies data to guard cells in non-uniform partitions
! guard cell on last processor is presumed already set.
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
      integer :: ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! special case for one processor
      if (nvp==1) return
      ks = kstrt - 1
      moff = nypmx*nvp + 2
! copy to guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kr < nvp) then
         call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      endif
      if (kl >= 0) then
         call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      endif
      if (kr < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNLAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! for scalar data.  no copying is done at the boundary edges.
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
      if (nvp==1) return
      ks = kstrt - 1
      moff = nypmx*nvp + 1
! add guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kl >= 0) then
         call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      endif
      if (kr < nvp) then
         call MPI_SEND(f(1,ks),nxv,mreal,kr,moff,lgrp,ierr)
         do j = 1, nx1
            f(j,ks) = 0.0
         enddo
      endif
      if (kl >= 0) then
         call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
         do j = 1, nx1
            f(j,1) = f(j,1) + scr(j)
         enddo
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNLACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! for vector data.  no copying is done at the boundary edges.
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
      if (nvp==1) return
      ks = kstrt - 1
      moff = nypmx*nvp + 1
      nnxv = ndim*nxv
! add guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kl >= 0) then
         call MPI_IRECV(scr,nnxv,mreal,kl,moff,lgrp,msid,ierr)
      endif
      if (kr < nvp) then
         call MPI_SEND(f(1,1,ks),nnxv,mreal,kr,moff,lgrp,ierr)
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,ks) = 0.0
            enddo
         enddo
      endif
      if (kl >= 0) then
         call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,1) = f(n,j,1) + scr(n,j)
            enddo
         enddo
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2C(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd&
     &,kyp2d,kshd)
! this subroutine creates a doubled vector array cu2 from a vector array
! cu, so that various 2d sine/cosine transforms can be performed with a
! 2d real to complex fft.  the x component is an odd function in y,
! and y component is an odd function in x.
! used for field solvers with dirichlet boundary conditions using ffts
! for distributed data
! input data cu must have a uniform partition
! cu array is temporarily modified and restored
! cus, cu2s = scratch arrays for shifts
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = second dimension of input array cu, must be >= nx+1
! kyp = number of data values per block in y
! kypd = third dimension of input array cu, must be >= kyp+1
! kyp2d = third dimension of output array cu2, must be >= 2*kyp
! kshd = third dimension scratch cu2s, kshd >= ksh, described below
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(2,nxv,kypd), intent(inout) :: cu
      real, dimension(2,2*nxv,kyp2d), intent(inout) :: cu2
      real, dimension(2,nxv), intent(inout) :: cus
      real, dimension(2,2*nxv,kshd), intent(inout) :: cu2s
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: i, j, k, nxs, nys, kyb, ny2, kyp2, kyb2, ks, moff, kyp1
      integer :: kyps, kyp2s, ksh, koff, l1, l2, m1, k1off, k2off, kypl
      integer :: kypr, k1, k2, joff, kypn, kypm, kk, nym, kl, kr, nnxv
      integer, dimension(lstat) :: istatus
      nxs = nx - 1
      nys = ny - 1
      nnxv = 2*nxv
! kyb = minimum number of processors used by input array
      kyb = (ny - 1)/kyp + 1
      ny2 = ny + ny
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny2 - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for received messages
      kyp1 = kyp + 1
! kyps = actual size used in input array on each node
      kyps = min(kyp,max(0,ny-kyp*ks))
! kyp2s = minimum positive actual size used in doubled array
      kyp2s = min(kyp2,max(0,ny2-kyp2*(kyb2-1)))
! ksh = final shift of received data
      ksh = kyp2 - kyp2s
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nys
            do j = 1, nxs
               cu2(1,j+1,k+1) = cu(1,j+1,k+1)
               cu2(2,j+1,k+1) = cu(2,j+1,k+1)
               cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
               cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
               cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
               cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
               cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
               cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
            enddo
            cu2(1,1,k+1) = cu(1,1,k+1)
            cu2(2,1,k+1) = 0.0
            cu2(1,nx+1,k+1) = cu(1,nx+1,k+1)
            cu2(2,nx+1,k+1) = 0.0
            cu2(1,1,k+ny+1) = -cu(1,1,ny-k+1)
            cu2(2,1,k+ny+1) = 0.0
            cu2(1,nx+1,k+ny+1) = -cu(1,nx+1,ny-k+1)
            cu2(2,nx+1,k+ny+1) = 0.0
         enddo
!$OMP END PARALLEL DO
         do j = 1, nxs
            cu2(1,j+1,1) = 0.0
            cu2(2,j+1,1) = cu(2,j+1,1)
            cu2(1,j+nx+1,1) = 0.0
            cu2(2,j+nx+1,1) = -cu(2,nx-j+1,1)
            cu2(1,j+1,ny+1) = 0.0
            cu2(2,j+1,ny+1) = cu(2,j+1,ny+1)
            cu2(1,j+nx+1,ny+1) = 0.0
            cu2(2,j+nx+1,ny+1) = -cu(2,nx-j+1,ny+1)
         enddo
         cu2(1,1,1) = 0.0
         cu2(2,1,1) = 0.0
         cu2(1,nx+1,1) = 0.0
         cu2(2,nx+1,1) = 0.0
         cu2(1,1,ny+1) = 0.0
         cu2(2,1,ny+1) = 0.0
         cu2(1,nx+1,ny+1) = 0.0
         cu2(2,nx+1,ny+1) = 0.0
         return
      endif
!
! copy to double array in x direction
!     cu2(1,j+1,k+1) = cu(1,j+1,k+1)
!     cu2(2,j+1,k+1) = cu(2,j+1,k+1)
      koff = kyp2*ks
! l1 = source location for data for first half
      l1 = koff/kyp
! l2 = source location for data for second half
      l2 = l1 + 1
      k1off = kyp*l1 - koff + 1
      k2off = kyp*l2 - koff + 1
      koff = kyp*ks
! m1 = destination location of original data
      m1 = koff/kyp2
      koff = koff - kyp2*m1
! kypl = room available at destination
      kypl = min(kyps,kyp2-koff)
! this segment is used for mpi computers
!
! receive first half
      if (l1 < nvp) then
         call MPI_IRECV(cu2,kyp*nnxv,mreal,l1,moff+1,lgrp,msid,ierr)
! receive second half
         if (l2 < nvp) then
            call MPI_IRECV(cu2(1,1,kyp+1),kyp*nnxv,mreal,l2,moff+1,lgrp,&
     &nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(cu,kypl*nnxv,mreal,m1,moff+1,lgrp,ierr)
      endif
! wait for data and unpack it
      if (l1 < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
         kypr = kypr/nnxv
         do k = 2, kypr
            k1 = kypr - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do j = 1, nxv
               do i = 1, 2
                  cu2(i,j,k1) = cu2(i,j+joff,k2)
               enddo
            enddo
         enddo
         if (l2 < nvp) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     cu2(i,j,k1+kyp) = cu2(i,j+joff,k2+kyp)
                  enddo
               enddo
            enddo
         endif
      endif
!
! copy to double array in y direction
!     cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
!     cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
      koff = kyp2*ks
! l1 = source location for data for second half
      l1 = (ny2 + ksh - koff - kyp2)/kyp
! l2 = source location for data for first half
      l2 = l1 + 1
! kypn = check if address for l1 receive is inbounds
      kypn = min(kyp,max(0,koff+ksh-ny+kyp1+kyp-1))
! kypm = check if address for l2 receive is inbounds
      kypm = min(kyp,max(0,koff+ksh-ny+kyp-1))
      koff = kyp*ks
! m1 = destination location of original data
      m1 = (ny2 + ksh - koff - 1)/kyp2
! amount of data available to send
      kypl = min(kyps,max(0,ny-koff-1))
! this segment is used for mpi computers
!
! first align data in input cu by left shifting
      do j = 1, nxv
         do i = 1, 2
            cus(i,j) = cu(i,j,1)
         enddo
      enddo
      do k = 2, kyp
         do j = 1, nxv
            do i = 1, 2
               cu(i,j,k-1) = cu(i,j,k)
            enddo
         enddo
      enddo
      if (ks < (nvp-1)) then
         call MPI_IRECV(cu(1,1,kyp),nnxv,mreal,ks+1,moff+2,lgrp,msid,   &
     &ierr)
      endif
      if (ks > 0) then
         call MPI_SEND(cus,nnxv,mreal,ks-1,moff+2,lgrp,ierr)
      endif
      if (ks < (nvp-1)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then continue to copy to double array in y direction
! post receive second half
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_IRECV(cu2(1,1,kyp1),kyp*nnxv,mreal,l1,moff+3,lgrp, &
     &msid,ierr)
         endif
      endif
! post receive first half
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_IRECV(cu2,kyp*nnxv,mreal,l2,moff+3,lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(cu,kypl*nnxv,mreal,m1,moff+3,lgrp,ierr)
      endif
! wait for data and unpack it
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     cu2(i,j,k1+kyp1-1) = cu2(i,j+joff,k2+kyp1-1)
                  enddo
               enddo
            enddo
         endif
      endif
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     cu2(i,j,k1) = cu2(i,j+joff,k2)
                  enddo
               enddo
            enddo
         endif
      endif
!
! restore original data in input cu
      do k = 2, kyp
         kk = kyp - k + 2
         do j = 1, nxv
            do i = 1, 2
               cu(i,j,kk) = cu(i,j,kk-1)
            enddo
         enddo
      enddo
      do j = 1, nxv
         do i = 1, 2
            cu(i,j,1) = cus(i,j)
         enddo
      enddo
!
! create mixed even/odd array
! first reverse local y index
      koff = kyp2*ks
      nym = ny2 - kyb*kyp
      do k = 1, kyp2
         kk = k + koff - ksh
         if (kk > nym) then
            if (k <= kyp) then
               k1 = kyp + 1
               do j = 1, nxs
                  cu2(1,nx+j+1,k1-k) = -cu2(1,nx-j+1,k)
                  cu2(2,nx+j+1,k1-k) = -cu2(2,nx-j+1,k)
               enddo
            else
               k1 = kyp2 + kyp + 1
               do j = 1, nxs
                  cu2(1,nx+j+1,k1-k) = -cu2(1,nx-j+1,k)
                  cu2(2,nx+j+1,k1-k) = -cu2(2,nx-j+1,k)
               enddo
            endif
            cu2(1,2*nx+1,k1-k) = -cu2(1,1,k)
            cu2(1,2*nx+2,k1-k) = -cu2(1,nx+1,k)
         endif
         cu2(1,1,k) = -cu2(1,1,k)
         cu2(2,1,k) = 0.0
         cu2(1,nx+1,k) = -cu2(1,nx+1,k)
         cu2(2,nx+1,k) = 0.0
      enddo
!
! ny+1 point is special
! l1 = source processor id for ny+1 point
      l1 = (ny - 1)/kyp
! m1 = destination processor id of ny+1 point
      m1 = ny/kyp2
! k2 = destination location for ny+1 point
      k2 = ny - kyp2*m1 + 1
      if (ks==m1) then
         call MPI_IRECV(cu2(1,1,k2),nnxv,mreal,l1,moff+6,lgrp,msid,ierr)
      endif
      if (ks==l1) then
         call MPI_SEND(cu(1,1,kyps+1),nnxv,mreal,m1,moff+6,lgrp,ierr)
      endif
      if (ks==m1) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then left shift rhs of output cu2 by ksh y indices
      if (ksh > 0) then
         nym = ny + 1
! kr = number of data points which arriving from right processor
         kr = min(ksh,max(0,koff+kyp2-nym))
! kl = number of data points which need to be moved to left processor
         kl = min(ksh,max(0,koff-nym))
! kk = lower bound for left shift
         kk = max(0,nym-koff)
         if (kl > 0) then
! copy data going to left processor
            do k = 1, kl
               do j = 1, 2*nxv
                  do i = 1, 2
                     cu2s(i,j,k) = cu2(i,j,k+ksh-kl)
                  enddo
               enddo
            enddo
         endif
! shift remaining data
         do k = ksh + kk + 1, kyp2
            do j = 1, 2*nxv
               do i = 1, 2
                  cu2(i,j,k-ksh) = cu2(i,j,k)
               enddo
            enddo
         enddo
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_IRECV(cu2(1,1,kyp2-kr+1),2*nnxv*ksh,mreal,ks+1, &
     &moff+4,lgrp,msid,ierr)
            endif
         endif
         if (kl > 0) then
            if (ks > 0) then
               call MPI_SEND(cu2s,2*nnxv*ksh,mreal,ks-1,moff+4,lgrp,ierr&
     &)
            endif
         endif
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_WAIT(msid,istatus,ierr)
            endif
         endif
      endif
!
! then create second half of mixed even/odd array
!     cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
!     cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if ((kk==1).or.(kk==(ny+1))) then
            do j = 1, nxs
               cu2(1,j+1,k) = 0.0
               cu2(2,j+1,k) = cu2(2,j+1,k)
               cu2(1,j+nx+1,k) = 0.0
               cu2(2,j+nx+1,k) = -cu2(2,nx-j+1,k)
            enddo
            cu2(1,1,k) = 0.0
            cu2(2,1,k) = 0.0
            cu2(1,nx+1,k) = 0.0
            cu2(2,nx+1,k) = 0.0
         else if (kk <= ny) then
            do j = 1, nxs
               cu2(1,nx+j+1,k) = cu2(1,nx-j+1,k)
               cu2(2,nx+j+1,k) = -cu2(2,nx-j+1,k)
            enddo
            cu2(1,1,k) = -cu2(1,1,k)
            cu2(2,1,k) = 0.0
            cu2(1,nx+1,k) = -cu2(1,nx+1,k)
            cu2(2,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
!
! finish mixed even/odd array
!     cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
!     cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if (kk > (ny+1)) then
            do j = 1, nxs
               cu2(1,nx-j+1,k) = cu2(1,nx+j+1,k)
               cu2(2,nx-j+1,k) = -cu2(2,nx+j+1,k)
            enddo
            cu2(1,1,k) = cu2(1,2*nx+1,k)
            cu2(1,nx+1,k) = cu2(1,2*nx+2,k)
            cu2(2,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2D(q,q2,qs,q2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,   &
     &kyp2d,kshd)
! this subroutine creates an odd array q2 from an array q, so that
! a 2d sine transform can be performed with a 2d real to complex fft.
! used for field solvers with dirichlet boundary conditions using ffts
! for distributed data
! input data q must have a uniform partition
! q array is temporarily modified and restored
! qs, q2s = scratch arrays for shifts
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = first dimension of input array q, must be >= nx+1
! kyp = number of data values per block in y
! kypd = second dimension of input array q, must be >= kyp+1
! kyp2d = second dimension of output array q2, must be >= 2*kyp
! kshd = second dimension scratch q2s, kshd >= ksh, described below
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(nxv,kypd), intent(inout) :: q
      real, dimension(2*nxv,kyp2d), intent(inout) :: q2
      real, dimension(nxv), intent(inout) :: qs
      real, dimension(2*nxv,kshd), intent(inout) :: q2s
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: j, k, nxs, nys, kyb, ny2, kyp2, kyb2, ks, moff, kyp1
      integer :: kyps, kyp2s, ksh, koff, l1, l2, m1, k1off, k2off, kypl
      integer :: kypr, k1, k2, joff, kypn, kypm, kk, nym, kl, kr
      integer, dimension(lstat) :: istatus
      nxs = nx - 1
      nys = ny - 1
! kyb = minimum number of processors used by input array
      kyb = (ny - 1)/kyp + 1
      ny2 = ny + ny
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny2 - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for received messages
      kyp1 = kyp + 1
! kyps = actual size used in input array on each node
      kyps = min(kyp,max(0,ny-kyp*ks))
! kyp2s = minimum positive actual size used in doubled array
      kyp2s = min(kyp2,max(0,ny2-kyp2*(kyb2-1)))
! ksh = final shift of received data
      ksh = kyp2 - kyp2s
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nys
            do j = 1, nxs
               q2(j+1,k+1) = q(j+1,k+1)
               q2(nx+j+1,k+1) = -q(nx-j+1,k+1)
               q2(j+1,ny+k+1) = -q(j+1,ny-k+1)
               q2(nx+j+1,ny+k+1) = q(nx-j+1,ny-k+1)
            enddo
            q2(1,k+1) = 0.0
            q2(nx+1,k+1) = 0.0
            q2(1,k+ny+1) = 0.0
            q2(nx+1,k+ny+1) = 0.0
         enddo
!$OMP END PARALLEL DO
         do j = 1, nx
            q2(j,1) = 0.0
            q2(j+nx,1) = 0.0
            q2(j,ny+1) = 0.0
            q2(j+nx,ny+1) = 0.0
         enddo
         return
      endif
!
! copy to double array in x direction
!     q2(j+1,k+1) = q(j+1,k+1)
      koff = kyp2*ks
! l1 = source location for data for first half
      l1 = koff/kyp
! l2 = source location for data for second half
      l2 = l1 + 1
      k1off = kyp*l1 - koff + 1
      k2off = kyp*l2 - koff + 1
      koff = kyp*ks
! m1 = destination location of original data
      m1 = koff/kyp2
      koff = koff - kyp2*m1
! kypl = room available at destination
      kypl = min(kyps,kyp2-koff)
! this segment is used for mpi computers
!
! receive first half
      if (l1 < nvp) then
         call MPI_IRECV(q2,kyp*nxv,mreal,l1,moff+1,lgrp,msid,ierr)
! receive second half
         if (l2 < nvp) then
            call MPI_IRECV(q2(1,kyp+1),kyp*nxv,mreal,l2,moff+1,lgrp,nsid&
     &,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(q,kypl*nxv,mreal,m1,moff+1,lgrp,ierr)
      endif
! wait for data and unpack it
      if (l1 < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
         kypr = kypr/nxv
         do k = 2, kypr
            k1 = kypr - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do j = 1, nxv
               q2(j,k1) = q2(j+joff,k2)
            enddo
         enddo
         if (l2 < nvp) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  q2(j,k1+kyp) = q2(j+joff,k2+kyp)
               enddo
            enddo
         endif
      endif
!
! copy to double array in y direction
!     q2(nx+j+1,ny+k+1) = q(nx-j+1,ny-k+1)
      koff = kyp2*ks
! l1 = source location for data for second half
      l1 = (ny2 + ksh - koff - kyp2)/kyp
! l2 = source location for data for first half
      l2 = l1 + 1
! kypn = check if address for l1 receive is inbounds
      kypn = min(kyp,max(0,koff+ksh-ny+kyp1+kyp-1))
! kypm = check if address for l2 receive is inbounds
      kypm = min(kyp,max(0,koff+ksh-ny+kyp-1))
      koff = kyp*ks
! m1 = destination location of original data
      m1 = (ny2 + ksh - koff - 1)/kyp2
! amount of data available to send
      kypl = min(kyps,max(0,ny-koff-1))
! this segment is used for mpi computers
!
! first align data in input q by left shifting
      do j = 1, nxv
         qs(j) = q(j,1)
      enddo
      do k = 2, kyp
         do j = 1, nxv
            q(j,k-1) = q(j,k)
         enddo
      enddo
      if (ks < (nvp-1)) then
         call MPI_IRECV(q(1,kyp),nxv,mreal,ks+1,moff+2,lgrp,msid,ierr)
      endif
      if (ks > 0) then
         call MPI_SEND(qs,nxv,mreal,ks-1,moff+2,lgrp,ierr)
      endif
      if (ks < (nvp-1)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then continue to copy to double array in y direction
! post receive second half
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_IRECV(q2(1,kyp1),kyp*nxv,mreal,l1,moff+3,lgrp,msid,&
     &ierr)
         endif
      endif
! post receive first half
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_IRECV(q2,kyp*nxv,mreal,l2,moff+3,lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(q,kypl*nxv,mreal,m1,moff+3,lgrp,ierr)
      endif
! wait for data and unpack it
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  q2(j,k1+kyp1-1) = q2(j+joff,k2+kyp1-1)
               enddo
            enddo
         endif
      endif
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  q2(j,k1) = q2(j+joff,k2)
               enddo
            enddo
         endif
      endif
!
! restore original data in input q
      do k = 2, kyp
         kk = kyp - k + 2
         do j = 1, nxv
            q(j,kk) = q(j,kk-1)
         enddo
      enddo
      do j = 1, nxv
         q(j,1) = qs(j)
      enddo
!
! create odd array
! first reverse local y index
      koff = kyp2*ks
      nym = ny2 - kyb*kyp
      do k = 1, kyp2
         kk = k + koff - ksh
         if (kk > nym) then
            if (k <= kyp) then
               k1 = kyp + 1
               do j = 1, nxs
                  q2(nx+j+1,k1-k) = q2(nx-j+1,k)
               enddo
            else
               k1 = kyp2 + kyp + 1
               do j = 1, nxs
                  q2(nx+j+1,k1-k) = q2(nx-j+1,k)
               enddo
            endif
         endif
         q2(1,k) = 0.0
         q2(nx+1,k) = 0.0
      enddo
!
! then left shift rhs of output q2 by ksh y indices
      if (ksh > 0) then
         nym = ny + 1
! kr = number of data points which arriving from right processor
         kr = min(ksh,max(0,koff+kyp2-nym))
! kl = number of data points which need to be moved to left processor
         kl = min(ksh,max(0,koff-nym))
! kk = lower bound for left shift
         kk = max(0,nym-koff)
         if (kl > 0) then
! copy data going to left processor
            do k = 1, kl
               do j = 1, 2*nxv
                  q2s(j,k) = q2(j,k+ksh-kl)
               enddo
            enddo
         endif
! shift remaining data
         do k = ksh + kk + 1, kyp2
            do j = 1, 2*nxv
               q2(j,k-ksh) = q2(j,k)
            enddo
         enddo
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_IRECV(q2(1,kyp2-kr+1),2*nxv*ksh,mreal,ks+1,     &
     &moff+4,lgrp,msid,ierr)
            endif
         endif
         if (kl > 0) then
            if (ks > 0) then
               call MPI_SEND(q2s,2*nxv*ksh,mreal,ks-1,moff+4,lgrp,ierr)
            endif
         endif
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_WAIT(msid,istatus,ierr)
            endif
         endif
      endif
!
! then create second half of odd array
!     q2(nx+j+1,k+1) = -q(nx-j+1,k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if ((kk==1).or.(kk==(ny+1))) then
            do j = 1, nx
               q2(j,k) = 0.0
               q2(j+nx,k) = 0.0
            enddo
         else if (kk <= ny) then
            do j = 1, nxs
               q2(nx+j+1,k) = -q2(nx-j+1,k)
            enddo
            q2(1,k) = 0.0
            q2(nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
!
! finish odd array
!     q2(j+1,ny+k+1) = -q(j+1,ny-k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if (kk > (ny+1)) then
            do j = 1, nxs
               q2(nx-j+1,k) = -q2(nx+j+1,k)
            enddo
            q2(nx+1,k) = -q2(nx+1,k)
         endif
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2B(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd&
     &,kyp2d,kshd)
! this subroutine creates a doubled vector array cu2 from a vector array
! cu, so that various 2d sine/cosine transforms can be performed with a
! 2d real to complex fft.  the x component is an odd function in y,
! y component is an odd function in x, and the z component is an odd
! function in both x and y.
! used for field solvers with dirichlet boundary conditions using ffts
! for distributed data
! input data cu must have a uniform partition
! cu array is temporarily modified and restored
! cus, cu2s = scratch arrays for shifts
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = second dimension of input array cu, must be >= nx+1
! kyp = number of data values per block in y
! kypd = third dimension of input array cu, must be >= kyp+1
! kyp2d = third dimension of output array cu2, must be >= 2*kyp
! kshd = third dimension scratch cu2s, kshd >= ksh, described below
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(3,nxv,kypd), intent(inout) :: cu
      real, dimension(3,2*nxv,kyp2d), intent(inout) :: cu2
      real, dimension(3,nxv), intent(inout) :: cus
      real, dimension(3,2*nxv,kshd), intent(inout) :: cu2s
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: i, j, k, nxs, nys, kyb, ny2, kyp2, kyb2, ks, moff, kyp1
      integer :: kyps, kyp2s, ksh, koff, l1, l2, m1, k1off, k2off, kypl
      integer :: kypr, k1, k2, joff, kypn, kypm, kk, nym, kl, kr, nnxv
      integer, dimension(lstat) :: istatus
      nxs = nx - 1
      nys = ny - 1
      nnxv = 3*nxv
! kyb = minimum number of processors used by input array
      kyb = (ny - 1)/kyp + 1
      ny2 = ny + ny
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny2 - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for received messages
      kyp1 = kyp + 1
! kyps = actual size used in input array on each node
      kyps = min(kyp,max(0,ny-kyp*ks))
! kyp2s = minimum positive actual size used in doubled array
      kyp2s = min(kyp2,max(0,ny2-kyp2*(kyb2-1)))
! ksh = final shift of received data
      ksh = kyp2 - kyp2s
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nys
            do j = 1, nxs
               cu2(1,j+1,k+1) = cu(1,j+1,k+1)
               cu2(2,j+1,k+1) = cu(2,j+1,k+1)
               cu2(3,j+1,k+1) = cu(3,j+1,k+1)
               cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
               cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
               cu2(3,nx+j+1,k+1) = -cu(3,nx-j+1,k+1)
               cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
               cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
               cu2(3,j+1,ny+k+1) = -cu(3,j+1,ny-k+1)
               cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
               cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
               cu2(3,nx+j+1,ny+k+1) = cu(3,nx-j+1,ny-k+1)
            enddo
            cu2(1,1,k+1) = cu(1,1,k+1)
            cu2(2,1,k+1) = 0.0
            cu2(3,1,k+1) = 0.0
            cu2(1,nx+1,k+1) = cu(1,nx+1,k+1)
            cu2(2,nx+1,k+1) = 0.0
            cu2(3,nx+1,k+1) = 0.0
            cu2(1,1,k+ny+1) = -cu(1,1,ny-k+1)
            cu2(2,1,k+ny+1) = 0.0
            cu2(3,1,k+ny+1) = 0.0
            cu2(1,nx+1,k+ny+1) = -cu(1,nx+1,ny-k+1)
            cu2(2,nx+1,k+ny+1) = 0.0
            cu2(3,nx+1,k+ny+1) = 0.0
         enddo
!$OMP END PARALLEL DO
         do j = 1, nxs
            cu2(1,j+1,1) = 0.0
            cu2(2,j+1,1) = cu(2,j+1,1)
            cu2(3,j+1,1) = 0.0
            cu2(1,j+nx+1,1) = 0.0
            cu2(2,j+nx+1,1) = -cu(2,nx-j+1,1)
            cu2(3,j+nx+1,1) = 0.0
            cu2(1,j+1,ny+1) = 0.0
            cu2(2,j+1,ny+1) = cu(2,j+1,ny+1)
            cu2(3,j+1,ny+1) = 0.0
            cu2(1,j+nx+1,ny+1) = 0.0
            cu2(2,j+nx+1,ny+1) = -cu(2,nx-j+1,ny+1)
            cu2(3,j+nx+1,ny+1) = 0.0
         enddo
         cu2(1,1,1) = 0.0
         cu2(2,1,1) = 0.0
         cu2(3,1,1) = 0.0
         cu2(1,nx+1,1) = 0.0
         cu2(2,nx+1,1) = 0.0
         cu2(3,nx+1,1) = 0.0
         cu2(1,1,ny+1) = 0.0
         cu2(2,1,ny+1) = 0.0
         cu2(3,1,ny+1) = 0.0
         cu2(1,nx+1,ny+1) = 0.0
         cu2(2,nx+1,ny+1) = 0.0
         cu2(3,nx+1,ny+1) = 0.0
         return
      endif
!
! copy to double array in x direction
!     cu2(1,j+1,k+1) = cu(1,j+1,k+1)
!     cu2(2,j+1,k+1) = cu(2,j+1,k+1)
!     cu2(3,j+1,k+1) = cu(3,j+1,k+1)
      koff = kyp2*ks
! l1 = source location for data for first half
      l1 = koff/kyp
! l2 = source location for data for second half
      l2 = l1 + 1
      k1off = kyp*l1 - koff + 1
      k2off = kyp*l2 - koff + 1
      koff = kyp*ks
! m1 = destination location of original data
      m1 = koff/kyp2
      koff = koff - kyp2*m1
! kypl = room available at destination
      kypl = min(kyps,kyp2-koff)
! this segment is used for mpi computers
!
! receive first half
      if (l1 < nvp) then
         call MPI_IRECV(cu2,kyp*nnxv,mreal,l1,moff+1,lgrp,msid,ierr)
! receive second half
         if (l2 < nvp) then
            call MPI_IRECV(cu2(1,1,kyp+1),kyp*nnxv,mreal,l2,moff+1,lgrp,&
     &nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(cu,kypl*nnxv,mreal,m1,moff+1,lgrp,ierr)
      endif
! wait for data and unpack it
      if (l1 < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
         kypr = kypr/nnxv
         do k = 2, kypr
            k1 = kypr - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do j = 1, nxv
               do i = 1, 3
                  cu2(i,j,k1) = cu2(i,j+joff,k2)
               enddo
            enddo
         enddo
         if (l2 < nvp) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 3
                     cu2(i,j,k1+kyp) = cu2(i,j+joff,k2+kyp)
                  enddo
               enddo
            enddo
         endif
      endif
!
! copy to double array in y direction
!     cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
!     cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
!     cu2(3,nx+j+1,ny+k+1) = cu(3,nx-j+1,ny-k+1)
      koff = kyp2*ks
! l1 = source location for data for second half
      l1 = (ny2 + ksh - koff - kyp2)/kyp
! l2 = source location for data for first half
      l2 = l1 + 1
! kypn = check if address for l1 receive is inbounds
      kypn = min(kyp,max(0,koff+ksh-ny+kyp1+kyp-1))
! kypm = check if address for l2 receive is inbounds
      kypm = min(kyp,max(0,koff+ksh-ny+kyp-1))
      koff = kyp*ks
! m1 = destination location of original data
      m1 = (ny2 + ksh - koff - 1)/kyp2
! amount of data available to send
      kypl = min(kyps,max(0,ny-koff-1))
! this segment is used for mpi computers
!
! first align data in input cu by left shifting
      do j = 1, nxv
         do i = 1, 3
            cus(i,j) = cu(i,j,1)
         enddo
      enddo
      do k = 2, kyp
         do j = 1, nxv
            do i = 1, 3
               cu(i,j,k-1) = cu(i,j,k)
            enddo
         enddo
      enddo
      if (ks < (nvp-1)) then
         call MPI_IRECV(cu(1,1,kyp),nnxv,mreal,ks+1,moff+2,lgrp,msid,   &
     &ierr)
      endif
      if (ks > 0) then
         call MPI_SEND(cus,nnxv,mreal,ks-1,moff+2,lgrp,ierr)
      endif
      if (ks < (nvp-1)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then continue to copy to double array in y direction
! post receive second half
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_IRECV(cu2(1,1,kyp1),kyp*nnxv,mreal,l1,moff+3,lgrp, &
     &msid,ierr)
         endif
      endif
! post receive first half
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_IRECV(cu2,kyp*nnxv,mreal,l2,moff+3,lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(cu,kypl*nnxv,mreal,m1,moff+3,lgrp,ierr)
      endif
! wait for data and unpack it
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 3
                     cu2(i,j,k1+kyp1-1) = cu2(i,j+joff,k2+kyp1-1)
                  enddo
               enddo
            enddo
         endif
      endif
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 3
                     cu2(i,j,k1) = cu2(i,j+joff,k2)
                  enddo
               enddo
            enddo
         endif
      endif
!
! restore original data in input cu
      do k = 2, kyp
         kk = kyp - k + 2
         do j = 1, nxv
            do i = 1, 3
               cu(i,j,kk) = cu(i,j,kk-1)
            enddo
         enddo
      enddo
      do j = 1, nxv
         do i = 1, 3
            cu(i,j,1) = cus(i,j)
         enddo
      enddo
!
! create mixed even/odd array
! first reverse local y index
      koff = kyp2*ks
      nym = ny2 - kyb*kyp
      do k = 1, kyp2
         kk = k + koff - ksh
         if (kk > nym) then
            if (k <= kyp) then
               k1 = kyp + 1
               do j = 1, nxs
                  cu2(1,nx+j+1,k1-k) = -cu2(1,nx-j+1,k)
                  cu2(2,nx+j+1,k1-k) = -cu2(2,nx-j+1,k)
                  cu2(3,nx+j+1,k1-k) = cu2(3,nx-j+1,k)
               enddo
            else
               k1 = kyp2 + kyp + 1
               do j = 1, nxs
                  cu2(1,nx+j+1,k1-k) = -cu2(1,nx-j+1,k)
                  cu2(2,nx+j+1,k1-k) = -cu2(2,nx-j+1,k)
                  cu2(3,nx+j+1,k1-k) = cu2(3,nx-j+1,k)
               enddo
            endif
            cu2(1,2*nx+1,k1-k) = -cu2(1,1,k)
            cu2(1,2*nx+2,k1-k) = -cu2(1,nx+1,k)
         endif
         cu2(1,1,k) = -cu2(1,1,k)
         cu2(2,1,k) = 0.0
         cu2(3,1,k) = 0.0
         cu2(1,nx+1,k) = -cu2(1,nx+1,k)
         cu2(2,nx+1,k) = 0.0
         cu2(3,nx+1,k) = 0.0
      enddo
!
! ny+1 point is special
! l1 = source processor id for ny+1 point
      l1 = (ny - 1)/kyp
! m1 = destination processor id of ny+1 point
      m1 = ny/kyp2
! k2 = destination location for ny+1 point
      k2 = ny - kyp2*m1 + 1
      if (ks==m1) then
         call MPI_IRECV(cu2(1,1,k2),nnxv,mreal,l1,moff+6,lgrp,msid,ierr)
      endif
      if (ks==l1) then
         call MPI_SEND(cu(1,1,kyps+1),nnxv,mreal,m1,moff+6,lgrp,ierr)
      endif
      if (ks==m1) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then left shift rhs of output cu2 by ksh y indices
      if (ksh > 0) then
         nym = ny + 1
! kr = number of data points which arriving from right processor
         kr = min(ksh,max(0,koff+kyp2-nym))
! kl = number of data points which need to be moved to left processor
         kl = min(ksh,max(0,koff-nym))
! kk = lower bound for left shift
         kk = max(0,nym-koff)
         if (kl > 0) then
! copy data going to left processor
            do k = 1, kl
               do j = 1, 2*nxv
                  do i = 1, 3
                     cu2s(i,j,k) = cu2(i,j,k+ksh-kl)
                  enddo
               enddo
            enddo
         endif
! shift remaining data
         do k = ksh + kk + 1, kyp2
            do j = 1, 2*nxv
               do i = 1, 3
                  cu2(i,j,k-ksh) = cu2(i,j,k)
               enddo
            enddo
         enddo
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_IRECV(cu2(1,1,kyp2-kr+1),2*nnxv*ksh,mreal,ks+1, &
     &moff+4,lgrp,msid,ierr)
            endif
         endif
         if (kl > 0) then
            if (ks > 0) then
               call MPI_SEND(cu2s,2*nnxv*ksh,mreal,ks-1,moff+4,lgrp,ierr&
     &)
            endif
         endif
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_WAIT(msid,istatus,ierr)
            endif
         endif
      endif
!
! then create second half of mixed even/odd array
!     cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
!     cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
!     cu2(3,nx+j+1,k+1) = -cu(3,nx-j+1,k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if ((kk==1).or.(kk==(ny+1))) then
            do j = 1, nxs
               cu2(1,j+1,k) = 0.0
               cu2(2,j+1,k) = cu2(2,j+1,k)
               cu2(3,j+1,k) = 0.0
               cu2(1,j+nx+1,k) = 0.0
               cu2(2,j+nx+1,k) = -cu2(2,nx-j+1,k)
               cu2(3,j+nx+1,k) = 0.0
            enddo
            cu2(1,1,k) = 0.0
            cu2(2,1,k) = 0.0
            cu2(3,1,k) = 0.0
            cu2(1,nx+1,k) = 0.0
            cu2(2,nx+1,k) = 0.0
            cu2(3,nx+1,k) = 0.0
         else if (kk <= ny) then
            do j = 1, nxs
               cu2(1,nx+j+1,k) = cu2(1,nx-j+1,k)
               cu2(2,nx+j+1,k) = -cu2(2,nx-j+1,k)
               cu2(3,nx+j+1,k) = -cu2(3,nx-j+1,k)
            enddo
            cu2(1,1,k) = -cu2(1,1,k)
            cu2(2,1,k) = 0.0
            cu2(3,1,k) = 0.0
            cu2(1,nx+1,k) = -cu2(1,nx+1,k)
            cu2(2,nx+1,k) = 0.0
            cu2(3,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
!
! finish mixed even/odd array
!     cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
!     cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
!     cu2(3,j+1,ny+k+1) = -cu(3,j+1,ny-k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if (kk > (ny+1)) then
            do j = 1, nxs
               cu2(1,nx-j+1,k) = cu2(1,nx+j+1,k)
               cu2(2,nx-j+1,k) = -cu2(2,nx+j+1,k)
               cu2(3,nx-j+1,k) = -cu2(3,nx+j+1,k)
            enddo
            cu2(1,1,k) = cu2(1,2*nx+1,k)
            cu2(1,nx+1,k) = cu2(1,2*nx+2,k)
            cu2(2,nx+1,k) = 0.0
            cu2(3,nx+1,k) = -cu2(3,nx+1,k)
         endif
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2M(amu,amu2,amus,amu2s,nx,ny,kstrt,nvp,nxv,kyp,&
     &kypd,kyp2d,kshd)
! this subroutine creates a doubled tensor array amu2 from a tensor array
! amu, so that various 2d sine/cosine transforms can be performed with a
! 2d real to complex fft.  the xx-yy component is an odd function in
! both x and y, the xy component is an even function in both x and y, the
! zx component is an even function in x and an odd function in y, the zy
! component is an odd function in x and and even function in y.
! used for field solvers with dirichlet boundary conditions using ffts
! for distributed data
! input data amu must have a uniform partition
! amu array is temporarily modified and restored
! amus, amu2s = scratch arrays for shifts
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = second dimension of input array cu, must be >= nx+1
! kyp = number of data values per block in y
! kypd = third dimension of input array cu, must be >= kyp+1
! kyp2d = third dimension of output array cu2, must be >= 2*kyp
! kshd = third dimension scratch cu2s, kshd >= ksh, described below
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(4,nxv,kypd), intent(inout) :: amu
      real, dimension(4,2*nxv,kyp2d), intent(inout) :: amu2
      real, dimension(4,nxv), intent(inout) :: amus
      real, dimension(4,2*nxv,kshd), intent(inout) :: amu2s
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: i, j, k, nxs, nys, kyb, ny2, kyp2, kyb2, ks, moff, kyp1
      integer :: kyps, kyp2s, ksh, koff, l1, l2, m1, k1off, k2off, kypl
      integer :: kypr, k1, k2, joff, kypn, kypm, kk, nym, kl, kr, nnxv
      integer, dimension(lstat) :: istatus
      nxs = nx - 1
      nys = ny - 1
      nnxv = 4*nxv
! kyb = minimum number of processors used by input array
      kyb = (ny - 1)/kyp + 1
      ny2 = ny + ny
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny2 - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for received messages
      kyp1 = kyp + 1
! kyps = actual size used in input array on each node
      kyps = min(kyp,max(0,ny-kyp*ks))
! kyp2s = minimum positive actual size used in doubled array
      kyp2s = min(kyp2,max(0,ny2-kyp2*(kyb2-1)))
! ksh = final shift of received data
      ksh = kyp2 - kyp2s
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nys
            do j = 1, nxs
               amu2(1,j+1,k+1) = amu(1,j+1,k+1)
               amu2(2,j+1,k+1) = amu(2,j+1,k+1)
               amu2(3,j+1,k+1) = amu(3,j+1,k+1)
               amu2(4,j+1,k+1) = amu(4,j+1,k+1)
               amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
               amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
               amu2(3,nx+j+1,k+1) = amu(3,nx-j+1,k+1)
               amu2(4,nx+j+1,k+1) = -amu(4,nx-j+1,k+1)
               amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
               amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
               amu2(3,j+1,ny+k+1) = -amu(3,j+1,ny-k+1)
               amu2(4,j+1,ny+k+1) = amu(4,j+1,ny-k+1)
               amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
               amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
               amu2(3,nx+j+1,ny+k+1) = -amu(3,nx-j+1,ny-k+1)
               amu2(4,nx+j+1,ny+k+1) = -amu(4,nx-j+1,ny-k+1)
            enddo
            amu2(1,1,k+1) = 0.0
            amu2(2,1,k+1) = amu(2,1,k+1)
            amu2(3,1,k+1) = amu(3,1,k+1)
            amu2(4,1,k+1) = 0.0
            amu2(1,nx+1,k+1) = 0.0
            amu2(2,nx+1,k+1) = amu(2,nx+1,k+1)
            amu2(3,nx+1,k+1) = amu(3,nx+1,k+1)
            amu2(4,nx+1,k+1) = 0.0
            amu2(1,1,k+ny+1) = 0.0
            amu2(2,1,k+ny+1) = amu(2,1,ny-k+1)
            amu2(3,1,k+ny+1) = -amu(3,1,ny-k+1)
            amu2(4,1,k+ny+1) = 0.0
            amu2(1,nx+1,k+ny+1) = 0.0
            amu2(2,nx+1,k+ny+1) = amu(2,nx+1,ny-k+1)
            amu2(3,nx+1,k+ny+1) = -amu(3,nx+1,ny-k+1)
            amu2(4,nx+1,k+ny+1) = 0.0
         enddo
!$OMP END PARALLEL DO
         do j = 1, nxs
            amu2(1,j+1,1) = 0.0
            amu2(2,j+1,1) = amu(2,j+1,1)
            amu2(3,j+1,1) = 0.0
            amu2(4,j+1,1) = amu(4,j+1,1)
            amu2(1,j+nx+1,1) = 0.0
            amu2(2,nx+j+1,1) = amu(2,nx-j+1,1)
            amu2(3,j+nx+1,1) = 0.0
            amu2(4,j+nx+1,1) = -amu(4,nx-j+1,1)
            amu2(1,j+1,ny+1) = 0.0
            amu2(2,j+1,ny+1) = amu(2,j+1,ny+1)
            amu2(3,j+1,ny+1) = 0.0
            amu2(4,j+1,ny+1) = amu(4,j+1,ny+1)
            amu2(1,j+nx+1,ny+1) = 0.0
            amu2(2,nx+j+1,ny+1) = amu(2,nx-j+1,ny+1)
            amu2(3,j+nx+1,ny+1) = 0.0
            amu2(4,j+nx+1,ny+1) = -amu(4,nx-j+1,ny+1)
         enddo
         amu2(1,1,1) = 0.0
         amu2(2,1,1) = amu(2,1,1)
         amu2(3,1,1) = 0.0
         amu2(4,1,1) = 0.0
         amu2(1,nx+1,1) = 0.0
         amu2(2,nx+1,1) = amu(2,nx+1,1)
         amu2(3,nx+1,1) = 0.0
         amu2(4,nx+1,1) = 0.0
         amu2(1,1,ny+1) = 0.0
         amu2(2,1,ny+1) = amu(2,1,ny+1)
         amu2(3,1,ny+1) = 0.0
         amu2(4,1,ny+1) = 0.0
         amu2(1,nx+1,ny+1) = 0.0
         amu2(2,nx+1,ny+1) = amu(2,nx+1,ny+1)
         amu2(3,nx+1,ny+1) = 0.0
         amu2(4,nx+1,ny+1) = 0.0
         return
      endif
!
! copy to double array in x direction
!     amu2(1,j+1,k+1) = amu(1,j+1,k+1)
!     amu2(2,j+1,k+1) = amu(2,j+1,k+1)
!     amu2(3,j+1,k+1) = amu(3,j+1,k+1)
!     amu2(4,j+1,k+1) = amu(4,j+1,k+1)
      koff = kyp2*ks
! l1 = source location for data for first half
      l1 = koff/kyp
! l2 = source location for data for second half
      l2 = l1 + 1
      k1off = kyp*l1 - koff + 1
      k2off = kyp*l2 - koff + 1
      koff = kyp*ks
! m1 = destination location of original data
      m1 = koff/kyp2
      koff = koff - kyp2*m1
! kypl = room available at destination
      kypl = min(kyps,kyp2-koff)
! this segment is used for mpi computers
!
! receive first half
      if (l1 < nvp) then
         call MPI_IRECV(amu2,kyp*nnxv,mreal,l1,moff+1,lgrp,msid,ierr)
! receive second half
         if (l2 < nvp) then
            call MPI_IRECV(amu2(1,1,kyp+1),kyp*nnxv,mreal,l2,moff+1,    &
     &lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(amu,kypl*nnxv,mreal,m1,moff+1,lgrp,ierr)
      endif
! wait for data and unpack it
      if (l1 < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
         kypr = kypr/nnxv
         do k = 2, kypr
            k1 = kypr - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do j = 1, nxv
               do i = 1, 4
                  amu2(i,j,k1) = amu2(i,j+joff,k2)
               enddo
            enddo
         enddo
         if (l2 < nvp) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 4
                     amu2(i,j,k1+kyp) = amu2(i,j+joff,k2+kyp)
                  enddo
               enddo
            enddo
         endif
      endif
!
! copy to double array in y direction
!     amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
!     amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
!     amu2(3,nx+j+1,ny+k+1) = -amu(3,nx-j+1,ny-k+1)
!     amu2(4,nx+j+1,ny+k+1) = -amu(4,nx-j+1,ny-k+1)
      koff = kyp2*ks
! l1 = source location for data for second half
      l1 = (ny2 + ksh - koff - kyp2)/kyp
! l2 = source location for data for first half
      l2 = l1 + 1
! kypn = check if address for l1 receive is inbounds
      kypn = min(kyp,max(0,koff+ksh-ny+kyp1+kyp-1))
! kypm = check if address for l2 receive is inbounds
      kypm = min(kyp,max(0,koff+ksh-ny+kyp-1))
      koff = kyp*ks
! m1 = destination location of original data
      m1 = (ny2 + ksh - koff - 1)/kyp2
! amount of data available to send
      kypl = min(kyps,max(0,ny-koff-1))
! this segment is used for mpi computers
!
! first align data in input cu by left shifting
      do j = 1, nxv
         do i = 1, 4
            amus(i,j) = amu(i,j,1)
         enddo
      enddo
      do k = 2, kyp
         do j = 1, nxv
            do i = 1, 4
               amu(i,j,k-1) = amu(i,j,k)
            enddo
         enddo
      enddo
      if (ks < (nvp-1)) then
         call MPI_IRECV(amu(1,1,kyp),nnxv,mreal,ks+1,moff+2,lgrp,msid,  &
     &ierr)
      endif
      if (ks > 0) then
         call MPI_SEND(amus,nnxv,mreal,ks-1,moff+2,lgrp,ierr)
      endif
      if (ks < (nvp-1)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then continue to copy to double array in y direction
! post receive second half
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_IRECV(amu2(1,1,kyp1),kyp*nnxv,mreal,l1,moff+3,lgrp,&
     &msid,ierr)
         endif
      endif
! post receive first half
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_IRECV(amu2,kyp*nnxv,mreal,l2,moff+3,lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(amu,kypl*nnxv,mreal,m1,moff+3,lgrp,ierr)
      endif
! wait for data and unpack it
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 4
                     amu2(i,j,k1+kyp1-1) = amu2(i,j+joff,k2+kyp1-1)
                  enddo
               enddo
            enddo
         endif
      endif
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 4
                     amu2(i,j,k1) = amu2(i,j+joff,k2)
                  enddo
               enddo
            enddo
         endif
      endif
!
! restore original data in input amu
      do k = 2, kyp
         kk = kyp - k + 2
         do j = 1, nxv
            do i = 1, 4
               amu(i,j,kk) = amu(i,j,kk-1)
            enddo
         enddo
      enddo
      do j = 1, nxv
         do i = 1, 4
            amu(i,j,1) = amus(i,j)
         enddo
      enddo
!
! create mixed even/odd array
! first reverse local y index
      koff = kyp2*ks
      nym = ny2 - kyb*kyp
      do k = 1, kyp2
         kk = k + koff - ksh
         if (kk > nym) then
            if (k <= kyp) then
               k1 = kyp + 1
               do j = 1, nxs
                  amu2(1,nx+j+1,k1-k) = amu2(1,nx-j+1,k)
                  amu2(2,nx+j+1,k1-k) = amu2(2,nx-j+1,k)
                  amu2(3,nx+j+1,k1-k) = -amu2(3,nx-j+1,k)
                  amu2(4,nx+j+1,k1-k) = -amu2(4,nx-j+1,k)
               enddo
            else
               k1 = kyp2 + kyp + 1
               do j = 1, nxs
                  amu2(1,nx+j+1,k1-k) = amu2(1,nx-j+1,k)
                  amu2(2,nx+j+1,k1-k) = amu2(2,nx-j+1,k)
                  amu2(3,nx+j+1,k1-k) = -amu2(3,nx-j+1,k)
                  amu2(4,nx+j+1,k1-k) = -amu2(4,nx-j+1,k)

               enddo
            endif
            amu2(2,2*nx+1,k1-k) = amu2(2,1,k)
            amu2(2,2*nx+2,k1-k) = amu2(2,nx+1,k)
            amu2(3,2*nx+1,k1-k) = -amu2(3,1,k)
            amu2(3,2*nx+2,k1-k) = -amu2(3,nx+1,k)
         endif
         amu2(1,1,k) = 0.0
         amu2(3,1,k) = -amu2(3,1,k)
         amu2(4,1,k) = 0.0
         amu2(1,nx+1,k) = 0.0
         amu2(3,nx+1,k) = -amu2(3,nx+1,k)
         amu2(4,nx+1,k) = 0.0
      enddo
!
! ny+1 point is special
! l1 = source processor id for ny+1 point
      l1 = (ny - 1)/kyp
! m1 = destination processor id of ny+1 point
      m1 = ny/kyp2
! k2 = destination location for ny+1 point
      k2 = ny - kyp2*m1 + 1
      if (ks==m1) then
         call MPI_IRECV(amu2(1,1,k2),nnxv,mreal,l1,moff+6,lgrp,msid,ierr&
     &)
      endif
      if (ks==l1) then
         call MPI_SEND(amu(1,1,kyps+1),nnxv,mreal,m1,moff+6,lgrp,ierr)
      endif
      if (ks==m1) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then left shift rhs of output cu2 by ksh y indices
      if (ksh > 0) then
         nym = ny + 1
! kr = number of data points which arriving from right processor
         kr = min(ksh,max(0,koff+kyp2-nym))
! kl = number of data points which need to be moved to left processor
         kl = min(ksh,max(0,koff-nym))
! kk = lower bound for left shift
         kk = max(0,nym-koff)
         if (kl > 0) then
! copy data going to left processor
            do k = 1, kl
               do j = 1, 2*nxv
                  do i = 1, 4
                     amu2s(i,j,k) = amu2(i,j,k+ksh-kl)
                  enddo
               enddo
            enddo
         endif
! shift remaining data
         do k = ksh + kk + 1, kyp2
            do j = 1, 2*nxv
               do i = 1, 4
                  amu2(i,j,k-ksh) = amu2(i,j,k)
               enddo
            enddo
         enddo
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_IRECV(amu2(1,1,kyp2-kr+1),2*nnxv*ksh,mreal,ks+1,&
     &moff+4,lgrp,msid,ierr)
            endif
         endif
         if (kl > 0) then
            if (ks > 0) then
               call MPI_SEND(amu2s,2*nnxv*ksh,mreal,ks-1,moff+4,lgrp,   &
     &ierr)
            endif
         endif
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_WAIT(msid,istatus,ierr)
            endif
         endif
      endif
!
! then create second half of mixed even/odd array
!     amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
!     amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
!     amu2(3,nx+j+1,k+1) = amu(3,nx-j+1,k+1)
!     amu2(4,nx+j+1,k+1) = -amu(4,nx-j+1,k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if ((kk==1).or.(kk==(ny+1))) then
            do j = 1, nxs
               amu2(1,j+1,k) = 0.0
               amu2(3,j+1,k) = 0.0
               amu2(4,j+1,k) = amu2(4,j+1,k)
               amu2(1,j+nx+1,k) = 0.0
               amu2(2,j+nx+1,k) = amu2(2,nx-j+1,k)
               amu2(3,j+nx+1,k) = 0.0
               amu2(4,j+nx+1,k) = -amu2(4,nx-j+1,k)
            enddo
            amu2(1,1,k) = 0.0
            amu2(3,1,k) = 0.0
            amu2(4,1,k) = 0.0
            amu2(1,nx+1,k) = 0.0
            amu2(3,nx+1,k) = 0.0
            amu2(4,nx+1,k) = 0.0
         else if (kk <= ny) then
            do j = 1, nxs
               amu2(1,nx+j+1,k) = -amu2(1,nx-j+1,k)
               amu2(2,nx+j+1,k) = amu2(2,nx-j+1,k)
               amu2(3,nx+j+1,k) = amu2(3,nx-j+1,k)
               amu2(4,nx+j+1,k) = -amu2(4,nx-j+1,k)
            enddo
            amu2(1,1,k) = 0.0
            amu2(3,1,k) = -amu2(3,1,k)
            amu2(4,1,k) = 0.0
            amu2(1,nx+1,k) = 0.0
            amu2(3,nx+1,k) = -amu2(3,nx+1,k)
            amu2(4,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
!
! finish mixed even/odd array
!     amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
!     amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
!     amu2(3,j+1,ny+k+1) = -amu(3,j+1,ny-k+1)
!     amu2(4,j+1,ny+k+1) = amu(4,j+1,ny-k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if (kk > (ny+1)) then
            do j = 1, nxs
               amu2(1,nx-j+1,k) = -amu2(1,nx+j+1,k)
               amu2(2,nx-j+1,k) = amu2(2,nx+j+1,k)
               amu2(3,nx-j+1,k) = amu2(3,nx+j+1,k)
               amu2(4,nx-j+1,k) = -amu2(4,nx+j+1,k)
            enddo
            amu2(1,nx+1,k) = -amu2(1,nx+1,k)
            amu2(2,1,k) = amu2(2,2*nx+1,k)
            amu2(2,nx+1,k) = amu2(2,2*nx+2,k)
            amu2(3,1,k) = amu2(3,2*nx+1,k)
            amu2(3,nx+1,k) = amu2(3,2*nx+2,k)
            amu2(4,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN22M(amu,amu2,amus,amu2s,nx,ny,kstrt,nvp,nxv,kyp&
     &,kypd,kyp2d,kshd)
! this subroutine creates a doubled tensor array amu2 from a tensor array
! amu, so that various 2d sine/cosine transforms can be performed with a
! 2d real to complex fft.  the xx-yy component is an odd function in
! both x and y, the xy component is an even function in both x and y.
! used for field solvers with dirichlet boundary conditions using ffts
! for distributed data
! input data amu must have a uniform partition
! amu array is temporarily modified and restored
! amus, amu2s = scratch arrays for shifts
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = second dimension of input array cu, must be >= nx+1
! kyp = number of data values per block in y
! kypd = third dimension of input array cu, must be >= kyp+1
! kyp2d = third dimension of output array cu2, must be >= 2*kyp
! kshd = third dimension scratch cu2s, kshd >= ksh, described below
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(2,nxv,kypd), intent(inout) :: amu
      real, dimension(2,2*nxv,kyp2d), intent(inout) :: amu2
      real, dimension(2,nxv), intent(inout) :: amus
      real, dimension(2,2*nxv,kshd), intent(inout) :: amu2s
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: i, j, k, nxs, nys, kyb, ny2, kyp2, kyb2, ks, moff, kyp1
      integer :: kyps, kyp2s, ksh, koff, l1, l2, m1, k1off, k2off, kypl
      integer :: kypr, k1, k2, joff, kypn, kypm, kk, nym, kl, kr, nnxv
      integer, dimension(lstat) :: istatus
      nxs = nx - 1
      nys = ny - 1
      nnxv = 2*nxv
! kyb = minimum number of processors used by input array
      kyb = (ny - 1)/kyp + 1
      ny2 = ny + ny
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny2 - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for received messages
      kyp1 = kyp + 1
! kyps = actual size used in input array on each node
      kyps = min(kyp,max(0,ny-kyp*ks))
! kyp2s = minimum positive actual size used in doubled array
      kyp2s = min(kyp2,max(0,ny2-kyp2*(kyb2-1)))
! ksh = final shift of received data
      ksh = kyp2 - kyp2s
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nys
            do j = 1, nxs
               amu2(1,j+1,k+1) = amu(1,j+1,k+1)
               amu2(2,j+1,k+1) = amu(2,j+1,k+1)
               amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
               amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
               amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
               amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
               amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
               amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
            enddo
            amu2(1,1,k+1) = 0.
            amu2(2,1,k+1) = amu(2,1,k+1)
            amu2(1,nx+1,k+1) = 0.
            amu2(2,nx+1,k+1) = amu(2,nx+1,k+1)
            amu2(1,1,k+ny+1) = 0.
            amu2(2,1,k+ny+1) = amu(2,1,ny-k+1)
            amu2(1,nx+1,k+ny+1) = 0.
            amu2(2,nx+1,k+ny+1) = amu(2,nx+1,ny-k+1)
         enddo
!$OMP END PARALLEL DO
         do j = 1, nxs
            amu2(1,j+1,1) = 0.
            amu2(2,j+1,1) = amu(2,j+1,1)
            amu2(1,j+nx+1,1) = 0.
            amu2(2,nx+j+1,1) = amu(2,nx-j+1,1)
            amu2(1,j+1,ny+1) = 0.
            amu2(2,j+1,ny+1) = amu(2,j+1,ny+1)
            amu2(1,j+nx+1,ny+1) = 0.
            amu2(2,nx+j+1,ny+1) = amu(2,nx-j+1,ny+1)
         enddo
         amu2(1,1,1) = 0.
         amu2(2,1,1) = amu(2,1,1)
         amu2(1,nx+1,1) = 0.
         amu2(2,nx+1,1) = amu(2,nx+1,1)
         amu2(1,1,ny+1) = 0.
         amu2(2,1,ny+1) = amu(2,1,ny+1)
         amu2(1,nx+1,ny+1) = 0.
         amu2(2,nx+1,ny+1) = amu(2,nx+1,ny+1)
         return
      endif
!
! copy to double array in x direction
!     amu2(1,j+1,k+1) = amu(1,j+1,k+1)
!     amu2(2,j+1,k+1) = amu(2,j+1,k+1)
      koff = kyp2*ks
! l1 = source location for data for first half
      l1 = koff/kyp
! l2 = source location for data for second half
      l2 = l1 + 1
      k1off = kyp*l1 - koff + 1
      k2off = kyp*l2 - koff + 1
      koff = kyp*ks
! m1 = destination location of original data
      m1 = koff/kyp2
      koff = koff - kyp2*m1
! kypl = room available at destination
      kypl = min(kyps,kyp2-koff)
! this segment is used for mpi computers
!
! receive first half
      if (l1 < nvp) then
         call MPI_IRECV(amu2,kyp*nnxv,mreal,l1,moff+1,lgrp,msid,ierr)
! receive second half
         if (l2 < nvp) then
            call MPI_IRECV(amu2(1,1,kyp+1),kyp*nnxv,mreal,l2,moff+1,    &
     &lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(amu,kypl*nnxv,mreal,m1,moff+1,lgrp,ierr)
      endif
! wait for data and unpack it
      if (l1 < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
         kypr = kypr/nnxv
         do k = 2, kypr
            k1 = kypr - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do j = 1, nxv
               do i = 1, 2
                  amu2(i,j,k1) = amu2(i,j+joff,k2)
               enddo
            enddo
         enddo
         if (l2 < nvp) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     amu2(i,j,k1+kyp) = amu2(i,j+joff,k2+kyp)
                  enddo
               enddo
            enddo
         endif
      endif
!
! copy to double array in y direction
!     amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
!     amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
      koff = kyp2*ks
! l1 = source location for data for second half
      l1 = (ny2 + ksh - koff - kyp2)/kyp
! l2 = source location for data for first half
      l2 = l1 + 1
! kypn = check if address for l1 receive is inbounds
      kypn = min(kyp,max(0,koff+ksh-ny+kyp1+kyp-1))
! kypm = check if address for l2 receive is inbounds
      kypm = min(kyp,max(0,koff+ksh-ny+kyp-1))
      koff = kyp*ks
! m1 = destination location of original data
      m1 = (ny2 + ksh - koff - 1)/kyp2
! amount of data available to send
      kypl = min(kyps,max(0,ny-koff-1))
! this segment is used for mpi computers
!
! first align data in input cu by left shifting
      do j = 1, nxv
         do i = 1, 2
            amus(i,j) = amu(i,j,1)
         enddo
      enddo
      do k = 2, kyp
         do j = 1, nxv
            do i = 1, 2
               amu(i,j,k-1) = amu(i,j,k)
            enddo
         enddo
      enddo
      if (ks < (nvp-1)) then
         call MPI_IRECV(amu(1,1,kyp),nnxv,mreal,ks+1,moff+2,lgrp,msid,  &
     &ierr)
      endif
      if (ks > 0) then
         call MPI_SEND(amus,nnxv,mreal,ks-1,moff+2,lgrp,ierr)
      endif
      if (ks < (nvp-1)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then continue to copy to double array in y direction
! post receive second half
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_IRECV(amu2(1,1,kyp1),kyp*nnxv,mreal,l1,moff+3,lgrp,&
     &msid,ierr)
         endif
      endif
! post receive first half
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_IRECV(amu2,kyp*nnxv,mreal,l2,moff+3,lgrp,nsid,ierr)
         endif
      endif
! send data
      if (m1 < nvp) then
         call MPI_SEND(amu,kypl*nnxv,mreal,m1,moff+3,lgrp,ierr)
      endif
! wait for data and unpack it
      if ((l1 >= 0).and.(l1 < nvp)) then
         if (kypn > 0) then
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     amu2(i,j,k1+kyp1-1) = amu2(i,j+joff,k2+kyp1-1)
                  enddo
               enddo
            enddo
         endif
      endif
      if ((l2 >= 0).and.(l2 < nvp)) then
         if (kypm > 0) then
            call MPI_WAIT(nsid,istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,kypr,ierr)
            kypr = kypr/nnxv
            do k = 2, kypr
               k1 = kypr - k + 2
               k2 = (k1 - 1)/2 + 1
               joff = nxv*(k1 - 2*k2 + 1)
               do j = 1, nxv
                  do i = 1, 2
                     amu2(i,j,k1) = amu2(i,j+joff,k2)
                  enddo
               enddo
            enddo
         endif
      endif
!
! restore original data in input amu
      do k = 2, kyp
         kk = kyp - k + 2
         do j = 1, nxv
            do i = 1, 2
               amu(i,j,kk) = amu(i,j,kk-1)
            enddo
         enddo
      enddo
      do j = 1, nxv
         do i = 1, 2
            amu(i,j,1) = amus(i,j)
         enddo
      enddo
!
! create mixed even/odd array
! first reverse local y index
      koff = kyp2*ks
      nym = ny2 - kyb*kyp
      do k = 1, kyp2
         kk = k + koff - ksh
         if (kk > nym) then
            if (k <= kyp) then
               k1 = kyp + 1
               do j = 1, nxs
                  amu2(1,nx+j+1,k1-k) = amu2(1,nx-j+1,k)
                  amu2(2,nx+j+1,k1-k) = amu2(2,nx-j+1,k)
               enddo
            else
               k1 = kyp2 + kyp + 1
               do j = 1, nxs
                  amu2(1,nx+j+1,k1-k) = amu2(1,nx-j+1,k)
                  amu2(2,nx+j+1,k1-k) = amu2(2,nx-j+1,k)

               enddo
            endif
            amu2(2,2*nx+1,k1-k) = amu2(2,1,k)
            amu2(2,2*nx+2,k1-k) = amu2(2,nx+1,k)
         endif
         amu2(1,1,k) = 0.0
         amu2(1,nx+1,k) = 0.0
      enddo
!
! ny+1 point is special
! l1 = source processor id for ny+1 point
      l1 = (ny - 1)/kyp
! m1 = destination processor id of ny+1 point
      m1 = ny/kyp2
! k2 = destination location for ny+1 point
      k2 = ny - kyp2*m1 + 1
      if (ks==m1) then
         call MPI_IRECV(amu2(1,1,k2),nnxv,mreal,l1,moff+6,lgrp,msid,ierr&
     &)
      endif
      if (ks==l1) then
         call MPI_SEND(amu(1,1,kyps+1),nnxv,mreal,m1,moff+6,lgrp,ierr)
      endif
      if (ks==m1) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
!
! then left shift rhs of output cu2 by ksh y indices
      if (ksh > 0) then
         nym = ny + 1
! kr = number of data points which arriving from right processor
         kr = min(ksh,max(0,koff+kyp2-nym))
! kl = number of data points which need to be moved to left processor
         kl = min(ksh,max(0,koff-nym))
! kk = lower bound for left shift
         kk = max(0,nym-koff)
         if (kl > 0) then
! copy data going to left processor
            do k = 1, kl
               do j = 1, 2*nxv
                  do i = 1, 2
                     amu2s(i,j,k) = amu2(i,j,k+ksh-kl)
                  enddo
               enddo
            enddo
         endif
! shift remaining data
         do k = ksh + kk + 1, kyp2
            do j = 1, 2*nxv
               do i = 1, 2
                  amu2(i,j,k-ksh) = amu2(i,j,k)
               enddo
            enddo
         enddo
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_IRECV(amu2(1,1,kyp2-kr+1),2*nnxv*ksh,mreal,ks+1,&
     &moff+4,lgrp,msid,ierr)
            endif
         endif
         if (kl > 0) then
            if (ks > 0) then
               call MPI_SEND(amu2s,2*nnxv*ksh,mreal,ks-1,moff+4,lgrp,   &
     &ierr)
            endif
         endif
         if (kr > 0) then
            if (ks < (nvp-1)) then
               call MPI_WAIT(msid,istatus,ierr)
            endif
         endif
      endif
!
! then create second half of mixed even/odd array
!     amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
!     amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if ((kk==1).or.(kk==(ny+1))) then
            do j = 1, nxs
               amu2(1,j+1,k) = 0.0
               amu2(1,j+nx+1,k) = 0.0
               amu2(2,j+nx+1,k) = amu2(2,nx-j+1,k)
            enddo
            amu2(1,1,k) = 0.0
            amu2(1,nx+1,k) = 0.0
         else if (kk <= ny) then
            do j = 1, nxs
               amu2(1,nx+j+1,k) = -amu2(1,nx-j+1,k)
               amu2(2,nx+j+1,k) = amu2(2,nx-j+1,k)
            enddo
            amu2(1,1,k) = 0.0
            amu2(1,nx+1,k) = 0.0
         endif
      enddo
!$OMP END PARALLEL DO
!
! finish mixed even/odd array
!     amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
!     amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
!$OMP PARALLEL DO PRIVATE(j,k,kk)
      do k = 1, kyp2
         kk = k + koff
         if (kk > (ny+1)) then
            do j = 1, nxs
               amu2(1,nx-j+1,k) = -amu2(1,nx+j+1,k)
               amu2(2,nx-j+1,k) = amu2(2,nx+j+1,k)
            enddo
            amu2(1,nx+1,k) = -amu2(1,nx+1,k)
            amu2(2,1,k) = amu2(2,2*nx+1,k)
            amu2(2,nx+1,k) = amu2(2,2*nx+2,k)
         endif
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPHAFDBL2D(q,q2,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d)
! this subroutine copies data from a double array q2 to regular array q,
! with guard cells for scalar field and linear interpolation
! q and q2 have uniform partitions for distributed data
! input q2 array may be modified
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! nxv = first dimension of input array q, must be >= nx+1
! kyp = number of data values per block in y
! kypd = second dimension of input array q, must be >= kyp+1
! kyp2d = second dimension of output array q2, must be >= 2*kyp
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      real, intent(inout), dimension(nxv,kypd) :: q
      real, intent(inout), dimension(2*nxv,kyp2d) :: q2
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: j, k, nx1, ny1, kyb, kyp2, kyb2, ks, kyp1
      integer :: joff, koff, moff, kk, m1, m2, l1, l2, lmx, kypl
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
! kyb = minimum number of processors used by output array
      kyb = (ny - 1)/kyp + 1
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny + ny - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for sent messages
      kyp1 = kyp + 1
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ny1
            do j = 1, nx1
               q(j,k) = q2(j,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
!
! copy from a double array
!     q(j,k) = q2(j,k)
      koff = kyp2*ks
! m1 = destination location for data for first half
      m1 = koff/kyp
! kypl = room available at destination
      kypl = min(kyp,max(0,ny+1-kyp*m1))
! m2 = destination location for data for second half
      m2 = m1 + 1
      koff = kyp*ks
! l1 = source location for data
      l1 = koff/kyp2
! l2 = source location for guard cell data
      l2 = ny/kyp2
! lmx = maximum source location for main data
      lmx = (ny - 1)/kyp2
! this segment is used for mpi computers
!
      if (ks < kyb) then
         if (l1 <= lmx) then
            call MPI_IRECV(q,nxv*kyp,mreal,l1,moff+5,lgrp,msid,ierr)
! receive guard cells if on another processor
            if (ks==(kyb-1)) then
               if (l2 > l1) then
                  call MPI_IRECV(q(1,kyp1),nxv,mreal,l2,moff+5,lgrp,nsid&
     &,ierr)
               endif
            endif
         endif
      endif
! pack data and send first half
      if (m1 < kyb) then
         do k = 2, kyp1
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do j = 1, nxv
               q2(j+joff,kk) = q2(j,k)
            enddo
         enddo
         call MPI_SEND(q2,nxv*kypl,mreal,m1,moff+5,lgrp,ierr)
! send guard cell if on another processor
      else if (kypl > 0) then
         call MPI_SEND(q2,nxv*kypl,mreal,m1-1,moff+5,lgrp,ierr)
      endif
! pack data and send second half
      if (m2 < kyb) then
         do k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do j = 1, nxv
               q2(j+joff,kk+kyp) = q2(j,k+kyp)
            enddo
         enddo
         call MPI_SEND(q2(1,kyp1),nxv*kypl,mreal,m2,moff+5,lgrp,ierr)
      endif
! wait for data
      if (ks < kyb) then
         if (l1 <= lmx) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
! wait for guard cells if on another processor
         if (ks==(kyb-1)) then
            if (l2 > l1) then
               call MPI_WAIT(nsid,istatus,ierr)
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPHAFDBL2C(fxy,fxy2,nx,ny,kstrt,nvp,ndim,nxv,kyp,kypd, &
     &kyp2d)
! this subroutine copies data from a double array fxy2 to regular array
! fxy, with guard cells for vector field and linear interpolation
! fxy and fxy2 have uniform partitions for distributed data
! input fxy2 array may be modified
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors requested
! ndim = first dimension of input array fxy
! nxv = second dimension of input array fxy, must be >= nx+1
! kyp = number of data values per block in y
! kypd = third dimension of input array fxy, must be >= kyp+1
! kyp2d = third dimension of output array fxy2, must be >= 2*kyp
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, ndim, nxv, kyp
      integer, intent(in) :: kypd, kyp2d
      real, intent(inout), dimension(ndim,nxv,kypd) :: fxy
      real, intent(inout), dimension(ndim,2*nxv,kyp2d) :: fxy2
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: msid, nsid, ierr
      integer :: i, j, k, nx1, ny1, nnxv, kyb, kyp2, kyb2, ks, kyp1
      integer :: joff, koff, moff, kk, m1, m2, l1, l2, lmx, kypl
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
      nnxv = ndim*nxv
! kyb = minimum number of processors used by output array
      kyb = (ny - 1)/kyp + 1
! kyp2 = doubled array partition size
      kyp2 = 2*kyp
! kyb2 = minimum number of processors used by doubled array
      kyb2 = (ny + ny - 1)/kyp2 + 1
! ks = processor id
      ks = kstrt - 1
! moff = initial message tag
      moff = kypd + kyb
! kyp1 = offset for sent messages
      kyp1 = kyp + 1
!
! return if too many processors
      if (kstrt > ny) return
! error: dimension kyp2d too small
      if (kyp2 > kyp2d) return
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ny1
            do j = 1, nx1
               do i = 1, ndim
                  fxy(i,j,k) = fxy2(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
!
! copy from a double array
!     fxy(i,j,k) = fxy2(i,j,k)
      koff = kyp2*ks
! m1 = destination location for data for first half
      m1 = koff/kyp
! kypl = room available at destination
      kypl = min(kyp,max(0,ny+1-kyp*m1))
! m2 = destination location for data for second half
      m2 = m1 + 1
      koff = kyp*ks
! l1 = source location for data
      l1 = koff/kyp2
! l2 = source location for guard cell data
      l2 = ny/kyp2
! lmx = maximum source location for main data
      lmx = (ny - 1)/kyp2
! this segment is used for mpi computers
!
      if (ks < kyb) then
         if (l1 <= lmx) then
            call MPI_IRECV(fxy,nnxv*kyp,mreal,l1,moff+5,lgrp,msid,ierr)
! receive guard cells if on another processor
            if (ks==(kyb-1)) then
               if (l2 > l1) then
                  call MPI_IRECV(fxy(1,1,kyp1),nnxv,mreal,l2,moff+5,lgrp&
     &,nsid,ierr)
               endif
            endif
         endif
      endif
! pack data and send first half
      if (m1 < kyb) then
         do k = 2, kyp1
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do j = 1, nxv
               do i = 1, ndim
                  fxy2(i,j+joff,kk) = fxy2(i,j,k)
               enddo
            enddo
         enddo
         call MPI_SEND(fxy2,nnxv*kypl,mreal,m1,moff+5,lgrp,ierr)
! send guard cell if on another processor
      else if (kypl > 0) then
         call MPI_SEND(fxy2,nnxv*kypl,mreal,m1-1,moff+5,lgrp,ierr)
      endif
! pack data and send second half
      if (m2 < kyb) then
         do k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do j = 1, nxv
               do i = 1, ndim
                  fxy2(i,j+joff,kk+kyp) = fxy2(i,j,k+kyp)
               enddo
            enddo
         enddo
         call MPI_SEND(fxy2(1,1,kyp1),nnxv*kypl,mreal,m2,moff+5,lgrp,   &
     &ierr)
      endif
! wait for data
      if (ks < kyb) then
         if (l1 <= lmx) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
! wait for guard cells if on another processor
         if (ks==(kyb-1)) then
            if (l2 > l1) then
               call MPI_WAIT(nsid,istatus,ierr)
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd, &
     &kypd)
! this subroutine performs a transpose of a real matrix f, distributed
! in y, to a real matrix g, distributed in x, that is,
! g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! includes an extra guard cell for last row and column
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = real input array
! g = real output array
! s, t = real scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, nxv >= nx+1
! nyv = first dimension of g, nyv >= ny+1
! kypd = second dimension of f, kypd >= kyp+1
! kxpd = second dimension of g, kxpd >= kxp+1
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      real, dimension(nxv,kypd), intent(in) :: f
      real, dimension(nyv,kxpd), intent(inout) :: g
      real, dimension((kxp+1)*(kyp+1)), intent(inout) :: s, t
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: nx1, ny1, kxb, kyb
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
! ks = processor id
      ks = kstrt - 1
! kxps = actual size used in x direction
      kxps = min(kxp,max(0,nx-kxp*ks))
! kyps = actual size used in y direction
      kyps = min(kyp,max(0,ny-kyp*ks))
! kxb = minimum number of processors needed in x direction
      kxb = (nx - 1)/kxp + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb-1)) kxps = kxps + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kyps = kyps + 1
! kxyp = maximum amount of data to be received
      kxyp = (kxp + 1)*(kyp + 1)
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ny1
            do j = 1, nx1
               g(k,j) = f(j,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        kyps = min(kyp,max(0,ny-koff))
!        if (m==kyb) kyps = kyps + 1
!        do k = 1, kyps
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              kxps = min(kxp,max(0,nx-joff))
!              if (l==kxb) kxps = kxps + 1
!              do j = 1, kxps
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
! add extra word for last processor in x
         if (id==(kxb-1)) ld = ld + 1
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyps
            do j = 1, ld
               s(j+ld*(k-1)) = f(j+joff,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         ld = ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mreal,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mreal,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
! add extra word for last processor in y
         if (id==(kyb-1)) ld = ld + 1
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
      subroutine PPRNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,&
     &kxpd,kypd)
! this subroutine performs a transpose of a real matrix f, distributed
! in y, to a real matrix g, distributed in x, that is,
! g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! includes an extra guard cell for last row and column
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = real input array
! g = real output array
! s, t = real scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! ndim = leading dimension of arrays f and g
! nxv = second dimension of f, nxv >= nx+1
! nyv = second dimension of g, nyv >= ny+1
! kypd = third dimension of f, kypd >= kyp+1
! kxpd = third dimension of g, kxpd >= kxp+1
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      real, dimension(ndim,nxv,kypd), intent(in) :: f
      real, dimension(ndim,nyv,kxpd), intent(inout) :: g
      real, dimension(ndim,(kxp+1)*(kyp+1)), intent(inout) :: s, t
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: nx1, ny1, kxb, kyb
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
! ks = processor id
      ks = kstrt - 1
! kxps = actual size used in x direction
      kxps = min(kxp,max(0,nx-kxp*ks))
! kyps = actual size used in y direction
      kyps = min(kyp,max(0,ny-kyp*ks))
! kxb = minimum number of processors needed in x direction
      kxb = (nx - 1)/kxp + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb-1)) kxps = kxps + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kyps = kyps + 1
! kxyp = maximum amount of data to be received
      kxyp = ndim*(kxp + 1)*(kyp + 1)
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ny1
            do j = 1, nx1
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
!        kyps = min(kyp,max(0,ny-koff))
!        if (m==kyb) kyps = kyps + 1
!        do k = 1, kyps
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              kxps = min(kxp,max(0,nx-joff))
!              if (l==kxb) kxps = kxps + 1
!              do j = 1, kxps
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
! add extra word for last processor in x
         if (id==(kxb-1)) ld = ld + 1
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
         call MPI_IRECV(t,kxyp,mreal,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mreal,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
! add extra word for last processor in y
         if (id==(kyb-1)) ld = ld + 1
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
      subroutine mplcguard2(f,nyp,tguard,kstrt,nvp)
! copies scalar data to interior guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnlcguard2(f,nyp,tguard,kstrt,nvp)
! copies vector data to interior guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1)*size(f,2); nypmx = size(f,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnlaguard2(f,nyp,tguard,nx,kstrt,nvp)
! adds scalar data from interior guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! local data
      integer :: nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(f,1); nypmx = size(f,2)
! check if required size of buffer has increased
      if (szscr < nxv) then
         if (szscr /= 0) deallocate(scr)
! allocate new buffer
         allocate(scr(nxv))
         szscr = nxv
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNLAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpnlacguard2(f,nyp,tguard,nx,kstrt,nvp)
! adds vector data from interior guard cells in non-uniform partitions
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! local data
      integer :: ndim, nxv, nypmx
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxv = size(f,2); nypmx = size(f,3)
! check if required size of buffer has increased
      if (szscr < ndim*nxv) then
         if (szscr /= 0) deallocate(scr)
! allocate new buffer
         allocate(scr(ndim*nxv))
         szscr = ndim*nxv
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPNLACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
! record time
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdblsin2d(q,q2,tfmov,nx,ny,kstrt,nvp,kyp)
! this subroutine creates an odd array q2 from an array q
! input data must have a uniform partition
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, kyp
      real, intent(inout) :: tfmov
      real, intent(inout), dimension(:,:) :: q, q2
! local data
      integer :: nxv, kypd, kyp2d, kyb2, ksh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(q,1); kypd = size(q,2); kyp2d = size(q2,2)
! calculate size of final shift of q2s
      kyb2 = (2*ny - 1)/(2*kyp) + 1
      ksh = 2*kyp - min(2*kyp,max(0,2*ny-2*kyp*(kyb2-1)))
! check if required size of qs buffer has increased
      if (szqs < nxv) then
         if (szqs > 0) deallocate(qs)
! allocate new buffer
         allocate(qs(nxv))
         szqs = nxv
      endif
! check if required size of q2s buffer has increased
      if (szq2s < 2*nxv*ksh) then
         if (szq2s >= 0) deallocate(q2s)
! allocate new buffer
         allocate(q2s(2*nxv,ksh))
         szq2s = 2*nxv*ksh
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPDBLSIN2D(q,q2,qs,q2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d,ksh&
     &)
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdblsin2c(cu,cu2,tfmov,nx,ny,kstrt,nvp,kyp)
! this subroutine creates an mixed even/odd array cu2 from an array cu
! input data must have a uniform partition
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, kyp
      real, intent(inout) :: tfmov
      real, intent(inout), dimension(:,:,:) :: cu, cu2
! local data
      integer :: ndim, nxv, kypd, kyp2d, kyb2, ksh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(cu,1); nxv = size(cu,2); kypd = size(cu,3)
      kyp2d = size(cu2,3)
! calculate size of final shift of cu2s
      kyb2 = (2*ny - 1)/(2*kyp) + 1
      ksh = 2*kyp - min(2*kyp,max(0,2*ny-2*kyp*(kyb2-1)))
! check if required size of cus buffer has increased
      if (szcus < ndim*nxv) then
         if (szcus > 0) deallocate(cus)
! allocate new buffer
         allocate(cus(ndim,nxv))
         szcus = ndim*nxv
      endif
! check if required size of cu2s buffer has increased
      if (szcu2s < 2*ndim*nxv*ksh) then
         if (szcu2s >= 0) deallocate(cu2s)
! allocate new buffer
         allocate(cu2s(ndim,2*nxv,ksh))
         szcu2s = 2*ndim*nxv*ksh
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (ndim==2) then
         call PPDBLSIN2C(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,  &
     &kyp2d,ksh)
      else if (ndim==3) then
         call PPDBLSIN2B(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,  &
     &kyp2d,ksh)
      else
         write (*,*) 'mpdblsin2c error: ndim=',ndim
      endif
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpdblsin2m(amu,amu2,tfmov,nx,ny,kstrt,nvp,kyp)
! this subroutine creates an mixed even/odd tensor array amu2 from an
! array amu
! input data must have a uniform partition
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, kyp
      real, intent(inout) :: tfmov
      real, intent(inout), dimension(:,:,:) :: amu, amu2
! local data
      integer :: mdim, nxv, kypd, kyp2d, kyb2, ksh
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      mdim = size(amu,1); nxv = size(amu,2); kypd = size(amu,3)
      kyp2d = size(amu2,3)
! calculate size of final shift of cu2s
      kyb2 = (2*ny - 1)/(2*kyp) + 1
      ksh = 2*kyp - min(2*kyp,max(0,2*ny-2*kyp*(kyb2-1)))
! check if required size of cus buffer has increased
      if (szamus < mdim*nxv) then
         if (szamus > 0) deallocate(amus)
! allocate new buffer
         allocate(amus(mdim,nxv))
         szamus = mdim*nxv
      endif
! check if required size of amu2s buffer has increased
      if (szamu2s < 2*mdim*nxv*ksh) then
         if (szamu2s >= 0) deallocate(amu2s)
! allocate new buffer
         allocate(amu2s(mdim,2*nxv,ksh))
         szamu2s = 2*mdim*nxv*ksh
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      if (mdim==2) then
         call PPDBLSIN22M(amu,amu2,amus,amu2s,nx,ny,kstrt,nvp,nxv,kyp,  &
     &kypd,kyp2d,ksh)
      else if (mdim==4) then
         call PPDBLSIN2M(amu,amu2,amus,amu2s,nx,ny,kstrt,nvp,nxv,kyp,   &
     &kypd,kyp2d,ksh)
      else
         write (*,*) 'mpdblsin2m error: mdim=',mdim
      endif
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mphafdbl2d(q,q2,tfmov,nx,ny,kstrt,nvp,kyp)
! this subroutine copies data from a double scalar array q2 to regular
! scalar array q
! input and output data must have a uniform partition
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, kyp
      real, intent(inout) :: tfmov
      real, intent(inout), dimension(:,:) :: q, q2
! local data
      integer :: nxv, kypd, kyp2d
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxv = size(q,1); kypd = size(q,2); kyp2d = size(q2,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPHAFDBL2D(q,q2,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d)
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mphafdbl2c(fxy,fxy2,tfmov,nx,ny,kstrt,nvp,kyp)
! this subroutine copies data from a double vector array fxy2 to regular
! vector array fxy
! input and output data must have a uniform partition
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, kyp
      real, intent(inout) :: tfmov
      real, intent(inout), dimension(:,:,:) :: fxy, fxy2
! local data
      integer :: ndim, nxv, kypd, kyp2d
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(fxy,1); nxv = size(fxy,2); kypd = size(fxy,3)
      kyp2d = size(fxy2,3)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPHAFDBL2C(fxy,fxy2,nx,ny,kstrt,nvp,ndim,nxv,kyp,kypd,kyp2d)
! record time
      call dtimer(dtime,itime,1)
      tfmov = tfmov + real(dtime)
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
      subroutine PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
      use mpdplib2, only: SUB => PPNLCGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      call SUB(f,nyp,kstrt,nvp,nxv,nypmx)
      end subroutine
!
      subroutine PPNLAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      use mpdplib2, only: SUB => PPNLAGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      end subroutine
!
      subroutine PPNLACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      use mpdplib2, only: SUB => PPNLACGUARD2L
      implicit none
      integer, intent(in) ::  nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2C(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd&
     &,kyp2d,kshd)
      use mpdplib2, only: SUB => PPDBLSIN2C
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(2,nxv,kypd), intent(inout) :: cu
      real, dimension(2,2*nxv,kyp2d), intent(inout) :: cu2
      real, dimension(2,nxv), intent(inout) :: cus
      real, dimension(2,2*nxv,kshd), intent(inout) :: cu2s
      call SUB(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d,kshd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2D(q,q2,qs,q2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,   &
     &kyp2d,kshd)
      use mpdplib2, only: SUB => PPDBLSIN2D
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(nxv,kypd), intent(inout) :: q
      real, dimension(2*nxv,kyp2d), intent(inout) :: q2
      real, dimension(nxv), intent(inout) :: qs
      real, dimension(2*nxv,kshd), intent(inout) :: q2s
      call SUB(q,q2,qs,q2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d,kshd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDBLSIN2B(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd&
     &,kyp2d,kshd)
      use mpdplib2, only: SUB => PPDBLSIN2B
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      integer, intent(in) :: kshd
      real, dimension(3,nxv,kypd), intent(inout) :: cu
      real, dimension(3,2*nxv,kyp2d), intent(inout) :: cu2
      real, dimension(3,nxv), intent(inout) :: cus
      real, dimension(3,2*nxv,kshd), intent(inout) :: cu2s
      call SUB(cu,cu2,cus,cu2s,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d,kshd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPHAFDBL2D(q,q2,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d)
      use mpdplib2, only: SUB => PPHAFDBL2D
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, nxv, kyp, kypd, kyp2d
      real, intent(inout), dimension(nxv,kypd) :: q
      real, intent(inout), dimension(2*nxv,kyp2d) :: q2
      call SUB(q,q2,nx,ny,kstrt,nvp,nxv,kyp,kypd,kyp2d)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPHAFDBL2C(fxy,fxy2,nx,ny,kstrt,nvp,ndim,nxv,kyp,kypd, &
     &kyp2d)
      use mpdplib2, only: SUB => PPHAFDBL2C
      implicit none
      integer, intent(in) :: nx, ny, kstrt, nvp, ndim, nxv, kyp
      integer, intent(in) :: kypd, kyp2d
      real, intent(inout), dimension(ndim,nxv,kypd) :: fxy
      real, intent(inout), dimension(ndim,2*nxv,kyp2d) :: fxy2
      call SUB(fxy,fxy2,nx,ny,kstrt,nvp,ndim,nxv,kyp,kypd,kyp2d)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd, &
     &kypd)
      use mpdplib2, only: SUB => PPRTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      real, dimension(nxv,kypd), intent(in) :: f
      real, dimension(nyv,kxpd), intent(inout) :: g
      real, dimension((kxp+1)*(kyp+1)), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,&
     &kxpd,kypd)
      use mpdplib2, only: SUB => PPRNTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      real, dimension(ndim,nxv,kypd), intent(in) :: f
      real, dimension(ndim,nyv,kxpd), intent(inout) :: g
      real, dimension(ndim,kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,kxpd,kypd)
      end subroutine


