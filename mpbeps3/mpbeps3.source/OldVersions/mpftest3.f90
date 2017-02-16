!-----------------------------------------------------------------------
! 3D Electrostatic MPI/OpenMP PIC code
! Field Partition Unit Test
! written by Viktor K. Decyk, UCLA
      program mpbeps3
      use mpplib3       ! use with mpplib3.f90
      use omplib
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
!     integer, parameter :: indx =   7, indy =   7, indz =   7
      integer, parameter :: indx =   4, indy =   4, indz =   4
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! declare scalars for standard code
      integer :: nx, ny, nz, nxe
      integer :: isign, ierr
!
! declare scalars for MPI code
      integer :: nvpy, nvpz, nvp, idproc, kstrt, kyp, kzp
      integer :: nypmx, nzpmx
      integer :: mterf, nterf
!
! declare scalars for OpenMP code
      integer :: irc
      integer :: nvpp
!
! declare arrays for standard code:
      double precision, dimension(4) :: wtot, work
!
! declare arrays for MPI code:
! edges(1:2) = lower:upper y boundaries of particle partition
! edges(3:4) = back:front z boundaries of particle partition
      real, dimension(:), pointer  :: edges
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z
      integer, dimension(:), pointer :: nyzp, noff
      integer, dimension(:), pointer :: myzp, moff
!
      integer :: i, j, k, l, js, ks, kyps, kzps, kk, ll, koff, loff
      integer :: kl, kr
      real :: eps, epsmax
      real, dimension(:,:), pointer :: g, h
      integer, dimension(:), pointer :: nyzps, noffs, nyzpd, noffd
      real, dimension(:,:,:), pointer :: fr, gr, fg
      real, dimension(:,:,:,:), pointer :: fvr, gvr, fgv
!
      irc = 0
! nvpp = number of shared memory nodes (0=default)
      nvpp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvpp
! initialize for shared memory parallel processing
      call INIT_OMP(nvpp)
!
! initialize scalars for standard code
! nx/ny/nz = number of grid points in x/y/z direction
      nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxe = nx + 2
!      
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
! obtain 2D partition (nvpy,nvpz) from nvp:
! nvpy/nvpz = number of processors in y/z
      call FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'FCOMP32 error: nvp,nvpy,nvpz=', nvp, nvpy, nvpz
         endif
         go to 3000
      endif
!
! initialize data for MPI code
      allocate(edges(idps),nyzp(idds),noff(idds))
      allocate(myzp(idds),moff(idds))
!
!     nypmx = 15;  nzpmx = 13
      nypmx = 30;  nzpmx = 26
!
!
! initialize additional scalars for MPI code
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvpy + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1
! mterf/nterf = number of shifts in y/z required by field manager
!               (0=initiate search)
      mterf = 0; nterf = 0
!
      allocate(g(ndim*nxe,nypmx*nzpmx),h(ndim*nxe,nypmx*nzpmx))
      allocate(nyzps(idds),noffs(idds),nyzpd(idds),noffd(idds))
      allocate(fr(nxe,nypmx,nzpmx),fvr(ndim,nxe,nypmx,nzpmx))
      allocate(gr(nxe,nypmx,nzpmx),gvr(ndim,nxe,nypmx,nzpmx))
      allocate(fg(nxe,ny+1,nz+1),fgv(ndim,nxe,ny+1,nz+1))
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kl = (ny - 1)/kyp
      kr = (nz - 1)/kzp
! set global data
      do l = 1, nz+1
      do k = 1, ny+1
      do j = 1, nxe
      fg(j,k,l) = j + nxe*(k - 1) + nxe*(ny+1)*(l - 1)
      if (kstrt==1) write (30,*) j,k,l,fg(j,k,l)
      enddo
      enddo
      enddo
! set local uniform data
      kk = kyps
      ll = kzps
      if (js==kl) kk = kk + 1
      if (ks==kr) ll = ll + 1
      do l = 1, ll
      loff = kzp*ks
      do k = 1, kk
      koff = kyp*js
      do j = 1, nxe
      fr(j,k,l) = fg(j,k+koff,l+loff)
      enddo
      enddo
      enddo
      gr = fr
!
! set non-uniform partition
      noff(1) = js; noff(2) = ks
      if (js.lt.(nvpy-1)) then
         nyzp(1) = 1
      else
         nyzp(1) = ny - nvpy + 1
      endif
      if (ks.lt.(nvpz-1)) then
         nyzp(2) = 1
      else
         nyzp(2) = nz - nvpz + 1
      endif
!
! set non-uniform partition
      if (js.gt.0) then
         myzp(1) = 1
         moff(1) = ny - nvpy + js
      else
         moff(1) = 0
         myzp(1) = ny - nvpy + 1
      endif
      if (ks.gt.0) then
         moff(2) = nz - nvpz + ks
         myzp(2) = 1
      else
         moff(2) = 0
         myzp(2) = nz - nvpz + 1
      endif
!
      isign = 1
      call PPFYMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,mterf,ierr)
      if (kstrt==1) write(*,*) 'mterf,ierr=',mterf,ierr
!
! check local non-uniform in y data
      kk = nyzp(1); if (js==(nvpy-1)) kk = kk + 1
      ll = kzps; if (ks==kr) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = kzp*ks
      do k = 1, kk
      koff = noff(1)
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (50+kstrt,*) j,k,l,gr(j,k,l),koff,loff,fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'y local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) write (*,*) 'y global epsmax=',epsmax
!
      isign = 1
      call PPFZMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,nterf,ierr)
      if (kstrt==1) write(*,*) 'nterf,ierr=',nterf,ierr
!
! check local non-uniform data
      kk = nyzp(1); if (js==(nvpy-1)) kk = kk + 1
      ll = nyzp(2); if (ks==(nvpz-1)) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = noff(2)
      do k = 1, kk
      koff = noff(1)
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (70+kstrt,*) j,k,l,gr(j,k,l),fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'non-uniform local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) write (*,*) 'non-uniform global epsmax=',epsmax
!
      isign = -1
      call PPFZMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,nterf,ierr)
      if (kstrt==1) write(*,*) 'nterf,ierr=',nterf,ierr
!
! check local non-uniform in y data
      kk = nyzp(1); if (js==(nvpy-1)) kk = kk + 1
      ll = kzps; if (ks==kr) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = kzp*ks
      do k = 1, kk
      koff = noff(1)
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (50+kstrt,*) j,k,l,gr(j,k,l),koff,loff,fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'y local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) write (*,*) 'y global epsmax=',epsmax
!
      isign = -1
      call PPFYMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,mterf,ierr)
      if (kstrt==1) write(*,*) 'mterf,ierr=',mterf,ierr
!
! check local uniform data
      kk = kyps
      ll = kzps
      if (js==kl) kk = kk + 1
      if (ks==kr) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = kzp*ks
      do k = 1, kk
      koff = kyp*js
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (50+kstrt,*) j,k,l,gr(j,k,l),koff,loff,fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'uniform local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) write (*,*) 'uniform global epsmax=',epsmax
!
! move data to first non-local fields
      isign = 1
      call PPFYMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,mterf,ierr)
      call PPFZMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,nterf,ierr)
!
! move data to second non-local fields
      isign = 0
      noffs = noff; noffd = moff
      nyzps = nyzp; nyzpd = myzp
      nterf = 0; mterf = 0.
      call PPFYMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,mterf,ierr)
      if (kstrt==1) write(*,*) 'mterf,ierr=',mterf,ierr
!
! check local non-uniform data
      kk = myzp(1); if (js==(nvpy-1)) kk = kk + 1
      ll = nyzp(2); if (ks==(nvpz-1)) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = noff(2)
      do k = 1, kk
      koff = moff(1)
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (70+kstrt,*) j,k,l,gr(j,k,l),fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'special non-uniform in y local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) then
         write (*,*) 'special non-uniform in y global epsmax=',epsmax
      endif
!
      isign = 0
      noffs(1) = moff(1); noffs(2) = noff(2)
      nyzps(1) = myzp(1); nyzps(2) = nyzp(2)
      noffd = moff; nyzpd = myzp
      call PPFZMOVE32(gr,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,isign,kyp&
     &,kzp,ny,nz,kstrt,nvpy,nvpz,nxe,nypmx,nzpmx,idds,nterf,ierr)
!
! check local non-uniform data
      kk = myzp(1); if (js==(nvpy-1)) kk = kk + 1
      ll = myzp(2); if (ks==(nvpz-1)) ll = ll + 1
      epsmax = 0.0
      do l = 1, ll
      loff = moff(2)
      do k = 1, kk
      koff = moff(1)
      do j = 1, nxe
      eps = abs(gr(j,k,l)-fg(j,k+koff,l+loff))
      if (eps > 0.0) then
         epsmax = eps
         write (70+kstrt,*) j,k,l,gr(j,k,l),fg(j,k+koff,l+loff)
      endif
      enddo
      enddo
      enddo  
!     write (*,*) kstrt,'special non-uniform local epsmax=',epsmax
      wtot(1) = epsmax
      call PPDSUM(wtot,work,1)
      epsmax = wtot(1)
      if (kstrt==1) write (*,*) 'special non-uniform global epsmax=',epsmax
!
 3000 continue
      call PPEXIT()
      end program
