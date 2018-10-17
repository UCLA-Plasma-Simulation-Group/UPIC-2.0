!-----------------------------------------------------------------------
! Fortran Library for initialization of domains and particles
! 2D MPI/OpenMP PIC Codes:
! NEXTRAN2 = skips over nextrand groups of random numbers
! PDICOMP2L determines uniform integer spatial decomposition for
!           uniform distribution of particles for 2d code
!           integer boundaries set, of equal size except for remainders
! PDNICOMP2L determines uniform integer spatial decomposition for
!            uniform distribution of particles for 2d code
!            integer boundaries set, but might be of unequal size
! PDCOMP2L determines uniform real spatial boundaries for uniform
!          distribution of particles for 2d code
!          real boundaries set, of equal size
! PUDISTR2 calculates initial particle co-ordinates with uniform density
!          for 2d or 2-1/2d code
! PLDISTR2 calculates initial particle co-ordinates with linear density
!          profile for 2d or 2-1/2d code
! PFDISTR2 calculates initial particle co-ordinates with general
!          distribution in space for 2d or 2-1/2d code
! PGFDISTR2 initializes x and co-ordinate for 2d or 2-1/2d code, with
!           general distribution in space with limited x and y
!           coordinate range
! PVDISTR2 calculates initial particle velocities with maxwellian
!          velocity with drift for 2d code
! PVDISTR2H calculates initial particle velocities with maxwellian
!           velocity with drift for 2-1/2d code
! PVRDISTR2 initial particle momenta with maxwell-juttner distribution
!           with drift for 2d code
! PVRDISTR2H calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 2-1/2d code
! PVBDISTR2H calculates initial particle velocities for a magnetized
!            plasma with maxwellian velocity with drift in direction
!            parallel to B, ring distribution in directions
!            perpendicular to B for 2-1/2d code
! PPDBLKP2L finds the maximum number of particles in each tile
! PFEDGES2 = finds new partitions boundaries (edges,noff,nyp)
!            from analytic general density profile.
! PGFEDGES2 = finds new partitions boundaries (edges,noff,nyp)
!             from analytic general density profile with limited y
!             coordinate range
! PFHOLES2 determines list of particles which are leaving this node
! ranorm gaussian random number generator
! randum uniform random number generator
! rd1rand 1d rayleigh random number generator
! The following functions are used by FDISTR1 and GFDISTR1:
! DLDISTR1 calculates either a density function or its integral
!          for a linear density profile with uniform background
! DSDISTR1 calculates either a density function or its integral
!          for a sinusoidal density profile with uniform background
! DGDISTR1 calculates either a density function or its integral
!          for a gaussian density profile with uniform background
! DHDISTR1 calculates either a density function or its integral
!          for a hyperbolic secant squared density profile
!          with uniform background
! DEDISTR1 calculates either a density function or its integral
!          for an exponential density profile with uniform background
! DGDISTR0 calculates either a density function or its integral
!          for a gaussian density profile with no background density
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 24, 2018
!-----------------------------------------------------------------------
      subroutine NEXTRAN2(nextrand,ndim,np)
! for 2d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
! nextrand = (0,N) = generate (default,Nth block) of random numbers
! ndim = number of velocity dimensions = 2 or 3
! np = number of particles in distribution
      implicit none
      integer nextrand, ndim, np
! local data
      integer j, n
      double precision d
      double precision ranorm
      n = ndim*np*nextrand
      do 10 j = 1, n
      d = ranorm()
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! integer boundaries are set, of equal size except for remainders
! nvp must be < ny.  some combinations of ny and nvp result in a zero
! value of nyp.  this is not supported.
! input: ny, kstrt, nvp, idps, output: edges, nyp, noff, nypmx, nypmn
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
! ny = system length in y direction
! kstrt = starting data block number (processor id + 1)
! nvp = number of real or virtual processors
! idps = number of partition boundaries
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
      real edges
      dimension edges(idps)
! local data
      integer kb, kyp
      real at1, any
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
      any = real(ny)
! determine decomposition
      kb = kstrt - 1
      kyp = (ny - 1)/nvp + 1
      at1 = real(kyp)
      edges(1) = at1*real(kb)
      if (edges(1).gt.any) edges(1) = any
      noff = edges(1)
      edges(2) = at1*real(kb + 1)
      if (edges(2).gt.any) edges(2) = any
      kb = edges(2)
      nyp = kb - noff
! find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps&
     &)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! integer boundaries are set, but might be of unequal size
! input: ny, kstrt, nvp, idps, output: edges, nyp, noff, nypmx, nypmn
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
! ny = system length in y direction
! kstrt = starting data block number (processor id + 1)
! nvp = number of real or virtual processors
! idps = number of partition boundaries
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
      real edges
      dimension edges(idps)
! local data
      integer kb
      real at1, at2
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
! determine decomposition
      kb = kstrt - 1
      at1 = real(ny)/real(nvp)
      at2 = at1*real(kb)
      noff = at2
      edges(1) = real(noff)
      at2 = at1*real(kb + 1)
      if (kstrt.eq.nvp) at2 = real(ny)
      kb = at2
      edges(2) = real(kb)
      nyp = kb - noff
! find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PDCOMP2L(edges,nyp,myp,lyp,noff,nypmx,ny,kstrt,nvp,idps&
     &)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! real boundaries set, of equal size
! input: ny, kstrt, nvp, idps, output: edges, nyp, noff, myp, lyp, nypmx
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! myp = number of full or partial grids in particle partition
! lyp = number of guard cells for processor on left
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! ny = system length in y direction
! kstrt = starting data block number (processor id + 1)
! nvp = number of real or virtual processors
! idps = number of partition boundaries
      implicit none
      integer nyp, noff, myp, lyp, nypmx, ny, kstrt, nvp, idps
      real edges
      dimension edges(idps)
! local data
      integer kb
      real at1, dt1
      integer mypm, iwork1
      dimension mypm(1), iwork1(1)
! determine decomposition
      kb = kstrt - 1
      at1 = real(ny)/real(nvp)
      edges(1) = at1*real(kb)
      noff = edges(1)
      dt1 = edges(1) - real(noff)
      if (dt1.eq.0.0) edges(1) = real(noff)
      edges(2) = at1*real(kb + 1)
      if (kstrt.eq.nvp) edges(2) = real(ny)
      kb = edges(2)
      dt1 = edges(2) - real(kb)
      if (dt1.eq.0.0) edges(2) = real(kb)
      nyp = kb - noff
      myp = nyp
      if (dt1.gt.0.0) myp = myp + 1
! find number of guard cells on the left
      mypm(1) = myp - nyp + 1
      call PPISHFTR(mypm,iwork1,1)
      lyp = mypm(1)
! find maximum partition size
      mypm(1) = myp
      call PPIMAX(mypm,iwork1,1)
      nypmx = mypm(1) + 1
      return
      end
!-----------------------------------------------------------------------
      subroutine PDISTR2(part,edges,npp,vtx,vty,vdx,vdy,npx,npy,nx,ny,  &
     &idimp,npmax,idps,ipbc,ierr)
! for 2d code, this subroutine calculates initial particle co-ordinates
! and velocities with uniform density and maxwellian velocity with drift
! for distributed data.
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, npp, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! npp = number of particles in partition, updated in this procedure
! vtx/vty = thermal velocity of particles in x/y direction
! vdx/vdy = drift velocity of particles in x/y direction
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, idimp, npmax, idps, ipbc, ierr
      real vtx, vty, vdx, vdy
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, nps, npt, npxyp
      real edgelx, edgely, at1, at2, xt, yt, vxt, vyt
      double precision dnpx, dnpxy, dt1
      integer ierr1, iwork1
      double precision sum3, work3
      dimension ierr1(1), iwork1(1), sum3(3), work3(3)
      double precision ranorm
      ierr = 0
      nps = npp + 1
! particle distribution constant
      dnpx = dble(npx)
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      xt = edgelx + at1*(real(j) - 0.5)
! maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         npt = npp + 1
         if (npt.le.npmax) then
            part(1,npt) = xt
            part(2,npt) = yt
            part(3,npt) = vxt
            part(4,npt) = vyt
            npp = npt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
   20 continue
      npxyp = 0
! add correct drift
      sum3(1) = 0.0d0
      sum3(2) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum3(1) = sum3(1) + part(3,j)
      sum3(2) = sum3(2) + part(4,j)
   30 continue
      sum3(3) = npxyp
      call PPDSUM(sum3,work3,3)
      dnpxy = sum3(3)
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      dt1 = 1.0d0/dnpxy
      sum3(1) = dt1*sum3(1) - vdx
      sum3(2) = dt1*sum3(2) - vdy
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum3(1)
      part(4,j) = part(4,j) - sum3(2)
   40 continue
! process errors
      dnpxy = dnpxy - dnpx*dble(npy)
      if (dnpxy.ne.0.0d0) ierr = dnpxy
      return
      end
!-----------------------------------------------------------------------
      subroutine PDISTR2H(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,nx,ny,idimp,npmax,idps,ipbc,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! co-ordinates and velocities with uniform density and maxwellian
! velocity with drift for distributed data.
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! npp = number of particles in partition, updated in this procedure
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift velocity of particles in x/y/z direction
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! idimp = size of phase space = 5
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, idimp, npmax, idps, ipbc, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, nps, npt, npxyp
      real edgelx, edgely, at1, at2, xt, yt, vxt, vyt, vzt
      double precision dnpx, dnpxy, dt1
      integer ierr1, iwork1
      double precision sum4, work4
      dimension ierr1(1), iwork1(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
      nps = npp + 1
! particle distribution constant
      dnpx = dble(npx)
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      xt = edgelx + at1*(real(j) - 0.5)
! maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         npt = npp + 1
         if (npt.le.npmax) then
            part(1,npt) = xt
            part(2,npt) = yt
            part(3,npt) = vxt
            part(4,npt) = vyt
            part(5,npt) = vzt
            npp = npt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
   20 continue
      npxyp = 0
! add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum4(1) = sum4(1) + part(3,j)
      sum4(2) = sum4(2) + part(4,j)
      sum4(3) = sum4(3) + part(5,j)
   30 continue
      sum4(4) = npxyp
      call PPDSUM(sum4,work4,4)
      dnpxy = sum4(4)
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      dt1 = 1.0d0/dnpxy
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum4(1)
      part(4,j) = part(4,j) - sum4(2)
      part(5,j) = part(5,j) - sum4(3)
   40 continue
! process errors
      dnpxy = dnpxy - dnpx*dble(npy)
      if (dnpxy.ne.0.0d0) ierr = dnpxy
      return
      end
!-----------------------------------------------------------------------
      subroutine PUDISTR2(part,edges,npp,npx,npy,nx,ny,idimp,npmax,idps,&
     &ipbc,ierr)
! for 2d code, this subroutine calculates initial particle co-ordinates
! with uniform density for distributed data.
! input: all except part, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! npp = number of particles in partition, updated in this procedure
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ierr = (0,1) = (no,yes) error condition exists
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, idimp, npmax, idps, ipbc, ierr
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, nps, npt
      real edgelx, edgely, at1, at2, xt, yt
      double precision dnpxy
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      ierr = 0
      nps = npp
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         do 10 j = 1, npx
         xt = edgelx + at1*(real(j) - 0.5)
         npt = npp + 1
         if (npt.le.npmax) then
            part(1,npt) = xt
            part(2,npt) = yt
            npp = npt
         else
            ierr = 1
         endif
   10    continue
      endif
   20 continue
! process errors
! first check if buffer overflow occurred
      if (ierr.eq.0) then
         ierr1(1) = 0  
      else
         ierr1(1) = npt
      endif
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.0) return
! then check if not all particles were distributed
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxy = sum1(1) - dble(npx)*dble(npy)
      if (dnpxy.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PLDISTR2(part,npp,anlx,anly,npx,npy,nx,ny,kstrt,nvp,   &
     &idimp,npmax,ipbc,ierr)
! for 2d code, this subroutine calculates initial particle co-ordinates
! with the following bi-linear density profile:
! n(x,y) = n(x)*n(y), where n(x) = n0x*(1. + anlx*(x/nx - .5)) and 
! n(y) = n0y*(1. + anly*(y/ny - .5)) and where
! n0x = npx/(nx - 2*edgelx) and n0y = npy/(ny - 2*edgely)
! for distributed data.
! the algorithm partitions the number of particles distributed in the
! y direction uniformly (with possible remainders).
! particles are not necessarily in the correct processor, and may need
! to be moved to another processor.
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! npp = number of particles in partition, updated in this procedure
! anlx/anly = initial linear density weight in x/y direction
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! idimp = size of phase space = 4
! kstrt = starting data block number
! nvp = number of real or virtual processors
! npmax = maximum number of particles in each partition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ierr = (0,1) = (no,yes) error condition exists
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, kstrt, nvp, idimp, npmax, ipbc
      integer ierr
      real anlx, anly
      real part
      dimension part(idimp,npmax)
! local data
      integer nps, npt, ks, j, k, kk, mpy, mpys, joff
      real edgelx, edgely, at1, at2, bt1, bt2, antx, anty, xt, yt
      double precision dnpxy
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      nps = npp
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! check if particle overflow will occurr
      npt = npp + npx*mpys
      ierr1(1) = npt
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.npmax) return
      ierr = 0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
      if (anlx.ne.0.0) then
         antx = anlx/real(nx)
         at1 = 2.0*antx*at1
         bt1 = 1.0 - 0.5*antx*(real(nx) - 2.0*edgelx)
      endif
      if (anly.ne.0.0) then
         anty = anly/real(ny)
         at2 = 2.0*anty*at2
         bt2 = 1.0 - 0.5*anty*(real(ny) - 2.0*edgely)
      endif
!$OMP PARALLEL DO PRIVATE(j,k,kk,joff,xt,yt)
      do 20 k = 1, mpys
      kk = k + mpy*ks
      joff = npx*(k - 1) + npp
! linear density in y
      if (anly.ne.0.0) then
         yt = edgely + (sqrt(bt2*bt2 + at2*(real(kk) - 0.5)) - bt2)/anty
! uniform density in y
      else
         yt = edgely + at2*(real(kk) - 0.5)
      endif
      do 10 j = 1, npx
! linear density in x
      if (anlx.ne.0.0) then
         xt = edgelx + (sqrt(bt1*bt1 + at1*(real(j) - 0.5)) - bt1)/antx
! uniform density in x
      else
         xt = edgelx + at1*(real(j) - 0.5)
      endif
      part(1,j+joff) = xt
      part(2,j+joff) = yt
   10 continue
   20 continue
!$OMP END PARALLEL DO
! update number of particles
      npp = npt
! check if not all particles were distributed
      dnpxy = dble(npx)*dble(npy)
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxy = sum1(1) - dnpxy
      if (dnpxy.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PFDISTR2(part,npp,fnx,argx1,argx2,argx3,fny,argy1,argy2&
     &,argy3,npx,npy,nx,ny,kstrt,nvp,idimp,npmax,ipbc,ierr)
! for 2d code, this subroutine calculates initial particle co-ordinates
! with general density profile n(x,y) = n(x)*n(y), 
! where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
! and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
! and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
! and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
! for distributed data.
! the algorithm partitions the number of particles distributed in the
! y direction uniformly (with possible remainders).
! particles are not necessarily in the correct processor, and may need
! to be moved to another processor.
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! npp = number of particles in partition, updated in this procedure
! fnx/fny = density and density integral function in x/y direction
! argx1,argx2,argx3 = arguments to fnx
! argy1,argy2,argy3 = arguments to fny
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ierr = (0,1) = (no,yes) error condition exists
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, kstrt, nvp, idimp, npmax, ipbc
      integer ierr
      double precision argx1, argx2, argx3, argy1, argy2, argy3
      real part
      dimension part(idimp,npmax)
      double precision fnx, fny
      external fnx, fny
! local data
      integer nps, npt, ks, i, j, k, kk, mpy, mpys, imax, joff, moff
      real edgelx, edgely, anx, any, bnx, bny, xt0, yt0, x0, y0
      real xn, yn, eps, big, f, fp
      double precision xt, yt, dnpxy
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      nps = npp
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! check if particle overflow will occurr
      npt = npp + npx*mpys
      ierr1(1) = npt
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.npmax) return
      ierr = 0
! eps = convergence criterion
      imax = max(nx,ny)
      eps = 0.0001
      big = 0.25
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
      else if (ipbc.eq.3) then
         edgelx = 1.0
      endif
! find normalization for function
      anx = real(nx) - edgelx
      any = real(ny) - edgely
      x0 = fnx(dble(edgelx),argx1,argx2,argx3,1)
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
      bnx = real(npx)/(fnx(dble(anx),argx1,argx2,argx3,1) - x0)
      bny = real(npy)/(fny(dble(any),argy1,argy2,argy3,1) - y0)
      x0 = bnx*x0 - 0.5
      y0 = bny*y0 - 0.5
! density profile in x
      do 20 j = 1, npx
      xn = real(j) + x0
! guess next value for xt
      if (j.eq.1) then
         xt0 = edgelx
         xt = xt0
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
         if (fp.gt.0.0) xt = xt + 0.5/fp
      else
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
         if (fp.eq.0.0) fp = 1.0
         xt = xt + 1.0/fp
      endif
      xt = max(edgelx,min(xt,anx))
      i = 0
   10 f = bnx*fnx(xt,argx1,argx2,argx3,1) - dble(xn)
! find improved value for xt
      if (abs(f).ge.eps) then
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
! newton's method
         if ((abs(f).lt.big).and.(fp.gt.0.0)) then
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
! bisection method
         else if (f.gt.0.0) then
            fp = 0.5*abs(xt0 - xt)
            xt = xt0 - fp
         else
            fp = abs(xt - xt0)
            xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
!        write (*,*) j,'newton iteration max exceeded, xt = ', xt
         ierr = ierr + 1
      endif
      part(1,j+npp) = xt
      xt0 = xt
   20 continue
! quit if error
      if (ierr.ne.0) return
! density profile in y
      moff = mpy*ks
      do 50 k = 1, mpys + moff
      kk = k - moff
      yn = real(k) + y0
! guess next value for yt
      if (k.eq.1) then
         yt0 = edgely
         yt = yt0
         fp = bny*fny(yt,argy1,argy2,argy3,0)
         if (fp.gt.0.0) yt = yt + 0.5/fp
      else
         fp = bny*fny(yt,argy1,argy2,argy3,0)
         if (fp.eq.0.0) fp = 1.0
         yt = yt + 1.0/fp
      endif
      yt = max(edgely,min(yt,any))
      i = 0
   30 f = bny*fny(yt,argy1,argy2,argy3,1) - dble(yn)
! find improved value for yt
      if (abs(f).ge.eps) then
         fp = bny*fny(yt,argy1,argy2,argy3,0)
! newton's method
         if ((abs(f).lt.big).and.(fp.gt.0.0)) then
            yt0 = yt
            yt = yt - f/fp
            yt = max(edgely,min(yt,any))
! bisection method
         else if (f.gt.0.0) then
            fp = 0.5*abs(yt0 - yt)
            yt = yt0 - fp
          else
            fp = abs(yt - yt0)
            yt0 = yt
            yt = yt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 30
!        write (*,*) k,'newton iteration max exceeded, yt = ', yt
         ierr = ierr + 1
      endif
! store co-ordinates
      if ((kk.ge.1).and.(kk.le.mpys)) then
         joff = npx*(kk-1) + npp
         do 40 j = 1, npx
         part(1,j+joff) = part(1,j+npp)
         part(2,j+joff) = yt
   40    continue
      endif
      yt0 = yt
   50 continue
! quit if error
      if (ierr.ne.0) return
! update number of particles
      npp = npt
! check if not all particles were distributed
      dnpxy = dble(npx)*dble(npy)
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxy = sum1(1) - dnpxy
      if (dnpxy.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PGFDISTR2(part,npp,fnx,argx1,argx2,argx3,fny,argy1,    &
     &argy2,argy3,xmin,xmax,ymin,ymax,npx,npy,nx,ny,kstrt,nvp,idimp,    &
     &npmax,ierr)
! for 2d code, this subroutine calculates initial particle co-ordinates
! with general density profile n(x,y) = n(x)*n(y), 
! where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
! and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
! and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
! and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
! for distributed data.
! x/y co-ordinates will be within range 0 < x < nx and 0 < y < ny
! the algorithm partitions the number of particles distributed in the
! y direction uniformly (with possible remainders).
! particles are not necessarily in the correct processor, and may need
! to be moved to another processor.
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! npp = number of particles in partition, updated in this procedure
! fnx/fny = density and density integral function in x/y direction
! argx1,argx2,argx3 = arguments to fnx
! argy1,argy2,argy3 = arguments to fny
! xmin/xmax = minimum/maximum range of particle coordinates in x
! ymin/ymax = minimum/maximum range of particle coordinates in y
! npx/npy = initial number of particles distributed in x/y direction
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! with spatial decomposition
      implicit none
      integer npp, npx, npy, nx, ny, kstrt, nvp, idimp, npmax, ierr
      double precision argx1, argx2, argx3, argy1, argy2, argy3
      real xmin, xmax, ymin, ymax
      real part
      dimension part(idimp,npmax)
      double precision fnx, fny
      external fnx, fny
! local data
      integer nps, npt, ks, i, j, k, kk, mpy, mpys, imax, joff, moff
      real edgelx, edgely, anx, any, bnx, bny, xt0, yt0, x0, y0
      real xn, yn, eps, big, f, fp
      double precision xt, yt, dnpxy
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      nps = npp
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! check if particle overflow will occurr
      npt = npp + npx*mpys
      ierr1(1) = npt
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.npmax) return
      ierr = 0
! eps = convergence criterion
      imax = max(nx,ny)
      eps = 0.0001
      big = 0.25
! set boundary value in x
      if ((xmin.lt.0.0).or.(xmin.ge.xmax).or.(xmax.gt.real(nx))) then
         ierr = -1
         return
      endif
      edgelx = xmin
! set boundary value in y
      if ((ymin.lt.0.0).or.(ymin.ge.ymax).or.(ymax.gt.real(ny))) then
         ierr = -2
         return
      endif
      edgely = ymin
! find normalization for function
      anx = xmax
      any = ymax
      x0 = fnx(dble(edgelx),argx1,argx2,argx3,1)
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
      bnx = real(npx)/(fnx(dble(anx),argx1,argx2,argx3,1) - x0)
      bny = real(npy)/(fny(dble(any),argy1,argy2,argy3,1) - y0)
      x0 = bnx*x0 - 0.5
      y0 = bny*y0 - 0.5
! density profile in x
      do 20 j = 1, npx
      xn = real(j) + x0
! guess next value for xt
      if (j.eq.1) then
         xt0 = edgelx
         xt = xt0
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
         if (fp.gt.0.0) xt = xt + 0.5/fp
      else
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
         if (fp.eq.0.0) fp = 1.0
         xt = xt + 1.0/fp
      endif
      xt = max(edgelx,min(xt,anx))
      i = 0
   10 f = bnx*fnx(xt,argx1,argx2,argx3,1) - dble(xn)
! find improved value for xt
      if (abs(f).ge.eps) then
         fp = bnx*fnx(xt,argx1,argx2,argx3,0)
! newton's method
         if ((abs(f).lt.big).and.(fp.gt.0.0)) then
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
! bisection method
         else if (f.gt.0.0) then
            fp = 0.5*abs(xt0 - xt)
            xt = xt0 - fp
         else
            fp = abs(xt - xt0)
            xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
!        write (*,*) j,'newton iteration max exceeded, xt = ', xt
         ierr = ierr + 1
      endif
      part(1,j+npp) = xt
      xt0 = xt
   20 continue
! quit if error
      if (ierr.ne.0) return
! density profile in y
      moff = mpy*ks
      do 50 k = 1, mpys + moff
      kk = k - moff
      yn = real(k) + y0
! guess next value for yt
      if (k.eq.1) then
         yt0 = edgely
         yt = yt0
         fp = bny*fny(yt,argy1,argy2,argy3,0)
         if (fp.gt.0.0) yt = yt + 0.5/fp
      else
         fp = bny*fny(yt,argy1,argy2,argy3,0)
         if (fp.eq.0.0) fp = 1.0
         yt = yt + 1.0/fp
      endif
      yt = max(edgely,min(yt,any))
      i = 0
   30 f = bny*fny(yt,argy1,argy2,argy3,1) - dble(yn)
! find improved value for yt
      if (abs(f).ge.eps) then
         fp = bny*fny(yt,argy1,argy2,argy3,0)
! newton's method
         if ((abs(f).lt.big).and.(fp.gt.0.0)) then
            yt0 = yt
            yt = yt - f/fp
            yt = max(edgely,min(yt,any))
! bisection method
         else if (f.gt.0.0) then
            fp = 0.5*abs(yt0 - yt)
            yt = yt0 - fp
          else
            fp = abs(yt - yt0)
            yt0 = yt
            yt = yt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 30
!        write (*,*) k,'newton iteration max exceeded, yt = ', yt
         ierr = ierr + 1
      endif
! store co-ordinates
      if ((kk.ge.1).and.(kk.le.mpys)) then
         joff = npx*(kk-1) + npp
         do 40 j = 1, npx
         part(1,j+joff) = part(1,j+npp)
         part(2,j+joff) = yt
   40    continue
      endif
      yt0 = yt
   50 continue
! quit if error
      if (ierr.ne.0) return
! update number of particles
      npp = npt
! check if not all particles were distributed
      dnpxy = dble(npx)*dble(npy)
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxy = sum1(1) - dnpxy
      if (dnpxy.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDISTR2(part,nps,npp,vtx,vty,vdx,vdy,npx,npy,kstrt,nvp&
     &,idimp,npmax,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! velocities with maxwellian velocity with drift for distributed data.
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty = thermal velocity of particles in x/y direction
! vdx/vdy = drift velocity of particles in x/y direction
! npx/npy = initial number of particles distributed in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer nps, npp, npx, npy, kstrt, nvp, idimp, npmax, ierr
      real vtx, vty, vdx, vdy
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyp, ks, j, k, kk, mpy, mpys
      real vxt, vyt
      double precision dnpxy, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum3, work3
      dimension sum3(3), work3(3)
      double precision ranorm
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! maxwellian velocity distribution
      npt = nps - 1
      do 20 kk = 1, npy
      k = kk - mpy*ks
      do 10 j = 1, npx
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      if ((k.ge.1).and.(k.le.mpys)) then
         npt = npt + 1
         if (npt.le.npp) then
            part(3,npt) = vxt
            part(4,npt) = vyt
         endif
      endif
   10 continue
   20 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyp = 0
      sum3(1) = 0.0d0
      sum3(2) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum3(1) = sum3(1) + part(3,j)
      sum3(2) = sum3(2) + part(4,j)
   30 continue
      sum3(3) = dble(npxyp)
      call PPDSUM(sum3,work3,3)
      dnpxy = sum3(3)
      dt1 = 1.0d0/dnpxy
      sum3(1) = dt1*sum3(1) - vdx
      sum3(2) = dt1*sum3(2) - vdy
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum3(1)
      part(4,j) = part(4,j) - sum3(2)
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &kstrt,nvp,idimp,npmax,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! velocities with maxwellian velocity with drift for distributed data.
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift velocity of particles in x/y/z direction
! npx/npy = initial number of particles distributed in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 5
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer nps, npp, npx, npy, kstrt, nvp, idimp, npmax, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyp, ks, j, k, kk, mpy, mpys
      real vxt, vyt, vzt
      double precision dnpxy, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision ranorm
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! maxwellian velocity distribution
      npt = nps - 1
      do 20 kk = 1, npy
      k = kk - mpy*ks
      do 10 j = 1, npx
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((k.ge.1).and.(k.le.mpys)) then
         npt = npt + 1
         if (npt.le.npp) then
            part(3,npt) = vxt
            part(4,npt) = vyt
            part(5,npt) = vzt
         endif
      endif
   10 continue
   20 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum4(1) = sum4(1) + part(3,j)
      sum4(2) = sum4(2) + part(4,j)
      sum4(3) = sum4(3) + part(5,j)
   30 continue
      sum4(4) = dble(npxyp)
      call PPDSUM(sum4,work4,4)
      dnpxy = sum4(4)
      dt1 = 1.0d0/dnpxy
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum4(1)
      part(4,j) = part(4,j) - sum4(2)
      part(5,j) = part(5,j) - sum4(3)
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVRDISTR2(part,nps,npp,vtx,vty,vdx,vdy,ci,npx,npy,kstrt&
     &,nvp,idimp,npmax,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! momentum with maxwell-juttner distribution with drift
! for relativistic particles and distributed data.
! f(p) = exp(-(gamma-1)*(m0c**2)/kT), where gamma = sqrt(1+(p/m0c)**2)
! since (gamma-1)*(m0c**2) = (p**2/m0)/(gamma+1), we can write
! f(p) = exp(-pt**2/2), where pt**2 = p**2/((gamma+1)/2)*m0*kT
! since pt is normally distributed, we can use a gaussian random number
! to calculate it.  We then solve the pt**2 equation to obtain:
! (p/m0)**2 = ((pt*vth)**2)*(1 + (pt/4c)**2),
! where vth = sqrt(kT/m0) is the thermal momentum/mass.  This equation
! is satisfied if we set the individual components j as follows:
! pj/m0 = ptj*vth*sqrt(1 + (pt/4c)**2)
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(3,n) = momentum px of particle n in partition
! part(4,n) = momentum py of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty = thermal velocity of particles in x/y direction
! vdx/vdy = drift momentum of particles in x/y direction
! ci = reciprocal of velocity of light
! npx/npy = initial number of particles distributed in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer nps, npp, npx, npy, kstrt, nvp, idimp, npmax, ierr
      real vtx, vty, vdx, vdy, ci
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyp, ks, j, k, kk, mpy, mpys
      real ci4, ptx, pty, pt2
      double precision dnpxy, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum3, work3
      dimension sum3(3), work3(3)
      double precision ranorm
      ci4 = 0.25*ci*ci
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! maxwell-juttner momentum distribution
      npt = nps - 1
      do 20 kk = 1, npy
      k = kk - mpy*ks
      do 10 j = 1, npx
      ptx = vtx*ranorm()
      pty = vty*ranorm()
      if ((k.ge.1).and.(k.le.mpys)) then
         pt2 = ptx*ptx + pty*pty
         pt2 = sqrt(1.0 + pt2*ci4)
         npt = npt + 1
         if (npt.le.npp) then
            part(3,npt) = ptx*pt2
            part(4,npt) = pty*pt2
         endif
      endif
   10 continue
   20 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyp = 0
      sum3(1) = 0.0d0
      sum3(2) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum3(1) = sum3(1) + part(3,j)
      sum3(2) = sum3(2) + part(4,j)
   30 continue
      sum3(3) = dble(npxyp)
      call PPDSUM(sum3,work3,3)
      dnpxy = sum3(3)
      dt1 = 1.0d0/dnpxy
      sum3(1) = dt1*sum3(1) - vdx
      sum3(2) = dt1*sum3(2) - vdy
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum3(1)
      part(4,j) = part(4,j) - sum3(2)
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVRDISTR2H(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,&
     &npy,kstrt,nvp,idimp,npmax,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! momentum with maxwell-juttner distribution with drift
! for relativistic particles and distributed data.
! f(p) = exp(-(gamma-1)*(m0c**2)/kT), where gamma = sqrt(1+(p/m0c)**2)
! since (gamma-1)*(m0c**2) = (p**2/m0)/(gamma+1), we can write
! f(p) = exp(-pt**2/2), where pt**2 = p**2/((gamma+1)/2)*m0*kT
! since pt is normally distributed, we can use a gaussian random number
! to calculate it.  We then solve the pt**2 equation to obtain:
! (p/m0)**2 = ((pt*vth)**2)*(1 + (pt/4c)**2),
! where vth = sqrt(kT/m0) is the thermal momentum/mass.  This equation
! is satisfied if we set the individual components j as follows:
! pj/m0 = ptj*vth*sqrt(1 + (pt/4c)**2)
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(3,n) = momentum px of particle n in partition
! part(4,n) = momentum py of particle n in partition
! part(5,n) = momentum pz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift momentum of particles in x/y/z direction
! ci = reciprocal of velocity of light
! npx/npy = initial number of particles distributed in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 5
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer nps, npp, npx, npy, idimp, kstrt, nvp, npmax, ierr
      real vtx, vty, vtz, vdx, vdy, vdz, ci
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyp, ks, j, k, kk, mpy, mpys
      real ci4, ptx, pty, ptz, pt2
      double precision dnpxy, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision ranorm
      ci4 = 0.25*ci*ci
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
! maxwell-juttner momentum distribution
      npt = nps - 1
      do 20 kk = 1, npy
      k = kk - mpy*ks
      do 10 j = 1, npx
      ptx = vtx*ranorm()
      pty = vty*ranorm()
      ptz = vtz*ranorm()
      if ((k.ge.1).and.(k.le.mpys)) then
         pt2 = ptx*ptx + pty*pty + ptz*ptz
         pt2 = sqrt(1.0 + pt2*ci4)
         npt = npt + 1
         if (npt.le.npp) then
            part(3,npt) = ptx*pt2
            part(4,npt) = pty*pt2
            part(5,npt) = ptz*pt2
         endif
      endif
   10 continue
   20 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum4(1) = sum4(1) + part(3,j)
      sum4(2) = sum4(2) + part(4,j)
      sum4(3) = sum4(3) + part(5,j)
   30 continue
      sum4(4) = dble(npxyp)
      call PPDSUM(sum4,work4,4)
      dnpxy = sum4(4)
      dt1 = 1.0d0/dnpxy
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum4(1)
      part(4,j) = part(4,j) - sum4(2)
      part(5,j) = part(5,j) - sum4(3)
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVBDISTR2H(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,npx&
     &,npy,kstrt,nvp,idimp,npmax,ierr)
! for 2-1/2d code, this subroutine calculates initial particle
! velocities for a magnetized plasma with maxwellian velocity with
! drift in direction parallel to B, ring distribution in directions
! perpendicular to B.  The algorithm due to phil pritchett proceeds by
! first assuming B is in the z direction, then rotating to the
! co-ordinates to match the actual B direction
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtr/vtz = thermal velocity of particles in azimuthal/B direction
! vdt/vdz = drift velocity of particles in azimuthal/B direction
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! npx/npy = initial number of particles distributed in x/y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 5
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! randum = uniform random number
! ranorm = gaussian random number with zero mean and unit variance
! with spatial decomposition
      implicit none
      integer nps, npp, npx, npy, kstrt, nvp, idimp, npmax, ierr
      real vtr, vtz, vdr, vdz, omx, omy, omz
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyp, ks, j, k, kk, mpy, mpys, ndir
      real twopi, at1, at2, vxt, vyt, vzt
      real ox, oy, oz, px, py, pz, qx, qy, qz
      double precision dnpxy, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision randum, ranorm
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvp + 1
      mpys = min(mpy,max(0,npy-mpy*ks))
      twopi = 6.28318530717959
! maxwell velocity distribution with ring distribution in x and y,
! maxwell velocity distribution in z
      npt = nps - 1
      do 20 kk = 1, npy
      k = kk - mpy*ks
      do 10 j = 1, npx
      at1 = twopi*randum()
      vxt = vdr*cos(at1) + vtr*ranorm()
      vyt = vdr*sin(at1) + vtr*ranorm()
      vzt = vtz*ranorm()
      if ((k.ge.1).and.(k.le.mpys)) then
         npt = npt + 1
         if (npt.le.npp) then
            part(3,npt) = vxt
            part(4,npt) = vyt
            part(5,npt) = vzt
         endif
      endif
   10 continue
   20 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum4(1) = sum4(1) + part(3,j)
      sum4(2) = sum4(2) + part(4,j)
      sum4(3) = sum4(3) + part(5,j)
   30 continue
      sum4(4) = dble(npxyp)
      call PPDSUM(sum4,work4,4)
      dnpxy = sum4(4)
      dt1 = 1.0d0/dnpxy
      sum4(1) = dt1*sum4(1)
      sum4(2) = dt1*sum4(2)
      sum4(3) = dt1*sum4(3) - vdz
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum4(1)
      part(4,j) = part(4,j) - sum4(2)
      part(5,j) = part(5,j) - sum4(3)
   40 continue
! rotate perpendicular and parallel component to align with actual B field
      at1 = sqrt(omx*omx + omy*omy + omz*omz)
! exit if zero B field
      if (at1.eq.0.0) return
! first create unit vector in B direction
      at1 = 1.0/at1
      ox = omx*at1
      oy = omy*at1
      oz = omz*at1
! then create unit vector in first perpendicular direction
! find direction with smallest component of B
      ndir = 1
      at1 = abs(omx)
      at2 = abs(omy)
      if (at2.le.at1) then
         ndir = 2
         at1 = at2
      endif
      if (abs(omz).lt.at1) ndir = 3
! take the cross product of that direction with B
! vpr1 = x cross B
      if (ndir.eq.1) then
         at1 = 1.0/sqrt(oy*oy + oz*oz)
         px = 0.0
         py = -oz*at1
         pz = oy*at1
! vpr1 = y cross B
      else if (ndir.eq.2) then
         at1 = 1.0/sqrt(ox*ox + oz*oz)
         px = oz*at1
         py = 0.0
         pz = -ox*at1
! vpr1 = z cross B
      else if (ndir.eq.3) then
         at1 = 1.0/sqrt(ox*ox + oy*oy)
         px = -oy*at1
         py = ox*at1
         pz = 0.0
      endif
! finally create unit vector in second perpendicular direction
! vpr2 = B cross vpr1
      qx = oy*pz - oz*py
      qy = oz*px - ox*pz
      qz = ox*py - oy*px
! rotate particle velocities
      do 50 j = nps, npp
      vxt = part(3,j)
      vyt = part(4,j)
      vzt = part(5,j)
! store components in appropriate direction
      part(3,j) = vzt*ox + vxt*px + vyt*qx
      part(4,j) = vzt*oy + vxt*py + vyt*qy
      part(5,j) = vzt*oz + vxt*pz + vyt*qz
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,  &
     &mx1,mxyp1,irc)
! this subroutine finds the maximum number of particles in each tile of
! mx, my to calculate size of segmented particle array ppart
! linear interpolation, spatial decomposition in y direction
! input: all except kpic, nppmx, output: kpic, nppmx
! part = input particle array
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! nppmx = return maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part
      dimension part(idimp,npmax)
      dimension kpic(mxyp1)
! local data
      integer j, k, n, m, mnoff, isum, ist, npx, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      if (m.le.mxyp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyp1)
      endif
   20 continue
! find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyp1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
! check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.npp) then
         irc = -1
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PFEDGES2(edges,nyp,noff,fny,argy1,argy2,argy3,nypmx,   &
     &nypmn,ny,kstrt,nvp,idps,ipbc)
! this subroutines finds new partitions boundaries (edges,noff,nyp)
! from density integral given by fny(y,argy1,argy2,argy3,1)
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! fny = density and density integral function
! argy1,argy2,argy3 = arguments to fny
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
! ny = system length in y direction
! kstrt = starting data block number (processor id + 1)
! nvp = number of real or virtual processors
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps, ipbc
      double precision argy1, argy2, argy3
      real edges
      dimension edges(idps)
      double precision fny
      external fny
! local data
      integer kb
      real edgely, any1, any, y0, y1, anpav, anpl, anpr, sum1, at1, at2
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
! particle distribution constants
      kb = kstrt - 1
! set boundary values
      edgely = 0.0
      if (ipbc.eq.2) then
         edgely = 1.0
      endif
! find normalization for function
      any = real(ny)
      any1 = any - edgely
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
! anpav = desired number of particles per processor
      anpav = (fny(dble(any1),argy1,argy2,argy3,1) - y0)/real(nvp)
! search for boundaries
      anpl = real(kb)*anpav
      anpr = real(kb+1)*anpav
      y1 = edgely
      sum1 = 0.0
! first find left boundary
   10 at1 = sum1
      sum1 = fny(dble(y1),argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpl).and.(y1.le.any)) go to 10 
      if (sum1.gt.at1) then
         at2 = (y1 - 2.0) + (anpl - at1)/(sum1 - at1)
      else
         at2 = y1 - 1.0
      endif
      edges(1) = at2
! set leftmost edge to zero
      if (kb.eq.0) edges(1) = 0.0
! then find right boundary
   20 at1 = sum1
      sum1 = fny(dble(y1),argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpr).and.(y1.le.any)) go to 20
      at2 = (y1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(2) = at2
! set rightmost edge to ny
      if ((kb+1).eq.nvp) edges(2) = any
! calculate number of grids and offsets in new partitions
      noff = edges(1) + 0.5
      kb = edges(2) + 0.5
      nyp = kb - noff
      edges(1) = real(noff)
      edges(2) = real(kb)
! find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PGFEDGES2(edges,nyp,noff,fny,argy1,argy2,argy3,ymin,   &
     &ymax,nypmx,nypmn,ny,kstrt,nvp,idps)
! this subroutines finds new partitions boundaries (edges,noff,nyp)
! from density integral given by fny(y,argy1,argy2,argy3,1)
! edges(1) = lower boundary of particle partition
! edges(2) = upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! fny = density and density integral function
! argy1,argy2,argy3 = arguments to fny
! ymin/ymax = minimum/maximum range of particle coordinates in y
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
! ny = system length in y direction
! kstrt = starting data block number (processor id + 1)
! nvp = number of real or virtual processors
! idps = number of partition boundaries
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
      double precision argy1, argy2, argy3
      real ymin, ymax
      real edges
      dimension edges(idps)
      double precision fny
      external fny
! local data
      integer kb
      real edgely, any, y0, y1, anpav, anpl, anpr, sum1, at1, at2
      real ymn, ymx
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
! particle distribution constants
      kb = kstrt - 1
! set boundary value in y
      ymn = ymin
      ymx = ymax
      if ((ymin.lt.0.0).or.(ymin.ge.ymax).or.(ymax.gt.real(ny))) then
         if (ymin.lt.0.0) ymn = 0.0
         if (ymax.gt.real(ny)) ymx = real(ny)
         if (kstrt.eq.1) then
            write (*,*) 'PGFEDGES2: invalid ymin,ymax=', ymin, ymax
         endif
      endif
      edgely = ymn
! find normalization for function
      any = ymx
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
! anpav = desired number of particles per processor
      anpav = (fny(dble(any),argy1,argy2,argy3,1) - y0)/real(nvp)
! search for boundaries
      anpl = real(kb)*anpav
      anpr = real(kb+1)*anpav
      y1 = edgely
      sum1 = 0.0
! first find left boundary
   10 at1 = sum1
      sum1 = fny(dble(y1),argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpl).and.(y1.le.any)) go to 10 
      if (sum1.gt.at1) then
         at2 = (y1 - 2.0) + (anpl - at1)/(sum1 - at1)
      else
         at2 = y1 - 1.0
      endif
      edges(1) = at2
! set leftmost edge to zero
      if (kb.eq.0) edges(1) = 0.0
! then find right boundary
   20 at1 = sum1
      sum1 = fny(dble(y1),argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpr).and.(y1.le.any)) go to 20
      at2 = (y1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(2) = at2
! set rightmost edge to ny
      if ((kb+1).eq.nvp) edges(2) = real(ny)
! calculate number of grids and offsets in new partitions
      noff = edges(1) + 0.5
      kb = edges(2) + 0.5
      nyp = kb - noff
      edges(1) = real(noff)
      edges(2) = real(kb)
! find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PFHOLES2(part,edges,npp,ihole,ndim,nc,idimp,npmax,idps,&
     &ntmax)
! for 2d code, this subroutine determines list of particles which are
! leaving this processor
! input: all except ihole, output: ihole
! part(2,n) = position y of particle n in partition
! edges(1:2) = lower:upper boundary of particle partition
! npp = number of particles in partition
! ihole = location of hole left in particle arrays
! ihole(1) = ih, number of holes left (error, if negative)
! ndim = number of velocity dimensions = 2 or 3
! nc = (1,2) = (normal,tagged) partitioned co-ordinate to sort by
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ntmax = size of hole array for particles leaving processors
      implicit none
      integer npp, ndim, nc, idimp, npmax, idps, ntmax
      real part, edges
      integer ihole
      dimension part(idimp,npmax)
      dimension edges(idps), ihole(ntmax+1)
! local data
      integer j, iy, ih, nh
      real dy
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      if ((nc.lt.1).or.(nc.gt.2)) return
! determine co-ordinate to sort by
      if (nc.eq.1) then
         iy = 2
      else if (idimp.gt.(ndim+2)) then
         iy = ndim + 3
      else
         return
      endif
      ih = 0
      nh = 0
      do 10 j = 1, npp
      dy = part(iy,j)
! find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
   10 continue
! set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
! set global error flag if needed
      ierr1(1) = -ih
      call PPIMAX(ierr1,iwork1,1)
      nh = -ierr1(1)
      if (nh.lt.0) ihole(1) = nh
      return
      end
!-----------------------------------------------------------------------
      function ranorm()
! this program calculates a random number y from a gaussian distribution
! with zero mean and unit variance, according to the method of
! mueller and box:
!    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
!    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
! where x is a random number uniformly distributed on (0,1).
! written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
!-----------------------------------------------------------------------
      function randum()
! this is a version of the random number generator dprandom due to
! c. bingham and the yale computer center, producing numbers
! in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
!-----------------------------------------------------------------------
      function rd1rand()
! random number generator for 1d rayleigh distribution:
! f(x) = x*exp(-x*x/2), using inverse transformation method
      implicit none
      double precision rd1rand, randum
      rd1rand = dsqrt(-2.0d0*dlog(randum()))
      return
      end
!-----------------------------------------------------------------------
      function DLDISTR1(x,anlx,anxi,shift,intg)
! this function calculates either a density function or its integral
! for a linear density profile.  Used in initializing particle
! coordinates.  The three parameters are redundant, and one can set one
! of them arbitrarily.  A convenient choice is to set  anxi = 1/Lx,
! anlx = NH - NL, shift = (1 - NL)/(NH - NL), where NL is the density
! at the left, and NH at the right compared to the average density
! if intg = 0, n(x) = 1. + anlx*(x*anxi - shift)
! if intg = 1, n(x) = x + .5*anlx*x*(x*anxi - 2.*shift)
      implicit none
      integer intg
      double precision x, anlx, anxi, shift
! local data
      double precision DLDISTR1, f
      if (intg.eq.0) then
         f = 1.0d0 + anlx*(x*anxi - shift)
      else if (intg.eq.1) then
         if (anxi.eq.0.0d0) then
            f = x
         else
            f = x + 0.5d0*anlx*x*(x*anxi - 2.0d0*shift)
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DLDISTR1 Error: f = ', f
      DLDISTR1 = f
      return
      end
!-----------------------------------------------------------------------
      function DSDISTR1(x,ans,dkx,phase,intg)
! this function calculates either a density function or its integral
! for a sinusoidal density profile.  Used in initializing particle
! coordinates.
! if intg = 0, n(x) = 1.0 + ans*sin(dkx*x - phase)
! if intg = 1, n(x) = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
      implicit none
      integer intg
      double precision x, ans, dkx, phase
! local data
      double precision DSDISTR1, f
      if (intg.eq.0) then
         f = 1.0d0 + ans*sin(dkx*x - phase)
      else if (intg.eq.1) then
         if (dkx.eq.0.0d0) then
            f = x - ans*sin(phase)*x
         else
            f = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DSDISTR1 Error: f = ', f
      DSDISTR1 = f
      return
      end
!-----------------------------------------------------------------------
      function DGDISTR1(x,ang,wi,x0,intg)
! this function calculates either a density function or its integral
! for a gaussian density profile.  Used in initializing particle
! coordinates.
! if intg = 0, n(x) = 1.0 + ang*exp(-((x-x0)*wi)**2/2.)
! if intg = 1, n(x) = x + (ang*sqrt(pi/2)/wi)*
!                         (erf((x-x0)*wi/sqrt(2)) + erf(x0*wi/sqrt(2)))
      implicit none
      integer intg
      double precision x, ang, x0, wi
! local data
      double precision DGDISTR1, f, sqrt2i, sqtpih, aw, t, derfn
      external derfn
      data sqrt2i, sqtpih /0.7071067811865476,1.253314137397325/
      save sqrt2i, sqtpih
      aw = wi*sqrt2i
      t = (x - x0)*aw
      if (intg.eq.0) then
         if (abs(t).lt.8.0d0) then
            f = 1.0d0 + ang*exp(-t**2)
         else
            f = 1.0d0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.0d0) then
            f = (1.0d0 + ang)*x
         else
            f = x + (ang*sqtpih/wi)*(derfn(t) + derfn(x0*aw))
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DGDISTR1 Error: f = ', f
      DGDISTR1 = f
      return
      end
!-----------------------------------------------------------------------
      function DHDISTR1(x,anh,wi,x0,intg)
! this function calculates either a density function or its integral
! for a hyperbolic secant squared density profile.  Used in initializing
! particle coordinates.
! if intg = 0, n(x) = 1.0 + anh*sech((x-x0)*wi)**2
! if intg = 1, n(x) = x + (anh/wi)*(tanh((x-x0)*wi) + tanh(x0*wi))
      implicit none
      integer intg
      double precision x, anh, x0, wi
! local data
      double precision DHDISTR1, f, g, t, u
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (abs(t).lt.32.0d0) then
            u = exp(-abs(t))
            f = 1.0d0 + anh*(2.0d0*u/(1.0d0 + u*u))**2
         else
            f = 1.0d0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.0d0) then
            f = (1.0d0 + anh)*x
         else
            if (abs(t).lt.32.0d0) then
               u = exp(-abs(t))**2
               f = (1.0d0 - u)/(1.0d0 + u)
            else
               f = 1.0d0
            endif
            if (t.lt.0.0d0) f = -f
            t = x0*wi
            if (abs(t).lt.32.0d0) then
               u = exp(-abs(t))**2
               g = (1.0d0 - u)/(1.0d0 + u)
            else
               g = 1.0d0
            endif
            if (t.lt.0.0d0) g = -g
            f = x + (anh/wi)*(f + g)
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DHDISTR1 Error: f = ', f
      DHDISTR1 = f
      return
      end
!-----------------------------------------------------------------------
      function DEDISTR1(x,ane,wi,x0,intg)
! this function calculates either a density function or its integral
! for an exponential density profile.  Used in initializing particle
! coordinates.
! if intg = 0, n(x) = 1.0 + ane*exp((x-x0)*wi)
! if intg = 1, n(x) = x + (ane/wi)*(exp((x-x0)*wi) - exp(-x0*wi)
      implicit none
      integer intg
      double precision x, ane, x0, wi
! local data
      double precision DEDISTR1, f, t
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (t.gt.(-64.0d0)) then
            f = 1.0d0 + ane*exp(t)
         else
            f = 1.0d0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.0d0) then
            f = (1.0d0 + ane)*x
         else
            f = x0*wi
            if (f.gt.64) then
               f = 0.0d0
            else
               f = exp(-f)
            endif
            f = x + (ane/wi)*(exp(t) - f)
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DEDISTR1 Error: f = ', f
      DEDISTR1 = f
      return
      end
!-----------------------------------------------------------------------
      function DGDISTR0(x,ang,wi,x0,intg)
! this function calculates either a density function or its integral
! for a gaussian density profile.  Used in initializing particle
! coordinates.  No background density
! if intg = 0, n(x) = ang*exp(-((x-x0)*wi)**2/2.)
! if intg = 1, n(x) = (ang*sqrt(pi/2)/wi)*
!                         (erf((x-x0)*wi/sqrt(2)) + erf(x0*wi/sqrt(2)))
      implicit none
      integer intg
      double precision x, ang, x0, wi
! local data
      double precision DGDISTR0, f, sqrt2i, sqtpih, aw, t, derfn
      external derfn
      data sqrt2i, sqtpih /0.7071067811865476,1.253314137397325/
      save sqrt2i, sqtpih
      aw = wi*sqrt2i
      t = (x - x0)*aw
      if (intg.eq.0) then
         if (abs(t).lt.8.0d0) then
            f = ang*exp(-t**2)
         else
            f = 0.0d0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.0d0) then
            f = ang*x
         else
            f = (ang*sqtpih/wi)*(derfn(t) + derfn(x0*aw))
         endif
      else
         f = -1.0d0
      endif
      if (f.lt.0.0d0) write (*,*) 'DGDISTR0 Error: f = ', f
      DGDISTR0 = f
      return
      end
!-----------------------------------------------------------------------
      function derfn(x)
! this function calculates the real error function, according to the
! formulae given in Abramowitz and Stegun, Handbook of Mathematical
! Functions, p. 299.  Error is < 1.5 x 10-7.
      implicit none
      double precision x
! local data
      double precision derfn, p, a1, a2, a3, a4, a5, t, f
      data p, a1, a2 /0.3275911,0.254829592,-0.284496736/
      data a3, a4, a5 /1.421413741,-1.453152027,1.061405429/
      save p, a1, a2, a3, a4, a5
      f = abs(x)
      t = 1.0d0/(1.0d0 + p*f)
      if (f.le.8.0d0) then
         derfn = 1.0d0 - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))       &
     &                    *exp(-x*x)
      else
         derfn = 1.0d0
      endif
      if (x.lt.0.0d0) derfn = -derfn
      return
      end
