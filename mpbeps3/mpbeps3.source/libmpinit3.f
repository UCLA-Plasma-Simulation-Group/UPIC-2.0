!-----------------------------------------------------------------------
! Fortran Library for initialization of domains and particles
! 3D MPI/OpenMP PIC Codes:
! NEXTRAN3 = skips over nextrand groups of random numbers
! PDICOMP32L determines uniform integer spatial decomposition for
!            uniform distribution of particles for 3d code
!            integer boundaries set, of equal size except for remainders
! PDNICOMP2L determines uniform integer spatial decomposition for
!            uniform distribution of particles for 3d code
!            integer boundaries set, but might be of unequal size
! PDCOMP32L determines uniform real spatial boundaries for uniform
!           distribution of particles for 3d code
!           real boundaries set, of equal size
! FCOMP32 determines optimal partition for nvp processors
! PDISTR32 calculates initial particle co-ordinates and velocities with
!          uniform density and maxwellian velocity with drift
!          for 3d code
! PRDISTR32 calculates initial particle co-ordinates and momenta with
!           uniform density and maxwell-juttner momentum with drift
!           for 3d code with relativistic particles
! PUDISTR32 calculates initial particle co-ordinates with uniform density
!           for 3d code
! PLDISTR32 calculates initial particle co-ordinates with linear density
!           profile for 3d code
! PFDISTR32 calculates initial particle co-ordinates with general
!           distribution in space for 3d code
! PVDISTR32 calculates initial particle velocities with maxwellian
!           velocity with drift for 3d code
! PVRDISTR32 calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 3d code with relativistic
!            particles
! PVBDISTR32 calculates initial particle velocities for a magnetized
!            plasma with maxwellian velocity with drift in direction
!            parallel to B, ring distribution in directions
!            perpendicular to B for 3d code
! PPDBLKP3L finds the maximum number of particles in each tile
! PFEDGES32 = finds new partitions boundaries (edges,noff,nyzp)
!             from analytic general density profiles.
! PFHOLES32 determines list of particles which are leaving this node
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
! update: november 1, 2020
!-----------------------------------------------------------------------
      subroutine NEXTRAN3(nextrand,ndim,np)
! for 2d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
! nextrand = (0,N) = generate (default,Nth block) of random numbers
! ndim = number of velocity dimensions = 3
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
      subroutine PDICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny, &
     &nz,kstrt,nvpy,nvpz,idps,idds)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! integer boundaries are set, of equal size except for remainders
! nvpy must be < ny and nvpz must be < nz.
! some combinations of ny and nvpy and nz and nvpz result in a zero
! value of nyzp.  this is not supported.
! input: ny, nz, kstrt, nvpy, nvpz, idps, idds
! output: edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
! for 2D spatial decomposition
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! nypmn = minimum value of nyzp(1)
! nzpmn = minimum value of nyzp(2)
! ny/nz = system length in y/z direction
! kstrt = starting data block number (processor id + 1)
! nvpy/nvpz = number of real or virtual processors in y/z
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nypmx, nzpmx, nypmn, nzpmn, ny, nz, kstrt, nvpy, nvpz
      integer idps, idds
      integer nyzp, noff
      real edges
      dimension nyzp(idds), noff(idds)
      dimension edges(idps)
! local data
      integer jb, kb, kyp, kzp
      real at1, at2, any, anz
      integer myzpm, iwork4
      dimension myzpm(4), iwork4(4)
      any = real(ny)
      anz = real(nz)
! determine decomposition
! find processor id in y/z
      kb = (kstrt - 1)/nvpy
      jb = kstrt - nvpy*kb - 1
! boundaries in y
      kyp = (ny - 1)/nvpy + 1
      at1 = real(kyp)
      edges(1) = at1*real(jb)
      if (edges(1).gt.any) edges(1) = any
      noff(1) = edges(1)
      edges(2) = at1*real(jb + 1)
      if (edges(2).gt.any) edges(2) = any
      jb = edges(2)
      nyzp(1) = jb - noff(1)
! boundaries in z
      kzp = (nz - 1)/nvpz + 1
      at2 = real(kzp)
      edges(3) = at2*real(kb)
      if (edges(3).gt.anz) edges(3) = anz
      noff(2) = edges(3)
      edges(4) = at2*real(kb + 1)
      if (edges(4).gt.anz) edges(4) = anz
      kb = edges(4)
      nyzp(2) = kb - noff(2)
! find maximum/minimum partition size in y and z
      myzpm(1) = nyzp(1)
      myzpm(2) = -nyzp(1)
      myzpm(3) = nyzp(2)
      myzpm(4) = -nyzp(2)
      call PPIMAX(myzpm,iwork4,4)
      nypmx = myzpm(1) + 1
      nypmn = -myzpm(2)
      nzpmx = myzpm(3) + 1
      nzpmn = -myzpm(4)
      return
      end
!-----------------------------------------------------------------------
      subroutine PDNICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,&
     &nz,kstrt,nvpy,nvpz,idps,idds)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! integer boundaries are set, of equal size except for remainders
! input: ny, nz, kstrt, nvpy, nvpz, idps, idds
! output: edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
! for 2D spatial decomposition
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! nypmn = minimum value of nyzp(1)
! nzpmn = minimum value of nyzp(2)
! ny/nz = system length in y/z direction
! kstrt = starting data block number (processor id + 1)
! nvpy/nvpz = number of real or virtual processors in y/z
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nypmx, nzpmx, nypmn, nzpmn, ny, nz, kstrt, nvpy, nvpz
      integer idps, idds
      integer nyzp, noff
      real edges
      dimension nyzp(idds), noff(idds)
      dimension edges(idps)
! local data
      integer jb, kb
      real at1, at2
      integer myzpm, iwork4
      dimension myzpm(4), iwork4(4)
! determine decomposition
! find processor id in y/z
      kb = (kstrt - 1)/nvpy
      jb = kstrt - nvpy*kb - 1
! boundaries in y
      at1 = real(ny)/real(nvpy)
      at2 = at1*real(jb)
      noff(1) = at2
      edges(1) = real(noff(1))
      at2 = at1*real(jb + 1)
      if (jb.eq.nvpy) at2 = real(ny)
      jb = at2
      edges(2) = real(jb)
      nyzp(1) = jb - noff(1)
! boundaries in z
      at1 = real(nz)/real(nvpz)
      at2 = at1*real(kb)
      noff(2) = at2
      edges(3) = real(noff(2))
      at2 = at1*real(kb + 1)
      if (kb.eq.nvpz) at2 = real(nz)
      kb = at2
      edges(4) = real(kb)
      nyzp(2) = kb - noff(2)
! find maximum/minimum partition size in y and z
      myzpm(1) = nyzp(1)
      myzpm(2) = -nyzp(1)
      myzpm(3) = nyzp(2)
      myzpm(4) = -nyzp(2)
      call PPIMAX(myzpm,iwork4,4)
      nypmx = myzpm(1) + 1
      nypmn = -myzpm(2)
      nzpmx = myzpm(3) + 1
      nzpmn = -myzpm(4)
      return
      end
!-----------------------------------------------------------------------
      subroutine PDCOMP32L(edges,nyzp,myzp,lyzp,noff,nypmx,nzpmx,ny,nz, &
     &kstrt,nvpy,nvpz,idps,idds)
! this subroutine determines spatial boundaries for uniform particle
! decomposition, calculates number of grid points in each spatial
! region, and the offset of these grid points from the global address
! input: ny, nz, kstrt, nvpy, nvpz, idps, idds
! output: edges, nyzp, myzp, lyzp, noff, nypmx, nzpmx
! for 2D spatial decomposition
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! myzp(1:2) = number of full or partial grids in y/z
! lyzp(1:2) = number of guard cells for processor on left in y/z
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! ny/nz = system length in y/z direction
! kstrt = starting data block number (processor id + 1)
! nvpy/nvpz = number of real or virtual processors in y/z
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nypmx, nzpmx, ny, nz, kstrt, nvpy, nvpz, idps, idds
      integer nyzp, myzp, lyzp, noff
      real edges
      dimension nyzp(idds), myzp(idds), lyzp(idds), noff(idds)
      dimension edges(idps)
! local data
      integer jb, kb
      real at1, at2, dt1, dt2
      integer myzpm, iwork2
      dimension myzpm(2), iwork2(2)
! determine decomposition
      kb = (kstrt - 1)/nvpy
      jb = kstrt - nvpy*kb - 1
      at1 = real(ny)/real(nvpy)
      at2 = real(nz)/real(nvpz)
! boundaries in y
      edges(1) = at1*real(jb)
      noff(1) = edges(1)
      dt1 = edges(1) - real(noff(1))
      if (dt1.eq.0.0) edges(1) = real(noff(1))
      edges(2) = at1*real(jb + 1)
      if (kstrt.eq.nvpy) edges(2) = real(ny)
      jb = edges(2)
      dt1 = edges(2) - real(jb)
      if (dt1.eq.0.0) edges(2) = real(jb)
      nyzp(1) = jb - noff(1)
      myzp(1) = nyzp(1)
      if (dt1.gt.0.0) myzp(1) = myzp(1) + 1
! boundaries in z
      edges(3) = at2*real(kb)
      noff(2) = edges(3)
      dt2 = edges(3) - real(noff(2))
      if (dt2.eq.0.0) edges(3) = real(noff(2))
      edges(4) = at2*real(kb + 1)
      if (kstrt.eq.nvpz) edges(4) = real(nz)
      kb = edges(4)
      dt2 = edges(4) - real(kb)
      if (dt2.eq.0.0) edges(4) = real(kb)
      nyzp(2) = kb - noff(2)
      myzp(2) = nyzp(2)
      if (dt1.gt.0.0) myzp(2) = myzp(2) + 1
! find number of guard cells on the left in y and z
      myzpm(1) = myzp(1) - nyzp(1) + 1
      myzpm(2) = myzp(2) - nyzp(2) + 1
      call PPISHFTR2(myzpm,iwork2,1,nvpy,nvpz)
      lyzp(1) = myzpm(1)
      lyzp(2) = myzpm(2)
! find maximum partition size in y and z
      myzpm(1) = myzp(1)
      myzpm(2) = myzp(2)
      call PPIMAX(myzpm,iwork2,2)
      nypmx = myzpm(1) + 1
      nzpmx = myzpm(2) + 1
      return
      end
!-----------------------------------------------------------------------
      subroutine FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
! determines optimal partition for nvp processors
! input: nvp, number of processors, nx, ny, nz = number of grids
! output: nvpy, nvpz, processors in y, z direction, ierr = error code
! nvp = number of real or virtual processors obtained
! nx/ny/nz = system length in x/y/z direction
! nvpy/nvpz = number of real or virtual processors in y/z
! ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer nvp, nx, ny, nz, nvpy, nvpz, ierr
! local data
      integer nxh, lvp
      double precision dt1
      nxh = nx/2
      ierr = 0
! algorithm 1: prefer equal number of grids in y and z partitions
      dt1 = sqrt(dble(nvp)*dble(ny)/dble(nz))
! algorithm 2: prefer equal number of grids in x*y and y*z partitions
!     dt1 = sqrt(nvp*sqrt(dble(nxh)/dble(nz)))
! return total number of processors in y and z
      nvpy = real(dt1)
      if (nvpy.lt.1) nvpy = 1
      nvpz = nvp/nvpy
      lvp = nvpy*nvpz
      if (lvp.gt.nvp) then
         write (*,*) 'invalid partition:nvpy,nvpz,nvp=', nvpy, nvpz, nvp
         ierr = 1
         return
      endif
   10 if (lvp.ne.nvp) then
         nvpy = nvpy - 1
         nvpz = nvp/nvpy
         lvp = nvpy*nvpz
         go to 10
      endif
      nvp = lvp
      return
      end
!-----------------------------------------------------------------------
      subroutine PDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
! for 3d code, this subroutine calculates initial particle co-ordinates
! and velocities with uniform density and maxwellian velocity with drift
! for distributed data
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, npp, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! part(4,n) = velocity vx of particle n in partition
! part(5,n) = velocity vy of particle n in partition
! part(6,n) = velocity vz of particle n in partition
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = back boundary in z of particle partition
! edges(4) = front boundary in z of particle partition
! npp = number of particles in partition, updated in this procedure
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift velocity of particles in x/y/z direction
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! nx/ny/nz = system length in x/y/z direction
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer npp, npx, npy, npz, nx, ny, nz, idimp, npmax, idps
      integer ipbc, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, l, nps, npt, npxyzp
      real edgelx, edgely, edgelz, at1, at2, at3
      real xt, yt, zt, vxt, vyt, vzt
      double precision dnpxyz, dt1
      integer ierr1, iwork1
      double precision sum4, work4
      dimension ierr1(1), iwork1(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
      nps = npp + 1
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 0.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz)/real(npz)
      endif
      do 30 l = 1, npz
      zt = edgelz + at3*(real(l) - 0.5)
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
! uniform density profile
      xt = edgelx + at1*(real(j) - 0.5)
! maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         if ((zt.ge.edges(3)).and.(zt.lt.edges(4))) then
            npt = npp + 1
            if (npt.le.npmax) then
               part(1,npt) = xt
               part(2,npt) = yt
               part(3,npt) = zt
               part(4,npt) = vxt
               part(5,npt) = vyt
               part(6,npt) = vzt
               npp = npt
            else
               ierr = 1
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
! process error: check if buffer overflow occurred
      if (ierr.eq.0) then
         ierr1(1) = 0  
      else
         ierr1(1) = npt
      endif
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.0) return
! add correct drift
      npxyzp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = npxyzp
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
! process errors
      dnpxyz = dnpxyz - dble(npx)*dble(npy)*dble(npz)
      if (dnpxyz.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PRDISTR32(part,edges,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx&
     &,npy,npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
! for 3d code, this subroutine calculates initial particle co-ordinates
! and momenta with uniform density and maxwell-juttner distribution with
! drift for relativistic particles and distributed data.
! algorithm is only approximate, valid when gamma of distribution < 4
! f(p) = exp(-(gamma-1)*(m0c**2)/kT), where gamma = sqrt(1+(p/m0c)**2)
! since (gamma-1)*(m0c**2) = (p**2/m0)/(gamma+1), we can write
! f(p) = exp(-pt**2/2), where pt**2 = p**2/((gamma+1)/2)*m0*kT
! since pt is normally distributed, we can use a gaussian random number
! to calculate it.  We then solve the pt**2 equation to obtain:
! (p/m0)**2 = ((pt*vth)**2)*(1 + (pt/4c)**2),
! where vth = sqrt(kT/m0) is the thermal velocity.  This equation is
! satisfied if we set the individual components j as follows:
! pj/m0 = ptj*vth*sqrt(1 + (pt/4c)**2)
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, npp, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! part(4,n) = momentum px of particle n in partition
! part(5,n) = momentum py of particle n in partition
! part(6,n) = momentum pz of particle n in partition
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = back boundary in z of particle partition
! edges(4) = front boundary in z of particle partition
! npp = number of particles in partition, updated in this procedure
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift momentum of particles in x/y/z direction
! ci = reciprocal of velocity of light
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! nx/ny/nz = system length in x/y/z direction
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer npp, npx, npy, npz, nx, ny, nz, idimp, npmax, idps
      integer ipbc, ierr
      real vtx, vty, vtz, vdx, vdy, vdz, ci
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, l, nps, npt, npxyzp
      real edgelx, edgely, edgelz, at1, at2, at3
      real ci4, xt, yt, zt, ptx, pty, ptz, pt2
      double precision dnpxyz, dt1
      integer ierr1, iwork1
      double precision sum4, work4
      dimension ierr1(1), iwork1(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
      ci4 = 0.25*ci*ci
      nps = npp + 1
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 0.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz)/real(npz)
      endif
      do 30 l = 1, npz
      zt = edgelz + at3*(real(l) - 0.5)
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
! uniform density profile
      xt = edgelx + at1*(real(j) - 0.5)
! maxwell-juttner momentum distribution
      ptx = vtx*ranorm()
      pty = vty*ranorm()
      ptz = vtz*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         if ((zt.ge.edges(3)).and.(zt.lt.edges(4))) then
            pt2 = ptx*ptx + pty*pty + ptz*ptz
            pt2 = sqrt(1.0 + pt2*ci4)
            npt = npp + 1
            if (npt.le.npmax) then
               part(1,npt) = xt
               part(2,npt) = yt
               part(3,npt) = zt
               part(4,npt) = ptx*pt2
               part(5,npt) = pty*pt2
               part(6,npt) = ptz*pt2
               npp = npt
            else
               ierr = 1
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
! process error: check if buffer overflow occurred
      if (ierr.eq.0) then
         ierr1(1) = 0  
      else
         ierr1(1) = npt
      endif
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.0) return
! add correct drift
      npxyzp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = npxyzp
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
! process errors
      dnpxyz = dnpxyz - dble(npx)*dble(npy)*dble(npz)
      if (dnpxyz.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PUDISTR32(part,edges,npp,npx,npy,npz,nx,ny,nz,idimp,   &
     &npmax,idps,ipbc,ierr)
! for 3d code, this subroutine calculates initial particle co-ordinates
! with uniform density for distributed data
! input: all except part, ierr, output: part, npp, ierr
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = back boundary in z of particle partition
! edges(4) = front boundary in z of particle partition
! npp = number of particles in partition, updated in this procedure
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! nx/ny/nz = system length in x/y/z direction
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! idps = number of partition boundaries
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
! ierr = (0,1) = (no,yes) error condition exists
! with 2D spatial decomposition
      implicit none
      integer npp, npx, npy, npz, nx, ny, nz, idimp, npmax, idps, ipbc
      integer ierr
      real part, edges
      dimension part(idimp,npmax), edges(idps)
! local data
      integer j, k, l, nps, npt
      real edgelx, edgely, edgelz, at1, at2, at3, xt, yt, zt
      double precision dnpxyz
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      ierr = 0
      nps = npp
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      endif
! uniform density profile
      do 30 l = 1, npz
      zt = edgelz + at3*(real(l) - 0.5)
      if ((zt.ge.edges(3)).and.(zt.lt.edges(4))) then
         do 20 k = 1, npy
         yt = edgely + at2*(real(k) - 0.5)
         if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
            do 10 j = 1, npx
            xt = edgelx + at1*(real(j) - 0.5)
            npt = npp + 1
            if (npt.le.npmax) then
               part(1,npt) = xt
               part(2,npt) = yt
               part(3,npt) = zt
               npp = npt
            else
               ierr = 1
            endif
   10       continue
         endif
   20    continue
      endif
   30 continue
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
      dnpxyz = sum1(1) - dble(npx)*dble(npy)*dble(npz)
      if (dnpxyz.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PLDISTR32(part,npp,anlx,anly,anlz,npx,npy,npz,nx,ny,nz,&
     &kstrt,nvpy,nvpz,idimp,npmax,ipbc,ierr)
! for 3d code, this subroutine calculates initial particle co-ordinates
! with the following tri-linear density profile:
! n(x,y,z) = n(x)*n(y)*n(z), 
! where n(x) = n0x*(1. + anlx*(x/nx - .5))
!   and n(y) = n0y*(1. + anly*(y/ny - .5))
!   and n(z) = n0z*(1. + anlz*(z/nz - .5))
! n0x = npx/(nx - 2*edgelx), n0y = npy/(ny - 2*edgely)
! and n0z = npz/(nz - 2*edgelz)
! for distributed data with 2D spatial decomposition
! the algorithm partitions the number of particles distributed in the
! y/z direction uniformly (with possible remainders).
! particles are not necessarily in the correct processor, and may need
! to be moved to another processor.
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! npp = number of particles in partition, updated in this procedure
! anlx/anly/anlz = initial linear density weight in x/y/z direction
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
! ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npp, npx, npy, npz, nx, ny, nz, kstrt, nvpy, nvpz
      integer idimp, npmax, ipbc, ierr
      real anlx, anly, anlz
      real part
      dimension part(idimp,npmax)
! local data
      integer nps, npt, js, ks, j, k, l, kk, ll, mpy, mpys, mpz, mpzs
      integer joff, koff
      real edgelx, edgely, edgelz, at1, at2, at3, bt1, bt2, bt3
      real antx, anty, antz, xt, yt, zt
      double precision dnpxyz
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      nps = npp
! particle distribution constants
      ks = kstrt - 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvpy + 1
      mpys = min(mpy,max(0,npy-mpy*js))
! mpz = number of particles per processor in z direction
      mpz = (npz - 1)/nvpz + 1
      mpzs = min(mpz,max(0,npz-mpz*ks))
! check if particle overflow will occurr
      npt = npp + npx*mpys*mpzs
      ierr1(1) = npt
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.npmax) return
      ierr = 0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      endif
      if (anlx.ne.0.0) then
         antx = anlx/real(nx)
         at1 = 2.0*antx*at1
         bt1 = 1.0 - 0.5*antx*(float(nx) - 2.0*edgelx)
      endif
      if (anly.ne.0.0) then
         anty = anly/real(ny)
         at2 = 2.0*anty*at2
         bt2 = 1.0 - 0.5*anty*(float(ny) - 2.0*edgely)
      endif
      if (anlz.ne.0.0) then
         antz = anlz/real(nz)
         at3 = 2.0*antz*at3
         bt3 = 1.0 - 0.5*antz*(float(nz) - 2.0*edgelz)
      endif
! assign particle spatial coordinates
!$OMP PARALLEL DO PRIVATE(j,k,l,kk,ll,joff,koff,xt,yt,zt)
      do 30 l = 1, mpzs
      ll = l + mpz*ks
      koff = npx*mpys*(l - 1) + npp
! linear density in z
      if (anlz.ne.0.0) then
         zt = edgelz + (sqrt(bt3*bt3 + at3*(real(ll) - 0.5)) - bt3)/antz
! uniform density in z
      else
         zt = edgelz + at3*(real(ll) - 0.5)
      endif
      do 20 k = 1, mpys
      kk = k + mpy*js
      joff = npx*(k - 1) + koff
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
      part(3,j+joff) = zt
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
! update number of particles
      npp = npt
! check if not all particles were distributed
      dnpxyz = dble(npx)*dble(npy)*dble(npz)
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxyz = sum1(1) - dnpxyz
      if (dnpxyz.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PFDISTR32(part,npp,fnx,argx1,argx2,argx3,fny,argy1,    &
     &argy2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,kstrt,nvpy,&
     &nvpz,idimp,npmax,ipbc,ierr)
! for 3d code, this subroutine calculates initial particle co-ordinates
! with general density profile n(x,y) = n(x)*n(y)*n(z), 
! where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
! and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
! and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
! and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
! and where density in z is given by n(z) = fnz(z,argz1,argz2,argz3,0)
! and integral of the density is given by = fnz(z,argz1,argz2,argz3,1)
! for distributed data with 2D spatial decomposition.
! the algorithm partitions the number of particles distributed in the
! y/z direction uniformly (with possible remainders).
! particles are not necessarily in the correct processor, and may need
! to be moved to another processor.
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! npp = number of particles in partition, updated in this procedure
! fnx/fny/fnz = density and density integral function in x/y/z direction
! argx1,argx2,argx3 = arguments to fnx
! argy1,argy2,argy3 = arguments to fny
! argz1,argz2,argz3 = arguments to fnz
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
! ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npp, npx, npy, npz, nx, ny, nz, kstrt, nvpy, nvpz
      integer idimp, npmax, ipbc, ierr
      double precision argx1, argx2, argx3, argy1, argy2, argy3
      double precision argz1, argz2, argz3
      real part
      dimension part(idimp,npmax)
      double precision fnx, fny, fnz
      external fnx, fny, fnz
! local data
      integer nps, npt, js, ks, i, j, k, l, kk, ll, mpy, mpys, mpz, mpzs
      integer imax, joff, koff, moff, loff
      real edgelx, edgely, edgelz, anx, any, anz, bnx, bny, bnz
      real xt0, yt0, zt0, x0, y0, z0, xn, yn, zn, eps, big, f, fp
      double precision xt, yt, zt, dnpxyz
      integer ierr1, iwork1
      double precision sum1, work1
      dimension ierr1(1), iwork1(1), sum1(1), work1(1)
      nps = npp
! particle distribution constants
      ks = kstrt - 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvpy + 1
      mpys = min(mpy,max(0,npy-mpy*js))
! mpz = number of particles per processor in z direction
      mpz = (npz - 1)/nvpz + 1
      mpzs = min(mpz,max(0,npz-mpz*ks))
! check if particle overflow will occurr
      npt = npp + npx*mpys*mpzs
      ierr1(1) = npt
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.gt.npmax) return
      ierr = 0
! eps = convergence criterion
      imax = max(nx,ny,nz)
      eps = 0.0001
      big = 0.25
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
      endif
! find normalization for function
      anx = real(nx) - edgelx
      any = real(ny) - edgely
      anz = real(nz) - edgelz
      x0 = fnx(dble(edgelx),argx1,argx2,argx3,1)
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
      z0 = fnz(dble(edgelz),argz1,argz2,argz3,1)
      bnx = real(npx)/(fnx(dble(anx),argx1,argx2,argx3,1) - x0)
      bny = real(npy)/(fny(dble(any),argy1,argy2,argy3,1) - y0)
      bnz = real(npz)/(fnz(dble(anz),argz1,argz2,argz3,1) - z0)
      x0 = bnx*x0 - 0.5
      y0 = bny*y0 - 0.5
      z0 = bnz*z0 - 0.5
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
      moff = mpy*js
      do 40 k = 1, mpys + moff
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
         part(2,joff+1) = yt
      endif
      yt0 = yt
   40 continue
! quit if error
      if (ierr.ne.0) return
! density profile in z
      loff = mpz*ks
      do 80 l = 1, mpzs + loff
      ll = l - loff
      zn = real(l) + z0
! guess next value for zt
      if (l.eq.1) then
         zt0 = edgelz
         zt = zt0
         fp = bnz*fnz(zt,argz1,argz2,argz3,0)
         if (fp.gt.0.0) zt = zt + 0.5/fp
      else
         fp = bnz*fnz(zt,argz1,argz2,argz3,0)
         if (fp.eq.0.0) fp = 1.0
         zt = zt + 1.0/fp
      endif
      zt = max(edgelz,min(zt,anz))
      i = 0
   50 f = bnz*fnz(zt,argz1,argz2,argz3,1) - dble(zn)
! find improved value for zt
      if (abs(f).ge.eps) then
         fp = bnz*fnz(zt,argz1,argz2,argz3,0)
! newton's method
         if ((abs(f).lt.big).and.(fp.gt.0.0)) then
            zt0 = zt
            zt = zt - f/fp
            zt = max(edgelz,min(zt,anz))
! bisection method
         else if (f.gt.0.0) then
            fp = 0.5*abs(zt0 - zt)
            zt = zt0 - fp
         else
            fp = abs(zt - zt0)
            zt0 = zt
            zt = zt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 50
!        write (*,*) l,'newton iteration max exceeded, zt = ', zt
         ierr = ierr + 1
      endif
! store co-ordinates
      if ((ll.ge.1).and.(ll.le.mpzs)) then
         koff = npx*mpys*(ll - 1)
         do 70 k = 1, mpys
         joff = npx*(k - 1) + npp
         do 60 j = 1, npx
         part(1,j+joff+koff) = part(1,j+npp)
         part(2,j+joff+koff) = part(2,joff+1)
         part(3,j+joff+koff) = zt
   60    continue
   70    continue
      endif
      zt0 = zt
   80 continue
! quit if error
      if (ierr.ne.0) return
! update number of particles
      npp = npt
! check if not all particles were distributed
      dnpxyz = dble(npx)*dble(npy)*dble(npz)
      sum1(1) = dble(npp-nps)
      call PPDSUM(sum1,work1,1)
      dnpxyz = sum1(1) - dnpxyz
      if (dnpxyz.ne.0.0d0) ierr = -1
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &npz,kstrt,nvpy,nvpz,idimp,npmax,ierr)
! for 3d code, this subroutine calculates initial particle velocities
! with maxwellian velocity with drift for distributed data
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(4,n) = velocity vx of particle n in partition
! part(5,n) = velocity vy of particle n in partition
! part(6,n) = velocity vz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift velocity of particles in x/y/z direction
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer nps, npp, npx, npy, npz, kstrt, nvpy, nvpz, idimp, npmax
      integer ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyzp, js, ks, j, k, l, kk, ll, mpy, mpys, mpz, mpzs
      real vxt, vyt, vzt
      double precision dnpxyz, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision ranorm
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvpy + 1
      mpys = min(mpy,max(0,npy-mpy*js))
! mpz = number of particles per processor in z direction
      mpz = (npz - 1)/nvpz + 1
      mpzs = min(mpz,max(0,npz-mpz*ks))
! maxwellian velocity distribution
      npt = nps - 1
      do 30 ll = 1, npz
      l = ll - mpz*ks
      do 20 kk = 1, npy
      k = kk - mpy*js
      do 10 j = 1, npx
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((l.ge.1).and.(l.le.mpzs)) then
         if ((k.ge.1).and.(k.le.mpys)) then
            npt = npt + 1
            if (npt.le.npmax) then
               part(4,npt) = vxt
               part(5,npt) = vyt
               part(6,npt) = vzt
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyzp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = dble(npxyzp)
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVRDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,ci,npx,&
     &npy,npz,kstrt,nvpy,nvpz,idimp,npmax,ierr)
! for 3d code, this subroutine calculates initial particle velocities
! momentum with maxwell-juttner distribution with drift
! for relativistic particles and distributed data.
! algorithm is only approximate, valid when gamma of distribution < 4
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
! part(4,n) = momentum px of particle n in partition
! part(5,n) = momentum py of particle n in partition
! part(6,n) = momentum pz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift momentum of particles in x/y/z direction
! ci = reciprocal of velocity of light
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer nps, npp, npx, npy, npz, kstrt, nvpy, nvpz, idimp, npmax
      integer ierr
      real vtx, vty, vtz, vdx, vdy, vdz, ci
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyzp, js, ks, j, k, l, kk, ll, mpy, mpys, mpz, mpzs
      real ci4, ptx, pty, ptz, pt2
      double precision dnpxyz, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision ranorm
      ci4 = 0.25*ci*ci
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvpy + 1
      mpys = min(mpy,max(0,npy-mpy*js))
! mpz = number of particles per processor in z direction
      mpz = (npz - 1)/nvpz + 1
      mpzs = min(mpz,max(0,npz-mpz*ks))
! maxwell-juttner momentum distribution
      npt = nps - 1
      do 30 ll = 1, npz
      l = ll - mpz*ks
      do 20 kk = 1, npy
      k = kk - mpy*js
      do 10 j = 1, npx
      ptx = vtx*ranorm()
      pty = vty*ranorm()
      ptz = vtz*ranorm()
      if ((l.ge.1).and.(l.le.mpzs)) then
         if ((k.ge.1).and.(k.le.mpys)) then
            pt2 = ptx*ptx + pty*pty + ptz*ptz
            pt2 = sqrt(1.0 + pt2*ci4)
            npt = npt + 1
            if (npt.le.npmax) then
               part(4,npt) = ptx*pt2
               part(5,npt) = pty*pt2
               part(6,npt) = ptz*pt2
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyzp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = dble(npxyzp)
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PVBDISTR32(part,nps,npp,vtr,vtz,vdr,vdz,omx,omy,omz,npx&
     &,npy,npz,kstrt,nvpy,nvpz,idimp,npmax,ierr)
! for 3d code, this subroutine calculates initial particle velocities
! velocities for a magnetized plasma with maxwellian velocity with
! drift in direction parallel to B, ring distribution in directions
! perpendicular to B.  The algorithm due to phil pritchett proceeds by
! first assuming B is in the z direction, then rotating to the
! co-ordinates to match the actual B direction
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, ierr
! part(4,n) = velocity vx of particle n in partition
! part(5,n) = velocity vy of particle n in partition
! part(6,n) = velocity vz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtr/vtz = thermal velocity of particles in azimuthal/B direction
! vdt/vdz = drift velocity of particles in azimuthal/B direction
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! randum = uniform random number
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer nps, npp, npx, npy, npz, kstrt, nvpy, nvpz, idimp, npmax
      integer ierr
      real vtr, vtz, vdr, vdz, omx, omy, omz
      real part
      dimension part(idimp,npmax)
! local data
      integer npt, npxyzp, js, ks, j, k, l, kk, ll, mpy, mpys, mpz, mpzs
      integer ndir
      real twopi, at1, at2, vxt, vyt, vzt
      real ox, oy, oz, px, py, pz, qx, qy, qz
      double precision dnpxyz, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision sum4, work4
      dimension sum4(4), work4(4)
      double precision randum, ranorm
      ierr = 0
! particle distribution constants
      ks = kstrt - 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
! mpy = number of particles per processor in y direction
      mpy = (npy - 1)/nvpy + 1
      mpys = min(mpy,max(0,npy-mpy*js))
! mpz = number of particles per processor in z direction
      mpz = (npz - 1)/nvpz + 1
      mpzs = min(mpz,max(0,npz-mpz*ks))
      twopi = 6.28318530717959
! maxwell velocity distribution with ring distribution in x and y,
! maxwell velocity distribution in z
      npt = nps - 1
      do 30 ll = 1, npz
      l = ll - mpz*ks
      do 20 kk = 1, npy
      k = kk - mpy*js
      do 10 j = 1, npx
      at1 = twopi*randum()
      vxt = vdr*cos(at1) + vtr*ranorm()
      vyt = vdr*sin(at1) + vtr*ranorm()
      vzt = vtz*ranorm()
      if ((l.ge.1).and.(l.le.mpzs)) then
         if ((k.ge.1).and.(k.le.mpys)) then
            npt = npt + 1
            if (npt.le.npmax) then
               part(4,npt) = vxt
               part(5,npt) = vyt
               part(6,npt) = vzt
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
! process particle number error
      ierr = npp - npt
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      if (ierr.ne.0) return
! add correct drift
      npxyzp = 0
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = dble(npxyzp)
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1)
      sum4(2) = dt1*sum4(2)
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
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
      do 60 j = nps, npp
      vxt = part(4,j)
      vyt = part(5,j)
      vzt = part(6,j)
! store components in appropriate direction
      part(4,j) = vzt*ox + vxt*px + vyt*qx
      part(5,j) = vzt*oy + vxt*py + vyt*qy
      part(6,j) = vzt*oz + vxt*pz + vyt*qz
   60 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPDBLKP3L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mz&
     &,mx1,myp1,mxyzp1,idds,irc)
! this subroutine finds the maximum number of particles in each tile of
! mx, my, mz to calculate size of segmented particle array ppart
! linear interpolation, spatial decomposition in y/z direction
! input: all except kpic, nppmx, output: kpic, nppmx
! part = input particle array
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! kpic = output number of particles per tile
! npp = number of particles in partition
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nppmx = return maximum number of particles in tile
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! mx/my/mz = number of grids in sorting cell in x, y and z
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, nppmx, idimp, npmax, mx, my, mz, mx1, myp1, mxyzp1
      integer idds, irc
      integer kpic, noff
      real part
      dimension part(idimp,npmax), kpic(mxyzp1)
      dimension noff(idds)
! local datal, 
      integer j, k, n, m, l, mnoff, lnoff, mxyp1, isum, ist, npx, ierr
      mnoff = noff(1)
      lnoff = noff(2)
      ierr = 0
      mxyp1 = mx1*myp1
! clear counter array
      do 10 k = 1, mxyzp1
      kpic(k) = 0
   10 continue
! find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      l = part(3,j)
      l = (l - lnoff)/mz
      m = n + mx1*m + mxyp1*l
      if (m.le.mxyzp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyzp1)
      endif
   20 continue
! find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyzp1
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
      subroutine PFEDGES32(edges,nyzp,noff,fny,argy1,argy2,argy3,fnz,   &
     &argz1,argz2,argz3,nypmx,nzpmx,nypmn,nzpmn,ny,nz,kstrt,nvpy,nvpz,  &
     &idps,idds,ipbc)
! this subroutines finds new partitions boundaries (edges,noff,nyzp)
! from density integrals given by = fny(z,argy1,argy2,argy3,1) and
! fnz(z,argz1,argz2,argz3,1)
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! fny/fnz = density and density integral function in y/z
! argy1,argy2,argy3 = arguments to fny
! argz1,argz2,argz3 = arguments to fnz
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! nypmn = minimum value of nyzp(1)
! nzpmn = minimum value of nyzp(2)
! ny/nz = system length in y/z direction
! kstrt = starting data block number (processor id + 1)
! nvpy/nvpz = number of real or virtual processors in y/z
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nypmx, nzpmx, nypmn, nzpmn, ny, nz, kstrt, nvpy, nvpz
      integer idps, idds, ipbc
      double precision argy1, argy2, argy3, argz1, argz2, argz3
      integer nyzp, noff
      real edges
      dimension nyzp(idds), noff(idds)
      dimension edges(idps)
      double precision fny, fnz
      external fny, fnz
! local data
      integer jb, kb
      real edgely, edgelz, any1, any, anz1, anz, y0, y1, z0, z1
      real anpav, anpl, anpr, sum1, at1, at2
      integer myzpm, iwork4
      dimension myzpm(4), iwork4(4)
! determine decomposition
! find processor id in y/z
      kb = (kstrt - 1)/nvpy
      jb = kstrt - nvpy*kb - 1
! set boundary values
      edgely = 0.0
      edgelz = 0.0
      if (ipbc.eq.2) then
         edgely = 1.0
         edgelz = 1.0
      else if (ipbc.eq.3) then
         edgely = 1.0
      endif
! find normalization for function
      any = real(ny)
      anz = real(nz)
      any1 = any - edgely
      anz1 = anz - edgelz
      y0 = fny(dble(edgely),argy1,argy2,argy3,1)
      z0 = fnz(dble(edgelz),argz1,argz2,argz3,1)
! anpav = desired number of particles per processor in y
      anpav = (fny(dble(any1),argy1,argy2,argy3,1) - y0)/real(nvpy)
! search for boundaries in y
      anpl = real(jb)*anpav
      anpr = real(jb+1)*anpav
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
      if (jb.eq.0) edges(1) = 0.0
! then find right boundary
   20 at1 = sum1
      sum1 = fny(dble(y1),argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpr).and.(y1.le.any)) go to 20
      at2 = (y1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(2) = at2
! set rightmost edge to ny
      if ((jb+1).eq.nvpy) edges(2) = any
! anpav = desired number of particles per processor in z
      anpav = (fnz(dble(anz1),argz1,argz2,argz3,1) - z0)/real(nvpz)
! search for boundaries in z
      anpl = real(kb)*anpav
      anpr = real(kb+1)*anpav
      z1 = edgelz
      sum1 = 0.0
! first find left boundary
   30 at1 = sum1
      sum1 = fnz(dble(z1),argz1,argz2,argz3,1) - z0
      z1 = z1 + 1.0
      if ((sum1.lt.anpl).and.(z1.le.anz)) go to 30 
      if (sum1.gt.at1) then
         at2 = (z1 - 2.0) + (anpl - at1)/(sum1 - at1)
      else
         at2 = z1 - 1.0
      endif
      edges(3) = at2
! set leftmost edge to zero
      if (kb.eq.0) edges(3) = 0.0
! then find right boundary
   40 at1 = sum1
      sum1 = fnz(dble(z1),argz1,argz2,argz3,1) - z0
      z1 = z1 + 1.0
      if ((sum1.lt.anpr).and.(z1.le.anz)) go to 40
      at2 = (z1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(4) = at2
! set rightmost edge to nz
      if ((kb+1).eq.nvpz) edges(4) = anz
! calculate number of grids and offsets in new partitions
      noff(1) = edges(1) + 0.5
      jb = edges(2) + 0.5
      nyzp(1) = jb - noff(1)
      edges(1) = real(noff(1))
      edges(2) = real(jb)
      noff(2) = edges(3) + 0.5
      kb = edges(4) + 0.5
      nyzp(2) = kb - noff(2)
      edges(3) = real(noff(2))
      edges(4) = real(kb)
! find maximum/minimum partition size in y and z
      myzpm(1) = nyzp(1)
      myzpm(2) = -nyzp(1)
      myzpm(3) = nyzp(2)
      myzpm(4) = -nyzp(2)
      call PPIMAX(myzpm,iwork4,4)
      nypmx = myzpm(1) + 1
      nypmn = -myzpm(2)
      nzpmx = myzpm(3) + 1
      nzpmn = -myzpm(4)
      return
      end
!-----------------------------------------------------------------------
      subroutine PFHOLES32(part,edges,npp,ihole,nc,idimp,npmax,idps,    &
     &ntmax)
! for 3d code, this subroutine determines list of particles which are
! leaving this processor, with 2D spatial decomposition
! input: all except ihole, output: ihole
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! npp = number of particles in partition
! ihole(:,2) = location of holes left in y/z in particle arrays
! ihole(1,:) = ih, number of holes left in y/z (error, if negative)
! nc = (1,2) = (normal,tagged) partitioned co-ordinate to sort by
! idimp = size of phase space = 6 or 7
! npmax = maximum number of particles in each partition
! idps = number of particle partition boundaries = 4
! ntmax = size of hole array for particles leaving processors
      implicit none
      integer npp, nc, idimp, npmax, idps, ntmax
      real part, edges
      integer ihole
      dimension part(idimp,npmax)
      dimension edges(idps), ihole(ntmax+1,2)
! local data
      integer j, iy, iz, ih1, ih2, nh
      real dy, dz
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      if ((nc.lt.1).or.(nc.gt.2)) return
! determine co-ordinate to sort by
      if (nc.eq.1) then
         iy = 2
         iz = 3
      else if (idimp.gt.6) then
         iy = 7
         iz = 7
      else
         return
      endif
      ih1 = 0
      ih2 = 0
      nh = 0
      do 10 j = 1, npp
      dy = part(iy,j)
      dz = part(iz,j)
! find particles out of bounds
! check particles leaving in y direction or y and z
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih1 = ih1 + 1
         if (ih1.le.ntmax) then
            ihole(ih1+1,1) = j
         else
            nh = 1
         endif
! check particles leaving in z direction only
      else if ((dz.lt.edges(3)).or.(dz.ge.edges(4))) then
         ih2 = ih2 + 1
         if (ih2.le.ntmax) then
            ihole(ih2+1,2) = j
         else
            nh = 1
         endif
      endif
   10 continue
! set end of file flag
      if (nh.gt.0) ih1 = -max(ih1,ih2)
      ihole(1,1) = ih1
      ihole(1,2) = ih2
! set global error flag if needed
      ierr1(1) = -ih1
      call PPIMAX(ierr1,iwork1,1)
      nh = -ierr1(1)
      if (nh.lt.0) ihole(1,1) = nh
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
