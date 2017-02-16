!-----------------------------------------------------------------------
! Fortran Library for initialization of domains and particles
! 3D MPI/OpenMP PIC Codes:
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
! PVDISTR32 calculates initial particle velocities with maxwellian
!           velocity with drift for 3d code
! PVRDISTR32 calculates initial particle momenta with maxwell-juttner
!            distribution with drift for 3d code with relativistic
!            particles
! PPDBLKP3L finds the maximum number of particles in each tile
! ranorm gaussian random number generator
! randum uniform random number generator
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 15, 2017
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
      subroutine PVDISTR32(part,nps,npp,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &npz,idimp,npmax,ierr)
! for 3d code, this subroutine calculates initial particle velocities
! with maxwellian velocity with drift for distributed data
! algorithm designed to assign the same random numbers to particle
! velocities, independent of the number of processors
! input: all except part, ierr, output: part, npp, ierr
! part(4,n) = velocity vx of particle n in partition
! part(5,n) = velocity vy of particle n in partition
! part(6,n) = velocity vz of particle n in partition
! nps = starting address of particles in partition
! npp = number of particles in partition
! vtx/vty/vtz = thermal velocity of particles in x/y/z direction
! vdx/vdy/vdz = drift velocity of particles in x/y/z direction
! npx/npy/npz = initial number of particles distributed in x/y/z
! direction
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer nps, npp, npx, npy, npz, idimp, npmax, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,npmax)
! local data
      integer j, k, l, npt, npxyzp
      real vxt, vyt, vzt
      double precision dps, dnpx, dnpxy, dkoff, djoff, dj, dnpxyz, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision dpp, dpt, sum4, work4
      dimension dpp(1), dpt(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
! determine offsets for particle number
! each random triplet is associated with a global particle index
! 1:sum(npp-nps+1)
      dps = dble(npp-nps+1)
      dnpx = dble(npx)
      dnpxy = dble(npy)*dnpx
      dpp(1) = dps
      call PPDSCAN(dpp,dpt,1)
      dps = dpp(1) - dps
! maxwellian velocity distribution
      npt = nps - 1
      do 30 l = 1, npz
      dkoff = dnpxy*dble(l - 1)
      do 20 k = 1, npy
      djoff = dnpx*dble(k - 1) + dkoff
      do 10 j = 1, npx
      dj = dble(j - 1) + djoff
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((dj.ge.dps).and.(dj.lt.dpp(1))) then
         npt = npt + 1
         if (npt.le.npmax) then
            part(4,npt) = vxt
            part(5,npt) = vyt
            part(6,npt) = vzt
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
     &npy,npz,idimp,npmax,ierr)
! for 3d code, this subroutine calculates initial particle velocities
! momentum with maxwell-juttner distribution with drift
! for relativistic particles and distributed data.
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
! input: all except part, ierr, output: part, npp, ierr
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
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! ierr = (0,1) = (no,yes) error condition exists
! ranorm = gaussian random number with zero mean and unit variance
! with 2D spatial decomposition
      implicit none
      integer nps, npp, npx, npy, npz, idimp, npmax, ierr
      real vtx, vty, vtz, vdx, vdy, vdz, ci
      real part
      dimension part(idimp,npmax)
! local data
      integer j, k, l, npt, npxyzp
      real ci4, ptx, pty, ptz, pt2
      double precision dps, dnpx, dnpxy, dkoff, djoff, dj, dnpxyz, dt1
      integer ierr1, iwork1
      dimension ierr1(1), iwork1(1)
      double precision dpp, dpt, sum4, work4
      dimension dpp(1), dpt(1), sum4(4), work4(4)
      double precision ranorm
      ci4 = 0.25*ci*ci
      ierr = 0
! determine offsets for particle number
! each random triplet is associated with a global particle index
! 1:sum(npp-nps+1)
      dps = dble(npp-nps+1)
      dnpx = dble(npx)
      dnpxy = dble(npy)*dnpx
      dpp(1) = dps
      call PPDSCAN(dpp,dpt,1)
      dps = dpp(1) - dps
! maxwell-juttner momentum distribution
      npt = nps - 1
      do 30 l = 1, npz
      dkoff = dnpxy*dble(l - 1)
      do 20 k = 1, npy
      djoff = dnpx*dble(k - 1) + dkoff
      do 10 j = 1, npx
      dj = dble(j - 1) + djoff
      ptx = vtx*ranorm()
      pty = vty*ranorm()
      ptz = vtz*ranorm()
      if ((dj.ge.dps).and.(dj.lt.dpp(1))) then
         pt2 = ptx*ptx + pty*pty + ptz*ptz
         pt2 = sqrt(1.0 + pt2*ci4)
         npt = npt + 1
         if (npt.le.npmax) then
            part(4,npt) = ptx*pt2
            part(5,npt) = pty*pt2
            part(6,npt) = ptz*pt2
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
