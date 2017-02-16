!-----------------------------------------------------------------------
! Fortran Library for depositing current density
! 3D MPI/OpenMP PIC Codes:
! PPGJPPOST32L calculates particle current density using linear
!              interpolation and advances particle positions half a
!              time-step with various particle boundary conditions
! PPGJPPOSTF32L calculates particle current density using linear
!               interpolation, advances particle positions half a
!               time-step with periodic boundary conditions, and
!               determines list of particles which are leaving each tile
! PPGRJPPOST32L calculates particle current density using linear
!               interpolation for relativistic particles and advances
!               particle positions half a time-step with with various
!               particle boundary conditions
! PPGRJPPOSTF32L calculates particle current density using linear
!                interpolation for relativistic particles, advances
!                particle positions half a time-step with periodic
!                boundary conditions, and determines list of particles
!                which are leaving each tile
! PPGMJPPOST32L calculates particle momentum flux using linear
!               interpolation
! PPGRMJPPOST32L calculates relativistic particle momentum flux using
!                linear interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 16, 2016
!-----------------------------------------------------------------------
      subroutine PPGJPPOST32L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny&
     &,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation
! in addition, particle positions are advanced a half time-step
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 69 flops/particle, 30 loads, 27 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qm, dt
      real ppart, cu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z
      real scu
      dimension scu(3,MXV,MYV,MZV)
!     dimension scu(3,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgelz = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      edgerz = real(nz-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,scu)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + vx*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + vy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + vx*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + vy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + vx*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + vy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + vx*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + vy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + vx*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + vy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + vx*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + vy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + vx*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + vy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + vx*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + vy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 40
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   40 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
   50 continue
   60 continue
   70 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,j+moffp,1+loffp) = cu(1,i+noffp,j+moffp,1+loffp)     &
     &+ scu(1,i,j,1)
!$OMP ATOMIC
      cu(2,i+noffp,j+moffp,1+loffp) = cu(2,i+noffp,j+moffp,1+loffp)     &
     &+ scu(2,i,j,1)
!$OMP ATOMIC
      cu(3,i+noffp,j+moffp,1+loffp) = cu(3,i+noffp,j+moffp,1+loffp)     &
     &+ scu(3,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp,j+moffp,lm+loffp) = cu(1,i+noffp,j+moffp,lm+loffp)&
     &   + scu(1,i,j,lm)
!$OMP ATOMIC
         cu(2,i+noffp,j+moffp,lm+loffp) = cu(2,i+noffp,j+moffp,lm+loffp)&
     &   + scu(2,i,j,lm)
!$OMP ATOMIC
         cu(3,i+noffp,j+moffp,lm+loffp) = cu(3,i+noffp,j+moffp,lm+loffp)&
     &   + scu(3,i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp,k+loffp) = cu(1,i+noffp,1+moffp,k+loffp)     &
     &+ scu(1,i,1,k)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp,k+loffp) = cu(2,i+noffp,1+moffp,k+loffp)     &
     &+ scu(2,i,1,k)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp,k+loffp) = cu(3,i+noffp,1+moffp,k+loffp)     &
     &+ scu(3,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp,k+loffp) = cu(1,i+noffp,mm+moffp,k+loffp)&
     &   + scu(1,i,mm,k)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp,k+loffp) = cu(2,i+noffp,mm+moffp,k+loffp)&
     &   + scu(2,i,mm,k)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp,k+loffp) = cu(3,i+noffp,mm+moffp,k+loffp)&
     &   + scu(3,i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp,k+loffp) = cu(1,1+noffp,j+moffp,k+loffp)     &
     &+ scu(1,1,j,k)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp,k+loffp) = cu(2,1+noffp,j+moffp,k+loffp)     &
     &+ scu(2,1,j,k)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp,k+loffp) = cu(3,1+noffp,j+moffp,k+loffp)     &
     &+ scu(3,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp,j+moffp,k+loffp) = cu(1,nm+noffp,j+moffp,k+loffp)&
     &   + scu(1,nm,j,k)
!$OMP ATOMIC
         cu(2,nm+noffp,j+moffp,k+loffp) = cu(2,nm+noffp,j+moffp,k+loffp)&
     &   + scu(2,nm,j,k)
!$OMP ATOMIC
         cu(3,nm+noffp,j+moffp,k+loffp) = cu(3,nm+noffp,j+moffp,k+loffp)&
     &   + scu(3,nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp,1+moffp,lm+loffp) = cu(1,i+noffp,1+moffp,lm+loffp)&
     &   + scu(1,i,1,lm)
!$OMP ATOMIC
         cu(2,i+noffp,1+moffp,lm+loffp) = cu(2,i+noffp,1+moffp,lm+loffp)&
     &   + scu(2,i,1,lm)
!$OMP ATOMIC
         cu(3,i+noffp,1+moffp,lm+loffp) = cu(3,i+noffp,1+moffp,lm+loffp)&
     &   + scu(3,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(1,i+noffp,mm+moffp,lm+loffp) + scu(1,i,mm,lm)
!$OMP ATOMIC
            cu(2,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(2,i+noffp,mm+moffp,lm+loffp) + scu(2,i,mm,lm)
!$OMP ATOMIC
            cu(3,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(3,i+noffp,mm+moffp,lm+loffp) + scu(3,i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp,j+moffp,lm+loffp) = cu(1,1+noffp,j+moffp,lm+loffp)&
     &   + scu(1,1,j,lm)
!$OMP ATOMIC
         cu(2,1+noffp,j+moffp,lm+loffp) = cu(2,1+noffp,j+moffp,lm+loffp)&
     &   + scu(2,1,j,lm)
!$OMP ATOMIC
         cu(3,1+noffp,j+moffp,lm+loffp) = cu(3,1+noffp,j+moffp,lm+loffp)&
     &   + scu(3,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(1,nm+noffp,j+moffp,lm+loffp) + scu(1,nm,j,lm)
!$OMP ATOMIC
            cu(2,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(2,nm+noffp,j+moffp,lm+loffp) + scu(2,nm,j,lm)
!$OMP ATOMIC
            cu(3,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(3,nm+noffp,j+moffp,lm+loffp) + scu(3,nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt, &
     &nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation.
! in addition, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 69 flops/particle, 30 loads, 27 stores
! input: all except ncl, ihole, irc, output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ntmax = size of hole array for particles leaving tiles
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, ih, nh, nn, mm, ll, nm, lm
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real scu
      dimension scu(3,MXV,MYV,MZV)
!     dimension scu(3,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,ih,nh
!$OMP& ,nm,lm,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx
!$OMP& ,edgely,edgelz,edgerx,edgery,edgerz,scu)
      do 160 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
      ih = 0
      nh = 0
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + vx*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + vy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + vx*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + vy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + vx*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + vy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + vx*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + vy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + vx*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + vy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + vx*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + vy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + vx*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + vy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + vx*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + vy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 50
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
   50 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 80 k = 2, ll
      do 70 j = 2, mm
      do 60 i = 2, nn
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
   60 continue
   70 continue
   80 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 100 j = 2, mm
      do 90 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,j+moffp,1+loffp) = cu(1,i+noffp,j+moffp,1+loffp)     &
     &+ scu(1,i,j,1)
!$OMP ATOMIC
      cu(2,i+noffp,j+moffp,1+loffp) = cu(2,i+noffp,j+moffp,1+loffp)     &
     &+ scu(2,i,j,1)
!$OMP ATOMIC
      cu(3,i+noffp,j+moffp,1+loffp) = cu(3,i+noffp,j+moffp,1+loffp)     &
     &+ scu(3,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp,j+moffp,lm+loffp) = cu(1,i+noffp,j+moffp,lm+loffp)&
     &   + scu(1,i,j,lm)
!$OMP ATOMIC
         cu(2,i+noffp,j+moffp,lm+loffp) = cu(2,i+noffp,j+moffp,lm+loffp)&
     &   + scu(2,i,j,lm)
!$OMP ATOMIC
         cu(3,i+noffp,j+moffp,lm+loffp) = cu(3,i+noffp,j+moffp,lm+loffp)&
     &   + scu(3,i,j,lm)
      endif
   90 continue
  100 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 130 k = 1, ll
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp,k+loffp) = cu(1,i+noffp,1+moffp,k+loffp)     &
     &+ scu(1,i,1,k)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp,k+loffp) = cu(2,i+noffp,1+moffp,k+loffp)     &
     &+ scu(2,i,1,k)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp,k+loffp) = cu(3,i+noffp,1+moffp,k+loffp)     &
     &+ scu(3,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp,k+loffp) = cu(1,i+noffp,mm+moffp,k+loffp)&
     &   + scu(1,i,mm,k)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp,k+loffp) = cu(2,i+noffp,mm+moffp,k+loffp)&
     &   + scu(2,i,mm,k)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp,k+loffp) = cu(3,i+noffp,mm+moffp,k+loffp)&
     &   + scu(3,i,mm,k)
      endif
  110 continue
      do 120 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp,k+loffp) = cu(1,1+noffp,j+moffp,k+loffp)     &
     &+ scu(1,1,j,k)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp,k+loffp) = cu(2,1+noffp,j+moffp,k+loffp)     &
     &+ scu(2,1,j,k)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp,k+loffp) = cu(3,1+noffp,j+moffp,k+loffp)     &
     &+ scu(3,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp,j+moffp,k+loffp) = cu(1,nm+noffp,j+moffp,k+loffp)&
     &   + scu(1,nm,j,k)
!$OMP ATOMIC
         cu(2,nm+noffp,j+moffp,k+loffp) = cu(2,nm+noffp,j+moffp,k+loffp)&
     &   + scu(2,nm,j,k)
!$OMP ATOMIC
         cu(3,nm+noffp,j+moffp,k+loffp) = cu(3,nm+noffp,j+moffp,k+loffp)&
     &   + scu(3,nm,j,k)
      endif
  120 continue
  130 continue
      if (lm > mz) then
         do 140 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp,1+moffp,lm+loffp) = cu(1,i+noffp,1+moffp,lm+loffp)&
     &   + scu(1,i,1,lm)
!$OMP ATOMIC
         cu(2,i+noffp,1+moffp,lm+loffp) = cu(2,i+noffp,1+moffp,lm+loffp)&
     &   + scu(2,i,1,lm)
!$OMP ATOMIC
         cu(3,i+noffp,1+moffp,lm+loffp) = cu(3,i+noffp,1+moffp,lm+loffp)&
     &   + scu(3,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(1,i+noffp,mm+moffp,lm+loffp) + scu(1,i,mm,lm)
!$OMP ATOMIC
            cu(2,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(2,i+noffp,mm+moffp,lm+loffp) + scu(2,i,mm,lm)
!$OMP ATOMIC
            cu(3,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(3,i+noffp,mm+moffp,lm+loffp) + scu(3,i,mm,lm)
         endif
  140    continue
         do 150 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp,j+moffp,lm+loffp) = cu(1,1+noffp,j+moffp,lm+loffp)&
     &   + scu(1,1,j,lm)
!$OMP ATOMIC
         cu(2,1+noffp,j+moffp,lm+loffp) = cu(2,1+noffp,j+moffp,lm+loffp)&
     &   + scu(2,1,j,lm)
!$OMP ATOMIC
         cu(3,1+noffp,j+moffp,lm+loffp) = cu(3,1+noffp,j+moffp,lm+loffp)&
     &   + scu(3,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(1,nm+noffp,j+moffp,lm+loffp) + scu(1,nm,j,lm)
!$OMP ATOMIC
            cu(2,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(2,nm+noffp,j+moffp,lm+loffp) + scu(2,nm,j,lm)
!$OMP ATOMIC
            cu(3,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(3,nm+noffp,j+moffp,lm+loffp) + scu(3,nm,j,lm)
         endif
  150    continue
      endif
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
  160 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRJPPOST32L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp, &
     &nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! in addition, particle positions are advanced a half time-step
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = x momentum of particle n in partition in tile m
! ppart(5,n,m) = y momentum of particle n in partition in tile m
! ppart(6,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z, p2, gami
      real scu
      dimension scu(3,MXV,MYV,MZV)
!     dimension scu(3,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      ci2 = ci*ci
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgelz = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      edgerz = real(nz-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,p2,gami,scu)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + vx*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + vy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + vx*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + vy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + vx*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + vy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + vx*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + vy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + vx*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + vy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + vx*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + vy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + vx*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + vy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + vx*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + vy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 40
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -ppart(6,j,l)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   40 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
   50 continue
   60 continue
   70 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,j+moffp,1+loffp) = cu(1,i+noffp,j+moffp,1+loffp)     &
     &+ scu(1,i,j,1)
!$OMP ATOMIC
      cu(2,i+noffp,j+moffp,1+loffp) = cu(2,i+noffp,j+moffp,1+loffp)     &
     &+ scu(2,i,j,1)
!$OMP ATOMIC
      cu(3,i+noffp,j+moffp,1+loffp) = cu(3,i+noffp,j+moffp,1+loffp)     &
     &+ scu(3,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp,j+moffp,lm+loffp) = cu(1,i+noffp,j+moffp,lm+loffp)&
     &   + scu(1,i,j,lm)
!$OMP ATOMIC
         cu(2,i+noffp,j+moffp,lm+loffp) = cu(2,i+noffp,j+moffp,lm+loffp)&
     &   + scu(2,i,j,lm)
!$OMP ATOMIC
         cu(3,i+noffp,j+moffp,lm+loffp) = cu(3,i+noffp,j+moffp,lm+loffp)&
     &   + scu(3,i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp,k+loffp) = cu(1,i+noffp,1+moffp,k+loffp)     &
     &+ scu(1,i,1,k)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp,k+loffp) = cu(2,i+noffp,1+moffp,k+loffp)     &
     &+ scu(2,i,1,k)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp,k+loffp) = cu(3,i+noffp,1+moffp,k+loffp)     &
     &+ scu(3,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp,k+loffp) = cu(1,i+noffp,mm+moffp,k+loffp)&
     &   + scu(1,i,mm,k)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp,k+loffp) = cu(2,i+noffp,mm+moffp,k+loffp)&
     &   + scu(2,i,mm,k)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp,k+loffp) = cu(3,i+noffp,mm+moffp,k+loffp)&
     &   + scu(3,i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp,k+loffp) = cu(1,1+noffp,j+moffp,k+loffp)     &
     &+ scu(1,1,j,k)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp,k+loffp) = cu(2,1+noffp,j+moffp,k+loffp)     &
     &+ scu(2,1,j,k)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp,k+loffp) = cu(3,1+noffp,j+moffp,k+loffp)     &
     &+ scu(3,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp,j+moffp,k+loffp) = cu(1,nm+noffp,j+moffp,k+loffp)&
     &   + scu(1,nm,j,k)
!$OMP ATOMIC
         cu(2,nm+noffp,j+moffp,k+loffp) = cu(2,nm+noffp,j+moffp,k+loffp)&
     &   + scu(2,nm,j,k)
!$OMP ATOMIC
         cu(3,nm+noffp,j+moffp,k+loffp) = cu(3,nm+noffp,j+moffp,k+loffp)&
     &   + scu(3,nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp,1+moffp,lm+loffp) = cu(1,i+noffp,1+moffp,lm+loffp)&
     &   + scu(1,i,1,lm)
!$OMP ATOMIC
         cu(2,i+noffp,1+moffp,lm+loffp) = cu(2,i+noffp,1+moffp,lm+loffp)&
     &   + scu(2,i,1,lm)
!$OMP ATOMIC
         cu(3,i+noffp,1+moffp,lm+loffp) = cu(3,i+noffp,1+moffp,lm+loffp)&
     &   + scu(3,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(1,i+noffp,mm+moffp,lm+loffp) + scu(1,i,mm,lm)
!$OMP ATOMIC
            cu(2,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(2,i+noffp,mm+moffp,lm+loffp) + scu(2,i,mm,lm)
!$OMP ATOMIC
            cu(3,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(3,i+noffp,mm+moffp,lm+loffp) + scu(3,i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp,j+moffp,lm+loffp) = cu(1,1+noffp,j+moffp,lm+loffp)&
     &   + scu(1,1,j,lm)
!$OMP ATOMIC
         cu(2,1+noffp,j+moffp,lm+loffp) = cu(2,1+noffp,j+moffp,lm+loffp)&
     &   + scu(2,1,j,lm)
!$OMP ATOMIC
         cu(3,1+noffp,j+moffp,lm+loffp) = cu(3,1+noffp,j+moffp,lm+loffp)&
     &   + scu(3,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(1,nm+noffp,j+moffp,lm+loffp) + scu(1,nm,j,lm)
!$OMP ATOMIC
            cu(2,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(2,nm+noffp,j+moffp,lm+loffp) + scu(2,nm,j,lm)
!$OMP ATOMIC
            cu(3,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(3,nm+noffp,j+moffp,lm+loffp) + scu(3,nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,&
     &ci,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1, &
     &ntmax,idds,irc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! in addition, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
! input: all except ncl, ihole, irc, output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = x momentum of particle n in partition in tile m
! ppart(5,n,m) = y momentum of particle n in partition in tile m
! ppart(6,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ntmax = size of hole array for particles leaving tiles
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, ih, nh, nn, mm, ll, nm, lm
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z, ci2, p2, gami
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real scu
      dimension scu(3,MXV,MYV,MZV)
!     dimension scu(3,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,ih,nh
!$OMP& ,nm,lm,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx
!$OMP& ,edgely,edgelz,edgerx,edgery,edgerz,p2,gami,scu)
      do 160 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
      ih = 0
      nh = 0
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + vx*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + vy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + vx*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + vy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + vx*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + vy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + vx*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + vy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + vx*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + vy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + vx*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + vy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + vx*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + vy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + vx*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + vy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 50
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
   50 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 80 k = 2, ll
      do 70 j = 2, mm
      do 60 i = 2, nn
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
   60 continue
   70 continue
   80 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 100 j = 2, mm
      do 90 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,j+moffp,1+loffp) = cu(1,i+noffp,j+moffp,1+loffp)     &
     &+ scu(1,i,j,1)
!$OMP ATOMIC
      cu(2,i+noffp,j+moffp,1+loffp) = cu(2,i+noffp,j+moffp,1+loffp)     &
     &+ scu(2,i,j,1)
!$OMP ATOMIC
      cu(3,i+noffp,j+moffp,1+loffp) = cu(3,i+noffp,j+moffp,1+loffp)     &
     &+ scu(3,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp,j+moffp,lm+loffp) = cu(1,i+noffp,j+moffp,lm+loffp)&
     &   + scu(1,i,j,lm)
!$OMP ATOMIC
         cu(2,i+noffp,j+moffp,lm+loffp) = cu(2,i+noffp,j+moffp,lm+loffp)&
     &   + scu(2,i,j,lm)
!$OMP ATOMIC
         cu(3,i+noffp,j+moffp,lm+loffp) = cu(3,i+noffp,j+moffp,lm+loffp)&
     &   + scu(3,i,j,lm)
      endif
   90 continue
  100 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 130 k = 1, ll
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp,k+loffp) = cu(1,i+noffp,1+moffp,k+loffp)     &
     &+ scu(1,i,1,k)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp,k+loffp) = cu(2,i+noffp,1+moffp,k+loffp)     &
     &+ scu(2,i,1,k)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp,k+loffp) = cu(3,i+noffp,1+moffp,k+loffp)     &
     &+ scu(3,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp,k+loffp) = cu(1,i+noffp,mm+moffp,k+loffp)&
     &   + scu(1,i,mm,k)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp,k+loffp) = cu(2,i+noffp,mm+moffp,k+loffp)&
     &   + scu(2,i,mm,k)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp,k+loffp) = cu(3,i+noffp,mm+moffp,k+loffp)&
     &   + scu(3,i,mm,k)
      endif
  110 continue
      do 120 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp,k+loffp) = cu(1,1+noffp,j+moffp,k+loffp)     &
     &+ scu(1,1,j,k)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp,k+loffp) = cu(2,1+noffp,j+moffp,k+loffp)     &
     &+ scu(2,1,j,k)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp,k+loffp) = cu(3,1+noffp,j+moffp,k+loffp)     &
     &+ scu(3,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp,j+moffp,k+loffp) = cu(1,nm+noffp,j+moffp,k+loffp)&
     &   + scu(1,nm,j,k)
!$OMP ATOMIC
         cu(2,nm+noffp,j+moffp,k+loffp) = cu(2,nm+noffp,j+moffp,k+loffp)&
     &   + scu(2,nm,j,k)
!$OMP ATOMIC
         cu(3,nm+noffp,j+moffp,k+loffp) = cu(3,nm+noffp,j+moffp,k+loffp)&
     &   + scu(3,nm,j,k)
      endif
  120 continue
  130 continue
      if (lm > mz) then
         do 140 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp,1+moffp,lm+loffp) = cu(1,i+noffp,1+moffp,lm+loffp)&
     &   + scu(1,i,1,lm)
!$OMP ATOMIC
         cu(2,i+noffp,1+moffp,lm+loffp) = cu(2,i+noffp,1+moffp,lm+loffp)&
     &   + scu(2,i,1,lm)
!$OMP ATOMIC
         cu(3,i+noffp,1+moffp,lm+loffp) = cu(3,i+noffp,1+moffp,lm+loffp)&
     &   + scu(3,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(1,i+noffp,mm+moffp,lm+loffp) + scu(1,i,mm,lm)
!$OMP ATOMIC
            cu(2,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(2,i+noffp,mm+moffp,lm+loffp) + scu(2,i,mm,lm)
!$OMP ATOMIC
            cu(3,i+noffp,mm+moffp,lm+loffp) =                           &
     &      cu(3,i+noffp,mm+moffp,lm+loffp) + scu(3,i,mm,lm)
         endif
  140    continue
         do 150 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp,j+moffp,lm+loffp) = cu(1,1+noffp,j+moffp,lm+loffp)&
     &   + scu(1,1,j,lm)
!$OMP ATOMIC
         cu(2,1+noffp,j+moffp,lm+loffp) = cu(2,1+noffp,j+moffp,lm+loffp)&
     &   + scu(2,1,j,lm)
!$OMP ATOMIC
         cu(3,1+noffp,j+moffp,lm+loffp) = cu(3,1+noffp,j+moffp,lm+loffp)&
     &   + scu(3,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(1,nm+noffp,j+moffp,lm+loffp) + scu(1,nm,j,lm)
!$OMP ATOMIC
            cu(2,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(2,nm+noffp,j+moffp,lm+loffp) + scu(2,nm,j,lm)
!$OMP ATOMIC
            cu(3,nm+noffp,j+moffp,lm+loffp) =                           &
     &      cu(3,nm+noffp,j+moffp,lm+loffp) + scu(3,nm,j,lm)
         endif
  150    continue
      endif
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
  160 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGMJPPOST32L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,&
     &mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle momentum flux
! using first-order linear interpolation
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 121 flops/particle, 52 loads, 48 stores
! input: all, output: ppart, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where vj = vj(t-dt/2) and vk = vk(t-dt/2)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm
      real ppart, amu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, vx, vy, vz
      real x, y, z, v1, v2, v3, v4, v5, v6
      real samu
      dimension samu(6,MXV,MYV,MZV)
!     dimension samu(6,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,vx,vy,vz,v1,v2,v3,v4,v5,
!$OMP& v6,samu)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit momentum flux
      dx = amx*amz
      dy = amy*amz
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
   40 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      amu(1,i+noffp,j+moffp,k+loffp) = amu(1,i+noffp,j+moffp,k+loffp)   &
     &+ samu(1,i,j,k)
      amu(2,i+noffp,j+moffp,k+loffp) = amu(2,i+noffp,j+moffp,k+loffp)   &
     &+ samu(2,i,j,k)
      amu(3,i+noffp,j+moffp,k+loffp) = amu(3,i+noffp,j+moffp,k+loffp)   &
     &+ samu(3,i,j,k)
      amu(4,i+noffp,j+moffp,k+loffp) = amu(4,i+noffp,j+moffp,k+loffp)   &
     &+ samu(4,i,j,k)
      amu(5,i+noffp,j+moffp,k+loffp) = amu(5,i+noffp,j+moffp,k+loffp)   &
     &+ samu(5,i,j,k)
      amu(6,i+noffp,j+moffp,k+loffp) = amu(6,i+noffp,j+moffp,k+loffp)   &
     &+ samu(6,i,j,k)
   50 continue
   60 continue
   70 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,j+moffp,1+loffp) = amu(1,i+noffp,j+moffp,1+loffp)   &
     &+ samu(1,i,j,1)
!$OMP ATOMIC
      amu(2,i+noffp,j+moffp,1+loffp) = amu(2,i+noffp,j+moffp,1+loffp)   &
     &+ samu(2,i,j,1)
!$OMP ATOMIC
      amu(3,i+noffp,j+moffp,1+loffp) = amu(3,i+noffp,j+moffp,1+loffp)   &
     &+ samu(3,i,j,1)
!$OMP ATOMIC
      amu(4,i+noffp,j+moffp,1+loffp) = amu(4,i+noffp,j+moffp,1+loffp)   &
     &+ samu(4,i,j,1)
!$OMP ATOMIC
      amu(5,i+noffp,j+moffp,1+loffp) = amu(5,i+noffp,j+moffp,1+loffp)   &
     &+ samu(5,i,j,1)
!$OMP ATOMIC
      amu(6,i+noffp,j+moffp,1+loffp) = amu(6,i+noffp,j+moffp,1+loffp)   &
     &+ samu(6,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(1,i+noffp,j+moffp,lm+loffp) + samu(1,i,j,lm)
!$OMP ATOMIC
         amu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(2,i+noffp,j+moffp,lm+loffp) + samu(2,i,j,lm)
!$OMP ATOMIC
         amu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(3,i+noffp,j+moffp,lm+loffp) + samu(3,i,j,lm)
!$OMP ATOMIC
         amu(4,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(4,i+noffp,j+moffp,lm+loffp) + samu(4,i,j,lm)
!$OMP ATOMIC
         amu(5,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(5,i+noffp,j+moffp,lm+loffp) + samu(5,i,j,lm)
!$OMP ATOMIC
         amu(6,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(6,i+noffp,j+moffp,lm+loffp) + samu(6,i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp,k+loffp) = amu(1,i+noffp,1+moffp,k+loffp)   &
     &+ samu(1,i,1,k)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp,k+loffp) = amu(2,i+noffp,1+moffp,k+loffp)   &
     &+ samu(2,i,1,k)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp,k+loffp) = amu(3,i+noffp,1+moffp,k+loffp)   &
     &+ samu(3,i,1,k)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp,k+loffp) = amu(4,i+noffp,1+moffp,k+loffp)   &
     &+ samu(4,i,1,k)
!$OMP ATOMIC
      amu(5,i+noffp,1+moffp,k+loffp) = amu(5,i+noffp,1+moffp,k+loffp)   &
     &+ samu(5,i,1,k)
!$OMP ATOMIC
      amu(6,i+noffp,1+moffp,k+loffp) = amu(6,i+noffp,1+moffp,k+loffp)   &
     &+ samu(6,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(1,i+noffp,mm+moffp,k+loffp) + samu(1,i,mm,k)
!$OMP ATOMIC
         amu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(2,i+noffp,mm+moffp,k+loffp) + samu(2,i,mm,k)
!$OMP ATOMIC
         amu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(3,i+noffp,mm+moffp,k+loffp) + samu(3,i,mm,k)
!$OMP ATOMIC
         amu(4,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(4,i+noffp,mm+moffp,k+loffp) + samu(4,i,mm,k)
!$OMP ATOMIC
         amu(5,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(5,i+noffp,mm+moffp,k+loffp) + samu(5,i,mm,k)
!$OMP ATOMIC
         amu(6,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(6,i+noffp,mm+moffp,k+loffp) + samu(6,i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp,k+loffp) = amu(1,1+noffp,j+moffp,k+loffp)   &
     &+ samu(1,1,j,k)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp,k+loffp) = amu(2,1+noffp,j+moffp,k+loffp)   &
     &+ samu(2,1,j,k)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp,k+loffp) = amu(3,1+noffp,j+moffp,k+loffp)   &
     &+ samu(3,1,j,k)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp,k+loffp) = amu(4,1+noffp,j+moffp,k+loffp)   &
     &+ samu(4,1,j,k)
!$OMP ATOMIC
      amu(5,1+noffp,j+moffp,k+loffp) = amu(5,1+noffp,j+moffp,k+loffp)   &
     &+ samu(5,1,j,k)
!$OMP ATOMIC
      amu(6,1+noffp,j+moffp,k+loffp) = amu(6,1+noffp,j+moffp,k+loffp)   &
     &+ samu(6,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(1,nm+noffp,j+moffp,k+loffp) + samu(1,nm,j,k)
!$OMP ATOMIC
         amu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(2,nm+noffp,j+moffp,k+loffp) + samu(2,nm,j,k)
!$OMP ATOMIC
         amu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(3,nm+noffp,j+moffp,k+loffp) + samu(3,nm,j,k)
!$OMP ATOMIC
         amu(4,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(4,nm+noffp,j+moffp,k+loffp) + samu(4,nm,j,k)
!$OMP ATOMIC
         amu(5,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(5,nm+noffp,j+moffp,k+loffp) + samu(5,nm,j,k)
!$OMP ATOMIC
         amu(6,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(6,nm+noffp,j+moffp,k+loffp) + samu(6,nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(1,i+noffp,1+moffp,lm+loffp) + samu(1,i,1,lm)
!$OMP ATOMIC
         amu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(2,i+noffp,1+moffp,lm+loffp) + samu(2,i,1,lm)
!$OMP ATOMIC
         amu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(3,i+noffp,1+moffp,lm+loffp) + samu(3,i,1,lm)
!$OMP ATOMIC
         amu(4,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(4,i+noffp,1+moffp,lm+loffp) + samu(4,i,1,lm)
!$OMP ATOMIC
         amu(5,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(5,i+noffp,1+moffp,lm+loffp) + samu(5,i,1,lm)
!$OMP ATOMIC
         amu(6,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(6,i+noffp,1+moffp,lm+loffp) + samu(6,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(1,i+noffp,mm+moffp,lm+loffp) + samu(1,i,mm,lm)
!$OMP ATOMIC
            amu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(2,i+noffp,mm+moffp,lm+loffp) + samu(2,i,mm,lm)
!$OMP ATOMIC
            amu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(3,i+noffp,mm+moffp,lm+loffp) + samu(3,i,mm,lm)
!$OMP ATOMIC
            amu(4,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(4,i+noffp,mm+moffp,lm+loffp) + samu(4,i,mm,lm)
!$OMP ATOMIC
            amu(5,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(5,i+noffp,mm+moffp,lm+loffp) + samu(5,i,mm,lm)
!$OMP ATOMIC
            amu(6,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(6,i+noffp,mm+moffp,lm+loffp) + samu(6,i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(1,1+noffp,j+moffp,lm+loffp) + samu(1,1,j,lm)
!$OMP ATOMIC
         amu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(2,1+noffp,j+moffp,lm+loffp) + samu(2,1,j,lm)
!$OMP ATOMIC
         amu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(3,1+noffp,j+moffp,lm+loffp) + samu(3,1,j,lm)
!$OMP ATOMIC
         amu(4,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(4,1+noffp,j+moffp,lm+loffp) + samu(4,1,j,lm)
!$OMP ATOMIC
         amu(5,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(5,1+noffp,j+moffp,lm+loffp) + samu(5,1,j,lm)
!$OMP ATOMIC
         amu(6,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(6,1+noffp,j+moffp,lm+loffp) + samu(6,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(1,nm+noffp,j+moffp,lm+loffp) + samu(1,nm,j,lm)
!$OMP ATOMIC
            amu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(2,nm+noffp,j+moffp,lm+loffp) + samu(2,nm,j,lm)
!$OMP ATOMIC
            amu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(3,nm+noffp,j+moffp,lm+loffp) + samu(3,nm,j,lm)
!$OMP ATOMIC
            amu(4,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(4,nm+noffp,j+moffp,lm+loffp) + samu(4,nm,j,lm)
!$OMP ATOMIC
            amu(5,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(5,nm+noffp,j+moffp,lm+loffp) + samu(5,nm,j,lm)
!$OMP ATOMIC
            amu(6,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(6,nm+noffp,j+moffp,lm+loffp) + samu(6,nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRMJPPOST32L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx&
     &,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle momentum flux
! using first-order linear interpolation for relativistic particles
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 134 flops/particle, 1 divide, 52 loads, 48 stores
! input: all, output: ppart, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pj*pk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where pj = pj(t-dt/2) and pk = pk(t-dt/2)
! where gami2 = 1./(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum vz of particle n in partition in tile m
! at t - dt/2
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm, ci
      real ppart, amu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real ci2, gami2, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy
      real vx, vy, vz, x, y, z, p2, v1, v2, v3, v4, v5, v6
      real samu
      dimension samu(6,MXV,MYV,MZV)
!     dimension samu(6,mx+1,my+1,mz+1)
      ci2 = ci*ci
      mxyp1 = mx1*myp1
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,vx,vy,vz,v1,v2,v3,v4,v5,
!$OMP& v6,p2,gami2,samu)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit momentum flux
      dx = amx*amz
      dy = amy*amz
      v1 = (vx*vx)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vx*vz)*gami2
      v4 = (vy*vy)*gami2
      v5 = (vy*vz)*gami2
      v6 = (vz*vz)*gami2
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
   40 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      amu(1,i+noffp,j+moffp,k+loffp) = amu(1,i+noffp,j+moffp,k+loffp)   &
     &+ samu(1,i,j,k)
      amu(2,i+noffp,j+moffp,k+loffp) = amu(2,i+noffp,j+moffp,k+loffp)   &
     &+ samu(2,i,j,k)
      amu(3,i+noffp,j+moffp,k+loffp) = amu(3,i+noffp,j+moffp,k+loffp)   &
     &+ samu(3,i,j,k)
      amu(4,i+noffp,j+moffp,k+loffp) = amu(4,i+noffp,j+moffp,k+loffp)   &
     &+ samu(4,i,j,k)
      amu(5,i+noffp,j+moffp,k+loffp) = amu(5,i+noffp,j+moffp,k+loffp)   &
     &+ samu(5,i,j,k)
      amu(6,i+noffp,j+moffp,k+loffp) = amu(6,i+noffp,j+moffp,k+loffp)   &
     &+ samu(6,i,j,k)
   50 continue
   60 continue
   70 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,j+moffp,1+loffp) = amu(1,i+noffp,j+moffp,1+loffp)   &
     &+ samu(1,i,j,1)
!$OMP ATOMIC
      amu(2,i+noffp,j+moffp,1+loffp) = amu(2,i+noffp,j+moffp,1+loffp)   &
     &+ samu(2,i,j,1)
!$OMP ATOMIC
      amu(3,i+noffp,j+moffp,1+loffp) = amu(3,i+noffp,j+moffp,1+loffp)   &
     &+ samu(3,i,j,1)
!$OMP ATOMIC
      amu(4,i+noffp,j+moffp,1+loffp) = amu(4,i+noffp,j+moffp,1+loffp)   &
     &+ samu(4,i,j,1)
!$OMP ATOMIC
      amu(5,i+noffp,j+moffp,1+loffp) = amu(5,i+noffp,j+moffp,1+loffp)   &
     &+ samu(5,i,j,1)
!$OMP ATOMIC
      amu(6,i+noffp,j+moffp,1+loffp) = amu(6,i+noffp,j+moffp,1+loffp)   &
     &+ samu(6,i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(1,i+noffp,j+moffp,lm+loffp) + samu(1,i,j,lm)
!$OMP ATOMIC
         amu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(2,i+noffp,j+moffp,lm+loffp) + samu(2,i,j,lm)
!$OMP ATOMIC
         amu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(3,i+noffp,j+moffp,lm+loffp) + samu(3,i,j,lm)
!$OMP ATOMIC
         amu(4,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(4,i+noffp,j+moffp,lm+loffp) + samu(4,i,j,lm)
!$OMP ATOMIC
         amu(5,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(5,i+noffp,j+moffp,lm+loffp) + samu(5,i,j,lm)
!$OMP ATOMIC
         amu(6,i+noffp,j+moffp,lm+loffp) =                              &
     &   amu(6,i+noffp,j+moffp,lm+loffp) + samu(6,i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp,k+loffp) = amu(1,i+noffp,1+moffp,k+loffp)   &
     &+ samu(1,i,1,k)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp,k+loffp) = amu(2,i+noffp,1+moffp,k+loffp)   &
     &+ samu(2,i,1,k)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp,k+loffp) = amu(3,i+noffp,1+moffp,k+loffp)   &
     &+ samu(3,i,1,k)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp,k+loffp) = amu(4,i+noffp,1+moffp,k+loffp)   &
     &+ samu(4,i,1,k)
!$OMP ATOMIC
      amu(5,i+noffp,1+moffp,k+loffp) = amu(5,i+noffp,1+moffp,k+loffp)   &
     &+ samu(5,i,1,k)
!$OMP ATOMIC
      amu(6,i+noffp,1+moffp,k+loffp) = amu(6,i+noffp,1+moffp,k+loffp)   &
     &+ samu(6,i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(1,i+noffp,mm+moffp,k+loffp) + samu(1,i,mm,k)
!$OMP ATOMIC
         amu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(2,i+noffp,mm+moffp,k+loffp) + samu(2,i,mm,k)
!$OMP ATOMIC
         amu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(3,i+noffp,mm+moffp,k+loffp) + samu(3,i,mm,k)
!$OMP ATOMIC
         amu(4,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(4,i+noffp,mm+moffp,k+loffp) + samu(4,i,mm,k)
!$OMP ATOMIC
         amu(5,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(5,i+noffp,mm+moffp,k+loffp) + samu(5,i,mm,k)
!$OMP ATOMIC
         amu(6,i+noffp,mm+moffp,k+loffp) =                              &
     &   amu(6,i+noffp,mm+moffp,k+loffp) + samu(6,i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp,k+loffp) = amu(1,1+noffp,j+moffp,k+loffp)   &
     &+ samu(1,1,j,k)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp,k+loffp) = amu(2,1+noffp,j+moffp,k+loffp)   &
     &+ samu(2,1,j,k)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp,k+loffp) = amu(3,1+noffp,j+moffp,k+loffp)   &
     &+ samu(3,1,j,k)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp,k+loffp) = amu(4,1+noffp,j+moffp,k+loffp)   &
     &+ samu(4,1,j,k)
!$OMP ATOMIC
      amu(5,1+noffp,j+moffp,k+loffp) = amu(5,1+noffp,j+moffp,k+loffp)   &
     &+ samu(5,1,j,k)
!$OMP ATOMIC
      amu(6,1+noffp,j+moffp,k+loffp) = amu(6,1+noffp,j+moffp,k+loffp)   &
     &+ samu(6,1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(1,nm+noffp,j+moffp,k+loffp) + samu(1,nm,j,k)
!$OMP ATOMIC
         amu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(2,nm+noffp,j+moffp,k+loffp) + samu(2,nm,j,k)
!$OMP ATOMIC
         amu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(3,nm+noffp,j+moffp,k+loffp) + samu(3,nm,j,k)
!$OMP ATOMIC
         amu(4,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(4,nm+noffp,j+moffp,k+loffp) + samu(4,nm,j,k)
!$OMP ATOMIC
         amu(5,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(5,nm+noffp,j+moffp,k+loffp) + samu(5,nm,j,k)
!$OMP ATOMIC
         amu(6,nm+noffp,j+moffp,k+loffp) =                              &
     &   amu(6,nm+noffp,j+moffp,k+loffp) + samu(6,nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(1,i+noffp,1+moffp,lm+loffp) + samu(1,i,1,lm)
!$OMP ATOMIC
         amu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(2,i+noffp,1+moffp,lm+loffp) + samu(2,i,1,lm)
!$OMP ATOMIC
         amu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(3,i+noffp,1+moffp,lm+loffp) + samu(3,i,1,lm)
!$OMP ATOMIC
         amu(4,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(4,i+noffp,1+moffp,lm+loffp) + samu(4,i,1,lm)
!$OMP ATOMIC
         amu(5,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(5,i+noffp,1+moffp,lm+loffp) + samu(5,i,1,lm)
!$OMP ATOMIC
         amu(6,i+noffp,1+moffp,lm+loffp) =                              &
     &   amu(6,i+noffp,1+moffp,lm+loffp) + samu(6,i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(1,i+noffp,mm+moffp,lm+loffp) + samu(1,i,mm,lm)
!$OMP ATOMIC
            amu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(2,i+noffp,mm+moffp,lm+loffp) + samu(2,i,mm,lm)
!$OMP ATOMIC
            amu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(3,i+noffp,mm+moffp,lm+loffp) + samu(3,i,mm,lm)
!$OMP ATOMIC
            amu(4,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(4,i+noffp,mm+moffp,lm+loffp) + samu(4,i,mm,lm)
!$OMP ATOMIC
            amu(5,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(5,i+noffp,mm+moffp,lm+loffp) + samu(5,i,mm,lm)
!$OMP ATOMIC
            amu(6,i+noffp,mm+moffp,lm+loffp) =                          &
     &      amu(6,i+noffp,mm+moffp,lm+loffp) + samu(6,i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(1,1+noffp,j+moffp,lm+loffp) + samu(1,1,j,lm)
!$OMP ATOMIC
         amu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(2,1+noffp,j+moffp,lm+loffp) + samu(2,1,j,lm)
!$OMP ATOMIC
         amu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(3,1+noffp,j+moffp,lm+loffp) + samu(3,1,j,lm)
!$OMP ATOMIC
         amu(4,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(4,1+noffp,j+moffp,lm+loffp) + samu(4,1,j,lm)
!$OMP ATOMIC
         amu(5,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(5,1+noffp,j+moffp,lm+loffp) + samu(5,1,j,lm)
!$OMP ATOMIC
         amu(6,1+noffp,j+moffp,lm+loffp) =                              &
     &   amu(6,1+noffp,j+moffp,lm+loffp) + samu(6,1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(1,nm+noffp,j+moffp,lm+loffp) + samu(1,nm,j,lm)
!$OMP ATOMIC
            amu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(2,nm+noffp,j+moffp,lm+loffp) + samu(2,nm,j,lm)
!$OMP ATOMIC
            amu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(3,nm+noffp,j+moffp,lm+loffp) + samu(3,nm,j,lm)
!$OMP ATOMIC
            amu(4,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(4,nm+noffp,j+moffp,lm+loffp) + samu(4,nm,j,lm)
!$OMP ATOMIC
            amu(5,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(5,nm+noffp,j+moffp,lm+loffp) + samu(5,nm,j,lm)
!$OMP ATOMIC
            amu(6,nm+noffp,j+moffp,lm+loffp) =                          &
     &      amu(6,nm+noffp,j+moffp,lm+loffp) + samu(6,nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
