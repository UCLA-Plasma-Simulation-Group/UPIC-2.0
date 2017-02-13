!-----------------------------------------------------------------------
! Fortran Library for depositing current density
! 2-1/2D MPI/OpenMP PIC Codes:
! PPGJPPOST2L calculates particle current density using linear
!             interpolation and advances particle positions half a
!             time-step with various particle boundary conditions
! PPGJPPOSTF2L calculates particle current density using linear
!              interpolation, advances particle positions half a
!              time-step with periodic boundary conditions, and
!              determines list of particles which are leaving each tile
! PPGRJPPOST2L calculates particle current density using linear
!              interpolation for relativistic particles and advances
!              particle positions half a time-step with with various
!              particle boundary conditions
! PPGRJPPOSTF2L calculates particle current density using linear
!               interpolation for relativistic particles, advances
!               particle positions half a time-step with periodic
!               boundary conditions, and determines list of particles
!               which are leaving each tile
! PPGMJPPOST2L calculates particle momentum flux using linear
!              interpolation
! PPGRMJPPOST2L calculates relativistic particle momentum flux using
!               linear interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 25, 2016
!-----------------------------------------------------------------------
      subroutine PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,&
     &mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 41 flops/particle, 17 loads, 14 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
      dimension scu(3,MXV,MYV)
!     dimension scu(3,mx+1,my+1)
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,dx
!$OMP& ,dy,vx,vy,vz,scu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 30
      dx = x + vx*dt
      dy = y + vy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   40 continue
   50 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,   &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 41 flops/particle, 17 loads, 14 stores
! input: all except ncl, ihole, irc,
! output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(3,MXV,MYV)
!     dimension scu(3,mx+1,my+1)
      anx = real(nx)
      any = real(ny)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 40
      dx = x + vx*dt
      dy = y + vy*dt
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
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   50 continue
   60 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 70 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   70 continue
      nn = min(mx+1,nxv-noffp)
      do 80 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   80 continue
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx&
     &,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, gami
      real scu
      dimension scu(3,MXV,MYV)
!     dimension scu(3,mx+1,my+1)
      ci2 = ci*ci
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,dx
!$OMP& ,dy,vx,vy,vz,p2,gami,scu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 30
      dx = x + vx*dt
      dy = y + vy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   40 continue
   50 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci&
     &,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
! input: all except ncl, ihole, irc,
! output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm, ih, nh
      real ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, gami
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(3,MXV,MYV)
!     dimension scu(3,mx+1,my+1)
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,p2,gami,scu)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 40
      dx = x + vx*dt
      dy = y + vy*dt
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
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   50 continue
   60 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 70 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   70 continue
      nn = min(mx+1,nxv-noffp)
      do 80 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   80 continue
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my, &
     &nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 51 flops/particle, 21 loads, 16 stores
! input: all, output: amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where vj = vj(t-dt/2) and vk = vk(t-dt/2)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! at t - dt/2
! amu(i,j,k) = ith component of momentum flux at grid point j,kk
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of flux array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, v1, v2, v3, v4
      real samu
      dimension samu(4,MXV,MYV)
!     dimension samu(4,mx+1,my+1)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,  
!$OMP& dx,dy,vx,vy,vz,v1,v2,v3,v4,samu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit momentum flux
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
   30 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
   40 continue
   50 continue
! deposit momentum flux to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp,mm+moffp) = amu(1,i+noffp,mm+moffp)              &
     & + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noffp,mm+moffp) = amu(2,i+noffp,mm+moffp)              &
     & + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noffp,mm+moffp) = amu(3,i+noffp,mm+moffp)              &
     & + samu(3,i,mm)
!$OMP ATOMIC
         amu(4,i+noffp,mm+moffp) = amu(4,i+noffp,mm+moffp)              &
     & + samu(4,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp,j+moffp) = amu(1,nn+noffp,j+moffp)              &
     & + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noffp,j+moffp) = amu(2,nn+noffp,j+moffp)              &
     & + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noffp,j+moffp) = amu(3,nn+noffp,j+moffp)              &
     & + samu(3,nn,j)
!$OMP ATOMIC
         amu(4,nn+noffp,j+moffp) = amu(4,nn+noffp,j+moffp)              &
     & + samu(4,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx,&
     &my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation for relativistic particles
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 62 flops/particle, 1 divide, 21 loads, 16 stores
! input: all, output: amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where pj = pj(t-dt/2) and pk = pk(t-dt/2)
! where gami2 = 1./(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! at t - dt/2
! amu(i,j,k) = ith component of momentum flux at grid point j,kk
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of flux array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, mx, my, nxv, nypmx, mx1, mxyp1
      real qm, ci
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real ci2, gami2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, v1, v2, v3, v4
      real samu
      dimension samu(4,MXV,MYV)
!     dimension samu(4,mx+1,my+1)
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,  
!$OMP& dx,dy,vx,vy,vz,v1,v2,v3,v4,p2,gami2,samu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit momentum flux
      dx = amx*amy
      dy = dxp*amy
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
   30 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
   40 continue
   50 continue
! deposit momentum flux to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp,mm+moffp) = amu(1,i+noffp,mm+moffp)              &
     & + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noffp,mm+moffp) = amu(2,i+noffp,mm+moffp)              &
     & + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noffp,mm+moffp) = amu(3,i+noffp,mm+moffp)              &
     & + samu(3,i,mm)
!$OMP ATOMIC
         amu(4,i+noffp,mm+moffp) = amu(4,i+noffp,mm+moffp)              &
     & + samu(4,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp,j+moffp) = amu(1,nn+noffp,j+moffp)              &
     & + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noffp,j+moffp) = amu(2,nn+noffp,j+moffp)              &
     & + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noffp,j+moffp) = amu(3,nn+noffp,j+moffp)              &
     & + samu(3,nn,j)
!$OMP ATOMIC
         amu(4,nn+noffp,j+moffp) = amu(4,nn+noffp,j+moffp)              &
     & + samu(4,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
