!-----------------------------------------------------------------------
! Fortran Library for pushing electromagnetic particles
! 1-2/2D OpenMP PIC Codes:
! GBPPUSH13L updates magnetized particle co-ordinates and velocities
!            using leap-frog scheme in time and linear interpolation
!            in space with various particle boundary conditions
! GBPPUSHF13L updates magnetized particle co-ordinates and velocities
!             using leap-frog scheme in time and linear interpolation
!             in space with periodic particle boundary conditions,
!             determines list of particles which are leaving each tile
! GRBPPUSH13L updates relativistic magnetized particle co-ordinates
!             and momenta using leap-frog scheme in time and linear
!             interpolation in space with various particle boundary
!             conditions
! GRBPPUSHF13L updates relativistic magnetized particle co-ordinates
!              and momenta using leap-frog scheme in time and linear
!              interpolation in space with periodic particle boundary
!              conditions, determines list of particles which are
!              leaving each tile
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 23, 2016
!-----------------------------------------------------------------------
      subroutine GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp,&
     &nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! OpenMP version using guard cells
! data read in tiles
! particles stored segmented array
! 78 flops/particle, 1 divide, 14 loads, 4 stores
! input: all, output: ppart, ek
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
! position equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real omx, qbm, dt, dtc, ek
      real ppart, fxyz, byz
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, 
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,
!$OMP& rot9,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
! new position
      dx = x + dx*dtc
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + .5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,  &
     &dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data read in tiles
! particles stored segmented array
! 78 flops/particle, 1 divide, 14 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
! position equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real omx, qbm, dt, dtc, ek
      real ppart, fxyz, byz
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtmh, anx, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,
!$OMP& rot8,rot9,edgelx,edgerx,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! load local fields from global array
      nn = nn + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
! clear counters
      do 30 j = 1, 2
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
! new position
      dx = x + dx*dtc
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
! set new position
      ppart(1,j,k) = dx
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   40 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + .5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GRBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,  &
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! OpenMP version using guard cells
! data read in tiles
! particles stored segmented array
! 90 flops/particle, 4 divides, 2 sqrts, 14 loads, 4 stores
! input: all, output: ppart, ek
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fx(x(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fy(x(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fz(x(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equation used is:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! ppart(3,n,m) = momentum py of particle n in tile m
! ppart(4,n,m) = momentum pz of particle n in tile m
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real omx, qbm, dt, dtc, ci, ek
      real ppart, fxyz, byz
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, ci2, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, 
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,
!$OMP& rot9,p2,gami,qtmg,dtg,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GRBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, &
     &dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data read in tiles
! particles stored segmented array
! 90 flops/particle, 4 divides, 2 sqrts, 14 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fx(x(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fy(x(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
!    .5*(q/m)*fz(x(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equation used is:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! ppart(3,n,m) = momentum py of particle n in tile m
! ppart(4,n,m) = momentum pz of particle n in tile m
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real omx, qbm, dt, dtc, ci, ek
      real ppart, fxyz, byz
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtmh, ci2, anx, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,
!$OMP& rot8,rot9,p2,gami,qtmg,dtg,edgelx,edgerx,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! load local fields from global array
      nn = nn + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
! clear counters
      do 30 j = 1, 2
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtg
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
! set new position
      ppart(1,j,k) = dx
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   40 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
