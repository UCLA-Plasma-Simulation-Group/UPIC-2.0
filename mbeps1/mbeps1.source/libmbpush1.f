!-----------------------------------------------------------------------
! Fortran Library for pushing electromagnetic particles
! 1-2/2D OpenMP PIC Codes:
! GBPPUSH13L updates magnetized particle co-ordinates and velocities
!            using leap-frog scheme in time and linear interpolation
!            in space with various particle boundary conditions with
!            Boris pusher
! GBPPUSHF13L updates magnetized particle co-ordinates and velocities
!             using leap-frog scheme in time and linear interpolation
!             in space with periodic particle boundary conditions,
!             with  Boris pusher and determines list of particles which
!             are leaving each tile
! GABPPUSH13L updates magnetized particle co-ordinates and velocities
!             using leap-frog scheme in time and linear interpolation
!             in space with various particle boundary conditions with
!             analytic Boris pusher
! GABPPUSHF13L updates magnetized particle co-ordinates and velocities
!              using leap-frog scheme in time and linear interpolation
!              in space with periodic particle boundary conditions with
!              analytic Boris pusher and determines list of particles
!              which are leaving each tile
! GRBPPUSH13L updates relativistic magnetized particle co-ordinates
!             and momenta using leap-frog scheme in time and linear
!             interpolation in space with various particle boundary
!             conditions with Boris pusher
! GRBPPUSHF13L updates relativistic magnetized particle co-ordinates
!              and momenta using leap-frog scheme in time and linear
!              interpolation in space with periodic particle boundary
!              conditions, determines list of particles which are
!              leaving each tile
! GARBPPUSH13L updates relativistic magnetized particle co-ordinates
!              and momenta using leap-frog scheme in time and linear
!              interpolation in space with various particle boundary
!              conditions with analytic Boris pusher
! GARBPPUSHF13L updates relativistic magnetized particle co-ordinates
!               and momenta using leap-frog scheme in time and linear
!               interpolation in space with periodic particle boundary
!               conditions with analytic Boris pusher and determines
!               list of particles which are leaving each tile
! GEARBPPUSH13L updates relativistic magnetized particle co-ordinates
!               and momenta using leap-frog scheme in time and linear
!               interpolation in space with various particle boundary
!               conditions with exact analytic pusher
! GEARBPPUSHF13L updates relativistic magnetized particle co-ordinates
!                and momenta using leap-frog scheme in time and linear
!                interpolation in space with various particle boundary
!                conditions with exact analytic pusher and determines
!                list of particles which are leaving each tile
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 4, 2021
!-----------------------------------------------------------------------
      subroutine GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp,&
     &nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,sum1,sfxyz,sbyz) REDUCTION(+:sum2) SCHEDULE(dynamic)
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
      ek = ek + 0.5*sum2
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
! particles stored in segmented array
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
!       returns new ntmax required
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy&
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,&
!$OMP& rot8,rot9,edgelx,edgerx,sum1,sfxyz,sbyz) REDUCTION(+:sum2)       &
!$OMP& SCHEDULE(dynamic)
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
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
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
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 60 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GABPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp&
     &,nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! Using the Analytic Boris Mover,
! assumes constant E, B fields during a time step
! OpenMP version, parallelization over particles
! data read in tiles, particles stored in segmented array
! 121 flops/particle, 1 divide, 1 sqrt, 1 tan, 14 loads, 4 stores
! input: all, output: ppart, ek
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(2)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(3)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(5)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(6)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(8)*(vy(t-dt/2) + (q/m)*gy(x(t))t) +
!    rot(9)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
! gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
! gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
! where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
! flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
! position equation used is not exact:
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
!     integer MXV
!     parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real dth, edgelx, edgerx
      real x, dxp, amx, dx, dy, dz, ox, oy, oz, acx, acy, acz
      real omxt, omyt, omzt, qtmh, omt, omti, epl
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      dth = 0.5*dt
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,qtmh,omt,omti,epl,anorm,rot1,rot2,rot3,rot4,rot5, &
!$OMP& rot6,rot7,rot8,rot9,sum1,sfxyz,sbyz) REDUCTION(+:sum2)           &
!$OMP& SCHEDULE(dynamic)
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
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration
      acx = ppart(2,j,k)
      acy = ppart(3,j,k)
      acz = ppart(4,j,k)
      omxt = acx + dx*dth
      omyt = acy + dy*dth
      omzt = acz + dz*dth
c time-centered kinetic energy
      sum1 = sum1 + (omxt*omxt + omyt*omyt + omzt*omzt)
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
      qtmh = dth
c correct the half-acceleration by decomposing E field into components
c parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      epl = dx*ox + dy*oy + dz*oz
      if (omt.gt.0.0) then
         omti = 1.0/omt
         qtmh = omti*tan(dth*omt)
      endif
      epl = epl*omti*omti
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      dx = omxt*dth + (dx - omxt)*qtmh
      dy = omyt*dth + (dy - omyt)*qtmh
      dz = omzt*dth + (dz - omzt)*qtmh
      acx = acx + dx
      acy = acy + dy
      acz = acz + dz
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
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GABPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, &
     &dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! Using the Analytic Boris Moverr, with periodic boundary conditions,
! assumes constant E, B fields during a time step
! also determines list of particles which are leaving this tile
! OpenMP version, parallelization over particles
! data read in tiles, particles stored in segmented array
! 121 flops/particle, 1 divide, 1 sqrt, 1 tan, 14 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(2)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(3)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(5)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(6)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(8)*(vy(t-dt/2) + (q/m)*gy(x(t))t) +
!    rot(9)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
! gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
! gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
! where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
! flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
! position equation used is not exact:
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
!       returns new ntmax required
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
!     integer MXV
!     parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real dth, anx, edgelx, edgerx
      real x, dxp, amx, dx, dy, dz, ox, oy, oz, acx, acy, acz
      real omxt, omyt, omzt, qtmh, omt, omti, epl
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      dth = 0.5*dt
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy&
!$OMP& ,acz,omxt,omyt,omzt,qtmh,omt,omti,epl,anorm,rot1,rot2,rot3,rot4, &
!$OMP& rot5,rot6,rot7,rot8,rot9,edgelx,edgerx,sum1,sfxyz,sbyz)          &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
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
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration
      acx = ppart(2,j,k)
      acy = ppart(3,j,k)
      acz = ppart(4,j,k)
      omxt = acx + dx*dth
      omyt = acy + dy*dth
      omzt = acz + dz*dth
c time-centered kinetic energy
      sum1 = sum1 + (omxt*omxt + omyt*omyt + omzt*omzt)
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
      qtmh = dth
c correct the half-acceleration by decomposing E field into components
c parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      epl = dx*ox + dy*oy + dz*oz
      if (omt.gt.0.0) then
         omti = 1.0/omt
         qtmh = omti*tan(dth*omt)
      endif
      epl = epl*omti*omti
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      dx = omxt*dth + (dx - omxt)*qtmh
      dy = omyt*dth + (dy - omyt)*qtmh
      dz = omzt*dth + (dz - omzt)*qtmh
      acx = acx + dx
      acy = acy + dy
      acz = acz + dz
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
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
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
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 60 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
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
! data read in tiles, particles stored in segmented array
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
!     integer MXV
!     parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, ci2, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      real sfxyz, sbyz
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,p2,gami,qtmg,dtg,sum1,sfxyz,sbyz) REDUCTION(+:sum2)         &
!$OMP& SCHEDULE(dynamic)
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
! particles stored in segmented array
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
!       returns new ntmax required
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy&
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,&
!$OMP& rot8,rot9,p2,gami,qtmg,dtg,edgelx,edgerx,sum1,sfxyz,sbyz)        &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
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
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
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
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 60 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek, &
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Analytic Boris Mover,
! assumes constant E, B fields, and gamma during a time step
! OpenMP version, parallelization over particles
! data read in tiles, particles stored in segmented array
! 291 flops/particle, 7 divides, 3 sqrts, 1 tan, 14 loads, 4 stores
! input: all, output: ppart, ek
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
! py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t))t) +
!    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
! gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
! gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
! where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
! flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami, where
! we approximate analytic gamma function with 3rd order taylor series:
! gam = gam(0) + dg0*tau+d2g0*tau**2/2+d3g0*tau**3/6
! where gam(0) = sqrt(1.0 + p2(0)*ci*ci)
! p2(0) = px(t-dt/2)**2+py(t-dt/2)**2+pz(t-dt/2)**2
! dg0 = ug0/(ci*ci) = dgamma/dtau
! d2g0 = u2g0/(ci*ci) = d2gamma/dtau2
! d3g0 = u3g0/(ci*ci) = -omt*omt*(dgamma_perp/dtau)
! then using the result t = integral of gamma(tau)*dtau, we can
! approximate tau(t) using a one pass Newton method and set gami = tau/t
! position equation used is not exact:
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
! ek = p2(0)/(1.+sqrt(1.+p2(0)*ci*ci))+ug0*th+u2g0*th**2/2+u3g0*th**3/6
! where th = tau/2
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
!     integer MXV
!     parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real sixth, ci2, dth, edgelx, edgerx
      real x, dxp, amx, dx, dy, dz, ox, oy, oz, px, py, pz, p2, gam0, wk
      real ug0, acx, acy, acz, u2g0, omt, omti, epl, omxt, omyt, omzt
      real ugt0, gami, tau, at1, f, fp, qtmg, dtg
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      sixth = 1.0/6.0
      dth = 0.5*dt
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,px,py,pz,p2, &
!$OMP& gam0,wk,ug0,acx,acy,acz,u2g0,omt,omti,epl,omxt,omyt,omzt,ugt0,   &
!$OMP& gami,tau,at1,f,fp,qtmg,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7, &
!$OMP& rot8,rot9,dtg,sum1,sfxyz,sbyz) REDUCTION(+:sum2)                 &
!$OMP& SCHEDULE(dynamic)
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
! normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
! normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
! read momentum
      px = ppart(2,j,k)
      py = ppart(3,j,k)
      pz = ppart(4,j,k)
! find initial gamma
      p2 = px*px + py*py + pz*pz
      gam0 = sqrt(1.0 + p2*ci2)
! initial kinetic energy
      wk = p2/(1.0 + gam0)
! dgamma/dtau
      ug0 = dx*px + dy*py + dz*pz
! d2gamma/dtau2
      acx = py*oz - pz*oy
      acy = pz*ox - px*oz
      acz = px*oy - py*ox
      u2g0 = gam0*(dx*dx + dy*dy + dz*dz) + (dx*acx + dy*acy + dz*acz)
! correct the half-acceleration by decomposing E field into components
! parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      if (omt.gt.0.0) omti = 1.0/omt
      epl = dx*ox + dy*oy + dz*oz
      epl = epl*omti*omti
! E parallel
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
! E perp
      dx = dx - omxt
      dy = dy - omyt
      dz = dz - omzt
! dgamma_perp/dtau
      ugt0 = dx*px + dy*py + dz*pz
! find gami with one pass Newton method and fourth order taylor series
      gami = 1.0/gam0
      tau = dt*gami
      at1 = omt*omt
      f = (0.5*ug0 + sixth*(u2g0 - 0.25*at1*ugt0*tau)*tau)*tau
      fp = (ug0 + (0.5*u2g0 - at1*sixth*ugt0*tau)*tau)*tau
      gami = gami*(1.0 - f*ci2/(gam0 + fp*ci2))
! time-centered kinetic energy
      at1 = -at1*ugt0
      wk = wk + 0.5*(ug0 + 0.25*(u2g0 + sixth*at1*tau)*tau)*tau
      sum1 = sum1 + wk
! set proper time
      qtmg = dth*gami
      if (omt.gt.0.0) qtmg = omti*tan(qtmg*omt)
      at1 = qtmg/gami
! modified half acceleration
      dx = omxt*dth + dx*at1
      dy = omyt*dth + dy*at1
      dz = omzt*dth + dz*at1
      acx = px + dx
      acy = py + dy
      acz = pz + dz
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
      px = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      py = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      pz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = px
      ppart(3,j,k) = py
      ppart(4,j,k) = pz
! update inverse gamma
      p2 = px*px + py*py + pz*pz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + px*dtg
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
      subroutine GARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,&
     &dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine updates particle co-ordinate and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Analytic Boris Mover, with periodic boundary conditions.
! assumes constant E, B fields, and gamma during a time step
! also determines list of particles which are leaving this tile
! OpenMP version, parallelization over particles
! data read in tiles, particles stored in segmented array
! 209 flops/particle, 7 divides, 3 sqrts, 1 tan, 14 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
! py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t))) +
!    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t))) +
!    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t))t) +
!    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
! gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
! gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
! where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
! flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
! the rotation matrix is determined by:
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami, where
! we approximate analytic gamma function with 3rd order taylor series:
! gam = gam(0) + dg0*tau+d2g0*tau**2/2+d3g0*tau**3/6
! where gam(0) = sqrt(1.0 + p2(0)*ci*ci)
! p2(0) = px(t-dt/2)**2+py(t-dt/2)**2+pz(t-dt/2)**2
! dg0 = ug0/(ci*ci) = dgamma/dtau
! d2g0 = u2g0/(ci*ci) = d2gamma/dtau2
! d3g0 = u3g0/(ci*ci) = -omt*omt*(dgamma_perp/dtau)
! then using the result t = integral of gamma(tau)*dtau, we can
! approximate tau(t) using a one pass Newton method and set gami = tau/t
! position equation used is not exact:
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
! ek = p2(0)/(1.+sqrt(1.+p2(0)*ci*ci))+ug0*th+u2g0*th**2/2+u3g0*th**3/6
! where th = tau/2
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
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
!     integer MXV
!     parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real sixth, ci2, dth, anx, edgelx, edgerx
      real x, dxp, amx, dx, dy, dz, ox, oy, oz, px, py, pz, p2, gam0, wk
      real ug0, acx, acy, acz, u2g0, omt, omti, epl, omxt, omyt, omzt
      real ugt0, gami, tau, at1, f, fp, qtmg, dtg
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      sixth = 1.0/6.0
      dth = 0.5*dt
      ci2 = ci*ci
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,px,py, &
!$OMP& pz,p2,gam0,wk,ug0,acx,acy,acz,u2g0,omt,omti,epl,omxt,omyt,omzt,  &
!$OMP& ugt0,gami,tau,at1,f,fp,qtmg,anorm,rot1,rot2,rot3,rot4,rot5,rot6, &
!$OMP& rot7,rot8,rot9,dtg,edgelx,edgerx,sum1,sfxyz,sbyz)                &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
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
! normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
! normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
! read momentum
      px = ppart(2,j,k)
      py = ppart(3,j,k)
      pz = ppart(4,j,k)
! find initial gamma
      p2 = px*px + py*py + pz*pz
      gam0 = sqrt(1.0 + p2*ci2)
! initial kinetic energy
      wk = p2/(1.0 + gam0)
! dgamma/dtau
      ug0 = dx*px + dy*py + dz*pz
! d2gamma/dtau2
      acx = py*oz - pz*oy
      acy = pz*ox - px*oz
      acz = px*oy - py*ox
      u2g0 = gam0*(dx*dx + dy*dy + dz*dz) + (dx*acx + dy*acy + dz*acz)
! correct the half-acceleration by decomposing E field into components
! parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      if (omt.gt.0.0) omti = 1.0/omt
      epl = dx*ox + dy*oy + dz*oz
      epl = epl*omti*omti
! E parallel
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
! E perp
      dx = dx - omxt
      dy = dy - omyt
      dz = dz - omzt
! dgamma_perp/dtau
      ugt0 = dx*px + dy*py + dz*pz
! find gami with one pass Newton method and third order taylor series
      gami = 1.0/gam0
      tau = dt*gami
      at1 = omt*omt
      f = (0.5*ug0 + sixth*(u2g0 - 0.25*at1*ugt0*tau)*tau)*tau
      fp = (ug0 + (0.5*u2g0 - at1*sixth*ugt0*tau)*tau)*tau
      gami = gami*(1.0 - f*ci2/(gam0 + fp*ci2))
! time-centered kinetic energy
      at1 = -at1*ugt0
      wk = wk + 0.5*(ug0 + 0.25*(u2g0 + sixth*at1*tau)*tau)*tau
      sum1 = sum1 + wk
! set proper time
      qtmg = dth*gami
      if (omt.gt.0.0) qtmg = omti*tan(qtmg*omt)
      at1 = qtmg/gami
! modified half acceleration
      dx = omxt*dth + dx*at1
      dy = omyt*dth + dy*at1
      dz = omzt*dth + dz*at1
      acx = px + dx
      acy = py + dy
      acz = pz + dz
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
      px = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      py = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      pz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = px
      ppart(3,j,k) = py
      ppart(4,j,k) = pz
! update inverse gamma
      p2 = px*px + py*py + pz*pz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + px*dtg
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
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
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 60 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GEARBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,&
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
! for 1-2/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Exact Analytic mover, a variant of the algorithm developed
! by Fei Li et al., where E and B fields are constant and gamma varies
! during a time step.
! OpenMP version using guard cells
! 254 FLOPs, 11 divides, 5 sqrts per particle
! plus 51 FLOPs, 2 divides, 1 tan, 1 tanh per Newton iteration/particle
! input: all, output: part, ek
! The particle momenta are advanced in natural coordinates, defined as:
! up is in the direction ep = E parallel to the B field,
! ut is in the direction of et = E perpendicular to the B field,
! ud is in the ExB direction.  Momenta are transformed from cartesian
! coordinates to natural coordinates before the advance, then
! transformed back to cartesian coordinates.
! momenta equations used are:
! up = up + ep*dt
! ut = ut + etp*dt + (om/w)*(F*(cos(w*tau)-1.0) + G*sin(w*tau))
! ud = ved*gam(tau) + F*sin(w*tau) - H*cos(w*tau), where
! F = (w/omega)*ut(0) - om2t*(omega*w)*ci2*ep*up(0)
! G = gam0*(et - al2*om2t) - ud(0)*omega
! H = om2t*omega*gam0 - ud(0)
! om = sqrt(qbm*omx)**2 + (qbm*by)**2 + (qbm*bz)**2)
! etp = om2t*al*al, ved = om*om2t, om2t = et/(om*om + al*al
! the gyration frequency w and al are given by:
! w = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))+om*om-e*e))
! al = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))-om*om+e*e))
! where e = sqrt((qbm*ex)**2+(qbm*ey)**2+(qbm*ez)**2)
! gam(tau) is the relativistic factor given by:
! gam = A*sin(w*tau) - B*cos(w*tau) + C*sinh(al*tau) + D*cos(al*tau),
! where the constants A, B, C, D are defined in terms of initial
! derivatives of the gamma function, as follows:
!  w*A = (w*w*dg0 - om*om*dgp0)/(w*w + al*al)
!  w*B = (w*w*d2g0 - om*om*d2gp0)/(w*w + al*al)
!  al*C = (al*al*dg0 + om*om*dgp0)/(w*w + al*al)
!  al*C = (al*al*d2g0 + om*om*d2gp0)/(w*w + al*al)
! and where the initial derivatives of the gamma function are:
!  dgp0 = ci*ci*ep*up(0)
!  d2gp0 = ci*ci*ep*ep*gam(0)
!  dg0 = dgp0 + ci*ci*et*ut(0)
!  d2g0 = d2gp0 + ci*ci*et*(et*gam(0) - ud(0)*om)
! the proper time tau can be determined from the relativistic factor
! using the equation t = integral gam(tau)*dtau.  This involves finding
! zeros of the function f via an iteration scheme, where
! f = gam(0)*tau - dt - (A*(cos(w*tau)-1.0) + B*(sin(w*tau)-w*tau))/w
!                + (C*(cosh(al*tau)-1.0) + C*(sinh(al*tau)-w*tau))/al
! once we know tau(dt), we can evaluate the momenta up, ut, and ud
! and transform back to cartesian coordinates
! position equation used is not exact:
! x(t+dt) = x(t) + px(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.0+(px(t+dt/2)**2+py(t+dt/2)**2+pz(t+dt/2)**2))*ci*ci)
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
! kinetic energy/mass at time t is also calculated, using the average
! value of energies at the beginning and end of the time step:
! ek = 0.5*((px(t-dt/2)*2+py(t-dt/2)**2+pz(t-dt/2)**2)/(1.0+gam(0))
!    +      (px(t+dt/2)*2+py(t+dt/2)**2+pz(t+dt/2)**2)/(1.0+gam(tau)))
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
      integer imax
      real half, sixth, s1, s2, s3, s4, c1, c2, c3, c4, weps, teps, deps
      real erps
! imax = maximum iteration count
      parameter (imax=25)
      parameter (half=0.5,sixth=1.0/6.0)
      parameter (s1=1.0/6.0,s2=1.0/20.0,s3=1.0/42.0,s4=1.0/72.0)
      parameter (c1=0.5,c2=1.0/12.0,c3=1.0/30.0,c4=1.0/56.0)
! weps = taylor series expansion criteria
      parameter (weps=1.0e-1)
! deps = taylor series expansion criteria for differences
      parameter (deps=1.0e-3)
      integer noff, npp
      integer i, j, k, nn
      real ci2, erm, edgelx, edgerx, tn, x, dxp, amx, small, prec
      real dx, dy, dz, ox, oy, oz, tx, ty, tz, om2, e2, omega, di, ep
      real et, px, py, pz, up, ut, ud, w2, wl2, al2, w, al, wi, ali, wi2
      real ali2, om2t, ws, p2, gam0, dgp0, dgt0, d2gp0, d2gt0, dg0, d2g0
      real wt2, ea, eb, ec, ed, fa, fb, tau, cs, sn, csh, snh, fpp, fp
      real ft, f, t2, t3, wt, wd, tn2, tn2i, csd, snd, cshd, snhd, csc 
      real snc, cshc, snhc, fc, fpi, gam, gami, wk, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2, dt1, err
      data small /1.0e-12/
      prec = 1.0 + small
! teps = tau precision; erps = orthogonality precision
! detect autodouble precision
      if (digits(prec) > 24) then
         teps = 1.0e-14
         erps = 1.0e-10
! default single precision
      else
         teps = 6.0e-7
         erps = 5.0e-3
      endif
      ci2 = ci*ci
      erm = 0.0
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noff,npp,nn,x,dxp,amx,tn,dx,dy,dz,ox,oy,oz,tx,ty,  &
!$OMP& tz,om2,omega,di,ep,e2,et,w2,wl2,al2,w,al,wi,ali,wi2,ali2,om2t,ws,&
!$OMP& px,py,pz,up,ut,ud,p2,gam0,wk,dgp0,dgt0,d2gp0,d2gt0,dg0,d2g0,wt2, &
!$OMP& ea,eb,ec,ed,fa,fb,fc,tau,cs,sn,csh,snh,fpp,fp,ft,t2,t3,wt,wd,tn2,&
!$OMP& tn2i,f,csc,snc,csd,snd,cshc,snhc,cshd,snhd,gam,fpi,gami,err,dtg, &
!$OMP& dt1,sum1,sfxyz,sbyz) REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 50 k = 1, mx1
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
! normalize magnetic field
      ox = qbm*omx
      oy = qbm*oy
      oz = qbm*oz
! normalize electric field
      tx = qbm*dx
      ty = qbm*dy
      tz = qbm*dz
! create direction cosines to translate from/to cartesian coordinates
! first find the direction along the magnetic field B
      om2 = ox*ox + oy*oy + oz*oz
      omega = sqrt(om2)
      if (om2.gt.0.0) then
         di = 1.0/omega
         ox = ox*di
         oy = oy*di
         oz = oz*di
! if omega = 0, then use the y direction
      else
         ox = 0.0
         oy = 1.0
         oz = 0.0
      endif
! then find the direction along the electric field Et perpendicular to B
      dt1 = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
      tx = tx - dt1*ox
      ty = ty - dt1*oy
      tz = tz - dt1*oz
      ep = dt1
      e2 = tx*tx + ty*ty + tz*tz
      et = sqrt(e2)
      if (et > 0.0) then
         di = 1.0/et
         tx = tx*di
         ty = ty*di
         tz = tz*di
! then find the direction along Et x B
         dx = ty*oz - tz*oy
         dy = tz*ox - tx*oz
         dz = tx*oy - ty*ox
! check for roundoff error
         err = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
         if (err > erps) then
            write (*,*) 'Error: Et not normal to omega=',et,err
            et = 0.0
         else
            err = dble(tx*dx) + dble(ty*dy) + dble(tz*dz)
            if (err > erps) then
               write (*,*) 'Error: Et d=',et,err
            endif
         endif
      endif
! special case Et = 0, or round off error detected
      if (et==0.0) then      
! first find direction with smallest component of B
         i = 1
         dx = abs(ox)
         dy = abs(oy)
         dz = abs(oz)
         di = dx
         if (dy <= dx) then
            i = 2
            di = dy
            if ((dx.gt.dy).and.(dy.eq.dz)) i = 3
         endif
         if (dz.lt.di) i = 3 
! then the cross product of that direction with B
         if (i.eq.1) then
            dz = 1.0/sqrt(oy*oy + oz*oz)
            dx = 0.0
            dy = -oz*dz
            dz = oy*dz
         else if (i.eq.2) then
            dz = 1.0/sqrt(ox*ox + oz*oz)
            dx = oz*dz
            dy = 0.0
            dz = -ox*dz
         else if (i.eq.3) then
            dz = 1.0/sqrt(ox*ox + oy*oy)
            dx = -oy*dz
            dy = ox*dz
            dz = 0.0
         endif
! then find the direction along minus d x B
         tx = dz*oy - dy*oz
         ty = dx*oz - dz*ox
         tz = dy*ox - dx*oy
      endif
! calculate frequencies
      e2 = (e2 + ep*ep)*ci2
      w2 = om2 - e2
      wl2 = sqrt(w2*w2 + 4.0*ci2*(ep*omega)**2)
      al2 = 0.5*(wl2 - w2)
      w2 = 0.5*(wl2 + w2)
      w = sqrt(w2)
      al = sqrt(al2)
      wi = 0.0
      if (w > 0.0) wi = 1.0/w
      ali = 0.0
      if (al > 0.0) ali = 1.0/al
      wi2 = wi*wi
      ali2 = ali*ali
! calculate weights
      om2t = om2 + al2
      if (om2t > 0.0) om2t = et/om2t
      ws = 0.0
      if (omega /= 0.0) ws = w/omega
! translate momenta from cartesian coordinates
      px = ppart(2,j,k)
      py = ppart(3,j,k)
      pz = ppart(4,j,k)
      up = px*ox + py*oy + pz*oz
      ut = px*tx + py*ty + pz*tz
      ud = px*dx + py*dy + pz*dz
! find initial gamma
      p2 = up*up + ut*ut + ud*ud
      gam0 = sqrt(1.0 + p2*ci2)
! calculate initial kinetic energy
      wk = p2/(1.0 + gam0)
! partial derivatives of gamma
      dgp0 = ci2*ep*up
      dgt0 = ci2*et*ut
      d2gp0 = ci2*ep*ep*gam0
      d2gt0 = ci2*et*(et*gam0 - ud*omega)
      dg0 = dgp0 + dgt0
      d2g0 = d2gp0 + d2gt0
! calculate trigonometric and hyperbolic coefficients
      wt2 = 0.0
      if (wl2 > 0.0) wt2 = 1.0/wl2
      ea = wt2*(w2*dg0 - om2*dgp0)
      eb = wt2*(w2*d2g0 - om2*d2gp0)
      ec = wt2*(al2*dg0 + om2*dgp0)
      ed = wt2*(al2*d2g0 + om2*d2gp0)
      if (wl2==0.0) then
         ea = dgt0
         eb = d2gt0
      endif
      fa = ws*ut - om2t*omega*wi*ci2*ep*up
      fb = gam0*(et - al2*om2t) - ud*omega
      fc = om2t*omega*gam0 - ud
! zeroth order guess for tau
      tau = dt/gam0
!
! iteration loop for finding tau(t)
      cs = 1.0
      sn = 0.0
      csh = 1.0
      snh = 0.0
      fpp = 0.0
      fp = 0.0
      ft = -tau
      do 30 i = 1, imax
      t2 = tau*tau; t3 = t2*tau
! calculate trigonometric functions
      wt = w*tau
      wd = w*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tan(0.5*wt)
            tn2 = tn*tn
            tn2i = 1.0/(1.0 + tn2)
            cs = (1.0 - tn2)*tn2i
            sn = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 - s1*wt2)
            wt2 = 1.0 - c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 - s1*wt2*(1.0 - s2*wt2))
!           wt2 = 1.0 - c1*wt2*(1.0 - c2*wt2)
            f = cs*wt2 - sn*wd
            sn = sn*wt2 + cs*wd
            cs = f
         endif
! calculate special functions
         csc = cs - 1.0
         snc = sn*wi
         csd = csc*wi2
         snd = (snc - tau)*wi2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2))*t2
         snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2))*t3
! fifth order taylor series approximation for wt
!        csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2*(1.0 - c4*wt2)))*t2
!        snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2*(1.0 - s4*wt2)))*t3
         csc = w2*csd
         snc = tau + w2*snd
         cs = csc + 1.0
         sn = w*snc
      endif
! calculate hyperbolic functions
      wt = al*tau
      wd = al*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tanh(0.5*wt)
            tn2 = tn*tn;
            tn2i = 1.0/(1.0 - tn2)
            csh = (1.0 + tn2)*tn2i
            snh = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 + s1*wt2)
            wt2 = 1.0 + c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 + s1*wt2*(1.0 + s2*wt2))
!           wt2 = 1.0 + c1*wt2*(1.0 + c2*wt2)
            f = csh*wt2 + snh*wd
            snh = snh*wt2 + csh*wd
            csh = f
         endif
! calculate special functions
         cshc = csh - 1.0
         snhc = snh*ali
         cshd = cshc*ali2
         snhd = (snhc - tau)*ali2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2))*t2
         snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2))*t3
! fifth order taylor series approximation for wt
!        cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2*(1.0 + c4*wt2)))*t2
!        snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2*(1.0 + s4*wt2)))*t3
         cshc = al2*cshd
         snhc = tau + al2*snhd
         csh = cshc + 1.0
         snh = al*snhc
      endif
! gam = gamma(tau)
      gam = gam0 + (ea*snc - eb*csd) + (ec*snhc + ed*cshd)
      fpi = 1.0/gam
! calculate time expression whose root we seek
      f = gam0*tau - (ea*csd + eb*snd) + (ec*cshd + ed*snhd) - dt
! newton's quadratic method
      ft = f*fpi
! either add Halley's optional cubic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
!     ft = ft/(1.0d0 - 0.5d0*ft*fpp*fpi)
! or add Householder's optional quartic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     wt2 = ft*fpp*fpi
!     ft = ft*(1.0d0 - 0.5d0*wt2)/(1.0d0 - wt2 + sixth*(ft*ft)*fp*fpi)
! update tau: ft = -delta tau
      tau = tau - ft
      if (abs(f) < teps) exit
   30 continue
!
! convergence failure
      if (i.gt.imax) write (*,*) i,'tau error:f,ft=',f,ft
! update gamma
! fpp = dgamma/dt
      if (fpp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! first order taylor series
      gam = gam - fpp*ft
! second order taylor series
! fpp = dgamma/dtau
!     if (fp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     if (fp.eq.0.0) fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     gam = gam - (fpp - half*fp*ft)*ft
!
! update sine/cosine functions
! first order taylor series
      wd = w*ft
      f = sn*wd
      csc = csc + f
      snc = snc - cs*ft
      f = cs + f
      sn = sn - cs*wd
      cs = f
! second order taylor series
!     wd = w*ft
!     wt2 = wd*wd
!     f = -ft*(1.0d0 - sixth*wt2)
!     wd = half*wt2
!     wt2 = 1.0d0 - wd
!     sc = csc*wt2 - wd
!     wd = w*f
!     sc = csc - sn*wd
!     snc = snc*wt2 + cs*f
!     f = cs*wt2 - sn*wd
!     sn = sn*wt2 + cs*wd
!     cs = f
! update momenta
      up = up + ep*dt
      ut = ut + om2t*al2*dt + (omega*wi*fa*csc + fb*snc)
      ud = om2t*omega*gam + fa*sn - fc*cs
! calculate inverse gamma and average kinetic energy
      gami = 1.0/gam
      p2 = up*up + ut*ut + ud*ud
      sum1 = sum1 + dble(0.5*(wk + gami*p2/(1.0 + gami)))
! sanity check
!     erm = max(erm,abs(gam-sqrt(1.0d0 + p2*ci2)))
! translate momenta to cartesian coordinates
      px = up*ox + ut*tx + ud*dx
      py = up*oy + ut*ty + ud*dy
      pz = up*oz + ut*tz + ud*dz
      ppart(2,j,k) = px
      ppart(3,j,k) = py
      ppart(4,j,k) = pz
! new position
      dtg = dtc*gami
      dx = x + px*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(2,j,k) = -px
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   40 continue
      sum2 = sum2 + sum1
   50 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
! sanity check
!     write (*,*) 'gamma sanity check=',erm
      return
      end
!-----------------------------------------------------------------------
      subroutine GEARBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt&
     &,dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Exact Analytic mover, a variant of the algorithm developed
! by Fei Li et al., where E and B fields are constant and gamma varies
! during a time step.
! also determines list of particles which are leaving this tile
! OpenMP version, parallelization over particles
! 254 FLOPs, 11 divides, 5 sqrts per particle
! plus 51 FLOPs, 2 divides, 1 tan, 1 tanh per Newton iteration/particle
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! The particle momenta are advanced in natural coordinates, defined as:
! up is in the direction ep = E parallel to the B field,
! ut is in the direction of et = E perpendicular to the B field,
! ud is in the ExB direction.  Momenta are transformed from cartesian
! coordinates to natural coordinates before the advance, then
! transformed back to cartesian coordinates.
! momenta equations used are:
! up = up + ep*dt
! ut = ut + etp*dt + (om/w)*(F*(cos(w*tau)-1.0) + G*sin(w*tau))
! ud = ved*gam(tau) + F*sin(w*tau) - H*cos(w*tau), where
! F = (w/omega)*ut(0) - om2t*(omega*w)*ci2*ep*up(0)
! G = gam0*(et - al2*om2t) - ud(0)*omega
! H = om2t*omega*gam0 - ud(0)
! om = sqrt(qbm*omx)**2 + (qbm*by)**2 + (qbm*bz)**2)
! etp = om2t*al*al, ved = om*om2t, om2t = et/(om*om + al*al
! the gyration frequency w and al are given by:
! w = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))+om*om-e*e))
! al = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))-om*om+e*e))
! where e = sqrt((qbm*ex)**2+(qbm*ey)**2+(qbm*ez)**2)
! gam(tau) is the relativistic factor given by:
! gam = A*sin(w*tau) - B*cos(w*tau) + C*sinh(al*tau) + D*cos(al*tau),
! where the constants A, B, C, D are defined in terms of initial
! derivatives of the gamma function, as follows:
!  w*A = (w*w*dg0 - om*om*dgp0)/(w*w + al*al)
!  w*B = (w*w*d2g0 - om*om*d2gp0)/(w*w + al*al)
!  al*C = (al*al*dg0 + om*om*dgp0)/(w*w + al*al)
!  al*C = (al*al*d2g0 + om*om*d2gp0)/(w*w + al*al)
! and where the initial derivatives of the gamma function are:
!  dgp0 = ci*ci*ep*up(0)
!  d2gp0 = ci*ci*ep*ep*gam(0)
!  dg0 = dgp0 + ci*ci*et*ut(0)
!  d2g0 = d2gp0 + ci*ci*et*(et*gam(0) - ud(0)*om)
! the proper time tau can be determined from the relativistic factor
! using the equation t = integral gam(tau)*dtau.  This involves finding
! zeros of the function f via an iteration scheme, where
! f = gam(0)*tau - dt - (A*(cos(w*tau)-1.0) + B*(sin(w*tau)-w*tau))/w
!                + (C*(cosh(al*tau)-1.0) + C*(sinh(al*tau)-w*tau))/al
! once we know tau(dt), we can evaluate the momenta up, ut, and ud
! and transform back to cartesian coordinates
! position equation used is not exact:
! x(t+dt) = x(t) + px(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.0+(px(t+dt/2)**2+py(t+dt/2)**2+pz(t+dt/2)**2))*ci*ci)
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
! kinetic energy/mass at time t is also calculated, using the average
! value of energies at the beginning and end of the time step:
! ek = 0.5*((px(t-dt/2)*2+py(t-dt/2)**2+pz(t-dt/2)**2)/(1.0+gam(0))
!    +      (px(t+dt/2)*2+py(t+dt/2)**2+pz(t+dt/2)**2)/(1.0+gam(tau)))
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
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
      integer imax
      real half, sixth, s1, s2, s3, s4, c1, c2, c3, c4, weps, teps, deps
      real erps
! imax = maximum iteration count
      parameter (imax=25)
      parameter (half=0.5,sixth=1.0/6.0)
      parameter (s1=1.0/6.0,s2=1.0/20.0,s3=1.0/42.0,s4=1.0/72.0)
      parameter (c1=0.5,c2=1.0/12.0,c3=1.0/30.0,c4=1.0/56.0)
! weps = taylor series expansion criteria
      parameter (weps=1.0e-1)
! deps = taylor series expansion criteria for differences
      parameter (deps=1.0e-3)
      integer noff, npp
      integer i, j, k, ih, nh, nn
      real ci2, erm, anx, edgelx, edgerx, tn, x, dxp, amx, small, prec
      real dx, dy, dz, ox, oy, oz, tx, ty, tz, om2, e2, omega, di, ep
      real et, px, py, pz, up, ut, ud, w2, wl2, al2, w, al, wi, ali, wi2
      real ali2, om2t, ws, p2, gam0, dgp0, dgt0, d2gp0, d2gt0, dg0, d2g0
      real wt2, ea, eb, ec, ed, fa, fb, tau, cs, sn, csh, snh, fpp, fp
      real ft, f, t2, t3, wt, wd, tn2, tn2i, csd, snd, cshd, snhd, csc 
      real snc, cshc, snhc, fc, fpi, gam, gami, wk, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2, dt1, err
      data small /1.0e-12/
      prec = 1.0 + small
! teps = tau precision; erps = orthogonality precision
! detect autodouble precision
      if (digits(prec) > 24) then
         teps = 1.0e-14
         erps = 1.0e-10
! default single precision
      else
         teps = 6.0e-7
         erps = 5.0e-3
      endif
      ci2 = ci*ci
      anx = real(nx)
      erm = 0.0
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noff,npp,nn,ih,nh,x,dxp,amx,tn,dx,dy,dz,ox,oy,oz,  &
!$OMP& tx,ty,tz,om2,omega,di,ep,e2,et,w2,wl2,al2,w,al,wi,ali,wi2,ali2,  &
!$OMP& om2t,ws,px,py,pz,up,ut,ud,p2,gam0,wk,dgp0,dgt0,d2gp0,d2gt0,dg0,  &
!$OMP& d2g0,wt2,ea,eb,ec,ed,fa,fb,fc,tau,cs,sn,csh,snh,fpp,fp,ft,t2,t3, &
!$OMP& wt,wd,tn2,tn2i,f,csc,snc,csd,snd,cshc,snhc,cshd,snhd,gam,fpi,    &
!$OMP& gami,err,dtg,dt1,edgelx,edgerx,sum1,sfxyz,sbyz) REDUCTION(+:sum2)&
!$OMP& SCHEDULE(dynamic)
      do 60 k = 1, mx1
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
      do 50 j = 1, npp
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
! normalize magnetic field
      ox = qbm*omx
      oy = qbm*oy
      oz = qbm*oz
! normalize electric field
      tx = qbm*dx
      ty = qbm*dy
      tz = qbm*dz
! create direction cosines to translate from/to cartesian coordinates
! first find the direction along the magnetic field B
      om2 = ox*ox + oy*oy + oz*oz
      omega = sqrt(om2)
      if (om2.gt.0.0) then
         di = 1.0/omega
         ox = ox*di
         oy = oy*di
         oz = oz*di
! if omega = 0, then use the y direction
      else
         ox = 0.0
         oy = 1.0
         oz = 0.0
      endif
! then find the direction along the electric field Et perpendicular to B
      dt1 = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
      tx = tx - dt1*ox
      ty = ty - dt1*oy
      tz = tz - dt1*oz
      ep = dt1
      e2 = tx*tx + ty*ty + tz*tz
      et = sqrt(e2)
      if (et > 0.0) then
         di = 1.0/et
         tx = tx*di
         ty = ty*di
         tz = tz*di
! then find the direction along Et x B
         dx = ty*oz - tz*oy
         dy = tz*ox - tx*oz
         dz = tx*oy - ty*ox
! check for roundoff error
         err = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
         if (err > erps) then
            write (*,*) 'Error: Et not normal to omega=',et,err
            et = 0.0
         else
            err = dble(tx*dx) + dble(ty*dy) + dble(tz*dz)
            if (err > erps) then
               write (*,*) 'Error: Et d=',et,err
            endif
         endif
      endif
! special case Et = 0, or round off error detected
      if (et==0.0) then      
! first find direction with smallest component of B
         i = 1
         dx = abs(ox)
         dy = abs(oy)
         dz = abs(oz)
         di = dx
         if (dy <= dx) then
            i = 2
            di = dy
            if ((dx.gt.dy).and.(dy.eq.dz)) i = 3
         endif
         if (dz.lt.di) i = 3 
! then the cross product of that direction with B
         if (i.eq.1) then
            dz = 1.0/sqrt(oy*oy + oz*oz)
            dx = 0.0
            dy = -oz*dz
            dz = oy*dz
         else if (i.eq.2) then
            dz = 1.0/sqrt(ox*ox + oz*oz)
            dx = oz*dz
            dy = 0.0
            dz = -ox*dz
         else if (i.eq.3) then
            dz = 1.0/sqrt(ox*ox + oy*oy)
            dx = -oy*dz
            dy = ox*dz
            dz = 0.0
         endif
! then find the direction along minus d x B
         tx = dz*oy - dy*oz
         ty = dx*oz - dz*ox
         tz = dy*ox - dx*oy
      endif
! calculate frequencies
      e2 = (e2 + ep*ep)*ci2
      w2 = om2 - e2
      wl2 = sqrt(w2*w2 + 4.0*ci2*(ep*omega)**2)
      al2 = 0.5*(wl2 - w2)
      w2 = 0.5*(wl2 + w2)
      w = sqrt(w2)
      al = sqrt(al2)
      wi = 0.0
      if (w > 0.0) wi = 1.0/w
      ali = 0.0
      if (al > 0.0) ali = 1.0/al
      wi2 = wi*wi
      ali2 = ali*ali
! calculate weights
      om2t = om2 + al2
      if (om2t > 0.0) om2t = et/om2t
      ws = 0.0
      if (omega /= 0.0) ws = w/omega
! translate momenta from cartesian coordinates
      px = ppart(2,j,k)
      py = ppart(3,j,k)
      pz = ppart(4,j,k)
      up = px*ox + py*oy + pz*oz
      ut = px*tx + py*ty + pz*tz
      ud = px*dx + py*dy + pz*dz
! find initial gamma
      p2 = up*up + ut*ut + ud*ud
      gam0 = sqrt(1.0 + p2*ci2)
! calculate initial kinetic energy
      wk = p2/(1.0 + gam0)
! partial derivatives of gamma
      dgp0 = ci2*ep*up
      dgt0 = ci2*et*ut
      d2gp0 = ci2*ep*ep*gam0
      d2gt0 = ci2*et*(et*gam0 - ud*omega)
      dg0 = dgp0 + dgt0
      d2g0 = d2gp0 + d2gt0
! calculate trigonometric and hyperbolic coefficients
      wt2 = 0.0
      if (wl2 > 0.0) wt2 = 1.0/wl2
      ea = wt2*(w2*dg0 - om2*dgp0)
      eb = wt2*(w2*d2g0 - om2*d2gp0)
      ec = wt2*(al2*dg0 + om2*dgp0)
      ed = wt2*(al2*d2g0 + om2*d2gp0)
      if (wl2==0.0) then
         ea = dgt0
         eb = d2gt0
      endif
      fa = ws*ut - om2t*omega*wi*ci2*ep*up
      fb = gam0*(et - al2*om2t) - ud*omega
      fc = om2t*omega*gam0 - ud
! zeroth order guess for tau
      tau = dt/gam0
!
! iteration loop for finding tau(t)
      cs = 1.0
      sn = 0.0
      csh = 1.0
      snh = 0.0
      fpp = 0.0
      fp = 0.0
      ft = -tau
      do 40 i = 1, imax
      t2 = tau*tau; t3 = t2*tau
! calculate trigonometric functions
      wt = w*tau
      wd = w*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tan(0.5*wt)
            tn2 = tn*tn
            tn2i = 1.0/(1.0 + tn2)
            cs = (1.0 - tn2)*tn2i
            sn = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 - s1*wt2)
            wt2 = 1.0 - c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 - s1*wt2*(1.0 - s2*wt2))
!           wt2 = 1.0 - c1*wt2*(1.0 - c2*wt2)
            f = cs*wt2 - sn*wd
            sn = sn*wt2 + cs*wd
            cs = f
         endif
! calculate special functions
         csc = cs - 1.0
         snc = sn*wi
         csd = csc*wi2
         snd = (snc - tau)*wi2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2))*t2
         snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2))*t3
! fifth order taylor series approximation for wt
!        csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2*(1.0 - c4*wt2)))*t2
!        snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2*(1.0 - s4*wt2)))*t3
         csc = w2*csd
         snc = tau + w2*snd
         cs = csc + 1.0
         sn = w*snc
      endif
! calculate hyperbolic functions
      wt = al*tau
      wd = al*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tanh(0.5*wt)
            tn2 = tn*tn;
            tn2i = 1.0/(1.0 - tn2)
            csh = (1.0 + tn2)*tn2i
            snh = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 + s1*wt2)
            wt2 = 1.0 + c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 + s1*wt2*(1.0 + s2*wt2))
!           wt2 = 1.0 + c1*wt2*(1.0 + c2*wt2)
            f = csh*wt2 + snh*wd
            snh = snh*wt2 + csh*wd
            csh = f
         endif
! calculate special functions
         cshc = csh - 1.0
         snhc = snh*ali
         cshd = cshc*ali2
         snhd = (snhc - tau)*ali2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2))*t2
         snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2))*t3
! fifth order taylor series approximation for wt
!        cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2*(1.0 + c4*wt2)))*t2
!        snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2*(1.0 + s4*wt2)))*t3
         cshc = al2*cshd
         snhc = tau + al2*snhd
         csh = cshc + 1.0
         snh = al*snhc
      endif
! gam = gamma(tau)
      gam = gam0 + (ea*snc - eb*csd) + (ec*snhc + ed*cshd)
      fpi = 1.0/gam
! calculate time expression whose root we seek
      f = gam0*tau - (ea*csd + eb*snd) + (ec*cshd + ed*snhd) - dt
! newton's quadratic method
      ft = f*fpi
! either add Halley's optional cubic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
!     ft = ft/(1.0d0 - 0.5d0*ft*fpp*fpi)
! or add Householder's optional quartic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     wt2 = ft*fpp*fpi
!     ft = ft*(1.0d0 - 0.5d0*wt2)/(1.0d0 - wt2 + sixth*(ft*ft)*fp*fpi)
! update tau: ft = -delta tau
      tau = tau - ft
      if (abs(f) < teps) exit
   40 continue
!
! convergence failure
      if (i.gt.imax) write (*,*) i,'tau error:f,ft=',f,ft
! update gamma
! fpp = dgamma/dt
      if (fpp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! first order taylor series
      gam = gam - fpp*ft
! second order taylor series
! fpp = dgamma/dtau
!     if (fp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     if (fp.eq.0.0) fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     gam = gam - (fpp - half*fp*ft)*ft
!
! update sine/cosine functions
! first order taylor series
      wd = w*ft
      f = sn*wd
      csc = csc + f
      snc = snc - cs*ft
      f = cs + f
      sn = sn - cs*wd
      cs = f
! second order taylor series
!     wd = w*ft
!     wt2 = wd*wd
!     f = -ft*(1.0d0 - sixth*wt2)
!     wd = half*wt2
!     wt2 = 1.0d0 - wd
!     sc = csc*wt2 - wd
!     wd = w*f
!     sc = csc - sn*wd
!     snc = snc*wt2 + cs*f
!     f = cs*wt2 - sn*wd
!     sn = sn*wt2 + cs*wd
!     cs = f
! update momenta
      up = up + ep*dt
      ut = ut + om2t*al2*dt + (omega*wi*fa*csc + fb*snc)
      ud = om2t*omega*gam + fa*sn - fc*cs
! calculate inverse gamma and average kinetic energy
      gami = 1.0/gam
      p2 = up*up + ut*ut + ud*ud
      sum1 = sum1 + dble(0.5*(wk + gami*p2/(1.0 + gami)))
! sanity check
!     erm = max(erm,abs(gam-sqrt(1.0d0 + p2*ci2)))
! translate momenta to cartesian coordinates
      px = up*ox + ut*tx + ud*dx
      py = up*oy + ut*ty + ud*dy
      pz = up*oz + ut*tz + ud*dz
      ppart(2,j,k) = px
      ppart(3,j,k) = py
      ppart(4,j,k) = pz
! new position
      dtg = dtc*gami
      dx = x + px*dtg
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
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
   50 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   60 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 70 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   70    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
! sanity check
!     write (*,*) 'gamma sanity check=',erm
      return
      end
