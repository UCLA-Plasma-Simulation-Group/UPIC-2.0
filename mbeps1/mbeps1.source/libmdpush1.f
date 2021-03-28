!-----------------------------------------------------------------------
! Fortran Library for depositing time derivative of current
! 1-2/2D OpenMP PIC Codes:
! FWPMINMX1 calculates maximum and minimum plasma frequency
! FWPTMINMX1 calculates maximum and minimum total plasma frequency
! GDJPPOST1L calculates particle momentum flux and acceleration
!            density using linear interpolation
! GDCJPPOST1L calculates particle momentum flux, acceleration density
!             and current density using linear interpolation.
! GRDJPPOST1L calculates relativistic particle momentum flux and
!             acceleration density using linear interpolation
! GRDCJPPOST1L calculates relativistic particle momentum flux,
!              acceleration density and current density using linear
!              interpolation.
! ASCFGUARD1L add scaled vector field to extended periodic field
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 4, 2021
!-----------------------------------------------------------------------
      subroutine FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
! calculates maximum and minimum plasma frequency.  assumes guard cells
! have already been added
! qe = charge density for electrons
! qbme = charge/mass ratio for electrons
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qbme, wpmax, wpmin
      integer nx, nxe
      dimension qe(nxe)
! local data
      integer j
      real at1
      wpmax = qbme*qe(1)
      wpmin = wpmax
      do 10 j = 1, nx
      at1 = qbme*qe(j)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FWPTMINMX1(qe,qi,qbme,qbmi,wpmax,wpmin,nx,nxe)
! calculates maximum and minimum total plasma frequency.  assumes guard
! cells have already been added
! qe/qi = charge density for electrons/ions
! qbme/qbmi = charge/mass ratio for electrons/ions
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qi, qbme, qbmi, wpmax, wpmin
      integer nx, nxe
      dimension qe(nxe), qi(nxe)
c local data
      integer j
      real at1
      wpmax = qbme*qe(1) + qbmi*qi(1)
      wpmin = wpmax
      do 10 j = 1, nx
      at1 = qbme*qe(j) + qbmi*qi(j)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine GDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt,  &
     &idimp,nppmx,nx,mx,nxv,mx1)
! for 1-2/2d code, this subroutine calculates particle momentum flux
! and acceleration density using first-order spline interpolation.
! OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 100 flops/particle, 1 divide, 12 loads, 8 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! velocity equations at t=t+dt/2 are calculated from:
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
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = velocity vx of particle n in tile m at t - dt/2
! ppart(3,n,m) = velocity vy of particle n in tile m at t - dt/2
! ppart(4,n,m) = velocity vz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! dcu(i,j) = ith component of acceleration density at grid point j
! for i = 1, 2
! amu(i,j) = ith component of momentum flux at grid point j
! for i = 1, 2
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1
      real omx, qm, qbm, dt
      real ppart, fxyz, byz, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension dcu(2,nxv), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, dti, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, vx, vy, vz, v1, v2
      real sfxyz, sbyz, sdcu, samu
      dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sdcu(2,MXV), samu(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
!     dimension sdcu(2,mx+1), samu(2,mx+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,vx,vy,vz,acx,&
!$OMP& acy,acz,v1,v2,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5, &
!$OMP& rot6,rot7,rot8,rot9,sfxyz,sbyz,samu,sdcu) SCHEDULE(dynamic)
      do 60 k = 1, mx1
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
! zero out local accumulators
      do 30 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
   30 continue
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
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
! deposit momentum flux and acceleration density to local accumulator
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
      sdcu(1,nn) = sdcu(1,nn) + vy*amx
      sdcu(2,nn) = sdcu(2,nn) + vz*amx
      sdcu(1,nn+1) = sdcu(1,nn+1) + vy*dxp
      sdcu(2,nn+1) = sdcu(2,nn+1) + vz*dxp
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 50 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
      dcu(1,j+noff) = dcu(1,j+noff) + sdcu(1,j)
      dcu(2,j+noff) = dcu(2,j+noff) + sdcu(2,j)
   50 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
!$OMP ATOMIC
      dcu(1,1+noff) = dcu(1,1+noff) + sdcu(1,1)
!$OMP ATOMIC
      dcu(2,1+noff) = dcu(2,1+noff) + sdcu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
!$OMP ATOMIC
         dcu(1,nn+noff) = dcu(1,nn+noff) + sdcu(1,nn)
!$OMP ATOMIC
         dcu(2,nn+noff) = dcu(2,nn+noff) + sdcu(2,nn)
      endif
   60 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,qbm, &
     &dt,idimp,nppmx,nx,mx,nxv,mx1)
! for 1-2/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation.
! OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 108 flops/particle, 1 divide, 16 loads, 12 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vj, where j = y,z, for i = 1, 2
! where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! velocity equations at t=t+dt/2 are calculated from:
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
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = velocity vx of particle n in tile m at t - dt/2
! ppart(3,n,m) = velocity vy of particle n in tile m at t - dt/2
! ppart(4,n,m) = velocity vz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! cu(i,j) = ith component of current density at grid point j
! for i = 1, 2
! dcu(i,j) = ith component of acceleration density at grid point j
! for i = 1, 2
! amu(i,j) = ith component of momentum flux at grid point j
! for i = 1, 2
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1
      real omx, qm, qbm, dt
      real ppart, fxyz, byz, cu, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, dti, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, vx, vy, vz, v1, v2
      real sfxyz, sbyz, scu, sdcu, samu
      dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension scu(2,MXV), sdcu(2,MXV), samu(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
!     dimension scu(2,mx+1), sdcu(2,mx+1), samu(2,mx+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,vx,vy,vz,acx,&
!$OMP& acy,acz,v1,v2,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5, &
!$OMP& rot6,rot7,rot8,rot9,sfxyz,sbyz,samu,sdcu,scu) SCHEDULE(dynamic)
      do 60 k = 1, mx1
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
! zero out local accumulators
      do 30 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   30 continue
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
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
! deposit momentum flux, acceleration density, and current density
! to local accumulator
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
      sdcu(1,nn) = sdcu(1,nn) + vy*amx
      sdcu(2,nn) = sdcu(2,nn) + vz*amx
      sdcu(1,nn+1) = sdcu(1,nn+1) + vy*dxp
      sdcu(2,nn+1) = sdcu(2,nn+1) + vz*dxp
      scu(1,nn) = scu(1,nn) + oy*amx
      scu(2,nn) = scu(2,nn) + oz*amx
      scu(1,nn+1) = scu(1,nn+1) + oy*dxp
      scu(2,nn+1) = scu(2,nn+1) + oz*dxp
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 50 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
      dcu(1,j+noff) = dcu(1,j+noff) + sdcu(1,j)
      dcu(2,j+noff) = dcu(2,j+noff) + sdcu(2,j)
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   50 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
!$OMP ATOMIC
      dcu(1,1+noff) = dcu(1,1+noff) + sdcu(1,1)
!$OMP ATOMIC
      dcu(2,1+noff) = dcu(2,1+noff) + sdcu(2,1)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
!$OMP ATOMIC
         dcu(1,nn+noff) = dcu(1,nn+noff) + sdcu(1,nn)
!$OMP ATOMIC
         dcu(2,nn+noff) = dcu(2,nn+noff) + sdcu(2,nn)
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   60 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,ci,qbm, &
     &dt,idimp,nppmx,nx,mx,nxv,mx1)
! for 1-2/2d code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation
! for relativistic particles.
! OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 120 flops/particle, 2 divide, 1 sqrt, 16 loads, 12 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*dvj*gami/dt, where j = y,z, for i = 1, 2
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pj*pk*gami**2, where jk = xy,xz, for i = 1, 2
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
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
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami.
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = momentum px of particle n in tile m at t - dt/2
! ppart(3,n,m) = momentum py of particle n in tile m at t - dt/2
! ppart(4,n,m) = momentum pz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! dcu(i,j) = ith component of acceleration density at grid point j
! for i = 1, 2
! amu(i,j) = ith component of momentum flux at grid point j
! for i = 1, 2
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1
      real omx, qm, qbm, dt, ci
      real ppart, fxyz, byz, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension dcu(2,nxv), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, vx, vy, vz, p2, v1, v2, v3
      real sfxyz, sbyz, sdcu, samu
      dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension sdcu(2,MXV), samu(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
!     dimension sdcu(2,mx+1), samu(2,mx+1)
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      dti = 1.0/dt
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,vx,vy,vz,acx,&
!$OMP& acy,acz,v1,v2,v3,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,   &
!$OMP& rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,gh,sfxyz,sbyz,sdcu,samu)   &
!$OMP& SCHEDULE(dynamic)
      do 60 k = 1, mx1
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
! zero out local accumulators
      do 30 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
   30 continue
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
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
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
! new momentum
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! deposit momentum flux and acceleration density to local accumulator
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      v1 = ox*oy
      v2 = ox*oz
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
      sdcu(1,nn) = sdcu(1,nn) + vy*amx
      sdcu(2,nn) = sdcu(2,nn) + vz*amx
      sdcu(1,nn+1) = sdcu(1,nn+1) + vy*dxp
      sdcu(2,nn+1) = sdcu(2,nn+1) + vz*dxp
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 50 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
      dcu(1,j+noff) = dcu(1,j+noff) + sdcu(1,j)
      dcu(2,j+noff) = dcu(2,j+noff) + sdcu(2,j)
   50 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
!$OMP ATOMIC
      dcu(1,1+noff) = dcu(1,1+noff) + sdcu(1,1)
!$OMP ATOMIC
      dcu(2,1+noff) = dcu(2,1+noff) + sdcu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
!$OMP ATOMIC
         dcu(1,nn+noff) = dcu(1,nn+noff) + sdcu(1,nn)
!$OMP ATOMIC
         dcu(2,nn+noff) = dcu(2,nn+noff) + sdcu(2,nn)
      endif
   60 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,ci, &
     &qbm,dt,idimp,nppmx,nx,mx,nxv,mx1)
! for 1-2/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 128 flops/particle, 2 divide, 1 sqrt, 16 loads, 12 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pj*gami, where j = y,z, for i = 1, 2
! where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*dvj*gami/dt, where j = y,z, for i = 1, 2
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pj*pk*gami**2, where jk = xy,xz, for i = 1, 2
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
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
! omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami.
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = momentum px of particle n in tile m at t - dt/2
! ppart(3,n,m) = momentum py of particle n in tile m at t - dt/2
! ppart(4,n,m) = momentum pz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! cu(i,j) = ith component of current density at grid point j
! for i = 1, 2
! dcu(i,j) = ith component of acceleration density at grid point j
! for i = 1, 2
! amu(i,j) = ith component of momentum flux at grid point j
! for i = 1, 2
! kpic = number of particles per tile
! omx = magnetic field electron cyclotron frequency in x
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1
      real omx, qm, qbm, dt, ci
      real ppart, fxyz, byz, cu, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, vx, vy, vz, p2, v1, v2, v3
      real sfxyz, sbyz, scu, sdcu, samu
      dimension sfxyz(3,MXV), sbyz(2,MXV)
      dimension scu(2,MXV), sdcu(2,MXV), samu(2,MXV)
!     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
!     dimension scu(2,mx+1), sdcu(2,mx+1), samu(2,mx+1)
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      dti = 1.0/dt
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,vx,vy,vz,acx,&
!$OMP& acy,acz,v1,v2,v3,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,   &
!$OMP& rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,gh,sfxyz,sbyz,samu,sdcu,   &
!$OMP& scu) SCHEDULE(dynamic)
      do 60 k = 1, mx1
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
! zero out local accumulators
      do 30 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   30 continue
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
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
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
! new momentum
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! deposit momentum flux, acceleration density, and current density
! to local accumulator
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      v1 = ox*oy
      v2 = ox*oz
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
      sdcu(1,nn) = sdcu(1,nn) + vy*amx
      sdcu(2,nn) = sdcu(2,nn) + vz*amx
      sdcu(1,nn+1) = sdcu(1,nn+1) + vy*dxp
      sdcu(2,nn+1) = sdcu(2,nn+1) + vz*dxp
      scu(1,nn) = scu(1,nn) + oy*amx
      scu(2,nn) = scu(2,nn) + oz*amx
      scu(1,nn+1) = scu(1,nn+1) + oy*dxp
      scu(2,nn+1) = scu(2,nn+1) + oz*dxp
   40 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 50 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
      dcu(1,j+noff) = dcu(1,j+noff) + sdcu(1,j)
      dcu(2,j+noff) = dcu(2,j+noff) + sdcu(2,j)
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   50 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
!$OMP ATOMIC
      dcu(1,1+noff) = dcu(1,1+noff) + sdcu(1,1)
!$OMP ATOMIC
      dcu(2,1+noff) = dcu(2,1+noff) + sdcu(2,1)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
!$OMP ATOMIC
         dcu(1,nn+noff) = dcu(1,nn+noff) + sdcu(1,nn)
!$OMP ATOMIC
         dcu(2,nn+noff) = dcu(2,nn+noff) + sdcu(2,nn)
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   60 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
! add scaled field to extended periodic field
! linear interpolation
      implicit none
      real dcu, cus, q2m0
      integer nx, nxe
      dimension dcu(2,nxe), cus(2,nxe)
! local data
      integer i, j
      do 20 j = 1, nx
      do 10 i = 1, 2
      dcu(i,j) = dcu(i,j) - q2m0*cus(i,j)
   10 continue
   20 continue
      return
      end
