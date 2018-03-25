!-----------------------------------------------------------------------
! Fortran Library for depositing time derivative of current
! 2-1/2D MPI/OpenMP PIC Codes:
! PPFWPMINMX2 calculates maximum and minimum plasma frequency
! PPFWPTMINMX2 calculates maximum and minimum total plasma frequency
! PPGDJPPOST2L calculates particle momentum flux and acceleration
!              density using linear interpolation
! PPGDCJPPOST2L calculates particle momentum flux, acceleration density
!               and current density using linear interpolation.
! PPGRDJPPOST2L calculates relativistic particle momentum flux and
!               acceleration density using linear interpolation
! PPGRDCJPPOST2L calculates relativistic particle momentum flux,
!                acceleration density and current density using linear
!                interpolation.
! PPASCFGUARD2L add scaled vector field to extended periodic field
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 26, 2018
!-----------------------------------------------------------------------
      subroutine PPFWPMINMX2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
! calculates maximum and minimum plasma frequency.  assumes guard cells
! have already been added
! qe = charge density for electrons
! nyp = number of primary gridpoints in particle partition
! qbme = charge/mass ratio for electrons
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer nyp, nx, nxe, nypmx
      real qbme, wpmax, wpmin
      real qe
      dimension qe(nxe,nypmx)
! local data
      integer j, k
      real at1
      wpmax = qbme*qe(1,1)
      wpmin = wpmax
!$OMP PARALLEL DO PRIVATE(j,k)                                          &
!$OMP& REDUCTION(max:wpmax), REDUCTION(min:wpmin)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      at1 = qbme*qe(j,k)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFWPTMINMX2(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,   &
     &nypmx)
! calculates maximum and minimum total plasma frequency.  assumes guard
! cells have already been added
! qe/qi = charge density for electrons/ions
! nyp = number of primary gridpoints in particle partition
! qbme/qbmi = charge/mass ratio for electrons/ions
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer nyp, nx, nxe, nypmx
      real qbme, qbmi, wpmax, wpmin
      real qe, qi
      dimension qe(nxe,nypmx), qi(nxe,nypmx)
! local data
      integer j, k
      real at1
      wpmax = qbme*qe(1,1) + qbmi*qi(1,1)
      wpmin = wpmax
!$OMP PARALLEL DO PRIVATE(j,k)                                          &
!$OMP& REDUCTION(max:wpmax), REDUCTION(min:wpmin)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      at1 = qbme*qe(j,k) + qbmi*qi(j,k)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,qbm&
     &,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation.
! OpenMP version using guard cells, for distributed data
! data read/written in tiles
! particles stored in segmented array
! 194 flops/particle, 1 divide, 57 loads, 28 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
! dcu(i,n+1,m)=qci*dx*(1.-dy)
! dcu(i,n,m+1)=qci*(1.-dx)*dy
! dcu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! where n,m = nearest grid points and dx = x-n, dy = y-m
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! dcu(i,j,k) = ith component of acceleration density
! at grid point j,kk for i = 1, 3
! amu(i,j,k) = ith component of momentum flux
! at grid point j,kk for i = 1, 4
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real qm, qbm, dt
      real ppart, fxy, bxy, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension dcu(3,nxv,nypmx), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension sdcu(3,MXV,MYV), samu(4,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,vz,v1,v2,v3,&
!$OMP& v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, &
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,&
!$OMP& sdcu,samu)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
   50 continue
   60 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
! deposit momentum flux and acceleration density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
   70 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
      dcu(1,i+noffp,j+moffp) = dcu(1,i+noffp,j+moffp) + sdcu(1,i,j)
      dcu(2,i+noffp,j+moffp) = dcu(2,i+noffp,j+moffp) + sdcu(2,i,j)
      dcu(3,i+noffp,j+moffp) = dcu(3,i+noffp,j+moffp) + sdcu(3,i,j)
   80 continue
   90 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp) = dcu(3,i+noffp,1+moffp) + sdcu(3,i,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp)              &
     & + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp)              &
     & + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp) = dcu(3,i+noffp,mm+moffp)              &
     & + sdcu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp) = dcu(3,1+noffp,j+moffp) + sdcu(3,1,j)
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
!$OMP ATOMIC
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp)              &
     & + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp)              &
     & + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noffp,j+moffp) = dcu(3,nn+noffp,j+moffp)              &
     & + sdcu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,qm&
     &,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation.
! OpenMP version using guard cells, for distributed data
! data read/written in tiles
! particles stored in segmented array
! 218 flops/particle, 1 divide, 69 loads, 40 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*vj, where j = x,y,z, for i = 1, 3
! where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
! dcu(i,n+1,m)=qci*dx*(1.-dy)
! dcu(i,n,m+1)=qci*(1.-dx)*dy
! dcu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! where n,m = nearest grid points and dx = x-n, dy = y-m
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! cu(i,j,k) = ith component of current density
! at grid point j,kk for i = 1, 3
! dcu(i,j,k) = ith component of acceleration density
! at grid point j,kk for i = 1, 3
! amu(i,j,k) = ith component of momentum flux
! at grid point j,kk for i = 1, 4
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real qm, qbm, dt
      real ppart, fxy, bxy, cu, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension cu(3,nxv,nypmx), dcu(3,nxv,nypmx), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension scu(3,MXV,MYV), sdcu(3,MXV,MYV), samu(4,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,vz,v1,v2,v3,&
!$OMP& v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, &
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,&
!$OMP& scu,sdcu,samu)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   50 continue
   60 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      scu(1,nn,mm) = scu(1,nn,mm) + ox*dx
      scu(2,nn,mm) = scu(2,nn,mm) + oy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + oz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + ox*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + oy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + oz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + ox*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + oy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + oz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + ox*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + oy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + oz*dy
   70 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
      dcu(1,i+noffp,j+moffp) = dcu(1,i+noffp,j+moffp) + sdcu(1,i,j)
      dcu(2,i+noffp,j+moffp) = dcu(2,i+noffp,j+moffp) + sdcu(2,i,j)
      dcu(3,i+noffp,j+moffp) = dcu(3,i+noffp,j+moffp) + sdcu(3,i,j)
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   80 continue
   90 continue
! deposit currents to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp) = dcu(3,i+noffp,1+moffp) + sdcu(3,i,1)
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp)              &
     & + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp)              &
     & + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp) = dcu(3,i+noffp,mm+moffp)              &
     & + sdcu(3,i,mm)
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp) = dcu(3,1+noffp,j+moffp) + sdcu(3,1,j)
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
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
!$OMP ATOMIC
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp)              &
     & + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp)              &
     & + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noffp,j+moffp) = dcu(3,nn+noffp,j+moffp)              &
     & + sdcu(3,nn,j)
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,  &
     &qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation
! for relativistic particles.
! OpenMP version using guard cells, for distributed data
! data read/written in tiles
! particles stored in segmented array
! 217 flops/particle, 2 divide, 1 sqrt, 57 loads, 28 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
! dcu(i,n+1,m)=qci*dx*(1.-dy)
! dcu(i,n,m+1)=qci*(1.-dx)*dy
! dcu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! and Ej = jth component of electric field
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
! where n,m = nearest grid points and dx = x-n, dy = y-m
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami.
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! dcu(i,j,k) = ith component of acceleration density
! at grid point j,kk for i = 1, 3
! amu(i,j,k) = ith component of momentum flux
! at grid point j,kk for i = 1, 4
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real qm, qbm, dt, ci
      real ppart, fxy, bxy, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension dcu(3,nxv,nypmx), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension sdcu(3,MXV,MYV), samu(4,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,vz,v1,v2,v3,&
!$OMP& v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, &
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,  &
!$OMP& qtmg,gh,sfxy,sbxy,sdcu,samu)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
   50 continue
   60 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
! deposit momentum flux and acceleration density
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
   70 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
      dcu(1,i+noffp,j+moffp) = dcu(1,i+noffp,j+moffp) + sdcu(1,i,j)
      dcu(2,i+noffp,j+moffp) = dcu(2,i+noffp,j+moffp) + sdcu(2,i,j)
      dcu(3,i+noffp,j+moffp) = dcu(3,i+noffp,j+moffp) + sdcu(3,i,j)
   80 continue
   90 continue
! deposit currents to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp) = dcu(3,i+noffp,1+moffp) + sdcu(3,i,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp)              &
     & + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp)              &
     & + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp) = dcu(3,i+noffp,mm+moffp)              &
     & + sdcu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp) = dcu(3,1+noffp,j+moffp) + sdcu(3,1,j)
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
!$OMP ATOMIC
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp)              &
     & + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp)              &
     & + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noffp,j+moffp) = dcu(3,nn+noffp,j+moffp)              &
     & + sdcu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp, &
     &qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! OpenMP version using guard cells, for distributed data
! data read/written in tiles
! particles stored in segmented array
! 241 flops/particle, 2 divide, 1 sqrt, 69 loads, 40 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
! where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
! dcu(i,n+1,m)=qci*dx*(1.-dy)
! dcu(i,n,m+1)=qci*(1.-dx)*dy
! dcu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! and Ej = jth component of electric field
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
! where n,m = nearest grid points and dx = x-n, dy = y-m
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami.
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! cu(i,j,k) = ith component of current density
! at grid point j,kk for i = 1, 3
! dcu(i,j,k) = ith component of acceleration density
! at grid point j,kk for i = 1, 3
! amu(i,j,k) = ith component of momentum flux
! at grid point j,kk for i = 1, 4
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real qm, qbm, dt, ci
      real ppart, fxy, bxy, cu, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension cu(3,nxv,nypmx), dcu(3,nxv,nypmx), amu(4,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension scu(3,MXV,MYV), sdcu(3,MXV,MYV), samu(4,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,vz,v1,v2,v3,&
!$OMP& v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, &
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,  &
!$OMP& qtmg,gh,sfxy,sbxy,scu,sdcu,samu)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   50 continue
   60 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      scu(1,nn,mm) = scu(1,nn,mm) + ox*dx
      scu(2,nn,mm) = scu(2,nn,mm) + oy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + oz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + ox*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + oy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + oz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + ox*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + oy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + oz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + ox*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + oy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + oz*dy
   70 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)
      amu(4,i+noffp,j+moffp) = amu(4,i+noffp,j+moffp) + samu(4,i,j)
      dcu(1,i+noffp,j+moffp) = dcu(1,i+noffp,j+moffp) + sdcu(1,i,j)
      dcu(2,i+noffp,j+moffp) = dcu(2,i+noffp,j+moffp) + sdcu(2,i,j)
      dcu(3,i+noffp,j+moffp) = dcu(3,i+noffp,j+moffp) + sdcu(3,i,j)
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   80 continue
   90 continue
! deposit currents to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noffp,1+moffp) = amu(4,i+noffp,1+moffp) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp) = dcu(3,i+noffp,1+moffp) + sdcu(3,i,1)
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp)              &
     & + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp)              &
     & + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp) = dcu(3,i+noffp,mm+moffp)              &
     & + sdcu(3,i,mm)
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noffp,j+moffp) = amu(4,1+noffp,j+moffp) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp) = dcu(3,1+noffp,j+moffp) + sdcu(3,1,j)
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
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
!$OMP ATOMIC
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp)              &
     & + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp)              &
     & + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noffp,j+moffp) = dcu(3,nn+noffp,j+moffp)              &
     & + sdcu(3,nn,j)
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
! add scaled vector field to extended periodic field
! linear interpolation, for distributed data
! nyp = number of primary (complete) gridpoints in particle partition
! q2m0 = wp0/affp, where
! wp0 = normalized total plasma frequency squared
! affp = normalization constant = nx*ny/np, where np=number of particles
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, nxe, nypmx
      real dcu, cus, q2m0
      dimension dcu(3,nxe,nypmx), cus(3,nxe,nypmx)
! local data
      integer i, j, k
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 30 k = 1, nyp
      do 20 j = 1, nx
      do 10 i = 1, 3
      dcu(i,j,k) = dcu(i,j,k) - q2m0*cus(i,j,k)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
