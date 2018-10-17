!-----------------------------------------------------------------------
! Fortran Library for depositing time derivative of current
! 2-1/2D Vector/MPI/OpenMP PIC Codes:
! VPPGDJPPOST2L calculates particle momentum flux and acceleration
!              density using linear interpolation
! VPPGDCJPPOST2L calculates particle momentum flux, acceleration density
!               and current density using linear interpolation.
! VPPGRDJPPOST2L calculates relativistic particle momentum flux and
!               acceleration density using linear interpolation
! VPPGRDCJPPOST2L calculates relativistic particle momentum flux,
!                acceleration density and current density using linear
!                interpolation.
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 25, 2018
!-----------------------------------------------------------------------
      subroutine VPPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,  &
     &qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation.
! vectorizable/OpenMP version using guard cells, for distributed data
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
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension dcu(3,nxv*nypmx), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, sdcu, samu
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
!     dimension sdcu(3,MXV*MYV), samu(4,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
      dimension sdcu(3,(mx+1)*(my+1)), samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,v1,v2,v3,v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,sfxy,sbxy,sdcu,samu,n,s,t)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm
!dir$ ivdep
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! zero out local accumulators
      do 50 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
   50 continue
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      t(j,1) = ppart(3,j+joff,k)
      t(j,2) = ppart(4,j+joff,k)
      t(j,3) = ppart(5,j+joff,k)
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      dz = sfxy(3,nn)*s(j,1) + sfxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      dz = dz + sfxy(3,mm)*s(j,3) + sfxy(3,mm+1)*s(j,4)
      ox = sbxy(1,nn)*s(j,1) + sbxy(1,nn+1)*s(j,2)
      oy = sbxy(2,nn)*s(j,1) + sbxy(2,nn+1)*s(j,2)
      oz = sbxy(3,nn)*s(j,1) + sbxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      ox = ox + sbxy(1,mm)*s(j,3) + sbxy(1,mm+1)*s(j,4)
      oy = oy + sbxy(2,mm)*s(j,3) + sbxy(2,mm+1)*s(j,4)
      oz = oz + sbxy(3,mm)*s(j,3) + sbxy(3,mm+1)*s(j,4)
      t(j,4) = dx
      t(j,5) = dy
      t(j,6) = dz
      t(j,7) = ox
      t(j,8) = oy
      t(j,9) = oz
   70 continue
! rescale weights for deposit
      do 90 i = 1, lvect
      do 80 j = 1, npblk
      s(j,i) = qm*s(j,i)
   80 continue
   90 continue
! new velocity
      do 100 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
! calculate half impulse
      dx = qtmh*t(j,4)
      dy = qtmh*t(j,5)
      dz = qtmh*t(j,6)
! half acceleration
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
! calculate cyclotron frequency
      omxt = qtmh*t(j,7)
      omyt = qtmh*t(j,8)
      omzt = qtmh*t(j,9)
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
! new velocity sums and differences
      t(j,1) = 0.5*(dx + vx)
      t(j,2) = 0.5*(dy + vy)
      t(j,3) = 0.5*(dz + vz)
      t(j,4) = dti*(dx - vx)
      t(j,5) = dti*(dy - vy)
      t(j,6) = dti*(dz - vz)
  100 continue
! deposit momentum flux and acceleration density within tile to local
! accumulator
      do 120 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
!dir$ ivdep
      do 110 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
  110 continue
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dz = amy*(amx*sfxy(3,nn) + dxp*sfxy(3,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
      dz = dz + dyp*(amx*sfxy(3,nn+lxv) + dxp*sfxy(3,nn+1+lxv))
! find magnetic field
      ox = amy*(amx*sbxy(1,nn) + dxp*sbxy(1,nn+1))
      oy = amy*(amx*sbxy(2,nn) + dxp*sbxy(2,nn+1))
      oz = amy*(amx*sbxy(3,nn) + dxp*sbxy(3,nn+1))
      ox = ox + dyp*(amx*sbxy(1,nn+lxv) + dxp*sbxy(1,nn+1+lxv)) 
      oy = oy + dyp*(amx*sbxy(2,nn+lxv) + dxp*sbxy(2,nn+1+lxv))
      oz = oz + dyp*(amx*sbxy(3,nn+lxv) + dxp*sbxy(3,nn+1+lxv))
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
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
  140 continue
! deposit momentum flux and acceleration density to interior points in
! global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 160 j = 2, mm
!dir$ ivdep
      do 150 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)) = dcu(1,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(1,i+lxv*(j-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)) = dcu(2,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(2,i+lxv*(j-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)) = dcu(3,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(3,i+lxv*(j-1))
  150 continue
  160 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      mm = min(my+1,nypmx-moffp)
      do 170 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp) = dcu(1,i+noffp+nxv*moffp) + sdcu(1,i)
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp) = dcu(2,i+noffp+nxv*moffp) + sdcu(2,i)
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp) = dcu(3,i+noffp+nxv*moffp) + sdcu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)) + sdcu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)) + sdcu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)) + sdcu(3,i+lxv*(mm-1))
      endif
  170 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)) = dcu(1,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(1,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)) = dcu(2,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(2,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)) = dcu(3,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(1,nn+noffp+nxv*(j+moffp-1)) + sdcu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(2,nn+noffp+nxv*(j+moffp-1)) + sdcu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(3,nn+noffp+nxv*(j+moffp-1)) + sdcu(3,nn+lxv*(j-1))
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp, &
     &qm,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation.
! vectorizable/OpenMP version using guard cells, for distributed data
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
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension cu(3,nxv*nypmx), dcu(3,nxv*nypmx), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
!     dimension scu(3,MXV*MYV), sdcu(3,MXV*MYV), samu(4,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
      dimension scu(3,(mx+1)*(my+1)), sdcu(3,(mx+1)*(my+1))
      dimension samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,v1,v2,v3,v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,sfxy,sbxy,scu,sdcu,samu,n,s,t)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm
!dir$ ivdep
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! zero out local accumulators
      do 50 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   50 continue
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      t(j,1) = ppart(3,j+joff,k)
      t(j,2) = ppart(4,j+joff,k)
      t(j,3) = ppart(5,j+joff,k)
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      dz = sfxy(3,nn)*s(j,1) + sfxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      dz = dz + sfxy(3,mm)*s(j,3) + sfxy(3,mm+1)*s(j,4)
      ox = sbxy(1,nn)*s(j,1) + sbxy(1,nn+1)*s(j,2)
      oy = sbxy(2,nn)*s(j,1) + sbxy(2,nn+1)*s(j,2)
      oz = sbxy(3,nn)*s(j,1) + sbxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      ox = ox + sbxy(1,mm)*s(j,3) + sbxy(1,mm+1)*s(j,4)
      oy = oy + sbxy(2,mm)*s(j,3) + sbxy(2,mm+1)*s(j,4)
      oz = oz + sbxy(3,mm)*s(j,3) + sbxy(3,mm+1)*s(j,4)
      t(j,4) = dx
      t(j,5) = dy
      t(j,6) = dz
      t(j,7) = ox
      t(j,8) = oy
      t(j,9) = oz
   70 continue
! rescale weights for deposit
      do 90 i = 1, lvect
      do 80 j = 1, npblk
      s(j,i) = qm*s(j,i)
   80 continue
   90 continue
! new velocity
      do 100 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
! calculate half impulse
      dx = qtmh*t(j,4)
      dy = qtmh*t(j,5)
      dz = qtmh*t(j,6)
! half acceleration
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
! calculate cyclotron frequency
      omxt = qtmh*t(j,7)
      omyt = qtmh*t(j,8)
      omzt = qtmh*t(j,9)
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
! new velocity sums and differences
      t(j,1) = 0.5*(dx + vx)
      t(j,2) = 0.5*(dy + vy)
      t(j,3) = 0.5*(dz + vz)
      t(j,4) = dti*(dx - vx)
      t(j,5) = dti*(dy - vy)
      t(j,6) = dti*(dz - vz)
  100 continue
! deposit momentum flux, acceleration density, and current density
! within tile to local accumulator
      do 120 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
!dir$ ivdep
      do 110 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + ox*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + oy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + oz*s(j,i)
  110 continue
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dz = amy*(amx*sfxy(3,nn) + dxp*sfxy(3,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
      dz = dz + dyp*(amx*sfxy(3,nn+lxv) + dxp*sfxy(3,nn+1+lxv))
! find magnetic field
      ox = amy*(amx*sbxy(1,nn) + dxp*sbxy(1,nn+1))
      oy = amy*(amx*sbxy(2,nn) + dxp*sbxy(2,nn+1))
      oz = amy*(amx*sbxy(3,nn) + dxp*sbxy(3,nn+1))
      ox = ox + dyp*(amx*sbxy(1,nn+lxv) + dxp*sbxy(1,nn+1+lxv)) 
      oy = oy + dyp*(amx*sbxy(2,nn+lxv) + dxp*sbxy(2,nn+1+lxv))
      oz = oz + dyp*(amx*sbxy(3,nn+lxv) + dxp*sbxy(3,nn+1+lxv))
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
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      scu(1,nn) = scu(1,nn) + ox*dx
      scu(2,nn) = scu(2,nn) + oy*dx
      scu(3,nn) = scu(3,nn) + oz*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      scu(1,nn+1) = scu(1,nn+1) + ox*dy
      scu(2,nn+1) = scu(2,nn+1) + oy*dy
      scu(3,nn+1) = scu(3,nn+1) + oz*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      scu(1,nn+lxv) = scu(1,nn+lxv) + ox*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + oy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + oz*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + ox*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + oy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + oz*dy
  140 continue
! deposit momentum flux, acceleration density, and current density to
! interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 160 j = 2, mm
!dir$ ivdep
      do 150 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)) = dcu(1,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(1,i+lxv*(j-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)) = dcu(2,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(2,i+lxv*(j-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)) = dcu(3,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(3,i+lxv*(j-1))
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
  150 continue
  160 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      mm = min(my+1,nypmx-moffp)
      do 170 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp) = dcu(1,i+noffp+nxv*moffp) + sdcu(1,i)
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp) = dcu(2,i+noffp+nxv*moffp) + sdcu(2,i)
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp) = dcu(3,i+noffp+nxv*moffp) + sdcu(3,i)
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)) + sdcu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)) + sdcu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)) + sdcu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(3,i+lxv*(mm-1))
      endif
  170 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)) = dcu(1,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(1,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)) = dcu(2,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(2,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)) = dcu(3,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(3,1+lxv*(j-1))
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(1,nn+noffp+nxv*(j+moffp-1)) + sdcu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(2,nn+noffp+nxv*(j+moffp-1)) + sdcu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(3,nn+noffp+nxv*(j+moffp-1)) + sdcu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm, &
     &qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation
! for relativistic particles.
! vectorizable/OpenMP version using guard cells, for distributed data
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
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension dcu(3,nxv*nypmx), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, sdcu, samu
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
!     dimension sdcu(3,MXV*MYV), samu(4,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
      dimension sdcu(3,(mx+1)*(my+1)), samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,v1,v2,v3,v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,p2,gami,qtmg,gh,sfxy,sbxy,sdcu,samu,n,s,t)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm
!dir$ ivdep
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! zero out local accumulators
      do 50 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
   50 continue
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      t(j,1) = ppart(3,j+joff,k)
      t(j,2) = ppart(4,j+joff,k)
      t(j,3) = ppart(5,j+joff,k)
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      dz = sfxy(3,nn)*s(j,1) + sfxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      dz = dz + sfxy(3,mm)*s(j,3) + sfxy(3,mm+1)*s(j,4)
      ox = sbxy(1,nn)*s(j,1) + sbxy(1,nn+1)*s(j,2)
      oy = sbxy(2,nn)*s(j,1) + sbxy(2,nn+1)*s(j,2)
      oz = sbxy(3,nn)*s(j,1) + sbxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      ox = ox + sbxy(1,mm)*s(j,3) + sbxy(1,mm+1)*s(j,4)
      oy = oy + sbxy(2,mm)*s(j,3) + sbxy(2,mm+1)*s(j,4)
      oz = oz + sbxy(3,mm)*s(j,3) + sbxy(3,mm+1)*s(j,4)
      t(j,4) = dx
      t(j,5) = dy
      t(j,6) = dz
      t(j,7) = ox
      t(j,8) = oy
      t(j,9) = oz
   70 continue
! rescale weights for deposit
      do 90 i = 1, lvect
      do 80 j = 1, npblk
      s(j,i) = qm*s(j,i)
   80 continue
   90 continue
! new momentum
      do 100 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
! calculate half impulse
      dx = qtmh*t(j,4)
      dy = qtmh*t(j,5)
      dz = qtmh*t(j,6)
! half acceleration
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
      omxt = qtmg*t(j,7)
      omyt = qtmg*t(j,8)
      omzt = qtmg*t(j,9)
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
! new velocity sums and differences
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      t(j,1) = ox
      t(j,2) = oy
      t(j,3) = oz
      t(j,4) = qtmg*(vx - ox*gh)
      t(j,5) = qtmg*(vy - oy*gh)
      t(j,6) = qtmg*(vz - oz*gh)
  100 continue
! deposit momentum flux and acceleration density within tile to local
! accumulator
      do 120 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
!dir$ ivdep
      do 110 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
  110 continue
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dz = amy*(amx*sfxy(3,nn) + dxp*sfxy(3,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
      dz = dz + dyp*(amx*sfxy(3,nn+lxv) + dxp*sfxy(3,nn+1+lxv))
! find magnetic field
      ox = amy*(amx*sbxy(1,nn) + dxp*sbxy(1,nn+1))
      oy = amy*(amx*sbxy(2,nn) + dxp*sbxy(2,nn+1))
      oz = amy*(amx*sbxy(3,nn) + dxp*sbxy(3,nn+1))
      ox = ox + dyp*(amx*sbxy(1,nn+lxv) + dxp*sbxy(1,nn+1+lxv)) 
      oy = oy + dyp*(amx*sbxy(2,nn+lxv) + dxp*sbxy(2,nn+1+lxv))
      oz = oz + dyp*(amx*sbxy(3,nn+lxv) + dxp*sbxy(3,nn+1+lxv))
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
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
  140 continue
! deposit momentum flux and acceleration density to interior points in
! global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 160 j = 2, mm
!dir$ ivdep
      do 150 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)) = dcu(1,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(1,i+lxv*(j-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)) = dcu(2,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(2,i+lxv*(j-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)) = dcu(3,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(3,i+lxv*(j-1))
  150 continue
  160 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      mm = min(my+1,nypmx-moffp)
      do 170 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp) = dcu(1,i+noffp+nxv*moffp) + sdcu(1,i)
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp) = dcu(2,i+noffp+nxv*moffp) + sdcu(2,i)
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp) = dcu(3,i+noffp+nxv*moffp) + sdcu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)) + sdcu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)) + sdcu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)) + sdcu(3,i+lxv*(mm-1))
      endif
  170 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)) = dcu(1,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(1,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)) = dcu(2,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(2,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)) = dcu(3,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(1,nn+noffp+nxv*(j+moffp-1)) + sdcu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(2,nn+noffp+nxv*(j+moffp-1)) + sdcu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(3,nn+noffp+nxv*(j+moffp-1)) + sdcu(3,nn+lxv*(j-1))
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp,&
     &qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! vectorizable/OpenMP version using guard cells, for distributed data
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
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension cu(3,nxv*nypmx), dcu(3,nxv*nypmx), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
!     dimension scu(3,MXV*MYV), sdcu(3,MXV*MYV), samu(4,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
      dimension scu(3,(mx+1)*(my+1)), sdcu(3,(mx+1)*(my+1))
      dimension samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,v1,v2,v3,v4,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz, &
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,&
!$OMP& rot9,p2,gami,qtmg,gh,sfxy,sbxy,scu,sdcu,samu,n,s,t)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm
!dir$ ivdep
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! zero out local accumulators
      do 50 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   50 continue
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      t(j,1) = ppart(3,j+joff,k)
      t(j,2) = ppart(4,j+joff,k)
      t(j,3) = ppart(5,j+joff,k)
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      dz = sfxy(3,nn)*s(j,1) + sfxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      dz = dz + sfxy(3,mm)*s(j,3) + sfxy(3,mm+1)*s(j,4)
      ox = sbxy(1,nn)*s(j,1) + sbxy(1,nn+1)*s(j,2)
      oy = sbxy(2,nn)*s(j,1) + sbxy(2,nn+1)*s(j,2)
      oz = sbxy(3,nn)*s(j,1) + sbxy(3,nn+1)*s(j,2)
      mm = nn + lxv
      ox = ox + sbxy(1,mm)*s(j,3) + sbxy(1,mm+1)*s(j,4)
      oy = oy + sbxy(2,mm)*s(j,3) + sbxy(2,mm+1)*s(j,4)
      oz = oz + sbxy(3,mm)*s(j,3) + sbxy(3,mm+1)*s(j,4)
      t(j,4) = dx
      t(j,5) = dy
      t(j,6) = dz
      t(j,7) = ox
      t(j,8) = oy
      t(j,9) = oz
   70 continue
! rescale weights for deposit
      do 90 i = 1, lvect
      do 80 j = 1, npblk
      s(j,i) = qm*s(j,i)
   80 continue
   90 continue
! new momentum
      do 100 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
! calculate half impulse
      dx = qtmh*t(j,4)
      dy = qtmh*t(j,5)
      dz = qtmh*t(j,6)
! half acceleration
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
      omxt = qtmg*t(j,7)
      omyt = qtmg*t(j,8)
      omzt = qtmg*t(j,9)
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
! new velocity sums and differences
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      t(j,1) = ox
      t(j,2) = oy
      t(j,3) = oz
      t(j,4) = qtmg*(vx - ox*gh)
      t(j,5) = qtmg*(vy - oy*gh)
      t(j,6) = qtmg*(vz - oz*gh)
  100 continue
! deposit momentum flux, acceleration density, and current density
! within tile to local accumulator
      do 120 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
!dir$ ivdep
      do 110 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + ox*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + oy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + oz*s(j,i)
  110 continue
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dz = amy*(amx*sfxy(3,nn) + dxp*sfxy(3,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
      dz = dz + dyp*(amx*sfxy(3,nn+lxv) + dxp*sfxy(3,nn+1+lxv))
! find magnetic field
      ox = amy*(amx*sbxy(1,nn) + dxp*sbxy(1,nn+1))
      oy = amy*(amx*sbxy(2,nn) + dxp*sbxy(2,nn+1))
      oz = amy*(amx*sbxy(3,nn) + dxp*sbxy(3,nn+1))
      ox = ox + dyp*(amx*sbxy(1,nn+lxv) + dxp*sbxy(1,nn+1+lxv)) 
      oy = oy + dyp*(amx*sbxy(2,nn+lxv) + dxp*sbxy(2,nn+1+lxv))
      oz = oz + dyp*(amx*sbxy(3,nn+lxv) + dxp*sbxy(3,nn+1+lxv))
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
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      scu(1,nn) = scu(1,nn) + ox*dx
      scu(2,nn) = scu(2,nn) + oy*dx
      scu(3,nn) = scu(3,nn) + oz*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      scu(1,nn+1) = scu(1,nn+1) + ox*dy
      scu(2,nn+1) = scu(2,nn+1) + oy*dy
      scu(3,nn+1) = scu(3,nn+1) + oz*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      scu(1,nn+lxv) = scu(1,nn+lxv) + ox*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + oy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + oz*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + ox*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + oy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + oz*dy
  140 continue
! deposit momentum flux, acceleration density, and current density to
! interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 160 j = 2, mm
!dir$ ivdep
      do 150 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)) = dcu(1,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(1,i+lxv*(j-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)) = dcu(2,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(2,i+lxv*(j-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)) = dcu(3,i+noffp+nxv*(j+moffp-1)) + &
     & sdcu(3,i+lxv*(j-1))
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
  150 continue
  160 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      mm = min(my+1,nypmx-moffp)
      do 170 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp) = dcu(1,i+noffp+nxv*moffp) + sdcu(1,i)
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp) = dcu(2,i+noffp+nxv*moffp) + sdcu(2,i)
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp) = dcu(3,i+noffp+nxv*moffp) + sdcu(3,i)
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)) + sdcu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)) + sdcu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)) + sdcu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     &   + scu(3,i+lxv*(mm-1))
      endif
  170 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)) = dcu(1,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(1,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)) = dcu(2,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(2,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)) = dcu(3,1+noffp+nxv*(j+moffp-1)) + &
     &sdcu(3,1+lxv*(j-1))
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(1,nn+noffp+nxv*(j+moffp-1)) + sdcu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(2,nn+noffp+nxv*(j+moffp-1)) + sdcu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         dcu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   dcu(3,nn+noffp+nxv*(j+moffp-1)) + sdcu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
