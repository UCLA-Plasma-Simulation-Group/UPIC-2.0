!-----------------------------------------------------------------------
! Fortran Library for depositing time derivative of current
! 3D Vector/MPI/OpenMP PIC Codes:
! VPPGDJPPOST32L calculates particle momentum flux and acceleration
!                density using linear interpolation
! VPPGDCJPPOST32L calculates particle momentum flux, acceleration
!                 density, and current density using linear
!                 interpolation.
! VPPGRDJPPOST32L calculates relativistic particle momentum flux and
!                 acceleration density using linear interpolation
! VPPGRDCJPPOST32L calculates relativistic particle momentum flux,
!                  acceleration density and current density using linear
!                  interpolation.
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: may 10, 2018
!-----------------------------------------------------------------------
      subroutine VPPGDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu, &
     &qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
! for 3 code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 350 flops/particle, 1 divide, 126 loads, 72 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
! omz = (q/m)*bz(x(t),y(t),z(t)).
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! dcu(i,j,k,l) = ith component of acceleration density
! at grid point j,k for i = 1, 3
! amu(i,j,k,l) = ith component of momentum flux
! at grid point j,k,l for i = 1, 6
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm, qbm, dt
      real ppart, fxyz, bxyz, dcu, amu
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension dcu(3,nxv*nypmx*nzpmx), amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, sdcu, samu
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
!     dimension sdcu(3,MXV*MYV*MZV), samu(6,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sdcu(3,(mx+1)*(my+1)*(mz+1))
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,mnoff,lnoff&
!$OMP& ,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,&
!$OMP& amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxyz,sbxyz,  &
!$OMP& sdcu,samu,n,s,t)
      do 270 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 70 j = 1, (mx+1)*(my+1)*(mz+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
   70 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 150 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 80 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      t(j,1) = ppart(4,j+joff,l)
      t(j,2) = ppart(5,j+joff,l)
      t(j,3) = ppart(6,j+joff,l)
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j)
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      ox = sbxyz(1,nn)*s(j,1) + sbxyz(1,nn+1)*s(j,2)
      oy = sbxyz(2,nn)*s(j,1) + sbxyz(2,nn+1)*s(j,2)
      oz = sbxyz(3,nn)*s(j,1) + sbxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ox = ox + sbxyz(1,mm)*s(j,3) + sbxyz(1,mm+1)*s(j,4)
      oy = oy + sbxyz(2,mm)*s(j,3) + sbxyz(2,mm+1)*s(j,4)
      oz = oz + sbxyz(3,mm)*s(j,3) + sbxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      acx = sbxyz(1,ll)*s(j,5) + sbxyz(1,ll+1)*s(j,6)
      acy = sbxyz(2,ll)*s(j,5) + sbxyz(2,ll+1)*s(j,6)
      acz = sbxyz(3,ll)*s(j,5) + sbxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      acx = acx + sbxyz(1,nn)*s(j,7) + sbxyz(1,nn+1)*s(j,8)
      acy = acy + sbxyz(2,nn)*s(j,7) + sbxyz(2,nn+1)*s(j,8)
      acz = acz + sbxyz(3,nn)*s(j,7) + sbxyz(3,nn+1)*s(j,8)
      t(j,4) = dx + vx
      t(j,5) = dy + vy
      t(j,6) = dz + vz
      t(j,7) = ox + acx
      t(j,8) = oy + acy
      t(j,9) = oz + acz
   90 continue
! rescale weights for deposit
      do 110 i = 1, lvect
      do 100 j = 1, npblk
      s(j,i) = qm*s(j,i)
  100 continue
  110 continue
! new velocity
      do 120 j = 1, npblk
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
  120 continue
! deposit momentum flux and acceleration density within tile to local
! accumulator
      do 140 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
!dir$ ivdep
      do 130 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
  130 continue
  140 continue
  150 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 160 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      ox = amx*sbxyz(1,nn) + amy*sbxyz(1,nn+1)
      oy = amx*sbxyz(2,nn) + amy*sbxyz(2,nn+1)
      oz = amx*sbxyz(3,nn) + amy*sbxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      ox = amz*(ox + dyp*sbxyz(1,nn+lxv) + dx1*sbxyz(1,nn+1+lxv))
      oy = amz*(oy + dyp*sbxyz(2,nn+lxv) + dx1*sbxyz(2,nn+1+lxv))
      oz = amz*(oz + dyp*sbxyz(3,nn+lxv) + dx1*sbxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      acx = amx*sbxyz(1,mm) + amy*sbxyz(1,mm+1)
      acy = amx*sbxyz(2,mm) + amy*sbxyz(2,mm+1)
      acz = amx*sbxyz(3,mm) + amy*sbxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
      ox = ox + dzp*(acx + dyp*sbxyz(1,mm+lxv) + dx1*sbxyz(1,mm+1+lxv))
      oy = oy + dzp*(acy + dyp*sbxyz(2,mm+lxv) + dx1*sbxyz(2,mm+1+lxv))
      oz = oz + dzp*(acz + dyp*sbxyz(3,mm+lxv) + dx1*sbxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
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
      amz = qm*amz
      dzp = qm*dzp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      sdcu(1,mm) = sdcu(1,mm) + vx*dx
      sdcu(2,mm) = sdcu(2,mm) + vy*dx
      sdcu(3,mm) = sdcu(3,mm) + vz*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      sdcu(1,mm+1) = sdcu(1,mm+1) + vx*dy
      sdcu(2,mm+1) = sdcu(2,mm+1) + vy*dy
      sdcu(3,mm+1) = sdcu(3,mm+1) + vz*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      sdcu(1,mm+lxv) = sdcu(1,mm+lxv) + vx*dx
      sdcu(2,mm+lxv) = sdcu(2,mm+lxv) + vy*dx
      sdcu(3,mm+lxv) = sdcu(3,mm+lxv) + vz*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
      sdcu(1,mm+1+lxv) = sdcu(1,mm+1+lxv) + vx*dy
      sdcu(2,mm+1+lxv) = sdcu(2,mm+1+lxv) + vy*dy
      sdcu(3,mm+1+lxv) = sdcu(3,mm+1+lxv) + vz*dy
  160 continue
! deposit momentum flux and acceleration density to interior points in
! global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 190 k = 2, ll
      do 180 j = 2, mm
!dir$ ivdep
      do 170 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(1,i+lxv*(j-1)+lxyv*(k-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(2,i+lxv*(j-1)+lxyv*(k-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(3,i+lxv*(j-1)+lxyv*(k-1))
  170 continue
  180 continue
  190 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      lm = min(mz+1,nzpmx-loffp)
      do 210 j = 2, mm
      do 200 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(1,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(2,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  200 continue
  210 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 220 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  220 continue
      do 230 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 250 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  250    continue
         do 260 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  260    continue
      endif
  270 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu, &
     &amu,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,   &
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 398 flops/particle, 1 divide, 150 loads, 96 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*vj, where j = x,y,z, for i = 1, 3
! where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
! where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
! and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
! omz = (q/m)*bz(x(t),y(t),z(t)).
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! cu(i,j,k,l) = ith component of current density
! at grid point j,k for i = 1, 3
! dcu(i,j,k,l) = ith component of acceleration density
! at grid point j,k for i = 1, 3
! amu(i,j,k,l) = ith component of momentum flux
! at grid point j,k,l for i = 1, 6
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm, qbm, dt
      real ppart, fxyz, bxyz, cu, dcu, amu
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension cu(3,nxv*nypmx*nzpmx), dcu(3,nxv*nypmx*nzpmx)
      dimension amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, scu, sdcu, samu
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
!     dimension scu(3,MXV*MYV*MZV), sdcu(3,MXV*MYV*MZV)
!     dimension samu(6,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
      dimension sdcu(3,(mx+1)*(my+1)*(mz+1))
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,mnoff,lnoff&
!$OMP& ,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,&
!$OMP& amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxyz,sbxyz,  &
!$OMP& scu,sdcu,samu,n,s,t)
      do 270 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 70 j = 1, (mx+1)*(my+1)*(mz+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   70 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 150 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 80 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      t(j,1) = ppart(4,j+joff,l)
      t(j,2) = ppart(5,j+joff,l)
      t(j,3) = ppart(6,j+joff,l)
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j)
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      ox = sbxyz(1,nn)*s(j,1) + sbxyz(1,nn+1)*s(j,2)
      oy = sbxyz(2,nn)*s(j,1) + sbxyz(2,nn+1)*s(j,2)
      oz = sbxyz(3,nn)*s(j,1) + sbxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ox = ox + sbxyz(1,mm)*s(j,3) + sbxyz(1,mm+1)*s(j,4)
      oy = oy + sbxyz(2,mm)*s(j,3) + sbxyz(2,mm+1)*s(j,4)
      oz = oz + sbxyz(3,mm)*s(j,3) + sbxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      acx = sbxyz(1,ll)*s(j,5) + sbxyz(1,ll+1)*s(j,6)
      acy = sbxyz(2,ll)*s(j,5) + sbxyz(2,ll+1)*s(j,6)
      acz = sbxyz(3,ll)*s(j,5) + sbxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      acx = acx + sbxyz(1,nn)*s(j,7) + sbxyz(1,nn+1)*s(j,8)
      acy = acy + sbxyz(2,nn)*s(j,7) + sbxyz(2,nn+1)*s(j,8)
      acz = acz + sbxyz(3,nn)*s(j,7) + sbxyz(3,nn+1)*s(j,8)
      t(j,4) = dx + vx
      t(j,5) = dy + vy
      t(j,6) = dz + vz
      t(j,7) = ox + acx
      t(j,8) = oy + acy
      t(j,9) = oz + acz
   90 continue
! rescale weights for deposit
      do 110 i = 1, lvect
      do 100 j = 1, npblk
      s(j,i) = qm*s(j,i)
  100 continue
  110 continue
! new velocity
      do 120 j = 1, npblk
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
  120 continue
! deposit momentum flux, acceleration density, and current density
! within tile to local accumulator
      do 140 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
!dir$ ivdep
      do 130 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + ox*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + oy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + oz*s(j,i)
  130 continue
  140 continue
  150 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 160 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      ox = amx*sbxyz(1,nn) + amy*sbxyz(1,nn+1)
      oy = amx*sbxyz(2,nn) + amy*sbxyz(2,nn+1)
      oz = amx*sbxyz(3,nn) + amy*sbxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      ox = amz*(ox + dyp*sbxyz(1,nn+lxv) + dx1*sbxyz(1,nn+1+lxv))
      oy = amz*(oy + dyp*sbxyz(2,nn+lxv) + dx1*sbxyz(2,nn+1+lxv))
      oz = amz*(oz + dyp*sbxyz(3,nn+lxv) + dx1*sbxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      acx = amx*sbxyz(1,mm) + amy*sbxyz(1,mm+1)
      acy = amx*sbxyz(2,mm) + amy*sbxyz(2,mm+1)
      acz = amx*sbxyz(3,mm) + amy*sbxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
      ox = ox + dzp*(acx + dyp*sbxyz(1,mm+lxv) + dx1*sbxyz(1,mm+1+lxv))
      oy = oy + dzp*(acy + dyp*sbxyz(2,mm+lxv) + dx1*sbxyz(2,mm+1+lxv))
      oz = oz + dzp*(acz + dyp*sbxyz(3,mm+lxv) + dx1*sbxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
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
      amz = qm*amz
      dzp = qm*dzp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      scu(1,nn) = scu(1,nn) + ox*dx
      scu(2,nn) = scu(2,nn) + oy*dx
      scu(3,nn) = scu(3,nn) + oz*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      scu(1,nn+1) = scu(1,nn+1) + ox*dy
      scu(2,nn+1) = scu(2,nn+1) + oy*dy
      scu(3,nn+1) = scu(3,nn+1) + oz*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      scu(1,nn+lxv) = scu(1,nn+lxv) + ox*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + oy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + oz*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + ox*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + oy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + oz*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      sdcu(1,mm) = sdcu(1,mm) + vx*dx
      sdcu(2,mm) = sdcu(2,mm) + vy*dx
      sdcu(3,mm) = sdcu(3,mm) + vz*dx
      scu(1,mm) = scu(1,mm) + ox*dx
      scu(2,mm) = scu(2,mm) + oy*dx
      scu(3,mm) = scu(3,mm) + oz*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      sdcu(1,mm+1) = sdcu(1,mm+1) + vx*dy
      sdcu(2,mm+1) = sdcu(2,mm+1) + vy*dy
      sdcu(3,mm+1) = sdcu(3,mm+1) + vz*dy
      scu(1,mm+1) = scu(1,mm+1) + ox*dy
      scu(2,mm+1) = scu(2,mm+1) + oy*dy
      scu(3,mm+1) = scu(3,mm+1) + oz*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      sdcu(1,mm+lxv) = sdcu(1,mm+lxv) + vx*dx
      sdcu(2,mm+lxv) = sdcu(2,mm+lxv) + vy*dx
      sdcu(3,mm+lxv) = sdcu(3,mm+lxv) + vz*dx
      scu(1,mm+lxv) = scu(1,mm+lxv) + ox*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + oy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + oz*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
      sdcu(1,mm+1+lxv) = sdcu(1,mm+1+lxv) + vx*dy
      sdcu(2,mm+1+lxv) = sdcu(2,mm+1+lxv) + vy*dy
      sdcu(3,mm+1+lxv) = sdcu(3,mm+1+lxv) + vz*dy
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + ox*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + oy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + oz*dy
  160 continue
! deposit momentum flux, acceleration density, and current density to
! interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 190 k = 2, ll
      do 180 j = 2, mm
!dir$ ivdep
      do 170 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(1,i+lxv*(j-1)+lxyv*(k-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(2,i+lxv*(j-1)+lxyv*(k-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(3,i+lxv*(j-1)+lxyv*(k-1))
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
  170 continue
  180 continue
  190 continue
! deposit momentum flux, acceleration density, and current density to
! edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 210 j = 2, mm
      do 200 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(1,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(2,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(3,i+lxv*(j-1))
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  200 continue
  210 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 220 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  220 continue
      do 230 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 250 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  250    continue
         do 260 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  260    continue
      endif
  270 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,&
     &qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,    &
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation
! for relativistic particles.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
! 384 flops/particle, 2 divides, 1 sqrt, 126 loads, 72 stores
! input: all, output: dcu, amu
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! and Ej = jth component of electric field
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*pj*pk*gami**2, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! dcu(i,j,k,l) = ith component of acceleration density
! at grid point j,k for i = 1, 3
! amu(i,j,k,l) = ith component of momentum flux
! at grid point j,k,l for i = 1, 6
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm, qbm, dt, ci
      real ppart, fxyz, bxyz, dcu, amu
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension dcu(3,nxv*nypmx*nzpmx), amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, dzp
      real amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, p2, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, sdcu, samu
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
!     dimension sdcu(3,MXV*MYV*MZV), samu(6,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sdcu(3,(mx+1)*(my+1)*(mz+1))
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,mnoff,lnoff&
!$OMP& ,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,&
!$OMP& amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg, &
!$OMP& gh,sfxyz,sbxyz,sdcu,samu,n,s,t)
      do 270 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 70 j = 1, (mx+1)*(my+1)*(mz+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
   70 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 150 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 80 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      t(j,1) = ppart(4,j+joff,l)
      t(j,2) = ppart(5,j+joff,l)
      t(j,3) = ppart(6,j+joff,l)
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j)
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      ox = sbxyz(1,nn)*s(j,1) + sbxyz(1,nn+1)*s(j,2)
      oy = sbxyz(2,nn)*s(j,1) + sbxyz(2,nn+1)*s(j,2)
      oz = sbxyz(3,nn)*s(j,1) + sbxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ox = ox + sbxyz(1,mm)*s(j,3) + sbxyz(1,mm+1)*s(j,4)
      oy = oy + sbxyz(2,mm)*s(j,3) + sbxyz(2,mm+1)*s(j,4)
      oz = oz + sbxyz(3,mm)*s(j,3) + sbxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      acx = sbxyz(1,ll)*s(j,5) + sbxyz(1,ll+1)*s(j,6)
      acy = sbxyz(2,ll)*s(j,5) + sbxyz(2,ll+1)*s(j,6)
      acz = sbxyz(3,ll)*s(j,5) + sbxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      acx = acx + sbxyz(1,nn)*s(j,7) + sbxyz(1,nn+1)*s(j,8)
      acy = acy + sbxyz(2,nn)*s(j,7) + sbxyz(2,nn+1)*s(j,8)
      acz = acz + sbxyz(3,nn)*s(j,7) + sbxyz(3,nn+1)*s(j,8)
      t(j,4) = dx + vx
      t(j,5) = dy + vy
      t(j,6) = dz + vz
      t(j,7) = ox + acx
      t(j,8) = oy + acy
      t(j,9) = oz + acz
   90 continue
! rescale weights for deposit
      do 110 i = 1, lvect
      do 100 j = 1, npblk
      s(j,i) = qm*s(j,i)
  100 continue
  110 continue
! new momentum
      do 120 j = 1, npblk
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
  120 continue
! deposit momentum flux and acceleration density within tile to local
! accumulator
      do 140 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
!dir$ ivdep
      do 130 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
  130 continue
  140 continue
  150 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 160 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      ox = amx*sbxyz(1,nn) + amy*sbxyz(1,nn+1)
      oy = amx*sbxyz(2,nn) + amy*sbxyz(2,nn+1)
      oz = amx*sbxyz(3,nn) + amy*sbxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      ox = amz*(ox + dyp*sbxyz(1,nn+lxv) + dx1*sbxyz(1,nn+1+lxv))
      oy = amz*(oy + dyp*sbxyz(2,nn+lxv) + dx1*sbxyz(2,nn+1+lxv))
      oz = amz*(oz + dyp*sbxyz(3,nn+lxv) + dx1*sbxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      acx = amx*sbxyz(1,mm) + amy*sbxyz(1,mm+1)
      acy = amx*sbxyz(2,mm) + amy*sbxyz(2,mm+1)
      acz = amx*sbxyz(3,mm) + amy*sbxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
      ox = ox + dzp*(acx + dyp*sbxyz(1,mm+lxv) + dx1*sbxyz(1,mm+1+lxv))
      oy = oy + dzp*(acy + dyp*sbxyz(2,mm+lxv) + dx1*sbxyz(2,mm+1+lxv))
      oz = oz + dzp*(acz + dyp*sbxyz(3,mm+lxv) + dx1*sbxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
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
      amz = qm*amz
      dzp = qm*dzp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      sdcu(1,mm) = sdcu(1,mm) + vx*dx
      sdcu(2,mm) = sdcu(2,mm) + vy*dx
      sdcu(3,mm) = sdcu(3,mm) + vz*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      sdcu(1,mm+1) = sdcu(1,mm+1) + vx*dy
      sdcu(2,mm+1) = sdcu(2,mm+1) + vy*dy
      sdcu(3,mm+1) = sdcu(3,mm+1) + vz*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      sdcu(1,mm+lxv) = sdcu(1,mm+lxv) + vx*dx
      sdcu(2,mm+lxv) = sdcu(2,mm+lxv) + vy*dx
      sdcu(3,mm+lxv) = sdcu(3,mm+lxv) + vz*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
      sdcu(1,mm+1+lxv) = sdcu(1,mm+1+lxv) + vx*dy
      sdcu(2,mm+1+lxv) = sdcu(2,mm+1+lxv) + vy*dy
      sdcu(3,mm+1+lxv) = sdcu(3,mm+1+lxv) + vz*dy
  160 continue
! deposit momentum flux and acceleration density to interior points in
! global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 190 k = 2, ll
      do 180 j = 2, mm
!dir$ ivdep
      do 170 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(1,i+lxv*(j-1)+lxyv*(k-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(2,i+lxv*(j-1)+lxyv*(k-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(3,i+lxv*(j-1)+lxyv*(k-1))
  170 continue
  180 continue
  190 continue
! deposit momentum flux and acceleration density to edge points in
! global array
      lm = min(mz+1,nzpmx-loffp)
      do 210 j = 2, mm
      do 200 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(1,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(2,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  200 continue
  210 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 220 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  220 continue
      do 230 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 250 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  250    continue
         do 260 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  260    continue
      endif
  270 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,&
     &amu,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,&
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented arra
! 432 flops/particle, 2 divides, 1 sqrt, 150 loads, 96 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
! where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! and Ej = jth component of electric field
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! and qci = qm*pj*pk*gami**2, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! cu(i,j,k,l) = ith component of current density
! at grid point j,k for i = 1, 3
! dcu(i,j,k,l) = ith component of acceleration density
! at grid point j,k for i = 1, 3
! amu(i,j,k,l) = ith component of momentum flux
! at grid point j,k,l for i = 1, 6
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm, qbm, dt, ci
      real ppart, fxyz, bxyz, cu, dcu, amu
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension cu(3,nxv*nypmx*nzpmx), dcu(3,nxv*nypmx*nzpmx)
      dimension amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, dzp
      real amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, p2, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, scu, sdcu, samu
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
!     dimension sdcu(3,MXV*MYV*MZV), samu(6,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
      dimension sdcu(3,(mx+1)*(my+1)*(mz+1))
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,mnoff,lnoff&
!$OMP& ,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,&
!$OMP& amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg, &
!$OMP& gh,sfxyz,sbxyz,scu,sdcu,samu,n,s,t)
      do 270 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
!dir$ ivdep
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 70 j = 1, (mx+1)*(my+1)*(mz+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
      sdcu(1,j) = 0.0
      sdcu(2,j) = 0.0
      sdcu(3,j) = 0.0
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   70 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 150 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 80 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      t(j,1) = ppart(4,j+joff,l)
      t(j,2) = ppart(5,j+joff,l)
      t(j,3) = ppart(6,j+joff,l)
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j)
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      ox = sbxyz(1,nn)*s(j,1) + sbxyz(1,nn+1)*s(j,2)
      oy = sbxyz(2,nn)*s(j,1) + sbxyz(2,nn+1)*s(j,2)
      oz = sbxyz(3,nn)*s(j,1) + sbxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ox = ox + sbxyz(1,mm)*s(j,3) + sbxyz(1,mm+1)*s(j,4)
      oy = oy + sbxyz(2,mm)*s(j,3) + sbxyz(2,mm+1)*s(j,4)
      oz = oz + sbxyz(3,mm)*s(j,3) + sbxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      acx = sbxyz(1,ll)*s(j,5) + sbxyz(1,ll+1)*s(j,6)
      acy = sbxyz(2,ll)*s(j,5) + sbxyz(2,ll+1)*s(j,6)
      acz = sbxyz(3,ll)*s(j,5) + sbxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      acx = acx + sbxyz(1,nn)*s(j,7) + sbxyz(1,nn+1)*s(j,8)
      acy = acy + sbxyz(2,nn)*s(j,7) + sbxyz(2,nn+1)*s(j,8)
      acz = acz + sbxyz(3,nn)*s(j,7) + sbxyz(3,nn+1)*s(j,8)
      t(j,4) = dx + vx
      t(j,5) = dy + vy
      t(j,6) = dz + vz
      t(j,7) = ox + acx
      t(j,8) = oy + acy
      t(j,9) = oz + acz
   90 continue
! rescale weights for deposit
      do 110 i = 1, lvect
      do 100 j = 1, npblk
      s(j,i) = qm*s(j,i)
  100 continue
  110 continue
! new momentum
      do 120 j = 1, npblk
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
  120 continue
! deposit momentum flux, acceleration density, and current density
! within tile to local accumulator
      do 140 j = 1, npblk
      ox = t(j,1)
      oy = t(j,2)
      oz = t(j,3)
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
!dir$ ivdep
      do 130 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
      sdcu(1,n(j)+mn(i)) = sdcu(1,n(j)+mn(i)) + vx*s(j,i)
      sdcu(2,n(j)+mn(i)) = sdcu(2,n(j)+mn(i)) + vy*s(j,i)
      sdcu(3,n(j)+mn(i)) = sdcu(3,n(j)+mn(i)) + vz*s(j,i)
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + ox*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + oy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + oz*s(j,i)
  130 continue
  140 continue
  150 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 160 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      ox = amx*sbxyz(1,nn) + amy*sbxyz(1,nn+1)
      oy = amx*sbxyz(2,nn) + amy*sbxyz(2,nn+1)
      oz = amx*sbxyz(3,nn) + amy*sbxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      ox = amz*(ox + dyp*sbxyz(1,nn+lxv) + dx1*sbxyz(1,nn+1+lxv))
      oy = amz*(oy + dyp*sbxyz(2,nn+lxv) + dx1*sbxyz(2,nn+1+lxv))
      oz = amz*(oz + dyp*sbxyz(3,nn+lxv) + dx1*sbxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      acx = amx*sbxyz(1,mm) + amy*sbxyz(1,mm+1)
      acy = amx*sbxyz(2,mm) + amy*sbxyz(2,mm+1)
      acz = amx*sbxyz(3,mm) + amy*sbxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
      ox = ox + dzp*(acx + dyp*sbxyz(1,mm+lxv) + dx1*sbxyz(1,mm+1+lxv))
      oy = oy + dzp*(acy + dyp*sbxyz(2,mm+lxv) + dx1*sbxyz(2,mm+1+lxv))
      oz = oz + dzp*(acz + dyp*sbxyz(3,mm+lxv) + dx1*sbxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
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
      amz = qm*amz
      dzp = qm*dzp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx
      sdcu(3,nn) = sdcu(3,nn) + vz*dx
      scu(1,nn) = scu(1,nn) + ox*dx
      scu(2,nn) = scu(2,nn) + oy*dx
      scu(3,nn) = scu(3,nn) + oz*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy
      sdcu(3,nn+1) = sdcu(3,nn+1) + vz*dy
      scu(1,nn+1) = scu(1,nn+1) + ox*dy
      scu(2,nn+1) = scu(2,nn+1) + oy*dy
      scu(3,nn+1) = scu(3,nn+1) + oz*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx
      sdcu(3,nn+lxv) = sdcu(3,nn+lxv) + vz*dx
      scu(1,nn+lxv) = scu(1,nn+lxv) + ox*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + oy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + oz*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy
      sdcu(3,nn+1+lxv) = sdcu(3,nn+1+lxv) + vz*dy
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + ox*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + oy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + oz*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      sdcu(1,mm) = sdcu(1,mm) + vx*dx
      sdcu(2,mm) = sdcu(2,mm) + vy*dx
      sdcu(3,mm) = sdcu(3,mm) + vz*dx
      scu(1,mm) = scu(1,mm) + ox*dx
      scu(2,mm) = scu(2,mm) + oy*dx
      scu(3,mm) = scu(3,mm) + oz*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      sdcu(1,mm+1) = sdcu(1,mm+1) + vx*dy
      sdcu(2,mm+1) = sdcu(2,mm+1) + vy*dy
      sdcu(3,mm+1) = sdcu(3,mm+1) + vz*dy
      scu(1,mm+1) = scu(1,mm+1) + ox*dy
      scu(2,mm+1) = scu(2,mm+1) + oy*dy
      scu(3,mm+1) = scu(3,mm+1) + oz*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      sdcu(1,mm+lxv) = sdcu(1,mm+lxv) + vx*dx
      sdcu(2,mm+lxv) = sdcu(2,mm+lxv) + vy*dx
      sdcu(3,mm+lxv) = sdcu(3,mm+lxv) + vz*dx
      scu(1,mm+lxv) = scu(1,mm+lxv) + ox*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + oy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + oz*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
      sdcu(1,mm+1+lxv) = sdcu(1,mm+1+lxv) + vx*dy
      sdcu(2,mm+1+lxv) = sdcu(2,mm+1+lxv) + vy*dy
      sdcu(3,mm+1+lxv) = sdcu(3,mm+1+lxv) + vz*dy
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + ox*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + oy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + oz*dy
  160 continue
! deposit momentum flux, acceleration density, and current density to
! interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 190 k = 2, ll
      do 180 j = 2, mm
!dir$ ivdep
      do 170 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(1,i+lxv*(j-1)+lxyv*(k-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(2,i+lxv*(j-1)+lxyv*(k-1))
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & sdcu(3,i+lxv*(j-1)+lxyv*(k-1))
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
  170 continue
  180 continue
  190 continue
! deposit momentum flux, acceleration density, and current density to
! edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 210 j = 2, mm
      do 200 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(1,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(2,i+lxv*(j-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sdcu(3,i+lxv*(j-1))
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  200 continue
  210 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 220 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &dcu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sdcu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  220 continue
      do 230 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ sdcu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + sdcu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 250 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   dcu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + sdcu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  250    continue
         do 260 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   dcu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + sdcu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      dcu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + sdcu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  260    continue
      endif
  270 continue
!$OMP END PARALLEL DO
      return
      end
