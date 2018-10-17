!-----------------------------------------------------------------------
! Fortran Library for depositing time derivative of current
! 3D MPI/OpenMP PIC Codes:
! PPFWPMINMX32 calculates maximum and minimum plasma frequency
! PPFWPTMINMX32 calculates maximum and minimum total plasma frequency
! PPGDJPPOST32L calculates particle momentum flux and acceleration
!               density using linear interpolation
! PPGDCJPPOST32L calculates particle momentum flux, acceleration density
!                and current density using linear interpolation.
! PPGRDJPPOST32L calculates relativistic particle momentum flux and
!                acceleration density using linear interpolation
! PPGRDCJPPOST32L calculates relativistic particle momentum flux,
!                 acceleration density and current density using linear
!                 interpolation.
! MPPASCFGUARD32L add scaled vector field to extended periodic field
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: jfebruary 15, 2018
!-----------------------------------------------------------------------
      subroutine PPFWPMINMX32(qe,nyzp,qbme,wpmax,wpmin,nx,nxe,nypmx,    &
     &nzpmx,idds)
! calculates maximum and minimum plasma frequency.  assumes guard cells
! have already been added
! for distributed data with 2D decomposition
! qe = charge density for electrons
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbme = charge/mass ratio for electrons
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real qbme, wpmax, wpmin
      integer nyzp
      real qe
      dimension qe(nxe,nypmx,nzpmx), nyzp(idds)
! local data
      integer j, k, l
      real at1
      wpmax = qbme*qe(1,1,1)
      wpmin = wpmax
!$OMP PARALLEL DO PRIVATE(j,k,l)                                        &
!$OMP& REDUCTION(max:wpmax), REDUCTION(min:wpmin)
      do 30 l = 1, nyzp(2)
      do 20 k = 1, nyzp(1)
      do 10 j = 1, nx
      at1 = qbme*qe(j,k,l)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFWPTMINMX32(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx,nxe, &
     &nypmx,nzpmx,idds)
! calculates maximum and minimum total plasma frequency.  assumes guard
! cells have already been added
! for distributed data with 2D decomposition
! qe/qi = charge density for electrons/ions
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbme/qbmi = charge/mass ratio for electrons/ions
! wpmax/wpmin = maximum/minimum plasma frequency
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real qbme, qbmi, wpmax, wpmin
      integer nyzp
      real qe, qi
      dimension qe(nxe,nypmx,nzpmx), qi(nxe,nypmx,nzpmx), nyzp(idds)
! local data
      integer j, k, l
      real at1
      wpmax = qbme*qe(1,1,1) + qbmi*qi(1,1,1)
      wpmin = wpmax
!$OMP PARALLEL DO PRIVATE(j,k,l)                                        &
!$OMP& REDUCTION(max:wpmax), REDUCTION(min:wpmin)
      do 30 l = 1, nyzp(2)
      do 20 k = 1, nyzp(1)
      do 10 j = 1, nx
      at1 = qbme*qe(j,k,l) + qbmi*qi(j,k,l)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,qm&
     &,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,  &
     &idds)
! for 3 code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
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
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension dcu(3,nxv,nypmx,nzpmx), amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, sdcu, samu
      dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
      dimension sdcu(3,MXV,MYV,MZV), samu(6,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
!     dimension sdcu(3,mx+1,my+1,mz+1), samu(6,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm&
!$OMP& ,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,&
!$OMP& dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,   &
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxyz,sbxyz,sdcu,samu)
      do 210 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 90 k = 1, mz+1
      do 80 j = 1, my+1
      do 70 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
      sdcu(1,i,j,k) = 0.0
      sdcu(2,i,j,k) = 0.0
      sdcu(3,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
! loop over particles in tile
      do 100 j = 1, nppp
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
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
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
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      sdcu(1,nn,mm,ll) = sdcu(1,nn,mm,ll) + vx*dx
      sdcu(2,nn,mm,ll) = sdcu(2,nn,mm,ll) + vy*dx
      sdcu(3,nn,mm,ll) = sdcu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      sdcu(1,nn+1,mm,ll) = sdcu(1,nn+1,mm,ll) + vx*dy
      sdcu(2,nn+1,mm,ll) = sdcu(2,nn+1,mm,ll) + vy*dy
      sdcu(3,nn+1,mm,ll) = sdcu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      sdcu(1,nn,mm+1,ll) = sdcu(1,nn,mm+1,ll) + vx*dx
      sdcu(2,nn,mm+1,ll) = sdcu(2,nn,mm+1,ll) + vy*dx
      sdcu(3,nn,mm+1,ll) = sdcu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      sdcu(1,nn+1,mm+1,ll) = sdcu(1,nn+1,mm+1,ll) + vx*dy
      sdcu(2,nn+1,mm+1,ll) = sdcu(2,nn+1,mm+1,ll) + vy*dy
      sdcu(3,nn+1,mm+1,ll) = sdcu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      sdcu(1,nn,mm,ll+1) = sdcu(1,nn,mm,ll+1) + vx*dx
      sdcu(2,nn,mm,ll+1) = sdcu(2,nn,mm,ll+1) + vy*dx
      sdcu(3,nn,mm,ll+1) = sdcu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      sdcu(1,nn+1,mm,ll+1) = sdcu(1,nn+1,mm,ll+1) + vx*dy
      sdcu(2,nn+1,mm,ll+1) = sdcu(2,nn+1,mm,ll+1) + vy*dy
      sdcu(3,nn+1,mm,ll+1) = sdcu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      sdcu(1,nn,mm+1,ll+1) = sdcu(1,nn,mm+1,ll+1) + vx*dx
      sdcu(2,nn,mm+1,ll+1) = sdcu(2,nn,mm+1,ll+1) + vy*dx
      sdcu(3,nn,mm+1,ll+1) = sdcu(3,nn,mm+1,ll+1) + vz*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
      sdcu(1,nn+1,mm+1,ll+1) = sdcu(1,nn+1,mm+1,ll+1) + vx*dy
      sdcu(2,nn+1,mm+1,ll+1) = sdcu(2,nn+1,mm+1,ll+1) + vy*dy
      sdcu(3,nn+1,mm+1,ll+1) = sdcu(3,nn+1,mm+1,ll+1) + vz*dy
  100 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
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
      dcu(1,i+noffp,j+moffp,k+loffp) = dcu(1,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,i,j,k)
      dcu(2,i+noffp,j+moffp,k+loffp) = dcu(2,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,i,j,k)
      dcu(3,i+noffp,j+moffp,k+loffp) = dcu(3,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,i,j,k)
  110 continue
  120 continue
  130 continue
! deposit currents to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 150 j = 2, mm
      do 140 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,j+moffp,1+loffp) = dcu(1,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(1,i,j,1)
!$OMP ATOMIC
      dcu(2,i+noffp,j+moffp,1+loffp) = dcu(2,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(2,i,j,1)
!$OMP ATOMIC
      dcu(3,i+noffp,j+moffp,1+loffp) = dcu(3,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(3,i,j,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,j+moffp,lm+loffp) + sdcu(1,i,j,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,j+moffp,lm+loffp) + sdcu(2,i,j,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,j+moffp,lm+loffp) + sdcu(3,i,j,lm)
      endif
  140 continue
  150 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 180 k = 1, ll
      do 160 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp,k+loffp) = dcu(1,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(1,i,1,k)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp,k+loffp) = dcu(2,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(2,i,1,k)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp,k+loffp) = dcu(3,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(3,i,1,k)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(1,i+noffp,mm+moffp,k+loffp) + sdcu(1,i,mm,k)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(2,i+noffp,mm+moffp,k+loffp) + sdcu(2,i,mm,k)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(3,i+noffp,mm+moffp,k+loffp) + sdcu(3,i,mm,k)
      endif
  160 continue
      do 170 j = 1, mm
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
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp,k+loffp) = dcu(1,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,1,j,k)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp,k+loffp) = dcu(2,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,1,j,k)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp,k+loffp) = dcu(3,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,1,j,k)
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
!$OMP ATOMIC
         dcu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(1,nm+noffp,j+moffp,k+loffp) + sdcu(1,nm,j,k)
!$OMP ATOMIC
         dcu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(2,nm+noffp,j+moffp,k+loffp) + sdcu(2,nm,j,k)
!$OMP ATOMIC
         dcu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(3,nm+noffp,j+moffp,k+loffp) + sdcu(3,nm,j,k)
      endif
  170 continue
  180 continue
      if (lm > mz) then
         do 190 i = 2, nn
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
!$OMP ATOMIC
         dcu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,1+moffp,lm+loffp) + sdcu(1,i,1,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,1+moffp,lm+loffp) + sdcu(2,i,1,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,1+moffp,lm+loffp) + sdcu(3,i,1,lm)
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
!$OMP ATOMIC
            dcu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(1,i+noffp,mm+moffp,lm+loffp) + sdcu(1,i,mm,lm)
!$OMP ATOMIC
            dcu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(2,i+noffp,mm+moffp,lm+loffp) + sdcu(2,i,mm,lm)
!$OMP ATOMIC
            dcu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(3,i+noffp,mm+moffp,lm+loffp) + sdcu(3,i,mm,lm)
         endif
  190    continue
         do 200 j = 1, mm
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
!$OMP ATOMIC
         dcu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,1+noffp,j+moffp,lm+loffp) + sdcu(1,1,j,lm)
!$OMP ATOMIC
         dcu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,1+noffp,j+moffp,lm+loffp)  + sdcu(2,1,j,lm)
!$OMP ATOMIC
         dcu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,1+noffp,j+moffp,lm+loffp) + sdcu(3,1,j,lm)
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
!$OMP ATOMIC
            dcu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(1,nm+noffp,j+moffp,lm+loffp) + sdcu(1,nm,j,lm)
!$OMP ATOMIC
            dcu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(2,nm+noffp,j+moffp,lm+loffp) + sdcu(2,nm,j,lm)
!$OMP ATOMIC
            dcu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(3,nm+noffp,j+moffp,lm+loffp) + sdcu(3,nm,j,lm)
         endif
  200    continue
      endif
  210 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,  &
     &amu,qm,qbm,dt,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,   &
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
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
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension cu(3,nxv,nypmx,nzpmx), dcu(3,nxv,nypmx,nzpmx)
      dimension amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, scu, sdcu, samu
      dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
      dimension scu(3,MXV,MYV,MZV), sdcu(3,MXV,MYV,MZV)
      dimension samu(6,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
!     dimension scu(3,mx+1,my+1,mz+1), sdcu(3,mx+1,my+1,mz+1)
!     dimension samu(6,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm&
!$OMP& ,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,&
!$OMP& dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,   &
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxyz,sbxyz,scu,sdcu,samu)
      do 210 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 90 k = 1, mz+1
      do 80 j = 1, my+1
      do 70 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
      sdcu(1,i,j,k) = 0.0
      sdcu(2,i,j,k) = 0.0
      sdcu(3,i,j,k) = 0.0
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
! loop over particles in tile
      do 100 j = 1, nppp
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
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
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
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      sdcu(1,nn,mm,ll) = sdcu(1,nn,mm,ll) + vx*dx
      sdcu(2,nn,mm,ll) = sdcu(2,nn,mm,ll) + vy*dx
      sdcu(3,nn,mm,ll) = sdcu(3,nn,mm,ll) + vz*dx
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + ox*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + oy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + oz*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      sdcu(1,nn+1,mm,ll) = sdcu(1,nn+1,mm,ll) + vx*dy
      sdcu(2,nn+1,mm,ll) = sdcu(2,nn+1,mm,ll) + vy*dy
      sdcu(3,nn+1,mm,ll) = sdcu(3,nn+1,mm,ll) + vz*dy
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + ox*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + oy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + oz*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      sdcu(1,nn,mm+1,ll) = sdcu(1,nn,mm+1,ll) + vx*dx
      sdcu(2,nn,mm+1,ll) = sdcu(2,nn,mm+1,ll) + vy*dx
      sdcu(3,nn,mm+1,ll) = sdcu(3,nn,mm+1,ll) + vz*dx
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + ox*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + oy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + oz*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      sdcu(1,nn+1,mm+1,ll) = sdcu(1,nn+1,mm+1,ll) + vx*dy
      sdcu(2,nn+1,mm+1,ll) = sdcu(2,nn+1,mm+1,ll) + vy*dy
      sdcu(3,nn+1,mm+1,ll) = sdcu(3,nn+1,mm+1,ll) + vz*dy
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + ox*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + oy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + oz*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      sdcu(1,nn,mm,ll+1) = sdcu(1,nn,mm,ll+1) + vx*dx
      sdcu(2,nn,mm,ll+1) = sdcu(2,nn,mm,ll+1) + vy*dx
      sdcu(3,nn,mm,ll+1) = sdcu(3,nn,mm,ll+1) + vz*dx
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + ox*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + oy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + oz*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      sdcu(1,nn+1,mm,ll+1) = sdcu(1,nn+1,mm,ll+1) + vx*dy
      sdcu(2,nn+1,mm,ll+1) = sdcu(2,nn+1,mm,ll+1) + vy*dy
      sdcu(3,nn+1,mm,ll+1) = sdcu(3,nn+1,mm,ll+1) + vz*dy
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + ox*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + oy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + oz*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      sdcu(1,nn,mm+1,ll+1) = sdcu(1,nn,mm+1,ll+1) + vx*dx
      sdcu(2,nn,mm+1,ll+1) = sdcu(2,nn,mm+1,ll+1) + vy*dx
      sdcu(3,nn,mm+1,ll+1) = sdcu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + ox*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + oy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + oz*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
      sdcu(1,nn+1,mm+1,ll+1) = sdcu(1,nn+1,mm+1,ll+1) + vx*dy
      sdcu(2,nn+1,mm+1,ll+1) = sdcu(2,nn+1,mm+1,ll+1) + vy*dy
      sdcu(3,nn+1,mm+1,ll+1) = sdcu(3,nn+1,mm+1,ll+1) + vz*dy
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + ox*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + oy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + oz*dy
  100 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
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
      dcu(1,i+noffp,j+moffp,k+loffp) = dcu(1,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,i,j,k)
      dcu(2,i+noffp,j+moffp,k+loffp) = dcu(2,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,i,j,k)
      dcu(3,i+noffp,j+moffp,k+loffp) = dcu(3,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,i,j,k)
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
  110 continue
  120 continue
  130 continue
! deposit currents to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 150 j = 2, mm
      do 140 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,j+moffp,1+loffp) = dcu(1,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(1,i,j,1)
!$OMP ATOMIC
      dcu(2,i+noffp,j+moffp,1+loffp) = dcu(2,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(2,i,j,1)
!$OMP ATOMIC
      dcu(3,i+noffp,j+moffp,1+loffp) = dcu(3,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(3,i,j,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,j+moffp,lm+loffp) + sdcu(1,i,j,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,j+moffp,lm+loffp) + sdcu(2,i,j,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,j+moffp,lm+loffp) + sdcu(3,i,j,lm)
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
  140 continue
  150 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 180 k = 1, ll
      do 160 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp,k+loffp) = dcu(1,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(1,i,1,k)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp,k+loffp) = dcu(2,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(2,i,1,k)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp,k+loffp) = dcu(3,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(3,i,1,k)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(1,i+noffp,mm+moffp,k+loffp) + sdcu(1,i,mm,k)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(2,i+noffp,mm+moffp,k+loffp) + sdcu(2,i,mm,k)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(3,i+noffp,mm+moffp,k+loffp) + sdcu(3,i,mm,k)
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
  160 continue
      do 170 j = 1, mm
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
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp,k+loffp) = dcu(1,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,1,j,k)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp,k+loffp) = dcu(2,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,1,j,k)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp,k+loffp) = dcu(3,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,1,j,k)
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
!$OMP ATOMIC
         dcu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(1,nm+noffp,j+moffp,k+loffp) + sdcu(1,nm,j,k)
!$OMP ATOMIC
         dcu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(2,nm+noffp,j+moffp,k+loffp) + sdcu(2,nm,j,k)
!$OMP ATOMIC
         dcu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(3,nm+noffp,j+moffp,k+loffp) + sdcu(3,nm,j,k)
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
  170 continue
  180 continue
      if (lm > mz) then
         do 190 i = 2, nn
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
!$OMP ATOMIC
         dcu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,1+moffp,lm+loffp) + sdcu(1,i,1,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,1+moffp,lm+loffp) + sdcu(2,i,1,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,1+moffp,lm+loffp) + sdcu(3,i,1,lm)
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
!$OMP ATOMIC
            dcu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(1,i+noffp,mm+moffp,lm+loffp) + sdcu(1,i,mm,lm)
!$OMP ATOMIC
            dcu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(2,i+noffp,mm+moffp,lm+loffp) + sdcu(2,i,mm,lm)
!$OMP ATOMIC
            dcu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(3,i+noffp,mm+moffp,lm+loffp) + sdcu(3,i,mm,lm)
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
  190    continue
         do 200 j = 1, mm
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
!$OMP ATOMIC
         dcu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,1+noffp,j+moffp,lm+loffp) + sdcu(1,1,j,lm)
!$OMP ATOMIC
         dcu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,1+noffp,j+moffp,lm+loffp)  + sdcu(2,1,j,lm)
!$OMP ATOMIC
         dcu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,1+noffp,j+moffp,lm+loffp) + sdcu(3,1,j,lm)
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
!$OMP ATOMIC
            dcu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(1,nm+noffp,j+moffp,lm+loffp) + sdcu(1,nm,j,lm)
!$OMP ATOMIC
            dcu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(2,nm+noffp,j+moffp,lm+loffp) + sdcu(2,nm,j,lm)
!$OMP ATOMIC
            dcu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(3,nm+noffp,j+moffp,lm+loffp) + sdcu(3,nm,j,lm)
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
  200    continue
      endif
  210 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu, &
     &qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,    &
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux and
! acceleration density using first-order spline interpolation
! for relativistic particles.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
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
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension dcu(3,nxv,nypmx,nzpmx), amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, dzp
      real amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, p2, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, sdcu, samu
      dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
      dimension sdcu(3,MXV,MYV,MZV), samu(6,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
!     dimension sdcu(3,mx+1,my+1,mz+1), samu(6,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm&
!$OMP& ,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,&
!$OMP& dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,   &
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,gh,sfxyz,sbxyz,  &
!$OMP& sdcu,samu)
      do 210 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 90 k = 1, mz+1
      do 80 j = 1, my+1
      do 70 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
      sdcu(1,i,j,k) = 0.0
      sdcu(2,i,j,k) = 0.0
      sdcu(3,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
! loop over particles in tile
      do 100 j = 1, nppp
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
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
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
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      sdcu(1,nn,mm,ll) = sdcu(1,nn,mm,ll) + vx*dx
      sdcu(2,nn,mm,ll) = sdcu(2,nn,mm,ll) + vy*dx
      sdcu(3,nn,mm,ll) = sdcu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      sdcu(1,nn+1,mm,ll) = sdcu(1,nn+1,mm,ll) + vx*dy
      sdcu(2,nn+1,mm,ll) = sdcu(2,nn+1,mm,ll) + vy*dy
      sdcu(3,nn+1,mm,ll) = sdcu(3,nn+1,mm,ll) + vz*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      sdcu(1,nn,mm+1,ll) = sdcu(1,nn,mm+1,ll) + vx*dx
      sdcu(2,nn,mm+1,ll) = sdcu(2,nn,mm+1,ll) + vy*dx
      sdcu(3,nn,mm+1,ll) = sdcu(3,nn,mm+1,ll) + vz*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      sdcu(1,nn+1,mm+1,ll) = sdcu(1,nn+1,mm+1,ll) + vx*dy
      sdcu(2,nn+1,mm+1,ll) = sdcu(2,nn+1,mm+1,ll) + vy*dy
      sdcu(3,nn+1,mm+1,ll) = sdcu(3,nn+1,mm+1,ll) + vz*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      sdcu(1,nn,mm,ll+1) = sdcu(1,nn,mm,ll+1) + vx*dx
      sdcu(2,nn,mm,ll+1) = sdcu(2,nn,mm,ll+1) + vy*dx
      sdcu(3,nn,mm,ll+1) = sdcu(3,nn,mm,ll+1) + vz*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      sdcu(1,nn+1,mm,ll+1) = sdcu(1,nn+1,mm,ll+1) + vx*dy
      sdcu(2,nn+1,mm,ll+1) = sdcu(2,nn+1,mm,ll+1) + vy*dy
      sdcu(3,nn+1,mm,ll+1) = sdcu(3,nn+1,mm,ll+1) + vz*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      sdcu(1,nn,mm+1,ll+1) = sdcu(1,nn,mm+1,ll+1) + vx*dx
      sdcu(2,nn,mm+1,ll+1) = sdcu(2,nn,mm+1,ll+1) + vy*dx
      sdcu(3,nn,mm+1,ll+1) = sdcu(3,nn,mm+1,ll+1) + vz*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
      sdcu(1,nn+1,mm+1,ll+1) = sdcu(1,nn+1,mm+1,ll+1) + vx*dy
      sdcu(2,nn+1,mm+1,ll+1) = sdcu(2,nn+1,mm+1,ll+1) + vy*dy
      sdcu(3,nn+1,mm+1,ll+1) = sdcu(3,nn+1,mm+1,ll+1) + vz*dy
  100 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
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
      dcu(1,i+noffp,j+moffp,k+loffp) = dcu(1,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,i,j,k)
      dcu(2,i+noffp,j+moffp,k+loffp) = dcu(2,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,i,j,k)
      dcu(3,i+noffp,j+moffp,k+loffp) = dcu(3,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,i,j,k)
  110 continue
  120 continue
  130 continue
! deposit currents to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 150 j = 2, mm
      do 140 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,j+moffp,1+loffp) = dcu(1,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(1,i,j,1)
!$OMP ATOMIC
      dcu(2,i+noffp,j+moffp,1+loffp) = dcu(2,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(2,i,j,1)
!$OMP ATOMIC
      dcu(3,i+noffp,j+moffp,1+loffp) = dcu(3,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(3,i,j,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,j+moffp,lm+loffp) + sdcu(1,i,j,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,j+moffp,lm+loffp) + sdcu(2,i,j,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,j+moffp,lm+loffp) + sdcu(3,i,j,lm)
      endif
  140 continue
  150 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 180 k = 1, ll
      do 160 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp,k+loffp) = dcu(1,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(1,i,1,k)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp,k+loffp) = dcu(2,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(2,i,1,k)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp,k+loffp) = dcu(3,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(3,i,1,k)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(1,i+noffp,mm+moffp,k+loffp) + sdcu(1,i,mm,k)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(2,i+noffp,mm+moffp,k+loffp) + sdcu(2,i,mm,k)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(3,i+noffp,mm+moffp,k+loffp) + sdcu(3,i,mm,k)
      endif
  160 continue
      do 170 j = 1, mm
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
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp,k+loffp) = dcu(1,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,1,j,k)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp,k+loffp) = dcu(2,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,1,j,k)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp,k+loffp) = dcu(3,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,1,j,k)
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
!$OMP ATOMIC
         dcu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(1,nm+noffp,j+moffp,k+loffp) + sdcu(1,nm,j,k)
!$OMP ATOMIC
         dcu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(2,nm+noffp,j+moffp,k+loffp) + sdcu(2,nm,j,k)
!$OMP ATOMIC
         dcu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(3,nm+noffp,j+moffp,k+loffp) + sdcu(3,nm,j,k)
      endif
  170 continue
  180 continue
      if (lm > mz) then
         do 190 i = 2, nn
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
!$OMP ATOMIC
         dcu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,1+moffp,lm+loffp) + sdcu(1,i,1,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,1+moffp,lm+loffp) + sdcu(2,i,1,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,1+moffp,lm+loffp) + sdcu(3,i,1,lm)
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
!$OMP ATOMIC
            dcu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(1,i+noffp,mm+moffp,lm+loffp) + sdcu(1,i,mm,lm)
!$OMP ATOMIC
            dcu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(2,i+noffp,mm+moffp,lm+loffp) + sdcu(2,i,mm,lm)
!$OMP ATOMIC
            dcu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(3,i+noffp,mm+moffp,lm+loffp) + sdcu(3,i,mm,lm)
         endif
  190    continue
         do 200 j = 1, mm
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
!$OMP ATOMIC
         dcu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,1+noffp,j+moffp,lm+loffp) + sdcu(1,1,j,lm)
!$OMP ATOMIC
         dcu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,1+noffp,j+moffp,lm+loffp)  + sdcu(2,1,j,lm)
!$OMP ATOMIC
         dcu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,1+noffp,j+moffp,lm+loffp) + sdcu(3,1,j,lm)
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
!$OMP ATOMIC
            dcu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(1,nm+noffp,j+moffp,lm+loffp) + sdcu(1,nm,j,lm)
!$OMP ATOMIC
            dcu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(2,nm+noffp,j+moffp,lm+loffp) + sdcu(2,nm,j,lm)
!$OMP ATOMIC
            dcu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(3,nm+noffp,j+moffp,lm+loffp) + sdcu(3,nm,j,lm)
         endif
  200    continue
      endif
  210 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDCJPPOST32L(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu, &
     &amu,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,&
     &mxyzp1,idds)
! for 3 code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data read/written in tiles
! particles stored in segmented array
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
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension cu(3,nxv,nypmx,nzpmx), dcu(3,nxv,nypmx,nzpmx)
      dimension amu(6,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, dzp
      real amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz, p2, v1, v2, v3, v4, v5, v6
      real sfxyz, sbxyz, scu, sdcu, samu
      dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
      dimension scu(3,MXV,MYV,MZV), sdcu(3,MXV,MYV,MZV)
      dimension samu(6,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
!     dimension scu(3,mx+1,my+1,mz+1), sdcu(3,mx+1,my+1,mz+1)
!     dimension samu(6,mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm&
!$OMP& ,x,y,z,vx,vy,vz,v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,&
!$OMP& dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,   &
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,gh,sfxyz,sbxyz,  &
!$OMP& scu,sdcu,samu)
      do 210 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 90 k = 1, mz+1
      do 80 j = 1, my+1
      do 70 i = 1, mx+1
      samu(1,i,j,k) = 0.0
      samu(2,i,j,k) = 0.0
      samu(3,i,j,k) = 0.0
      samu(4,i,j,k) = 0.0
      samu(5,i,j,k) = 0.0
      samu(6,i,j,k) = 0.0
      sdcu(1,i,j,k) = 0.0
      sdcu(2,i,j,k) = 0.0
      sdcu(3,i,j,k) = 0.0
      scu(1,i,j,k) = 0.0
      scu(2,i,j,k) = 0.0
      scu(3,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
! loop over particles in tile
      do 100 j = 1, nppp
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
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
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
      samu(1,nn,mm,ll) = samu(1,nn,mm,ll) + v1*dx
      samu(2,nn,mm,ll) = samu(2,nn,mm,ll) + v2*dx
      samu(3,nn,mm,ll) = samu(3,nn,mm,ll) + v3*dx
      samu(4,nn,mm,ll) = samu(4,nn,mm,ll) + v4*dx
      samu(5,nn,mm,ll) = samu(5,nn,mm,ll) + v5*dx
      samu(6,nn,mm,ll) = samu(6,nn,mm,ll) + v6*dx
      sdcu(1,nn,mm,ll) = sdcu(1,nn,mm,ll) + vx*dx
      sdcu(2,nn,mm,ll) = sdcu(2,nn,mm,ll) + vy*dx
      sdcu(3,nn,mm,ll) = sdcu(3,nn,mm,ll) + vz*dx
      scu(1,nn,mm,ll) = scu(1,nn,mm,ll) + ox*dx
      scu(2,nn,mm,ll) = scu(2,nn,mm,ll) + oy*dx
      scu(3,nn,mm,ll) = scu(3,nn,mm,ll) + oz*dx
      dx = dyp*amz
      samu(1,nn+1,mm,ll) = samu(1,nn+1,mm,ll) + v1*dy
      samu(2,nn+1,mm,ll) = samu(2,nn+1,mm,ll) + v2*dy
      samu(3,nn+1,mm,ll) = samu(3,nn+1,mm,ll) + v3*dy
      samu(4,nn+1,mm,ll) = samu(4,nn+1,mm,ll) + v4*dy
      samu(5,nn+1,mm,ll) = samu(5,nn+1,mm,ll) + v5*dy
      samu(6,nn+1,mm,ll) = samu(6,nn+1,mm,ll) + v6*dy
      sdcu(1,nn+1,mm,ll) = sdcu(1,nn+1,mm,ll) + vx*dy
      sdcu(2,nn+1,mm,ll) = sdcu(2,nn+1,mm,ll) + vy*dy
      sdcu(3,nn+1,mm,ll) = sdcu(3,nn+1,mm,ll) + vz*dy
      scu(1,nn+1,mm,ll) = scu(1,nn+1,mm,ll) + ox*dy
      scu(2,nn+1,mm,ll) = scu(2,nn+1,mm,ll) + oy*dy
      scu(3,nn+1,mm,ll) = scu(3,nn+1,mm,ll) + oz*dy
      dy = dx1*amz
      samu(1,nn,mm+1,ll) = samu(1,nn,mm+1,ll) + v1*dx
      samu(2,nn,mm+1,ll) = samu(2,nn,mm+1,ll) + v2*dx
      samu(3,nn,mm+1,ll) = samu(3,nn,mm+1,ll) + v3*dx
      samu(4,nn,mm+1,ll) = samu(4,nn,mm+1,ll) + v4*dx
      samu(5,nn,mm+1,ll) = samu(5,nn,mm+1,ll) + v5*dx
      samu(6,nn,mm+1,ll) = samu(6,nn,mm+1,ll) + v6*dx
      sdcu(1,nn,mm+1,ll) = sdcu(1,nn,mm+1,ll) + vx*dx
      sdcu(2,nn,mm+1,ll) = sdcu(2,nn,mm+1,ll) + vy*dx
      sdcu(3,nn,mm+1,ll) = sdcu(3,nn,mm+1,ll) + vz*dx
      scu(1,nn,mm+1,ll) = scu(1,nn,mm+1,ll) + ox*dx
      scu(2,nn,mm+1,ll) = scu(2,nn,mm+1,ll) + oy*dx
      scu(3,nn,mm+1,ll) = scu(3,nn,mm+1,ll) + oz*dx
      dx = amx*dzp
      samu(1,nn+1,mm+1,ll) = samu(1,nn+1,mm+1,ll) + v1*dy
      samu(2,nn+1,mm+1,ll) = samu(2,nn+1,mm+1,ll) + v2*dy
      samu(3,nn+1,mm+1,ll) = samu(3,nn+1,mm+1,ll) + v3*dy
      samu(4,nn+1,mm+1,ll) = samu(4,nn+1,mm+1,ll) + v4*dy
      samu(5,nn+1,mm+1,ll) = samu(5,nn+1,mm+1,ll) + v5*dy
      samu(6,nn+1,mm+1,ll) = samu(6,nn+1,mm+1,ll) + v6*dy
      sdcu(1,nn+1,mm+1,ll) = sdcu(1,nn+1,mm+1,ll) + vx*dy
      sdcu(2,nn+1,mm+1,ll) = sdcu(2,nn+1,mm+1,ll) + vy*dy
      sdcu(3,nn+1,mm+1,ll) = sdcu(3,nn+1,mm+1,ll) + vz*dy
      scu(1,nn+1,mm+1,ll) = scu(1,nn+1,mm+1,ll) + ox*dy
      scu(2,nn+1,mm+1,ll) = scu(2,nn+1,mm+1,ll) + oy*dy
      scu(3,nn+1,mm+1,ll) = scu(3,nn+1,mm+1,ll) + oz*dy
      dy = amy*dzp
      samu(1,nn,mm,ll+1) = samu(1,nn,mm,ll+1) + v1*dx
      samu(2,nn,mm,ll+1) = samu(2,nn,mm,ll+1) + v2*dx
      samu(3,nn,mm,ll+1) = samu(3,nn,mm,ll+1) + v3*dx
      samu(4,nn,mm,ll+1) = samu(4,nn,mm,ll+1) + v4*dx
      samu(5,nn,mm,ll+1) = samu(5,nn,mm,ll+1) + v5*dx
      samu(6,nn,mm,ll+1) = samu(6,nn,mm,ll+1) + v6*dx
      sdcu(1,nn,mm,ll+1) = sdcu(1,nn,mm,ll+1) + vx*dx
      sdcu(2,nn,mm,ll+1) = sdcu(2,nn,mm,ll+1) + vy*dx
      sdcu(3,nn,mm,ll+1) = sdcu(3,nn,mm,ll+1) + vz*dx
      scu(1,nn,mm,ll+1) = scu(1,nn,mm,ll+1) + ox*dx
      scu(2,nn,mm,ll+1) = scu(2,nn,mm,ll+1) + oy*dx
      scu(3,nn,mm,ll+1) = scu(3,nn,mm,ll+1) + oz*dx
      dx = dyp*dzp
      samu(1,nn+1,mm,ll+1) = samu(1,nn+1,mm,ll+1) + v1*dy
      samu(2,nn+1,mm,ll+1) = samu(2,nn+1,mm,ll+1) + v2*dy
      samu(3,nn+1,mm,ll+1) = samu(3,nn+1,mm,ll+1) + v3*dy
      samu(4,nn+1,mm,ll+1) = samu(4,nn+1,mm,ll+1) + v4*dy
      samu(5,nn+1,mm,ll+1) = samu(5,nn+1,mm,ll+1) + v5*dy
      samu(6,nn+1,mm,ll+1) = samu(6,nn+1,mm,ll+1) + v6*dy
      sdcu(1,nn+1,mm,ll+1) = sdcu(1,nn+1,mm,ll+1) + vx*dy
      sdcu(2,nn+1,mm,ll+1) = sdcu(2,nn+1,mm,ll+1) + vy*dy
      sdcu(3,nn+1,mm,ll+1) = sdcu(3,nn+1,mm,ll+1) + vz*dy
      scu(1,nn+1,mm,ll+1) = scu(1,nn+1,mm,ll+1) + ox*dy
      scu(2,nn+1,mm,ll+1) = scu(2,nn+1,mm,ll+1) + oy*dy
      scu(3,nn+1,mm,ll+1) = scu(3,nn+1,mm,ll+1) + oz*dy
      dy = dx1*dzp
      samu(1,nn,mm+1,ll+1) = samu(1,nn,mm+1,ll+1) + v1*dx
      samu(2,nn,mm+1,ll+1) = samu(2,nn,mm+1,ll+1) + v2*dx
      samu(3,nn,mm+1,ll+1) = samu(3,nn,mm+1,ll+1) + v3*dx
      samu(4,nn,mm+1,ll+1) = samu(4,nn,mm+1,ll+1) + v4*dx
      samu(5,nn,mm+1,ll+1) = samu(5,nn,mm+1,ll+1) + v5*dx
      samu(6,nn,mm+1,ll+1) = samu(6,nn,mm+1,ll+1) + v6*dx
      sdcu(1,nn,mm+1,ll+1) = sdcu(1,nn,mm+1,ll+1) + vx*dx
      sdcu(2,nn,mm+1,ll+1) = sdcu(2,nn,mm+1,ll+1) + vy*dx
      sdcu(3,nn,mm+1,ll+1) = sdcu(3,nn,mm+1,ll+1) + vz*dx
      scu(1,nn,mm+1,ll+1) = scu(1,nn,mm+1,ll+1) + ox*dx
      scu(2,nn,mm+1,ll+1) = scu(2,nn,mm+1,ll+1) + oy*dx
      scu(3,nn,mm+1,ll+1) = scu(3,nn,mm+1,ll+1) + oz*dx
      samu(1,nn+1,mm+1,ll+1) = samu(1,nn+1,mm+1,ll+1) + v1*dy
      samu(2,nn+1,mm+1,ll+1) = samu(2,nn+1,mm+1,ll+1) + v2*dy
      samu(3,nn+1,mm+1,ll+1) = samu(3,nn+1,mm+1,ll+1) + v3*dy
      samu(4,nn+1,mm+1,ll+1) = samu(4,nn+1,mm+1,ll+1) + v4*dy
      samu(5,nn+1,mm+1,ll+1) = samu(5,nn+1,mm+1,ll+1) + v5*dy
      samu(6,nn+1,mm+1,ll+1) = samu(6,nn+1,mm+1,ll+1) + v6*dy
      sdcu(1,nn+1,mm+1,ll+1) = sdcu(1,nn+1,mm+1,ll+1) + vx*dy
      sdcu(2,nn+1,mm+1,ll+1) = sdcu(2,nn+1,mm+1,ll+1) + vy*dy
      sdcu(3,nn+1,mm+1,ll+1) = sdcu(3,nn+1,mm+1,ll+1) + vz*dy
      scu(1,nn+1,mm+1,ll+1) = scu(1,nn+1,mm+1,ll+1) + ox*dy
      scu(2,nn+1,mm+1,ll+1) = scu(2,nn+1,mm+1,ll+1) + oy*dy
      scu(3,nn+1,mm+1,ll+1) = scu(3,nn+1,mm+1,ll+1) + oz*dy
  100 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
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
      dcu(1,i+noffp,j+moffp,k+loffp) = dcu(1,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,i,j,k)
      dcu(2,i+noffp,j+moffp,k+loffp) = dcu(2,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,i,j,k)
      dcu(3,i+noffp,j+moffp,k+loffp) = dcu(3,i+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,i,j,k)
      cu(1,i+noffp,j+moffp,k+loffp) = cu(1,i+noffp,j+moffp,k+loffp)     &
     &+ scu(1,i,j,k)
      cu(2,i+noffp,j+moffp,k+loffp) = cu(2,i+noffp,j+moffp,k+loffp)     &
     &+ scu(2,i,j,k)
      cu(3,i+noffp,j+moffp,k+loffp) = cu(3,i+noffp,j+moffp,k+loffp)     &
     &+ scu(3,i,j,k)
  110 continue
  120 continue
  130 continue
! deposit currents to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 150 j = 2, mm
      do 140 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,j+moffp,1+loffp) = dcu(1,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(1,i,j,1)
!$OMP ATOMIC
      dcu(2,i+noffp,j+moffp,1+loffp) = dcu(2,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(2,i,j,1)
!$OMP ATOMIC
      dcu(3,i+noffp,j+moffp,1+loffp) = dcu(3,i+noffp,j+moffp,1+loffp)   &
     &+ sdcu(3,i,j,1)
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
!$OMP ATOMIC
         dcu(1,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,j+moffp,lm+loffp) + sdcu(1,i,j,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,j+moffp,lm+loffp) + sdcu(2,i,j,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,j+moffp,lm+loffp) + sdcu(3,i,j,lm)
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
  140 continue
  150 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 180 k = 1, ll
      do 160 i = 2, nn
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
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp,k+loffp) = dcu(1,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(1,i,1,k)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp,k+loffp) = dcu(2,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(2,i,1,k)
!$OMP ATOMIC
      dcu(3,i+noffp,1+moffp,k+loffp) = dcu(3,i+noffp,1+moffp,k+loffp)   &
     &+ sdcu(3,i,1,k)
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
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(1,i+noffp,mm+moffp,k+loffp) + sdcu(1,i,mm,k)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(2,i+noffp,mm+moffp,k+loffp) + sdcu(2,i,mm,k)
!$OMP ATOMIC
         dcu(3,i+noffp,mm+moffp,k+loffp) =                              &
     &   dcu(3,i+noffp,mm+moffp,k+loffp) + sdcu(3,i,mm,k)
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
  160 continue
      do 170 j = 1, mm
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
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp,k+loffp) = dcu(1,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(1,1,j,k)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp,k+loffp) = dcu(2,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(2,1,j,k)
!$OMP ATOMIC
      dcu(3,1+noffp,j+moffp,k+loffp) = dcu(3,1+noffp,j+moffp,k+loffp)   &
     &+ sdcu(3,1,j,k)
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
!$OMP ATOMIC
         dcu(1,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(1,nm+noffp,j+moffp,k+loffp) + sdcu(1,nm,j,k)
!$OMP ATOMIC
         dcu(2,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(2,nm+noffp,j+moffp,k+loffp) + sdcu(2,nm,j,k)
!$OMP ATOMIC
         dcu(3,nm+noffp,j+moffp,k+loffp) =                              &
     &   dcu(3,nm+noffp,j+moffp,k+loffp) + sdcu(3,nm,j,k)
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
  170 continue
  180 continue
      if (lm > mz) then
         do 190 i = 2, nn
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
!$OMP ATOMIC
         dcu(1,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(1,i+noffp,1+moffp,lm+loffp) + sdcu(1,i,1,lm)
!$OMP ATOMIC
         dcu(2,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(2,i+noffp,1+moffp,lm+loffp) + sdcu(2,i,1,lm)
!$OMP ATOMIC
         dcu(3,i+noffp,1+moffp,lm+loffp) =                              &
     &   dcu(3,i+noffp,1+moffp,lm+loffp) + sdcu(3,i,1,lm)
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
!$OMP ATOMIC
            dcu(1,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(1,i+noffp,mm+moffp,lm+loffp) + sdcu(1,i,mm,lm)
!$OMP ATOMIC
            dcu(2,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(2,i+noffp,mm+moffp,lm+loffp) + sdcu(2,i,mm,lm)
!$OMP ATOMIC
            dcu(3,i+noffp,mm+moffp,lm+loffp) =                          &
     &      dcu(3,i+noffp,mm+moffp,lm+loffp) + sdcu(3,i,mm,lm)
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
  190    continue
         do 200 j = 1, mm
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
!$OMP ATOMIC
         dcu(1,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(1,1+noffp,j+moffp,lm+loffp) + sdcu(1,1,j,lm)
!$OMP ATOMIC
         dcu(2,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(2,1+noffp,j+moffp,lm+loffp)  + sdcu(2,1,j,lm)
!$OMP ATOMIC
         dcu(3,1+noffp,j+moffp,lm+loffp) =                              &
     &   dcu(3,1+noffp,j+moffp,lm+loffp) + sdcu(3,1,j,lm)
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
!$OMP ATOMIC
            dcu(1,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(1,nm+noffp,j+moffp,lm+loffp) + sdcu(1,nm,j,lm)
!$OMP ATOMIC
            dcu(2,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(2,nm+noffp,j+moffp,lm+loffp) + sdcu(2,nm,j,lm)
!$OMP ATOMIC
            dcu(3,nm+noffp,j+moffp,lm+loffp) =                          &
     &      dcu(3,nm+noffp,j+moffp,lm+loffp) + sdcu(3,nm,j,lm)
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
  200    continue
      endif
  210 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPASCFGUARD32L(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,  &
     &idds)
! add scaled field to extended periodic field
! linear interpolation, for distributed data with 2D decomposition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! q2m0 = wp0/affp, where
! wp0 = normalized total plasma frequency squared
! affp = normalization constant = nx*ny/np, where np=number of particles
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      integer nyzp
      real dcu, cus, q2m0
      dimension dcu(3,nxe,nypmx,nzpmx), cus(3,nxe,nypmx,nzpmx)
      dimension nyzp(idds)
! local data
      integer i, j, k, l, ll
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
      do 30 ll = 1, nyzp(1)*nyzp(2)
      l = (ll - 1)/nyzp(1)
      k = ll - nyzp(1)*l
      l = l + 1
      do 20 j = 1, nx
      do 10 i = 1, 3
      dcu(i,j,k,l) = dcu(i,j,k,l) - q2m0*cus(i,j,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end