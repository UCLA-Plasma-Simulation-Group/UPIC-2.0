!-----------------------------------------------------------------------
! Fortran Library for pushing electromagnetic particles
! 3D Vector/MPI/OpenMP PIC Codes:
! VPPGBPPUSH32L updates magnetized particle co-ordinates and velocities
!               using leap-frog scheme in time and linear interpolation
!               in space with various particle boundary conditions
! VPPGBPPUSHF32L updates magnetized particle co-ordinates and velocities
!                using leap-frog scheme in time and linear interpolation
!                in space with periodic particle boundary conditions,
!                determines list of particles which are leaving each
!                tile
! VPPGRBPPUSH32L updates relativistic magnetized particle co-ordinates
!                and momenta using leap-frog scheme in time and linear
!                interpolation in space with various particle boundary
!                conditions
! VPPGRBPPUSHF32L updates relativistic magnetized particle co-ordinates
!                 and momenta using leap-frog scheme in time and linear
!                 interpolation in space with periodic particle boundary
!                 conditions, determines list of particles which are
!                 leaving each tile
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: april 27, 2018
!-----------------------------------------------------------------------
      subroutine VPPGBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,dtc&
     &,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field.  Using the Boris Mover.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 190 flops/particle, 1 divide, 54 loads, 6 stores
! input: all, output: ppart, ek
! velocity equations used are:
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
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! z(t+dt)=z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
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
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qbm, dt, dtc, ek
      real ppart, fxyz, bxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=8,loopv=1)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, lxv, lxyv, nxyv
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz
      real sfxyz, sbxyz
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,   &
!$OMP& nps,nn,mm,ll,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,&
!$OMP& ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4&
!$OMP& ,rot5,rot6,rot7,rot8,rot9,sum1,sfxyz,sbxyz,n,s,t)                &
!$OMP& REDUCTION(+:sum2)
      do 130 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 70 j = 1, npblk
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
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   70 continue
! find acceleration
      do 80 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
      s(j,4) = ox + acx
      s(j,5) = oy + acy
      s(j,6) = oz + acz
   80 continue
! new velocity
! !dir$ vector aligned
      do 90 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      omxt = qtmh*s(j,4)
      omyt = qtmh*s(j,5)
      omzt = qtmh*s(j,6)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! new position and velocity
      s(j,1) = x + vx*dtc
      s(j,2) = y + vy*dtc
      s(j,3) = z + vz*dtc
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   90 continue
! check boundary conditions
! !dir$ vector aligned
      do 100 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
      vx = s(j,4)
      vy = s(j,5)
      vz = s(j,6)
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
      endif
! set new position
      ppart(1,j+joff,l) = dx
      ppart(2,j+joff,l) = dy
      ppart(3,j+joff,l) = dz
! set new velocity
      ppart(4,j+joff,l) = vx
      ppart(5,j+joff,l) = vy
      ppart(6,j+joff,l) = vz
  100 continue
  110 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 120 j = nps, nppp
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
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! new position
      dx = x + vx*dtc
      dy = y + vy*dtc
      dz = z + vz*dtc
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! set new velocity
      ppart(4,j,l) = vx
      ppart(5,j,l) = vy
      ppart(6,j,l) = vz
  120 continue
      sum2 = sum2 + sum1
  130 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp&
     &,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1, &
     &myp1,mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field.  Using the Boris Mover.
! with periodic boundary conditions
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 190 flops/particle, 1 divide, 54 loads, 6 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! velocity equations used are:
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
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! z(t+dt)=z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
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
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ntmax = size of hole array for particles leaving tiles
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qbm, dt, dtc, ek
      real ppart, fxyz, bxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, ih, nh, nn, mm, ll
      integer lxv, lxyv, nxyv
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real qtmh, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz
      real sfxyz, sbxyz
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,nps&
!$OMP& ,ih,nh,nn,mm,ll,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,&
!$OMP& dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3, &
!$OMP& rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgelz,edgerx,edgery,&
!$OMP& edgerz,sum1,sfxyz,sbxyz,n,s,t) REDUCTION(+:sum2)
      do 150 l = 1, mxyzp1
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
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, ll+1
      do 50 j = 1, mm+1
!dir$ ivdep
      do 40 i = 1, nn+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! clear counters
      do 70 j = 1, 26
      ncl(j,l) = 0
   70 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
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
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
      s(j,4) = ox + acx
      s(j,5) = oy + acy
      s(j,6) = oz + acz
   90 continue
! new velocity
! !dir$ vector aligned
      do 100 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      omxt = qtmh*s(j,4)
      omyt = qtmh*s(j,5)
      omzt = qtmh*s(j,6)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! new position and velocity
      s(j,1) = x + vx*dtc
      s(j,2) = y + vy*dtc
      s(j,3) = z + vz*dtc
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
  100 continue
! check boundary conditions
! !dir$ vector aligned
      do 110 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         mm = 2
      else if (dx.lt.edgelx) then
         mm = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               s(j,1) = 0.0
               mm = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         mm = mm + 6
      else if (dy.lt.edgely) then
         mm = mm + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               s(j,2) = 0.0
               mm = mm - 3
            endif
         endif
      endif
      if (dz.ge.edgerz) then
         mm = mm + 18
      else if (dz.lt.edgelz) then
         mm = mm + 9
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.ge.anz) then
               s(j,3) = 0.0
               mm = mm - 9
            endif
         endif
      endif
! set new position
      ppart(1,j+joff,l) = s(j,1)
      ppart(2,j+joff,l) = s(j,2)
      ppart(3,j+joff,l) = s(j,3)
! set new velocity
      ppart(4,j+joff,l) = s(j,4)
      ppart(5,j+joff,l) = s(j,5)
      ppart(6,j+joff,l) = s(j,6)
      n(j) = mm
  110 continue
! increment counters
      do 120 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j + joff
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
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
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! new position
      dx = x + vx*dtc
      dy = y + vy*dtc
      dz = z + vz*dtc
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! set new velocity
      ppart(4,j,l) = vx
      ppart(5,j,l) = vy
      ppart(6,j,l) = vz
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         mm = 2
      else if (dx.lt.edgelx) then
         mm = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,l) = 0.0
               mm = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         mm = mm + 6
      else if (dy.lt.edgely) then
         mm = mm + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               ppart(2,j,l) = 0.0
               mm = mm - 3
            endif
         endif
      endif
      if (dz.ge.edgerz) then
         mm = mm + 18
      else if (dz.lt.edgelz) then
         mm = mm + 9
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.ge.anz) then
               ppart(3,j,l) = 0.0
               mm = mm - 9
            endif
         endif
      endif
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
  140 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  150 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 160 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  160    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRBPPUSH32L(ppart,fxyz,bxyz,kpic,noff,nyzp,qbm,dt,  &
     &dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1, &
     &mxyzp1,idds,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
! input: all, output: ppart, ek
! momentum equations used are:
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
! omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
! omy = (q/m)*by(x(t),y(t),z(t))*gami,
! omz = (q/m)*bz(x(t),y(t),z(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
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
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qbm, dt, dtc, ci, ek
      real ppart, fxyz, bxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=8,loopv=1)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, lxv, lxyv, nxyv
      real qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg
      real x, y, z, vx, vy, vz
      real sfxyz, sbxyz
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,nps&
!$OMP& ,nn,mm,ll,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,&
!$OMP& oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,  &
!$OMP& rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxyz,sbxyz,n,s,t)&
!$OMP& REDUCTION(+:sum2)
      do 130 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 70 j = 1, npblk
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
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   70 continue
! find acceleration
      do 80 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
      s(j,4) = ox + acx
      s(j,5) = oy + acy
      s(j,6) = oz + acz
   80 continue
! new momentum
! !dir$ vector aligned
      do 90 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmg*s(j,4)
      omyt = qtmg*s(j,5)
      omzt = qtmg*s(j,6)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = z + vz*dtg
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   90 continue
! check boundary conditions
! !dir$ vector aligned
      do 100 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
      vx = s(j,4)
      vy = s(j,5)
      vz = s(j,6)
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
      endif
! set new position
      ppart(1,j+joff,l) = dx
      ppart(2,j+joff,l) = dy
      ppart(3,j+joff,l) = dz
! set new momentum
      ppart(4,j+joff,l) = vx
      ppart(5,j+joff,l) = vy
      ppart(6,j+joff,l) = vz
  100 continue
  110 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 120 j = nps, nppp
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
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
! new momentum
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + vx*dtg
      dy = y + vy*dtg
      dz = z + vz*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! set new momentum
      ppart(4,j,l) = vx
      ppart(5,j,l) = vy
      ppart(6,j,l) = vz
  120 continue
      sum2 = sum2 + sum1
  130 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRBPPUSHF32L(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,   &
     &nyzp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,    &
     &nzpmx,mx1,myp1,mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! momentum equations used are:
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
! omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
! omy = (q/m)*by(x(t),y(t),z(t))*gami,
! omz = (q/m)*bz(x(t),y(t),z(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
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
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ntmax = size of hole array for particles leaving tiles
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qbm, dt, dtc, ci, ek
      real ppart, fxyz, bxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv*nypmx*nzpmx), bxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, ih, nh, nn, mm, ll
      integer lxv, lxyv, nxyv
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg
      real qtmh, ci2, x, y, z, vx, vy, vz
      real sfxyz, sbxyz
!     dimension sfxyz(3,MXV*MYV*MZV), sbxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
      dimension sbxyz(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nypmx
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      ci2 = ci*ci
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,nps&
!$OMP& ,ih,nh,nn,mm,ll,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,&
!$OMP& dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3, &
!$OMP& rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,edgelx,edgely,    &
!$OMP& edgelz,edgerx,edgery,edgerz,sum1,sfxyz,sbxyz,n,s,t)              &
!$OMP& REDUCTION(+:sum2)
      do 150 l = 1, mxyzp1
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
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & fxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   10 continue
   20 continue
   30 continue
      do 60 k = 1, ll+1
      do 50 j = 1, mm+1
!dir$ ivdep
      do 40 i = 1, nn+1
      sbxyz(1,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(2,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
      sbxyz(3,i+lxv*(j-1)+lxyv*(k-1)) =                                 &
     & bxyz(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))
   40 continue
   50 continue
   60 continue
! clear counters
      do 70 j = 1, 26
      ncl(j,l) = 0
   70 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 130 m = 1, ipp
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
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   80 continue
! find acceleration
      do 90 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
      s(j,4) = ox + acx
      s(j,5) = oy + acy
      s(j,6) = oz + acz
   90 continue
! new momentum
! !dir$ vector aligned
      do 100 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmg*s(j,4)
      omyt = qtmg*s(j,5)
      omzt = qtmg*s(j,6)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = z + vz*dtg
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
  100 continue
! check boundary conditions
! !dir$ vector aligned
      do 110 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         mm = 2
      else if (dx.lt.edgelx) then
         mm = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               s(j,1) = 0.0
               mm = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         mm = mm + 6
      else if (dy.lt.edgely) then
         mm = mm + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               s(j,2) = 0.0
               mm = mm - 3
            endif
         endif
      endif
      if (dz.ge.edgerz) then
         mm = mm + 18
      else if (dz.lt.edgelz) then
         mm = mm + 9
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.ge.anz) then
               s(j,3) = 0.0
               mm = mm - 9
            endif
         endif
      endif
! set new position
      ppart(1,j+joff,l) = s(j,1)
      ppart(2,j+joff,l) = s(j,2)
      ppart(3,j+joff,l) = s(j,3)
! set new velocity
      ppart(4,j+joff,l) = s(j,4)
      ppart(5,j+joff,l) = s(j,5)
      ppart(6,j+joff,l) = s(j,6)
      n(j) = mm
  110 continue
! increment counters
      do 120 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j + joff
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  120 continue
  130 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 140 j = nps, nppp
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
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
! new momentum
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + vx*dtg
      dy = y + vy*dtg
      dz = z + vz*dtg
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
! set new momentum
      ppart(4,j,l) = vx
      ppart(5,j,l) = vy
      ppart(6,j,l) = vz
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         mm = 2
      else if (dx.lt.edgelx) then
         mm = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,l) = 0.0
               mm = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         mm = mm + 6
      else if (dy.lt.edgely) then
         mm = mm + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               ppart(2,j,l) = 0.0
               mm = mm - 3
            endif
         endif
      endif
      if (dz.ge.edgerz) then
         mm = mm + 18
      else if (dz.lt.edgelz) then
         mm = mm + 9
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.ge.anz) then
               ppart(3,j,l) = 0.0
               mm = mm - 9
            endif
         endif
      endif
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
  140 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  150 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 160 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  160    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
