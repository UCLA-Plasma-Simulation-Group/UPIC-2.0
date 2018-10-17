!-----------------------------------------------------------------------
! Fortran Library for pushing electromagnetic particles
! 2-1/2D Vector/MPI/OpenMP PIC Codes:
! VPPGBPPUSH23L updates magnetized particle co-ordinates and velocities
!               using leap-frog scheme in time and linear interpolation
! VPPGBPPUSHF23L updates magnetized particle co-ordinates and velocities
!                using leap-frog scheme in time and linear interpolation
!                in space with periodic particle boundary conditions,
!                determines list of particles which are leaving each
!                tile
! VPPGRBPPUSH23L updates relativistic magnetized particle co-ordinates
!                and momenta using leap-frog scheme in time and linear
!                interpolation in space with various particle boundary
!                conditions
! VPPGRBPPUSHF23L updates relativistic magnetized particle co-ordinates
!                 and momenta using leap-frog scheme in time and linear
!                 interpolation in space with periodic particle boundary
!                 conditions, determines list of particles which are
!                 leaving each tile
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 12, 2018
!-----------------------------------------------------------------------
      subroutine VPPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek&
     &,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 119 flops/particle, 1 divide, 29 loads, 5 stores
! input: all, output: ppart, ek
! velocity equations used are:
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
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
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
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,2*lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,   &
!$OMP& omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,&
!$OMP& sfxy,sbxy,n,s,t) REDUCTION(+:sum2)
      do 110 k = 1, mxyp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 50 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   50 continue
! find acceleration
      do 60 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
      s(j,4) = ox
      s(j,5) = oy
      s(j,6) = oz
   60 continue
! new velocity
! !dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
      acz = ppart(5,j+joff,k) + dz
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
      s(j,3) = vx
      s(j,4) = vy
      s(j,5) = vz
   70 continue
! check boundary conditions
! !dir$ vector aligned
      do 80 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      vx = s(j,3)
      vy = s(j,4)
      vz = s(j,5)
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
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
      endif
! set new position
      ppart(1,j+joff,k) = dx
      ppart(2,j+joff,k) = dy
! set new velocity
      ppart(3,j+joff,k) = vx
      ppart(4,j+joff,k) = vy
      ppart(5,j+joff,k) = vz
   80 continue
   90 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 100 j = nps, nppp
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
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
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
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new velocity
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
      ppart(5,j,k) = vz
  100 continue
      sum2 = sum2 + sum1
  110 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,  &
     &qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,  &
     &irc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 119 flops/particle, 1 divide, 29 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! velocity equations used are:
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
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
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
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, ih, nh, nn, mm, lxv
      real qtmh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,2*lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,ipp,joff,nps, &
!$OMP& x,y,vx,vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt, &
!$OMP& omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,&
!$OMP& edgelx,edgely,edgerx,edgery,sum1,sfxy,sbxy,n,s,t)                &
!$OMP& REDUCTION(+:sum2)
      do 130 k = 1, mxyp1
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
      mnoff = moffp + noff
! load local fields from global arrays
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
!dir$ ivdep
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 110 m = 1, ipp
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
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
      s(j,4) = ox
      s(j,5) = oy
      s(j,6) = oz
   70 continue
! new velocity
! !dir$ vector aligned
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
      acz = ppart(5,j+joff,k) + dz
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
      s(j,3) = vx
      s(j,4) = vy
      s(j,5) = vz
   80 continue
! check boundary conditions
! !dir$ vector aligned
      do 90 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
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
! set new position
      ppart(1,j+joff,k) = s(j,1)
      ppart(2,j+joff,k) = s(j,2)
! set new velocity
      ppart(3,j+joff,k) = s(j,3)
      ppart(4,j+joff,k) = s(j,4)
      ppart(5,j+joff,k) = s(j,5)
      n(j) = mm
   90 continue
! increment counters
      do 100 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 120 j = nps, nppp
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
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
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
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new velocity
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
      ppart(5,j,k) = vz
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
               ppart(1,j,k) = 0.0
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
               ppart(2,j,k) = 0.0
               mm = mm - 3
            endif
         endif
      endif
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
  120 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 140 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  140    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc, &
     &ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all, output: ppart, ek
! momentum equations used are:
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
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum vx of particle n in partition in tile m
! ppart(4,n,m) = momentum vy of particle n in partition in tile m
! ppart(5,n,m) = momentum vz of particle n in partition in tile m
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
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,2*lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
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
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,mnoff,ipp,joff,nps,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,   &
!$OMP& omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,  &
!$OMP& gami,qtmg,dtg,sum1,sfxy,sbxy,n,s,t) REDUCTION(+:sum2)
      do 110 k = 1, mxyp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 50 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   50 continue
! find acceleration
      do 60 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
      s(j,4) = ox
      s(j,5) = oy
      s(j,6) = oz
   60 continue
! new momentum
! !dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
      acz = ppart(5,j+joff,k) + dz
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
      s(j,3) = vx
      s(j,4) = vy
      s(j,5) = vz
   70 continue
! check boundary conditions
! !dir$ vector aligned
      do 80 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      vx = s(j,3)
      vy = s(j,4)
      vz = s(j,5)
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
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
      endif
! set new position
      ppart(1,j+joff,k) = dx
      ppart(2,j+joff,k) = dy
! set new momentum
      ppart(3,j+joff,k) = vx
      ppart(4,j+joff,k) = vy
      ppart(5,j+joff,k) = vz
   80 continue
   90 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 100 j = nps, nppp
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
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
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
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new momentum
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
      ppart(5,j,k) = vz
  100 continue
      sum2 = sum2 + sum1
  110 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp, &
     &qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax&
     &,irc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! momentum equations used are:
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
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum vx of particle n in partition in tile m
! ppart(4,n,m) = momentum vy of particle n in partition in tile m
! ppart(5,n,m) = momentum vz of particle n in partition in tile m
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
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv*nypmx), bxy(3,nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, ih, nh, nn, mm, lxv
      real qtmh, ci2, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
!     dimension sfxy(3,MXV*MYV), sbxy(3,MXV*MYV)
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,2*lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,ipp,joff,nps, &
!$OMP& x,y,vx,vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt, &
!$OMP& omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,&
!$OMP& edgelx,edgely,edgerx,edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy,n,s,t&
!$OMP& ) REDUCTION(+:sum2)
      do 130 k = 1, mxyp1
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
      mnoff = moffp + noff
! load local fields from global arrays
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
!dir$ ivdep
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 110 m = 1, ipp
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
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   60 continue
! find acceleration
      do 70 j = 1, npblk
      nn = n(j) + 1
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
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
      s(j,4) = ox
      s(j,5) = oy
      s(j,6) = oz
   70 continue
! new momentum
! !dir$ vector aligned
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
! calculate half impulse
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
      acz = ppart(5,j+joff,k) + dz
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
! new position and velocity
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = vx
      s(j,4) = vy
      s(j,5) = vz
   80 continue
! check boundary conditions
! !dir$ vector aligned
      do 90 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
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
! set new position
      ppart(1,j+joff,k) = s(j,1)
      ppart(2,j+joff,k) = s(j,2)
! set new momentum
      ppart(3,j+joff,k) = s(j,3)
      ppart(4,j+joff,k) = s(j,4)
      ppart(5,j+joff,k) = s(j,5)
      n(j) = mm
   90 continue
! increment counters
      do 100 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 120 j = nps, nppp
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
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
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
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new momenta
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
      ppart(5,j,k) = vz
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
               ppart(1,j,k) = 0.0
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
               ppart(2,j,k) = 0.0
               mm = mm - 3
            endif
         endif
      endif
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
  120 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 140 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  140    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
