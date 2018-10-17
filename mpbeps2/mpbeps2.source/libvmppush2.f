!-----------------------------------------------------------------------
! Fortran Library for pushing electrostatic particles and depositing
! charge
! 2D Vector/MPI/OpenMP PIC Codes:
! VPPGPPUSH2L updates particle co-ordinates and velocities using
!             electric field only, with linear interpolation and various
!             particle boundary conditions
! VPPGPPUSHF2L updates particle co-ordinates and velocities using
!              electric field only, with linear interpolation and
!              periodic particle boundary conditions.  also determines
!              list of particles which are leaving each tile
! VPPGRPPUSH2L updates relativistic particle co-ordinates and momenta
!              using electric field only, with linear interpolation and
!              various particle boundary conditions
! VPPGRPPUSHF2L updates relativistic particle co-ordinates and momenta
!               using electric field only, with linear interpolation and
!               periodic particle boundary conditions.  also determines
!               list of particles which are leaving each tile
! VPPGPPUSH2ZF update particle co-ordinate for particles with fixed
!              velocities
! VPPGPPUSHF2ZF update particle co-ordinate for particles with fixed
!               velocities with periodic particle boundary conditions.
!               also determines list of particles which are leaving each
!               tile
! VPPGRPPUSH2ZF update particle co-ordinates for particles with fixed
!               velocities, for 2d code, and relativistic particles.
! VPPGRPPUSHF2ZF update particle co-ordinates for particles with fixed
!                velocities with periodic particle boundary conditions,
!                for 2d code, and relativistic particles.  also
!                determines list of particles which are leaving each
!                tile
! VPPGPPOST2L calculates particle charge density using linear
!             interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 18, 2018
!-----------------------------------------------------------------------
      subroutine VPPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,mx,&
     &my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with various boundary conditions
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 42 flops/particle, 12 loads, 4 stores
! input: all, output: ppart, ek
! equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, ek
      real ppart, fxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=4,loopv=1)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real sfxy
!     dimension sfxy(2,MXV*MYV)
      dimension sfxy(2,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtm = qbm*dt
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,x,y,   &
!$OMP& dxp,dyp,amx,amy,dx,dy,vx,vy,sum1,sfxy,n,s,t) REDUCTION(+:sum2)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 70 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 30 j = 1, npblk
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
   30 continue
! find acceleration
      do 40 j = 1, npblk
      nn = n(j) + 1
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      s(j,1) = dx
      s(j,2) = dy
   40 continue
! new velocity
! !dir$ vector aligned
      do 50 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      dxp = ppart(3,j+joff,k)
      dyp = ppart(4,j+joff,k)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      sum1 = sum1 + (dxp*dxp + dyp*dyp)
! new position and velocity
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = vx
      s(j,4) = vy
   50 continue
! check boundary conditions
! !dir$ vector aligned
      do 60 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      vx = s(j,3)
      vy = s(j,4)
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
   60 continue
   70 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 80 j = nps, nppp
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
! find acceleration
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
! new velocity
      dxp = ppart(3,j,k)
      dyp = ppart(4,j,k)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      sum1 = sum1 + (dxp*dxp + dyp*dyp)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
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
   80 continue
      sum2 = sum2 + sum1
   90 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
     &ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with periodic boundary conditions
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 42 flops/particle, 12 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, ek
      real ppart, fxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, ih, nh, nn, mm, lxv
      real qtm, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real anx, any, edgelx, edgely, edgerx, edgery
      real sfxy
!     dimension sfxy(2,MXV*MYV)
      dimension sfxy(2,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,ih,nh, &
!$OMP& x,y,dxp,dyp,amx,amy,dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,&
!$OMP& sfxy,n,s,t) REDUCTION(+:sum2)
      do 110 k = 1, mxyp1
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
! load local fields from global array
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 40 j = 1, npblk
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
   40 continue
! find acceleration
      do 50 j = 1, npblk
      nn = n(j) + 1
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      s(j,1) = dx
      s(j,2) = dy
   50 continue
! new velocity
! !dir$ vector aligned
      do 60 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      dxp = ppart(3,j+joff,k)
      dyp = ppart(4,j+joff,k)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      sum1 = sum1 + (dxp*dxp + dyp*dyp)
! new position and velocity
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = vx
      s(j,4) = vy
   60 continue
! check boundary conditions
! !dir$ vector aligned
      do 70 j = 1, npblk
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
      n(j) = mm
   70 continue
! increment counters
      do 80 j = 1, npblk
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
! find acceleration
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
! new velocity
      dxp = ppart(3,j,k)
      dyp = ppart(4,j,k)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      sum1 = sum1 + (dxp*dxp + dyp*dyp)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new velocity
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
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
  100 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  110 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 120 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  120    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,ny&
     &,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with various boundary conditions
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 50 flops/particle, 2 divides, 2 sqrts, 12 loads, 4 stores
! input: all, output: ppart, ek
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, ci, ek
      real ppart, fxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv*nypmx)
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
      real x, y, dx, dy, vx, vy, acx, acy, p2, dtg
      real sfxy
!     dimension sfxy(2,MXV*MYV)
      dimension sfxy(2,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,2)
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,x,y,   &
!$OMP& dxp,dyp,amx,amy,dx,dy,vx,vy,acx,acy,p2,dtg,sum1,sfxy,n,s,t)      &
!$OMP& REDUCTION(+:sum2)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 70 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 30 j = 1, npblk
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
   30 continue
! find acceleration
      do 40 j = 1, npblk
      nn = n(j) + 1
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      s(j,1) = dx
      s(j,2) = dy
   40 continue
! calculate half impulse
! !dir$ vector aligned
      do 50 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
! time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
! update inverse gamma
      p2 = vx*vx + vy*vy
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = vx
      s(j,4) = vy
   50 continue
! check boundary conditions
! !dir$ vector aligned
      do 60 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      vx = s(j,3)
      vy = s(j,4)
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
   60 continue
   70 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 80 j = nps, nppp
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
! find acceleration
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
! time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
! update inverse gamma
      p2 = vx*vx + vy*vy
      dtg = dt/sqrt(1.0 + p2*ci2)
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
   80 continue
      sum2 = sum2 + sum1
   90 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &ci,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with periodic boundary conditions
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 50 flops/particle, 2 divides, 2 sqrts, 12 loads, 4 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, ci, ek
      real ppart, fxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv*nypmx)
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
      real x, y, dx, dy, vx, vy, acx, acy, p2, dtg
      real anx, any, edgelx, edgely, edgerx, edgery
      real sfxy
!     dimension sfxy(2,MXV*MYV)
      dimension sfxy(2,(mx+1)*(my+1))
! scratch arrays
      integer n
      real s, t
!dir$ attributes align : 64 :: n, s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,2)
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,ih,nh, &
!$OMP& x,y,dxp,dyp,amx,amy,dx,dy,vx,vy,acx,acy,p2,dtg,edgelx,edgely,    &
!$OMP& edgerx,edgery,sum1,sfxy,n,s,t) REDUCTION(+:sum2)
      do 110 k = 1, mxyp1
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
! load local fields from global array
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 40 j = 1, npblk
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
   40 continue
! find acceleration
      do 50 j = 1, npblk
      nn = n(j) + 1
      dx = sfxy(1,nn)*s(j,1) + sfxy(1,nn+1)*s(j,2)
      dy = sfxy(2,nn)*s(j,1) + sfxy(2,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxy(1,mm)*s(j,3) + sfxy(1,mm+1)*s(j,4)
      dy = dy + sfxy(2,mm)*s(j,3) + sfxy(2,mm+1)*s(j,4)
      s(j,1) = dx
      s(j,2) = dy
   50 continue
! calculate half impulse
! !dir$ vector aligned
      do 60 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
! half acceleration
      acx = ppart(3,j+joff,k) + dx
      acy = ppart(4,j+joff,k) + dy
! time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
! update inverse gamma
      p2 = vx*vx + vy*vy
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = vx
      s(j,4) = vy
   60 continue
! check boundary conditions
! !dir$ vector aligned
      do 70 j = 1, npblk
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
      n(j) = mm
   70 continue
! increment counters
      do 80 j = 1, npblk
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
! find acceleration
      dx = amy*(amx*sfxy(1,nn) + dxp*sfxy(1,nn+1))
      dy = amy*(amx*sfxy(2,nn) + dxp*sfxy(2,nn+1))
      dx = dx + dyp*(amx*sfxy(1,nn+lxv) + dxp*sfxy(1,nn+1+lxv)) 
      dy = dy + dyp*(amx*sfxy(2,nn+lxv) + dxp*sfxy(2,nn+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
! time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
! update inverse gamma
      p2 = vx*vx + vy*vy
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = x + vx*dtg
      dy = y + vy*dtg
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! set new momentum
      ppart(3,j,k) = vx
      ppart(4,j,k) = vy
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
  100 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  110 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 120 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  120    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1, &
     &ipbc)
! for 2d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with various boundary conditions.
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 9 flops/particle, 4 loads, 2 stores
! input: all, output: ppart, ek
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t-dt/2)**2+vy(t-dt/2)**2)
! nx/ny = system length in x/y direction
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nx, ny, idimp, nppmx, mxyp1, ipbc
      real dt, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer nppp
      integer j, k
      real edgelx, edgely, edgerx, edgery, x, y, dx, dy, vx, vy
      double precision sum1, sum2
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
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,x,y,dx,dy,vx,vy,sum1)                &
!$OMP& REDUCTION(+:sum2)
      do 20 k = 1, mxyp1
      nppp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
! position
      x = ppart(1,j,k)
      y = ppart(2,j,k)
! velocity
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
! average kinetic energy
      sum1 = sum1 + (vx*vx + vy*vy)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(4,j,k) = -vy
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -vx
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,ny&
     &,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 9 flops/particle, 4 loads, 2 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t-dt/2)**2+vy(t-dt/2)**2)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, mx1, mxyp1, ntmax
      integer irc
      real dt, ek
      real ppart
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=2)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer j, k, m, ih, nh, nn, mm
      real x, y, dx, dy
      real anx, any, edgelx, edgely, edgerx, edgery
! scratch arrays
      integer n
      real s
!dir$ attributes align : 64 :: n, s
      dimension n(npblk), s(npblk,lvect)
      double precision sum1, sum2
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,ih,nh,x,y,dx,dy&
!$OMP& ,edgelx,edgely,edgerx,edgery,sum1,n,s) REDUCTION(+:sum2)
      do 70 k = 1, mxyp1
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
! clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
! position
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
! velocity
      dx = ppart(3,j+joff,k)
      dy = ppart(4,j+joff,k)
! average kinetic energy
      sum1 = sum1 + (dx*dx + dy*dy)
! new position
      s(j,1) = x + dx*dt
      s(j,2) = y + dy*dt
   20 continue
! check boundary conditions
! !dir$ vector aligned
      do 30 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
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
      n(j) = mm
   30 continue
! increment counters
      do 40 j = 1, npblk
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
   40 continue
   50 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 60 j = nps, nppp
! position
      x = ppart(1,j,k)
      y = ppart(2,j,k)
! velocity
      dx = ppart(3,j,k)
      dy = ppart(4,j,k)
! average kinetic energy
      sum1 = sum1 + (dx*dx + dy*dy)
! new position
      dx = ppart(1,j,k) + dx*dt
      dy = ppart(2,j,k) + dy*dt
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
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
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 80 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   80    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,   &
     &mxyp1,ipbc)
! for 2d code, this subroutine updates particle co-ordinates for
! relativistic particles with fixed velocities, with various boundary
! conditions.
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 12 flops/particle, 1 divides, 1 sqrts, 4 loads, 2 stores
! input: all, output: ppart, ek
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum(px(t-dt/2)**2 + py(t-dt/2)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
! nx/ny = system length in x/y direction
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer  nx, ny, idimp, nppmx, mxyp1, ipbc
      real dt, ci, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer nppp
      integer j, k
      real ci2, edgelx, edgely, edgerx, edgery
      real x, y, dx, dy, vx, vy, p2, gam, dtg
      double precision sum1, sum2
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
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,x,y,dx,dy,vx,vy,p2,gam,dtg,sum1)     &
!$OMP& REDUCTION(+:sum2)
      do 20 k = 1, mxyp1
      nppp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
! position
      x = ppart(1,j,k)
      y = ppart(2,j,k)
! momentum
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
! average kinetic energy
      p2 = vx*vx + vy*vy
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = x + vx*dtg
      dy = y + vy*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(4,j,k) = -vy
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -vx
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek, &
     &nx,ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates for
! relativistic particles with fixed velocities, with periodic boundary
! conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored in segmented array
! 9 flops/particle, 1 divides, 1 sqrts, 4 loads, 2 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum(px(t-dt/2)**2 + py(t-dt/2)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nx, ny, mx, my, idimp, nppmx, mx1, mxyp1, ntmax
      integer irc
      real dt, ci, ek
      real ppart
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=2)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer j, k, m, ih, nh, nn, mm
      real ci2, x, y, dx, dy, p2, gam, dtg
      real anx, any, edgelx, edgely, edgerx, edgery
! scratch arrays
      integer n
      real s
!dir$ attributes align : 64 :: n, s
      dimension n(npblk), s(npblk,lvect)
      double precision sum1, sum2
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,ih,nh,x,y,dx,dy&
!$OMP& ,p2,gam,dtg,edgelx,edgely,edgerx,edgery,sum1,n,s)                &
!$OMP& REDUCTION(+:sum2)
      do 70 k = 1, mxyp1
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
! clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
! position
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
! momentum
      dx = ppart(3,j+joff,k)
      dy = ppart(4,j+joff,k)
! average kinetic energy
      p2 = dx*dx + dy*dy
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      s(j,1) = x + dx*dtg
      s(j,2) = y + dy*dtg
   20 continue
! check boundary conditions
! !dir$ vector aligned
      do 30 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
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
      n(j) = mm
   30 continue
! increment counters
      do 40 j = 1, npblk
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
   40 continue
   50 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 60 j = nps, nppp
! position
      x = ppart(1,j,k)
      y = ppart(2,j,k)
! momentum
      dx = ppart(3,j,k)
      dy = ppart(4,j,k)
! average kinetic energy
      p2 = dx*dx + dy*dy
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = x + dx*dtg
      dy = y + dy*dtg
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
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
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 80 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   80    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,&
     &nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 17 flops/particle, 6 loads, 4 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m)=qm*(1.-dx)*(1.-dy)
! q(n+1,m)=qm*dx*(1.-dy)
! q(n,m+1)=qm*(1.-dx)*dy
! q(n+1,m+1)=qm*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! q(j,k) = charge density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), q(nxv*nypmx), kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real x, y, z, w, dxp, dyp, amx, amy
      real sq
!     dimension sq(MXV*MYV)
      dimension sq((mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s
!dir$ attributes align : 64 :: n, mn, s
      dimension n(npblk), mn(lvect), s(npblk,lvect)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,x,y,z,w&
!$OMP& ,dxp,dyp,amx,amy,sq,n,s)
      do 110 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)
      sq(j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   20 continue
! deposit charge within tile to local accumulator
      do 40 j = 1, npblk
!dir$ ivdep
      do 30 i = 1, lvect
      sq(n(j)+mn(i)) = sq(n(j)+mn(i)) + s(j,i)
   30 continue
   40 continue
   50 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 60 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit charge within tile to local accumulator
      x = sq(nn) + amx*amy
      y = sq(nn+1) + dxp*amy
      z = sq(nn+lxv) + amx*dyp
      w = sq(nn+1+lxv) + dxp*dyp
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
   60 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      q(i+noffp+nxv*(j+moffp-1)) = q(i+noffp+nxv*(j+moffp-1)) +         &
     & sq(i+lxv*(j-1))
   70 continue
   80 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 90 i = 2, nn
!$OMP ATOMIC
      q(i+noffp+nxv*moffp) = q(i+noffp+nxv*moffp) + sq(i)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp+nxv*(mm+moffp-1)) = q(i+noffp+nxv*(mm+moffp-1)) +    &
     &sq(i+lxv*(mm-1))
      endif
   90 continue
      nn = min(mx+1,nxv-noffp)
      do 100 j = 1, mm
!$OMP ATOMIC
      q(1+noffp+nxv*(j+moffp-1)) = q(1+noffp+nxv*(j+moffp-1)) +         &
     &sq(1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp+nxv*(j+moffp-1)) = q(nn+noffp+nxv*(j+moffp-1)) +    &
     &sq(nn+lxv*(j-1))
      endif
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine V2PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv&
     &,nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 17 flops/particle, 6 loads, 4 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m)=qm*(1.-dx)*(1.-dy)
! q(n+1,m)=qm*dx*(1.-dy)
! q(n,m+1)=qm*(1.-dx)*dy
! q(n+1,m+1)=qm*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! q(j,k) = charge density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), q(nxv*nypmx), kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real x, y, z, w, dxp, dyp, amx, amy
      real sq
!     dimension sq(MXV*MYV)
      dimension sq((mx+1)*(my+1))
! scratch arrays
      integer n
      real s
!dir$ attributes align : 64 :: n, s
      dimension n(npblk), s(npblk,lvect)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
      lxv = mx + 1
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,mnoff,ipp,joff,nps,nn,mm,x,y,z,w&
!$OMP& ,dxp,dyp,amx,amy,sq,n,s)
      do 100 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)
      sq(j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 40 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   20 continue
! deposit charge within tile to local accumulator
      do 30 j = 1, npblk
      nn = n(j)
      x = sq(nn) + s(j,1)
      y = sq(nn+1) + s(j,2)
      z = sq(nn+lxv) + s(j,3)
      w = sq(nn+1+lxv) + s(j,4)
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
   30 continue
   40 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 50 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit charge within tile to local accumulator
      x = sq(nn) + amx*amy
      y = sq(nn+1) + dxp*amy
      z = sq(nn+lxv) + amx*dyp
      w = sq(nn+1+lxv) + dxp*dyp
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
   50 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
!dir$ ivdep
      do 70 j = 2, mm
      do 60 i = 2, nn
      q(i+noffp+nxv*(j+moffp-1)) = q(i+noffp+nxv*(j+moffp-1)) +         &
     & sq(i+lxv*(j-1))
   60 continue
   70 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 80 i = 2, nn
!$OMP ATOMIC
      q(i+noffp+nxv*moffp) = q(i+noffp+nxv*moffp) + sq(i)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp+nxv*(mm+moffp-1)) = q(i+noffp+nxv*(mm+moffp-1)) +    &
     &sq(i+lxv*(mm-1))
      endif
   80 continue
      nn = min(mx+1,nxv-noffp)
      do 90 j = 1, mm
!$OMP ATOMIC
      q(1+noffp+nxv*(j+moffp-1)) = q(1+noffp+nxv*(j+moffp-1)) +         &
     &sq(1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp+nxv*(j+moffp-1)) = q(nn+noffp+nxv*(j+moffp-1)) +    &
     &sq(nn+lxv*(j-1))
      endif
   90 continue
  100 continue
!$OMP END PARALLEL DO
      return
      end
