!-----------------------------------------------------------------------
! Fortran Library for pushing electrostatic particles and depositing
! charge
! 3D Vector/MPI/OpenMP PIC Codes:
! VPPGPPUSH32L updates particle co-ordinates and velocities using
!              electric field only, with linear interpolation and
!              various particle boundary conditions
! VPPGPPUSHF32L updates particle co-ordinates and velocities using
!               electric field only, with linear interpolation and
!               periodic particle boundary conditions.  also determines
!               list of particles which are leaving each tile
! VPPGRPPUSH32L updates relativistic particle co-ordinates and momenta
!               using electric field only, with linear interpolation and
!               various particle boundary conditions
! VPPGRPPUSHF32L updates relativistic particle co-ordinates and
!                momenta using electric field only, with linear
!                interpolation and periodic particle boundary conditions
!                also determines list of particles which are leaving each
!                tile
! VPPGPPUSH32ZF update particle co-ordinate for particles with fixed
!               velocities
! VPPGPPUSHF32ZF update particle co-ordinate for particles with fixed
!                velocities with periodic particle boundary conditions.
!                also determines list of particles which are leaving
!                each tile
! VPPGRPPUSH32ZF update particle co-ordinates for particles with fixed
!                velocities, for 3d code, and relativistic particles.
! VPPGRPPUSHF32ZF update particle co-ordinates for particles with fixed
!                 velocities with periodic particle boundary conditions,
!                 for 3d code, and relativistic particles.  also
!                 determines list of particles which are leaving each
!                 tile
! VPPGPPOST32L calculates particle charge density using linear
!              interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: june 28, 2018
!-----------------------------------------------------------------------
      subroutine VPPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,idimp,&
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with various boundary conditions.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 90 flops/particle, 30 loads, 6 stores
! input: all, output: part, ek
! equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
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
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
! (vz(t+dt/2)+vz(t-dt/2))**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field array, must be >= nx+1
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
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
! loopv = (0,1) = execute (scalar,vector) loop
      integer npblk, lvect, loopv
      parameter(npblk=32,lvect=8,loopv=1)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, lxv, lxyv, nxyv
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
!     dimension sfxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
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
      qtm = qbm*dt
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
!$OMP& nps,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,&
!$OMP& sum1,sfxyz,n,s,t) REDUCTION(+:sum2)
      do 100 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 80 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 40 j = 1, npblk
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
   40 continue
! find acceleration
      do 50 j = 1, npblk
      nn = n(j) + 1
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
   50 continue
! new velocity
! !dir$ vector aligned
      do 60 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(4,j+joff,l)
      dyp = ppart(5,j+joff,l)
      dzp = ppart(6,j+joff,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position and velocity
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   60 continue
! check boundary conditions
! !dir$ vector aligned
      do 70 j = 1, npblk
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
   70 continue
   80 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 90 j = nps, nppp
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
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
! new velocity
      dxp = ppart(4,j,l)
      dyp = ppart(5,j,l)
      dzp = ppart(6,j,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
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
   90 continue
      sum2 = sum2 + sum1
  100 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm, &
     &dt,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,     &
     &mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 90 flops/particle, 30 loads, 6 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
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
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
! (vz(t+dt/2)+vz(t-dt/2))**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field array, must be >= nx+1
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
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv*nypmx*nzpmx)
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
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real sfxyz
!     dimension sfxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
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
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,   &
!$OMP& nps,ih,nh,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,&
!$OMP& vy,vz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,sum1,sfxyz,n,s,t)&
!$OMP& REDUCTION(+:sum2)
      do 120 l = 1, mxyzp1
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
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 100 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 50 j = 1, npblk
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
   50 continue
! find acceleration
      do 60 j = 1, npblk
      nn = n(j) + 1
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
   60 continue
! new velocity
! !dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(4,j+joff,l)
      dyp = ppart(5,j+joff,l)
      dzp = ppart(6,j+joff,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position and velocity
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   70 continue
! check boundary conditions
! !dir$ vector aligned
      do 80 j = 1, npblk
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
   80 continue
! increment counters
      do 90 j = 1, npblk
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
   90 continue
  100 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 110 j = nps, nppp
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
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
! new velocity
      dxp = ppart(4,j,l)
      dyp = ppart(5,j,l)
      dzp = ppart(6,j,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
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
! use periodic boundary conditions and check for roundoff error
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
  110 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  120 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 130 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  130    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,  &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &idds,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with various boundary conditions.
! for distributed data with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells,
! data read in tiles
! particles stored in segmented array
! 100 flops/particle, 1 divide, 1 sqrt, 30 loads, 6 stores
! input: all, output: part, ek
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
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
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field array, must be >= nx+1
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
      real qbm, dt, ci, ek
      real ppart, fxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv*nypmx*nzpmx)
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
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz, acx, acy, acz, p2, dtg
      real sfxyz
!     dimension sfxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
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
!$OMP& ,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,acx&
!$OMP& ,acy,acz,p2,dtg,sum1,sfxyz,n,s,t) REDUCTION(+:sum2)
      do 100 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      ipp = 0
      if (loopv==1) ipp = nppp/npblk
! outer loop over number of full blocks
      do 80 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 40 j = 1, npblk
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
   40 continue
! find acceleration
      do 50 j = 1, npblk
      nn = n(j) + 1
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
   50 continue
! calculate half impulse
! !dir$ vector aligned
      do 60 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
      vz = acz + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = z + vz*dtg
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   60 continue
! check boundary conditions
! !dir$ vector aligned
      do 70 j = 1, npblk
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
   70 continue
   80 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 90 j = nps, nppp
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
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
! time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
      vz = acz + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dt/sqrt(1.0 + p2*ci2)
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
   90 continue
      sum2 = sum2 + sum1
  100 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,&
     &dt,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,  &
     &mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells,
! data read in tiles
! particles stored in segmented array
! 100 flops/particle, 1 divide, 1 sqrt, 30 loads, 6 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
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
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field array, must be >= nx+1
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
      real qbm, dt, ci, ek
      real ppart, fxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv*nypmx*nzpmx)
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
      real qtmh, ci2, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1
      real dx, dy, dz, vx, vy, vz, acx, acy, acz, p2, dtg
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real sfxyz
!     dimension sfxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
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
!$OMP& ,ih,nh,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,&
!$OMP& vz,acx,acy,acz,p2,dtg,edgelx,edgely,edgelz,edgerx,edgery,edgerz, &
!$OMP& sum1,sfxyz,n,s,t) REDUCTION(+:sum2)
      do 120 l = 1, mxyzp1
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
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 100 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 50 j = 1, npblk
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
   50 continue
! find acceleration
      do 60 j = 1, npblk
      nn = n(j) + 1
      dx = sfxyz(1,nn)*s(j,1) + sfxyz(1,nn+1)*s(j,2)
      dy = sfxyz(2,nn)*s(j,1) + sfxyz(2,nn+1)*s(j,2)
      dz = sfxyz(3,nn)*s(j,1) + sfxyz(3,nn+1)*s(j,2)
      mm = nn + lxv
      dx = dx + sfxyz(1,mm)*s(j,3) + sfxyz(1,mm+1)*s(j,4)
      dy = dy + sfxyz(2,mm)*s(j,3) + sfxyz(2,mm+1)*s(j,4)
      dz = dz + sfxyz(3,mm)*s(j,3) + sfxyz(3,mm+1)*s(j,4)
      ll = nn + lxyv
      vx = sfxyz(1,ll)*s(j,5) + sfxyz(1,ll+1)*s(j,6)
      vy = sfxyz(2,ll)*s(j,5) + sfxyz(2,ll+1)*s(j,6)
      vz = sfxyz(3,ll)*s(j,5) + sfxyz(3,ll+1)*s(j,6)
      nn = ll + lxv
      vx = vx + sfxyz(1,nn)*s(j,7) + sfxyz(1,nn+1)*s(j,8)
      vy = vy + sfxyz(2,nn)*s(j,7) + sfxyz(2,nn+1)*s(j,8)
      vz = vz + sfxyz(3,nn)*s(j,7) + sfxyz(3,nn+1)*s(j,8)
      s(j,1) = dx + vx
      s(j,2) = dy + vy
      s(j,3) = dz + vz
   60 continue
! calculate half impulse
! !dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dx = qtmh*s(j,1)
      dy = qtmh*s(j,2)
      dz = qtmh*s(j,3)
! half acceleration
      acx = ppart(4,j+joff,l) + dx
      acy = ppart(5,j+joff,l) + dy
      acz = ppart(6,j+joff,l) + dz
! time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
      vz = acz + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position and momentum
      s(j,1) = x + vx*dtg
      s(j,2) = y + vy*dtg
      s(j,3) = z + vz*dtg
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   70 continue
! check boundary conditions
! !dir$ vector aligned
      do 80 j = 1, npblk
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
! set new momentum
      ppart(4,j+joff,l) = s(j,4)
      ppart(5,j+joff,l) = s(j,5)
      ppart(6,j+joff,l) = s(j,6)
      n(j) = mm
   80 continue
! increment counters
      do 90 j = 1, npblk
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
   90 continue
  100 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 110 j = nps, nppp
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
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(4,j,l) + dx
      acy = ppart(5,j,l) + dy
      acz = ppart(6,j,l) + dz
! time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      vx = acx + dx
      vy = acy + dy
      vz = acz + dz
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dt/sqrt(1.0 + p2*ci2)
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
! use periodic boundary conditions and check for roundoff error
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
  110 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  120 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 130 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  130    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine V2PPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,idimp&
     &,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc&
     &)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with various boundary conditions.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 90 flops/particle, 30 loads, 6 stores
! input: all, output: part, ek
! equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
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
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
! (vz(t+dt/2)+vz(t-dt/2))**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of field array, must be >= nx+1
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
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, lxv, lxyv, nxyv
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
!     dimension sfxyz(3,MXV*MYV*MZV)
      dimension sfxyz(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
      double precision sum1, sum2
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
      qtm = qbm*dt
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
!$OMP& nps,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,&
!$OMP& sum1,sfxyz,n,s,t) REDUCTION(+:sum2)
      do 110 l = 1, mxyzp1
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
   40 continue
! find acceleration
      do 60 j = 1, npblk
      dx = 0.0
      dy = 0.0
      dz = 0.0
!dir$ ivdep
      do 50 i = 1, lvect
      dx = dx + sfxyz(1,n(j)+mn(i))*s(j,i)
      dy = dy + sfxyz(2,n(j)+mn(i))*s(j,i)
      dz = dz + sfxyz(3,n(j)+mn(i))*s(j,i)
   50 continue
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
   60 continue
! new velocity
! !dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(4,j+joff,l)
      dyp = ppart(5,j+joff,l)
      dzp = ppart(6,j+joff,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position and velocity
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   70 continue
! check boundary conditions
! !dir$ vector aligned
      do 80 j = 1, npblk
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
   80 continue
   90 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 100 j = nps, nppp
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
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
! new velocity
      dxp = ppart(4,j,l)
      dyp = ppart(5,j,l)
      dzp = ppart(6,j,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
! average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
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
  100 continue
      sum2 = sum2 + sum1
  110 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSH32ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,ny,nz,   &
     &mxyzp1,ipbc)
! for 3d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with various boundary conditions.
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version
! data read in tiles
! particles stored in segmented array
! 13 flops/particle, 6 loads, 3 stores
! input: all, output: part, ek
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t-dt/2)**2+(vy(t-dt/2)**2+(vz(t-dt/2)**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mxyzp1 = mx1*myp1*mzp1,
! where mx1 = (system length in x direction - 1)/mx + 1
! where myp1 = (partition length in y direction - 1)/my + 1
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mxyzp1, ipbc
      real dt, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1)
! local data
      integer nppp
      integer j, l
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dx, dy, dz, vx, vy, vz
      double precision sum1, sum2
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
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,l,nppp,x,y,z,dx,dy,dz,vx,vy,vz,sum1)        &
!$OMP& REDUCTION(+:sum2)
      do 20 l = 1, mxyzp1
      nppp = kpic(l)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
! position
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
! velocity
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
! average kinetic energy
      sum1 = sum1 + (vx*vx + vy*vy + vz*vz)
! new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ek,   &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version
! data read in tiles
! particles stored in segmented array
! 13 flops/particle, 6 loads, 3 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! kpic(l) = number of particles in tile l
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t-dt/2)**2+(vy(t-dt/2)**2+(vz(t-dt/2)**2)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
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
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, myp1, mxyzp1
      integer ntmax, idds, irc
      real dt, ek
      real ppart
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer j, k, l, m, ih, nh, nn, mm, ll
      real x, y, z, dx, dy, dz
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
! scratch arrays
      integer n
      real s
!dir$ attributes align : 64 :: n, s
      dimension n(npblk), s(npblk,lvect)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,ih,nh,nn,mm, &
!$OMP& ll,x,y,z,dx,dy,dz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,sum1,&
!$OMP& n,s) REDUCTION(+:sum2)
      do 70 l = 1, mxyzp1
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
! clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
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
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
! velocity
      dx = ppart(4,j+joff,l)
      dy = ppart(5,j+joff,l)
      dz = ppart(6,j+joff,l)
! average kinetic energy
      sum1 = sum1 + (dx*dx + dy*dy+ dz*dz)
! new position
      s(j,1) = x + dx*dt
      s(j,2) = y + dy*dt
      s(j,3) = z + dz*dt
   20 continue
! check boundary conditions
! !dir$ vector aligned
      do 30 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
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
      n(j) = mm
   30 continue
! increment counters
      do 40 j = 1, npblk
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
   40 continue
   50 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 60 j = nps, nppp
! position
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
! velocity
      dx = ppart(4,j,l)
      dy = ppart(5,j,l)
      dz = ppart(6,j,l)
! average kinetic energy
      sum1 = sum1 + (dx*dx + dy*dy+ dz*dz)
! new position
      dx = x + dx*dt
      dy = y + dy*dt
      dz = z + dz*dt
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
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
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
   70 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 80 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
   80    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSH32ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,ny,nz&
     &,mxyzp1,ipbc)
! for 3d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with various boundary conditions.
! with various boundary conditions.
! for distributed data with 2D spatial decomposition
! vectorizable/OpenMP version
! data read in tiles
! particles stored in segmented array
! 16 flops/particle, 1 divide, 1 sqrt, 6 loads, 3 stores
! input: all, output: part, ek
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum(px(t-dt/2)**2 + py(t-dt/2)**2 + pz(t-dt/2))**2/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mxyzp1 = mx1*myp1*mzp1,
! where mx1 = (system length in x direction - 1)/mx + 1
! where myp1 = (partition length in y direction - 1)/my + 1
! where mzp1 = (partition length in z direction - 1)/mz + 1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mxyzp1, ipbc
      real dt, ci, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1)
! local data
      integer nppp
      integer j, l
      real ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dx, dy, dz, vx, vy, vz, p2, gam, dtg
      double precision sum1, sum2
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
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,l,nppp,x,y,z,dx,dy,dz,vx,vy,vz,p2,gam,dtg,  &
!$OMP& sum1) REDUCTION(+:sum2)
      do 20 l = 1, mxyzp1
      nppp = kpic(l)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
! position
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
! momentum
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
! average kinetic energy
      p2 = vx*vx + vy*vy + vz*vz
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = x + vx*dtg
      dy = y + vy*dtg
      dz = z + vz*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -vz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -vy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRPPUSHF32ZF(ppart,kpic,ncl,ihole,noff,nyzp,dt,ci,  &
     &ek,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates for
! relativistic particles with fixed velocities, with periodic boundary
! conditions.
! also determines list of particles which are leaving this tile
! for distributed data with 2D spatial decomposition
! vectorizable/OpenMP version
! data read in tiles
! particles stored in segmented array
! 16 flops/particle, 1 divide, 1 sqrt, 6 loads, 3 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
! equations used are: 
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! kpic(l) = number of particles in tile l
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum(px(t-dt/2)**2 + py(t-dt/2)**2 + pz(t-dt/2))**2/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
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
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, myp1, mxyzp1
      integer ntmax, idds, irc
      real dt, ci, ek
      real ppart
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer j, k, l, m, ih, nh, nn, mm, ll
      real x, y, z, dx, dy, dz, ci2, p2, gam, dtg
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
! scratch arrays
      integer n
      real s
!dir$ attributes align : 64 :: n, s
      dimension n(npblk), s(npblk,lvect)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,ih,nh,nn,mm, &
!$OMP& ll,x,y,z,dx,dy,dz,p2,gam,dtg,edgelx,edgely,edgelz,edgerx,edgery, &
!$OMP& edgerz,sum1,n,s) REDUCTION(+:sum2)
      do 70 l = 1, mxyzp1
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
! clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
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
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
! momentum
      dx = ppart(4,j+joff,l)
      dy = ppart(5,j+joff,l)
      dz = ppart(6,j+joff,l)
! average kinetic energy
      p2 = dx*dx + dy*dy + dz*dz
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      s(j,1) = x + dx*dtg
      s(j,2) = y + dy*dtg
      s(j,3) = z + dz*dtg
   20 continue
! check boundary conditions
! !dir$ vector aligned
      do 30 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
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
      n(j) = mm
   30 continue
! increment counters
      do 40 j = 1, npblk
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
   40 continue
   50 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 60 j = nps, nppp
! position
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
! momentum
      dx = ppart(4,j,l)
      dy = ppart(5,j,l)
      dz = ppart(6,j,l)
! average kinetic energy
      p2 = dx*dx + dy*dy + dz*dz
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = x + dx*dtg
      dy = y + dy*dtg
      dz = z + dz*dtg
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
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
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
   70 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 80 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
   80    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my,mz,&
     &nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 33 flops/particle, 11 loads, 8 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
! q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
! q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
! q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
! q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
! q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
! q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
! q(n+1,m+1,l+1)=qm*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! q(j,k,l) = charge density at grid point (j,kk,ll),
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm
      real ppart, q
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), q(nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real x, y, z, w, dxp, dyp, dzp, amx, amy, amz, dx1
      real sq
!     dimension sq(MXV*MYV*MZV)
      dimension sq((mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s
!dir$ attributes align : 64 :: n, mn, s
      dimension n(npblk), mn(lvect), s(npblk,lvect)
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
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,mnoff,lnoff,ipp,joff,nps&
!$OMP& ,nn,mm,ll,nm,lm,x,y,z,w,dxp,dyp,dzp,amx,amy,amz,dx1,sq,n,s)
      do 170 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1)
      lnoff = loffp + noff(2)
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)*(mz+1)
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
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = qm - dxp
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
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit charge within tile to local accumulator
      x = sq(nn) + amx*amz
      y = sq(nn+1) + amy*amz
      z = sq(nn+lxv) + dyp*amz
      w = sq(nn+1+lxv) + dx1*amz
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
      mm = nn + lxyv
      x = sq(mm) + amx*dzp
      y = sq(mm+1) + amy*dzp
      z = sq(mm+lxv) + dyp*dzp
      w = sq(mm+1+lxv) + dx1*dzp
      sq(mm) = x
      sq(mm+1) = y
      sq(mm+lxv) = z
      sq(mm+1+lxv) = w
   60 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 90 k = 2, ll
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      q(i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                     &
     &q(i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                     &
     & sq(i+lxv*(j-1)+lxyv*(k-1))
   70 continue
   80 continue
   90 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 110 j = 2, mm
      do 100 i = 2, nn
!$OMP ATOMIC
      q(i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                           &
     &q(i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + sq(i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         q(i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =                 &
     &   q(i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                   &
     &   + sq(i+lxv*(j-1)+lxyv*(lm-1))
      endif
  100 continue
  110 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 140 k = 1, ll
      do 120 i = 2, nn
!$OMP ATOMIC
      q(i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                           &
     &q(i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + sq(i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &   q(i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                   &
     &   + sq(i+lxv*(mm-1)+lxyv*(k-1))
      endif
  120 continue
      do 130 j = 1, mm
!$OMP ATOMIC
      q(1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                     &
     &q(1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                       &
     &+ sq(1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         q(nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &   q(nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &   + sq(nm+lxv*(j-1)+lxyv*(k-1))
      endif
  130 continue
  140 continue
      if (lm > mz) then
         do 150 i = 2, nn
!$OMP ATOMIC
         q(i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                       &
     &   q(i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) + sq(i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            q(i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &      q(i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))               &
     &      + sq(i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  150    continue
         do 160 j = 1, mm
!$OMP ATOMIC
         q(1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =                 &
     &   q(1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                   &
     &   + sq(1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            q(nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &      q(nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &      + sq(nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  160    continue
      endif
  170 continue
!$OMP END PARALLEL DO
      return
      end
