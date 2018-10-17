!-----------------------------------------------------------------------
! Fortran Library for pushing electrostatic particles, depositing charge
! and copying particles
! 2D MPI/OpenMP PIC Codes:
! PPPMOVIN2L sorts particles by x,y grid in tiles of mx, my and
!            copies to segmented array ppart
! PPPMOVIN2LP sorts particles by x,ygrid in tiles of mx, my and
!             copies to segmented array ppart for NUMA architectures
! PPPCOPYOUT2 copies segmented particle data ppart to the array part
! PPPCOPYIN2 copies linear array part to segmented particle data ppart
! PPPCHECK2L performs a sanity check to make sure particles sorted
!            by x,y grid in tiles of mx, my, are all within bounds.
! PPGPPUSH2L updates particle co-ordinates and velocities using electric
!            field only, with linear interpolation and various particle
!            boundary conditions
! PPGPPUSHF2L updates particle co-ordinates and velocities using
!             electric field only, with linear interpolation and
!             periodic particle boundary conditions.  also determines
!             list of particles which are leaving each tile
! PPGRPPUSH2L updates relativistic particle co-ordinates and momenta
!             using electric field only, with linear interpolation and
!             various particle boundary conditions
! PPGRPPUSHF2L updates relativistic particle co-ordinates and momenta
!              using electric field only, with linear interpolation and
!              periodic particle boundary conditions.  also determines
!              list of particles which are leaving each tile
! PPGPPUSH2ZF update particle co-ordinate for particles with fixed
!             velocities
! PPGPPUSHF2ZF update particle co-ordinate for particles with fixed
!              velocities with periodic particle boundary conditions.
!              also determines list of particles which are leaving each
!              tile
! PPGRPPUSH2ZF update particle co-ordinates for particles with fixed
!              velocities, for 2d code, and relativistic particles.
! PPGRPPUSHF2ZF update particle co-ordinates for particles with fixed
!               velocities with periodic particle boundary conditions,
!               for 2d code, and relativistic particles.  also
!               determines list of particles which are leaving each tile
! PPGPPOST2L calculates particle charge density using linear
!            interpolation
! SET_SZERO2 zeros out charge density array.
! SET_PVZERO2 zeros out current density array.
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 31, 2018
!-----------------------------------------------------------------------
      subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, &
     &mx,my,mx1,mxyp1,irc)
! this subroutine sorts particles by x,y grid in tiles of
! mx, my and copies to segmented array ppart
! linear interpolation, spatial decomposition in y direction
! input: all except ppart, kpic, output: ppart, kpic
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! nppmx = maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, n, m, mnoff, ip, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,m) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPMOVIN2LP(part,ppart,kpic,kp,npp,noff,nppmx,idimp,   &
     &npmax,mx,my,mx1,mxyp1,irc)
! this subroutine sorts particles by x,y grid in tiles of
! mx, my and copies to segmented array ppart
! designed for NUMA architectures, where memory is associated with the
! processor which first writes a memory location.
! linear interpolation, spatial decomposition in y direction
! input: all except ppart, kpic, kp, output: ppart, kpic, kp
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! kp = original location of reordered particle
! nppmx = maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, kp, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1), kp(nppmx,mxyp1)
! local data
      integer i, j, k, n, m, mnoff, ip, nppp, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile to reorder particles
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         kp(ip,m) = j
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   20 continue
! check for overflow
      if (ierr.gt.0) then
         irc = ierr
         return
      endif
! copy reordered particles
!$OMP PARALLEL DO PRIVATE(i,j,k,m,nppp)
      do 50 k = 1, mxyp1
      nppp = kpic(k)
      do 40 j = 1, nppp
      m = kp(j,k)
      do 30 i = 1, idimp
      ppart(i,j,k) = part(i,m)
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1&
     &,irc)
! for 2d code, this subroutine copies segmented particle data ppart to
! the linear array part
! spatial decomposition in y direction
! input: all except part, npp, irc, output: part, npp, irc
! part(i,j) = i-th coordinate for particle j in partition
! ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
! kpic = number of particles per tile
! npp = number of particles in partition
! npmax = maximum number of particles in each partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mxyp1 = total number of tiles in partition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, npoff, nppp
! check for overflow
      nppp = 0
      do 10 k = 1, mxyp1
      nppp = nppp + kpic(k)
   10 continue
      if (nppp.gt.npmax) then
         irc = nppp
         return
      endif
      npp = nppp
! loop over tiles
      npoff = 0
      do 40 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 30 j = 1, nppp
      do 20 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   20 continue
   30 continue
      npoff = npoff + nppp
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCOPYIN2(part,ppart,kpic,npmax,nppmx,idimp,mxyp1,irc)
! for 2d code, this subroutine copies linear array part to
! segmented particle data ppart, assuming kpic values are known
! spatial decomposition in y direction
! used in resizing segmented particle array ppart if overflow occurs
! input: all except part, irc, output: part, irc
! part(i,j) = i-th coordinate for particle j in partition
! ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
! kpic = number of particles per tile
! npmax = maximum number of particles in each partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mxyp1 = total number of tiles in partition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npmax, nppmx, idimp, mxyp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, npoff, nppp
! check for overflow
      nppp = 0
      do 10 k = 1, mxyp1
      nppp = nppp + kpic(k)
   10 continue
      if (nppp.gt.npmax) then
         irc = nppp
         return
      endif
! loop over tiles
      npoff = 0
      do 40 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 30 j = 1, nppp
      do 20 i = 1, idimp
      ppart(i,j,k) = part(i,j+npoff)
   20 continue
   30 continue
      npoff = npoff + nppp
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1&
     &,myp1,irc)
! this subroutine performs a sanity check to make sure particles sorted
! by x,y grid in tiles of mx, my, are all within bounds.
! tiles are assumed to be arranged in 2D linear memory
! input: all except irc
! output: irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! kpic(k) = number of reordered output particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension kpic(mx1*myp1)
! local data
      integer mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ist,edgelx,edgely,edgerx,     &
!$OMP& edgery,dx,dy)
      do 20 k = 1, mxyp1
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
! loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
! find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny,mx, &
     &my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with various boundary conditions
! OpenMP version using guard cells, for distributed data
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
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real sfxy
      dimension sfxy(2,MXV,MYV)
!     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
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
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,dx&
!$OMP& ,dy,vx,vy,sum1,sfxy) REDUCTION(+:sum2)
      do 40 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, nppp
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
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
! new velocity
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
! average kinetic energy
      vx = vx + dx
      vy = vy + dy
      sum1 = sum1 + (vx*vx + vy*vy)
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
! new position
      dx = x + dx*dt
      dy = y + dy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt,ek&
     &,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with periodic boundary conditions
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
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
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtm, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real anx, any, edgelx, edgely, edgerx, edgery
      real sfxy
      dimension sfxy(2,MXV,MYV)
!     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,&
!$OMP& amy,dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,sfxy)           &
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mxyp1
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
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, nppp
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
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
! new velocity
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
! average kinetic energy
      vx = vx + dx
      vy = vy + dy
      sum1 = sum1 + (vx*vx + vy*vy)
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
! new position
      dx = x + dx*dt
      dy = y + dy*dt
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
         do 60 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ci,ek,nx,ny,&
     &mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
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
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, acx, acy, p2, dtg
      real sfxy
      dimension sfxy(2,MXV,MYV)
!     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
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
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,dx&
!$OMP& ,dy,vx,vy,acx,acy,p2,dtg,sum1,sfxy) REDUCTION(+:sum2)
      do 40 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, nppp
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
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
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
      dx = acx + dx
      dy = acy + dy
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
! update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtg
      dy = y + dy*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
     &ci,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates and
! momenta using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with periodic boundary conditions
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
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
      dimension ppart(idimp,nppmx,mxyp1), fxy(2,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, acx, acy, p2, dtg
      real anx, any, edgelx, edgely, edgerx, edgery
      real sfxy
      dimension sfxy(2,MXV,MYV)
!     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,ih,nh,x,y,dxp,dyp,amx,&
!$OMP& amy,dx,dy,vx,vy,acx,acy,p2,dtg,edgelx,edgely,edgerx,edgery,sum1, &
!$OMP& sfxy) REDUCTION(+:sum2)
      do 50 k = 1, mxyp1
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
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
! clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, nppp
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
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
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
      dx = acx + dx
      dy = acy + dy
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
! update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
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
         do 60 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   60    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPUSH2ZF(ppart,kpic,dt,ek,nx,ny,idimp,nppmx,mxyp1,  &
     &ipbc)
! for 2d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with various boundary conditions.
! OpenMP version using guard cells, for distributed data
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
      real edgelx, edgely, edgerx, edgery, dx, dy
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
!$OMP PARALLEL DO PRIVATE(j,k,nppp,dx,dy,sum1) REDUCTION(+:sum2)
      do 20 k = 1, mxyp1
      nppp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
! velocity
      dx = ppart(3,j,k)
      dy = ppart(4,j,k)
! average kinetic energy
      sum1 = sum1 + (dx*dx + dy*dy)
! new position
      dx = ppart(1,j,k) + dx*dt
      dy = ppart(2,j,k) + dy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
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
      subroutine PPGPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ek,nx,ny,&
     &mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates for
! particles with fixed velocities, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
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
      integer noffp, moffp, nppp
      integer j, k, ih, nh, nn, mm
      real dx, dy
      real anx, any, edgelx, edgely, edgerx, edgery
      double precision sum1, sum2
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ih,nh,dx,dy,edgelx,edgely,    &
!$OMP& edgerx,edgery,sum1) REDUCTION(+:sum2)
      do 30 k = 1, mxyp1
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
      do 20 j = 1, nppp
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
   20 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 40 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRPPUSH2ZF(ppart,kpic,dt,ci,ek,nx,ny,idimp,nppmx,    &
     &mxyp1,ipbc)
! for 2d code, this subroutine updates particle co-ordinates for
! relativistic particles with fixed velocities, with various boundary
! conditions.
! OpenMP version using guard cells, for distributed data
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
      real ci2, edgelx, edgely, edgerx, edgery, dx, dy, p2, gam, dtg
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
!$OMP PARALLEL DO PRIVATE(j,k,nppp,dx,dy,p2,gam,dtg,sum1)               &
!$OMP& REDUCTION(+:sum2)
      do 20 k = 1, mxyp1
      nppp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, nppp
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
      dx = ppart(1,j,k) + dx*dtg
      dy = ppart(2,j,k) + dy*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
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
      subroutine PPGRPPUSHF2ZF(ppart,kpic,ncl,ihole,noff,nyp,dt,ci,ek,nx&
     &,ny,mx,my,idimp,nppmx,mx1,mxyp1,ntmax,irc)
! for 2d code, this subroutine updates particle co-ordinates for
! relativistic particles with fixed velocities, with periodic boundary
! conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
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
      integer noffp, moffp, nppp
      integer j, k, ih, nh, nn, mm
      real ci2, dx, dy, p2, gam, dtg
      real anx, any, edgelx, edgely, edgerx, edgery
      double precision sum1, sum2
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ih,nh,dx,dy,p2,gam,dtg,edgelx,&
!$OMP& edgely,edgerx,edgery,sum1) REDUCTION(+:sum2)
      do 30 k = 1, mxyp1
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
      do 20 j = 1, nppp
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
      dx = ppart(1,j,k) + dx*dtg
      dy = ppart(2,j,k) + dy*dtg
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
   20 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 40 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv, &
     &nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! OpenMP version using guard cells, for distributed data
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
      dimension ppart(idimp,nppmx,mxyp1), q(nxv,nypmx), kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
!     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,  &
!$OMP& sq)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit charge within tile to local accumulator
      x = sq(nn,mm) + amx*amy
      y = sq(nn+1,mm) + dxp*amy
      sq(nn,mm) = x
      sq(nn+1,mm) = y
      x = sq(nn,mm+1) + amx*dyp
      y = sq(nn+1,mm+1) + dxp*dyp
      sq(nn,mm+1) = x
      sq(nn+1,mm+1) = y
   30 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noffp,j+moffp) = q(i+noffp,j+moffp) + sq(i,j)
   40 continue
   50 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp) = q(i+noffp,1+moffp) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp) = q(i+noffp,mm+moffp) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp) = q(1+noffp,j+moffp) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp,j+moffp) = q(nn+noffp,j+moffp) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine SET_PSZERO2(q,mx,my,nxv,nypmx,mx1,myp1)
! for 2d code, this subroutine zeros out charge density array.
! for Intel NUMA architecture with first touch policy, this associates
! array segments with appropriate threads
! OpenMP version
! input: all, output: q
! q(j,k) = charge density at grid point j,k
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = second dimension of charge array
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction-1)/my+1
      implicit none
      integer mx, my, nxv, nypmx, mx1, myp1
      real q
      dimension q(nxv,nypmx)
! local data
      integer mxyp1, noffp, moffp
      integer i, j, k, nn, mm
      mxyp1 = mx1*myp1
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,noffp,moffp,nn,mm)
      do 30 k = 1, mxyp1
      j = (k - 1)/mx1
      moffp = my*j
      mm = my
      if ((j+1).eq.myp1) mm = nypmx - moffp
      i = k - mx1*j
      noffp = mx*(i - 1)
      nn = mx
      if (i.eq.mx1) nn = nxv - noffp
! zero charge in global array
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      q(i+noffp,j+moffp) = 0.0
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine SET_PVZERO2(cu,mx,my,ndim,nxv,nypmx,mx1,myp1)
! for 2d code, this subroutine zeros out current density array.
! for Intel NUMA architecture with first touch policy, this associates
! array segments with appropriate threads
! OpenMP version
! input: all, output: cu
! cu(m,j,k) = current density at grid point m,j,k
! mx/my = number of grids in sorting cell in x/y
! ndim = first dimension of current array
! nxv = second dimension of current array, must be >= nx+1
! nypmx = third dimension of current array
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
      implicit none
      integer mx, my, ndim, nxv, nypmx, mx1, myp1
      real cu
      dimension cu(ndim,nxv,nypmx)
! local data
      integer mxyp1, noffp, moffp
      integer i, j, k, m, nn, mm
      mxyp1 = mx1*myp1
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,m,noffp,moffp,nn,mm)
      do 40 k = 1, mxyp1
      j = (k - 1)/mx1
      moffp = my*j
      mm = my
      if ((j+1).eq.myp1) mm = nypmx - moffp
      i = k - mx1*j
      noffp = mx*(i - 1)
      nn = mx
      if (i.eq.mx1) nn = nxv - noffp
! zero current in global array
      do 30 j = 1, mm
      do 20 i = 1, nn
!dir$ ivdep
      do 10 m = 1, ndim
      cu(m,i+noffp,j+moffp) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
!$OMP END PARALLEL DO
      return
      end
