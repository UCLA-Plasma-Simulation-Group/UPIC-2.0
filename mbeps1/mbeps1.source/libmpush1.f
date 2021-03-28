!-----------------------------------------------------------------------
! Fortran Library for pushing electrostatic particles, depositing charge
! and copying particles
! 1D OpenMP PIC Codes:
! PPMOVIN1L sorts particles by x grid in tiles of mx and copies to
!           segmented array ppart
! PPCOPYOUT1 copies segmented particle data ppart to linear array part
! PPCOPYIN1 copies linear array part to segmented particle data ppart
! PPCHECK1L performs a sanity check to make sure particles sorted
!           by x grid in tiles of mx, are all within bounds.
! GPPUSH1L updates particle co-ordinates and velocities using electric
!          field only, with linear interpolation and various particle
!          boundary conditions
! GPPUSHF1L updates particle co-ordinates and velocities using
!           electric field only, with linear interpolation and
!           periodic particle boundary conditions.  also determines
!           list of particles which are leaving each tile
! GRPPUSH1L updates relativistic particle co-ordinates and momenta
!           using electric field only, with linear interpolation and
!           various particle boundary conditions
! GRPPUSHF1L updates relativistic particle co-ordinates and velocities
!            using electric field only, with linear interpolation and
!            periodic particle boundary conditions.  also determines
!            list of particles which are leaving each tile
! PPUSH1ZF update particle co-ordinate for particles with fixed
!          velocities
! PPUSHF1ZF update particle co-ordinate for particles with fixed
!           velocities with periodic particle boundary conditions.  also
!           determines list of particles which are leaving each tile
! RPPUSH1ZF update particle co-ordinates for particles with fixed
!           velocities, for 1d code, and relativistic particles.
! RPPUSHF1ZF update particle co-ordinates for particles with fixed
!            velocities with periodic particle boundary conditions, for 
!            1d code, and relativistic particles.  also determines list
!            of particles which are leaving each tile
! GPPOST1L calculates particle charge density using linear interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 4, 2021
!-----------------------------------------------------------------------
      subroutine PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc)
! this subroutine sorts particles by x grid in tiles of mx and copies
! to segmented array ppart
! linear interpolation
! input: all except ppart, kpic, output: ppart, kpic, irc
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic = output number of particles per tile
! nppmx = rmaximum number of particles in tile
! idimp = size of phase space = 2
! nop = number of particles
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, mx1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer i, j, k, n, ip, ierr
      ierr = 0
! clear counter array
      do 10 k = 1, mx1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile and reorder particles
      do 30 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      ip = kpic(n) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,n) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(n) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPCOPYOUT1(part,ppart,kpic,np,nop,nppmx,idimp,mx1,irc)
! for 1d code, this subroutine copies segmented particle data ppart to
! the linear array part
! input: all except part, irc, output: part, np, irc
! part(i,j) = i-th coordinate for particle j
! ppart(i,j,k) = i-th coordinate for particle j in tile k
! kpic = number of particles per tile
! np = number of particles copied
! nop = size of linear particle array
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! mx1 = total number of tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer np, nop, nppmx, idimp, mx1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,nop), ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer i, j, k, npoff, npp, ne, ierr
      npoff = 0
      ierr = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
      ne = npp + npoff
      if (ne.gt.nop) ierr = max(ierr,ne-nop)
      if (ierr.gt.0) npp = 0
! loop over particles in tile
      do 20 j = 1, npp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + npp
   30 continue
      np = ne
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPCOPYIN1(part,ppart,kpic,nop,nppmx,idimp,mx1,irc)
! for 1d code, this subroutine copies linear array part to
! segmented particle data ppart, assuming kpic values are known
! used in resizing segmented particle array ppart if overflow occurs
! input: all except ppart, irc, output: ppart, irc
! part(i,j) = i-th coordinate for particle j
! ppart(i,j,k) = i-th coordinate for particle j in tile k
! kpic = number of particles per tile
! nop = size of linear particle array
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! mx1 = total number of tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nop, nppmx, idimp, mx1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,nop), ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer i, j, k, npoff, npp, ne, ierr
      npoff = 0
      ierr = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
      ne = npp + npoff
      if (ne.gt.nop) ierr = max(ierr,ne-nop)
      if (ierr.gt.0) npp = 0
! loop over particles in tile
      do 20 j = 1, npp
      do 10 i = 1, idimp
      ppart(i,j,k) = part(i,j+npoff)
   10 continue
   20 continue
      npoff = npoff + npp
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
! this subroutine performs a sanity check to make sure particles sorted
! by x grid in tiles of mx, are all within bounds.
! tiles are assumed to be arranged in 1D linear memory
! input: all except irc
! output: irc
! ppart(1,n,k) = position x of particle n in tile k
! kpic(k) = number of reordered output particles in tile k
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, mx1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer noff, npp, j, k, ist, nn
      real edgelx, edgerx, dx
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,noff,npp,nn,ist,edgelx,edgerx,dx)         &
!$OMP& SCHEDULE(dynamic)
      do 20 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
! loop over particles in tile
      do 10 j = 1, npp
      dx = ppart(1,j,k)
! find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GPPUSH1L(ppart,fx,kpic,qbm,dt,ek,idimp,nppmx,nx,mx,nxv,&
     &mx1,ipbc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, with various boundary conditions.
! OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 16 flops/particle, 4 loads, 2 stores
! input: all, output: ppart, ek
! equations used are:
! v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
! and x(t+dt) = x(t) + v(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! kpic = number of particles per tile
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real qbm, dt, ek
      real ppart, fx
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fx(nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtm, edgelx, edgerx, x, dx, vx
      real sfx
      dimension sfx(MXV)
!     dimension sfx(mx+1)
      double precision sum1, sum2
      qtm = qbm*dt
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dx,vx,sum1,sfx) REDUCTION(+:sum2)      &
!$OMP& SCHEDULE(dynamic)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      do 10 j = 1, min(mx,nx-noff)+1
      sfx(j) = fx(j+noff)
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = x - real(nn)
      nn = nn - noff + 1
! find acceleration
      dx = (1.0 - dx)*sfx(nn) + dx*sfx(nn+1)
! new velocity
      vx = ppart(2,j,k)
      dx = vx + qtm*dx
! average kinetic energy
      vx = vx + dx
      sum1 = sum1 + vx*vx
      ppart(2,j,k) = dx
! new position
      dx = x + dx*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   20 continue
      sum2 = sum2 + sum1
   30 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + .125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,idimp,nppmx&
     &,nx,mx,nxv,mx1,ntmax,irc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 16 flops/particle, 4 loads, 2 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
! and x(t+dt) = x(t) + v(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real qbm, dt, ek
      real ppart, fx
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), fx(nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtm, anx, edgelx, edgerx, x, dx, vx
      real sfx
      dimension sfx(MXV)
!     dimension sfx(mx+1)
      double precision sum1, sum2
      qtm = qbm*dt
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dx,vx,edgelx,edgerx,sum1,sfx)    &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! load local fields from global array
      do 10 j = 1, nn+1
      sfx(j) = fx(j+noff)
   10 continue
! clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = x - real(nn)
      nn = nn - noff + 1
! find acceleration
      dx = (1.0 - dx)*sfx(nn) + dx*sfx(nn+1)
! new velocity
      vx = ppart(2,j,k)
      dx = vx + qtm*dx
! average kinetic energy
      vx = vx + dx
      sum1 = sum1 + vx*vx
      ppart(2,j,k) = dx
! new position
      dx = x + dx*dt
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   30 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   40 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 50 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   50    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + .125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GRPPUSH1L(ppart,fx,kpic,qbm,dt,ci,ek,idimp,nppmx,nx,mx,&
     &nxv,mx1,ipbc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, for relativistic particles
! with various boundary conditions.
! OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 21 flops/particle, 2 divides, 2 sqrts, 4 loads, 2 stores
! input: all, output: ppart, ek
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2))*ci*ci)
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! kpic = number of particles per tile
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t))*ci*ci)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real qbm, dt, ci, ek
      real ppart, fx
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fx(nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, ci2, edgelx, edgerx, x, dx, acx, p2, dtg
      real sfx
      dimension sfx(MXV)
!     dimension sfx(mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dx,acx,p2,dtg,sum1,sfx)                &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      do 10 j = 1, min(mx,nx-noff)+1
      sfx(j) = fx(j+noff)
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = x - real(nn)
      nn = nn - noff + 1
! find acceleration
      dx = (1.0 - dx)*sfx(nn) + dx*sfx(nn+1)
! calculate half impulse
      dx = qtmh*dx
! half acceleration
      acx = ppart(2,j,k) + dx
! time-centered kinetic energy
      p2 = acx*acx
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      dx = acx + dx
      ppart(2,j,k) = dx
! update inverse gamma
      p2 = dx*dx
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = ppart(1,j,k) + dx*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   20 continue
      sum2 = sum2 + sum1
   30 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GRPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ci,ek,idimp, &
     &nppmx,nx,mx,nxv,mx1,ntmax,irc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, for relativistic particles
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data read in tiles
! particles stored in segmented array
! 21 flops/particle, 2 divides, 2 sqrts, 4 loads, 2 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2))*ci*ci)
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t)*px(t))*ci*ci)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real qbm, dt, ci, ek
      real ppart, fx
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), fx(nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtmh, ci2, anx, edgelx, edgerx, x, dx, acx, p2, dtg
      real sfx
      dimension sfx(MXV)
!     dimension sfx(mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      sum2 = 0.0d0
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dx,edgelx,edgerx,acx,p2,dtg,sum1,&
!$OMP& sfx) REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! load local fields from global array
      do 10 j = 1, nn+1
      sfx(j) = fx(j+noff)
   10 continue
! clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = x - real(nn)
      nn = nn - noff + 1
! find acceleration
      dx = (1.0 - dx)*sfx(nn) + dx*sfx(nn+1)
! calculate half impulse
      dx = qtmh*dx
! half acceleration
      acx = ppart(2,j,k) + dx
! time-centered kinetic energy
      p2 = acx*acx
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
! new momentum
      dx = acx + dx
      ppart(2,j,k) = dx
! update inverse gamma
      p2 = dx*dx
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = ppart(1,j,k) + dx*dtg
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   30 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
   40 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 50 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   50    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPUSH1ZF(ppart,kpic,dt,ek,idimp,nppmx,nx,mx1,ipbc)
! for 1d code, this subroutine updates particle co-ordinate for
! particles with fixed velocities, with various boundary conditions.
! OpenMP version
! 4 flops/particle, 2 loads, 1 stores
! input: all, output: ppart, ek
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t+dt/2)**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nppmx, nx, mx1, ipbc
      real dt, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer npp
      integer j, k
      real edgelx, edgerx, dx
      double precision sum1, sum2
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,dx,sum1) REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 20 k = 1, mx1
      npp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, npp
! velocity
      dx = ppart(2,j,k)
! average kinetic energy
      sum1 = sum1 + dx**2
! new position
      dx = ppart(1,j,k) + dx*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ek,idimp,nppmx,nx,mx,&
     &mx1,ntmax,irc)
! for 1d code, this subroutine updates particle co-ordinate for
! particles with fixed velocities, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version
! 4 flops/particle, 2 loads, 1 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum(vx(t+dt/2)**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, mx1, ntmax, irc
      real dt, ek
      real ppart
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer noff, npp
      integer j, k, nn, ih, nh
      real anx, edgelx, edgerx, dx
      double precision sum1, sum2
      anx = real(nx)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,dx,edgelx,edgerx,sum1)             &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! clear counters
      do 10 j = 1, 2
      ncl(j,k) = 0
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 20 j = 1, npp
! velocity
      dx = ppart(2,j,k)
! average kinetic energy
      sum1 = sum1 + dx**2
! new position
      dx = ppart(1,j,k) + dx*dt
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
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
         do 40 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine RPPUSH1ZF(ppart,kpic,dt,ci,ek,idimp,nppmx,nx,mx1,ipbc)
! for 1d code, this subroutine updates particle co-ordinate for
! relativistic particles with fixed velocities, with various boundary
! conditions. OpenMP version
! 7 flops/particle, 2 divides, 1 sqrts, 2 loads, 1 stores
! input: all, output: ppart, ek
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2))
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum vx of particle n in tile m
! kpic = number of particles per tile
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2))**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t+dt/2)**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nppmx, nx, mx1, ipbc
      real dt, ci, ek
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
! local data
      integer npp
      integer j, k
      real ci2, edgelx, edgerx, dx, p2, gam, dtg
      double precision sum1, sum2
      ci2 = ci*ci
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,dx,p2,gam,dtg,sum1) REDUCTION(+:sum2)            &
!$OMP& SCHEDULE(dynamic)
      do 20 k = 1, mx1
      npp = kpic(k)
      sum1 = 0.0d0
! loop over particles in tile
      do 10 j = 1, npp
! momentum
      dx = ppart(2,j,k)
! average kinetic energy
      p2 = dx*dx
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = ppart(1,j,k) + dx*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine RPPUSHF1ZF(ppart,kpic,ncl,ihole,dt,ci,ek,idimp,nppmx,nx&
     &,mx,mx1,ntmax,irc)
! for 1d code, this subroutine updates particle co-ordinate for
! relativistic particles with fixed velocities, with periodic boundary
! conditions.
! also determines list of particles which are leaving this tile
! OpenMP version
! 7 flops/particle, 2 divides, 1 sqrts, 2 loads, 1 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2))
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momemtum vx of particle n in tile m
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2))**2)/(1. + gamma)
! where gamma = sqrt(1.+(px(t+dt/2)**2)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer idimp, nppmx, nx, mx, mx1, ntmax, irc
      real dt, ci, ek
      real ppart
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer noff, npp
      integer j, k, nn, ih, nh
      real ci2, anx, edgelx, edgerx, dx, p2, gam, dtg
      double precision sum1, sum2
      ci2 = ci*ci
      anx = real(nx)
      sum2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,dx,p2,gam,dtg,edgelx,edgerx,sum1)  &
!$OMP& REDUCTION(+:sum2) SCHEDULE(dynamic)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! clear counters
      do 10 j = 1, 2
      ncl(j,k) = 0
   10 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 20 j = 1, npp
! momentum
      dx = ppart(2,j,k)
! average kinetic energy
      p2 = dx*dx
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
      dtg = dt/gam
! new position
      dx = ppart(1,j,k) + dx*dtg
! set new position
      ppart(1,j,k) = dx
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         nn = 2
      else if (dx.lt.edgelx) then
         nn = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               nn = 0
            endif
         endif
      endif
! increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
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
         do 40 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc = ih
      endif
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
! for 1d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 7 flops/particle, 3 loads, 3 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n)=qm*(1.-dx) and q(n+1)=qm*dx
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in tile m
! q(j) = charge density at grid point j
! kpic = number of particles per tile
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 2
! mx = number of grids in sorting cell in x
! nxv = first dimension of charge array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, mx, nxv, mx1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mx1), q(nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real x, dx
      real sq
!     dimension sq(MXV)
      dimension sq(mx+1)
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,noff,npp,nn,x,dx,sq) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 10 j = 1, mx+1
      sq(j) = 0.0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = qm*(x - real(nn))
      nn = nn - noff + 1
! deposit charge within tile to local accumulator
      sq(nn) = sq(nn) + (qm - dx)
      sq(nn+1) = sq(nn+1) + dx
   20 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      q(j+noff) = q(j+noff) + sq(j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      q(1+noff) = q(1+noff) + sq(1)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noff) = q(nn+noff) + sq(nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
