!-----------------------------------------------------------------------
! Fortran Library for pushing electrostatic particles, depositing charge
! and copying particles
! 3D MPI/OpenMP PIC Codes:
! PPPMOVIN3L sorts particles by x,y,z grid in tiles of mx, my, mz and
!            copies to segmented array ppart
! PPPCOPYOUT3 copies segmented particle data ppart to the array part
! PPPCHECK3L performs a sanity check to make sure particles sorted
!            by x,y grid in tiles of mx, my, are all within bounds.
! PPGPPUSH32L updates particle co-ordinates and velocities using
!             electric field only, with linear interpolation and various
!             particle boundary conditions
! PPGPPUSHF32L updates particle co-ordinates and velocities using
!              electric field only, with linear interpolation and
!              periodic particle boundary conditions.  also determines
!              list of particles which are leaving each tile
! PPGRPPUSH32L updates relativistic particle co-ordinates and momenta
!              using electric field only, with linear interpolation and
!              various particle boundary conditions
! PPGRPPUSHF32L updates relativistic particle co-ordinates and
!               velocities using electric field only, with linear
!               interpolation and periodic particle boundary conditions.
!               also determines list of particles which are leaving each
!               tile
! PPGPPOST32L calculates particle charge density using linear
!             interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 12, 2016
!-----------------------------------------------------------------------
      subroutine PPPMOVIN3L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, &
     &mx,my,mz,mx1,myp1,mxyzp1,idds,irc)
! this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
! and copies to segmented array ppart
! linear interpolation, spatial decomposition in y/z direction
! input: all except ppart, kpic, output: ppart, kpic
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = position z of particle n in tile m
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = velocity vy of particle n in tile m
! ppart(6,n,m) = velocity vz of particle n in tile m
! kpic = output number of particles per tile
! npp = number of particles in partition
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! mx/my/mz = number of grids in sorting cell in x, y and z
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, nppmx, idimp, npmax, mx, my, mz, mx1, myp1, mxyzp1
      integer idds, irc
      integer kpic, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer i, j, k, n, m, l, mnoff, lnoff, mxyp1, ip, ierr
      mnoff = noff(1)
      lnoff = noff(2)
      ierr = 0
      mxyp1 = mx1*myp1
! clear counter array
      do 10 k = 1, mxyzp1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      l = part(3,j)
      l = (l - lnoff)/mz
      m = n + mx1*m + mxyp1*l
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
      subroutine PPPCOPYOUT3(part,ppart,kpic,npp,npmax,nppmx,idimp,     &
     &mxyzp1,irc)
! for 3d code, this subroutine copies segmented particle data ppart to
! the array part with original tiled layout
! spatial decomposition in y/z direction
! input: all except part, npp, output: part, npp
! part(i,j) = i-th coordinate for particle j in partition
! ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
! kpic = number of particles per tile
! npp = number of particles in partition
! npmax = maximum number of particles in each partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mxyzp1 = total number of tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyzp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1)
! local data
      integer i, j, k, npoff, nppp, ne, ierr
      npoff = 0
      ierr = 0
! loop over tiles
      do 30 k = 1, mxyzp1
      nppp = kpic(k)
      ne = nppp + npoff
      if (ne.gt.npmax) ierr = max(ierr,ne-npmax)
      if (ierr.gt.0) nppp = 0
! loop over particles in tile
      do 20 j = 1, nppp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + nppp
   30 continue
      npp = npoff
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCHECK3L(ppart,kpic,noff,nyzp,idimp,nppmx,nx,mx,my,mz&
     &,mx1,myp1,mzp1,idds,irc)
! this subroutine performs a sanity check to make sure particles sorted
! by x,y,z grid in tiles of mx, my, mz, are all within bounds.
! tiles are assumed to be arranged in 3D linear memory
! input: all except irc
! output: irc
! ppart(1,n,l) = position x of particle n in tile l
! ppart(2,n,l) = position y of particle n in tile l
! ppart(3,n,l) = position a of particle n in tile l
! kpic(l) = number of reordered output particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx = system length in x
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
! irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, mx1, myp1, mzp1, idds, irc
      real ppart
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension kpic(mx1*myp1*mzp1), noff(idds), nyzp(idds)
! local data
      integer mxyp1, mxyzp1, noffp, moffp, loffp, nppp
      integer j, k, l, nn, mm, ll, ist
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz, dx, dy, dz
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,ist,edgelx,edgely, 
!$OMP& edgelz,edgerx,edgery,edgerz,dx,dy,dz)
      do 20 l = 1, mxyzp1
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
! loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(1,j,l)
      dy = ppart(2,j,l)
      dz = ppart(3,j,l)
! find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (dz.lt.edgelz) ist = ist + 9
      if (dz.ge.edgerz) ist = ist + 18
      if (ist.gt.0) irc = l
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,idimp, &
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with various boundary conditions.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data read in tiles
! particles stored segmented array
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
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
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
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,x,y,z
!$OMP& ,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 50 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, nppp
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
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! new velocity
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
      dz = vz + qtm*dz
! average kinetic energy
      vx = vx + dx
      vy = vy + dy
      vz = vz + dz
      sum1 = sum1 + (vx*vx + vy*vy+ vz*vz)
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
! new position
      dx = x + dx*dt
      dy = y + dy*dt
      dz = z + dz*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -ppart(6,j,l)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   40 continue
      sum2 = sum2 + sum1
   50 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt&
     &,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! for distributed data, with 2D spatial decomposition
! data read in tiles
! particles stored segmented array
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
! optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, ih, nh, nn, mm, ll
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,ih,nh,nn,mm,ll
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,     
!$OMP& edgely,edgelz,edgerx,edgery,edgerz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 60 l = 1, mxyzp1
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
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 50 j = 1, nppp
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
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! new velocity
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
      dz = vz + qtm*dz
! average kinetic energy
      vx = vx + dx
      vy = vy + dy
      vz = vz + dz
      sum1 = sum1 + (vx*vx + vy*vy+ vz*vz)
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
! new position
      dx = x + dx*dt
      dy = y + dy*dt
      dz = z + dz*dt
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
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
   50 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   60 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ci,ek,   &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds&
     &,ipbc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with various boundary conditions.
! with 2D spatial decomposition
! OpenMP version using guard cells,
! for distributed data with 2D spatial decomposition
! data read in tiles
! particles stored segmented array
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
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll
      real qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz, acx, acy, acz, p2, dtg
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
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
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,x,y,z
!$OMP& ,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,acx,acy,acz,p2,dtg
!$OMP& ,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 50 l = 1, mxyzp1
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
      sum1 = 0.0d0
! loop over particles in tile
      do 40 j = 1, nppp
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
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
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
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtg
      dy = y + dy*dtg
      dz = z + dz*dtg
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -ppart(6,j,l)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   40 continue
      sum2 = sum2 + sum1
   50 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm, &
     &dt,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,  &
     &mxyzp1,ntmax,idds,irc)
! for 3d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells,
! for distributed data with 2D spatial decomposition
! data read in tiles
! particles stored segmented array
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
! optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qbm, dt, ci, ek
      real ppart, fxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, ih, nh, nn, mm, ll
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real qtmh, ci2, vx, vy, vz, acx, acy, acz, p2, dtg
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
!     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,ih, nh, nn,mm,
!$OMP& ll,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,acx,acy,  
!$OMP& acz,p2,dtg,edgelx,edgely,edgelz,edgerx,edgery,edgerz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 60 l = 1, mxyzp1
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
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
! clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 50 j = 1, nppp
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
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
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
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtg
      dy = y + dy*dtg
      dz = z + dz*dtg
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
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
   50 continue
      sum2 = sum2 + sum1
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   60 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my,mz, &
     &nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
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
      dimension ppart(idimp,nppmx,mxyzp1), q(nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1
      real sq
!     dimension sq(MXV,MYV,MZV)
      dimension sq(mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,sq)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 40 j = 1, nppp
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
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit charge within tile to local accumulator
      x = sq(nn,mm,ll) + amx*amz
      y = sq(nn+1,mm,ll) + amy*amz
      sq(nn,mm,ll) = x
      sq(nn+1,mm,ll) = y
      x = sq(nn,mm+1,ll) + dyp*amz
      y = sq(nn+1,mm+1,ll) + dx1*amz
      sq(nn,mm+1,ll) = x
      sq(nn+1,mm+1,ll) = y
      x = sq(nn,mm,ll+1) + amx*dzp
      y = sq(nn+1,mm,ll+1) + amy*dzp
      sq(nn,mm,ll+1) = x
      sq(nn+1,mm,ll+1) = y
      x = sq(nn,mm+1,ll+1) + dyp*dzp
      y = sq(nn+1,mm+1,ll+1) + dx1*dzp
      sq(nn,mm+1,ll+1) = x
      sq(nn+1,mm+1,ll+1) = y
   40 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      q(i+noffp,j+moffp,k+loffp) = q(i+noffp,j+moffp,k+loffp)           &
     & + sq(i,j,k)
   50 continue
   60 continue
   70 continue
! deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,j+moffp,1+loffp) = q(i+noffp,j+moffp,1+loffp)           &
     & + sq(i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         q(i+noffp,j+moffp,lm+loffp) = q(i+noffp,j+moffp,lm+loffp)      &
     &   + sq(i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp,k+loffp) = q(i+noffp,1+moffp,k+loffp)           &
     & + sq(i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp,k+loffp) = q(i+noffp,mm+moffp,k+loffp)      &
     &   + sq(i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp,k+loffp) = q(1+noffp,j+moffp,k+loffp)           &
     & + sq(1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         q(nm+noffp,j+moffp,k+loffp) = q(nm+noffp,j+moffp,k+loffp)      &
     &   + sq(nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         q(i+noffp,1+moffp,lm+loffp) = q(i+noffp,1+moffp,lm+loffp)      &
     &   + sq(i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            q(i+noffp,mm+moffp,lm+loffp) = q(i+noffp,mm+moffp,lm+loffp) &
     &      + sq(i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         q(1+noffp,j+moffp,lm+loffp) = q(1+noffp,j+moffp,lm+loffp)      &
     &   + sq(1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            q(nm+noffp,j+moffp,lm+loffp) = q(nm+noffp,j+moffp,lm+loffp) &
     &      + sq(nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
