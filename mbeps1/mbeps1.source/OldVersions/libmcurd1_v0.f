!-----------------------------------------------------------------------
! Fortran Library for depositing current density
! 1-2/2D OpenMP PIC Codes:
! GJPPOST1L calculates particle current density using linear
!           interpolation and advances particle positions half a
!           time-step with various particle boundary conditions
! GJPPOSTF1L calculates particle current density using linear
!            interpolation, advances particle positions half a
!            time-step with periodic boundary conditions, and
!            determines list of particles which are leaving each tile
! GRJPPOST1L calculates particle current density using linear
!            interpolation for relativistic particles and advances
!            particle positions half a time-step with with various
!            particle boundary conditions
! GRJPPOSTF1L calculates particle current density using linear
!             interpolation for relativistic particles, advances
!             particle positions half a time-step with periodic
!             boundary conditions, and determines list of particles
!             which are leaving each tile
! GMJPPOST1L calculates particle momentum flux using linear interpolation
!            interpolation
! GRMJPPOST1L calculates relativistic particle momentum flux using
!             linear interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 28, 2016
!-----------------------------------------------------------------------
      subroutine GJPPOST1L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,mx,nxv,mx1&
     &,ipbc)
! for 1-2/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 14 flops/particle, 8 loads, 5 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vi, where i = y,z
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! cu(i,j) = ith component of current density at grid point j
! kpic = number of particles per tile
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of current array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real edgelx, edgerx, dxp, amx, vy, vz, x, dx
      real scu
      dimension scu(2,MXV)
!     dimension scu(2,mx+1)
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
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,vy,vz,scu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
! deposit current within tile to local accumulator
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
! advance position half a time-step
      if (dt.eq.0.0) go to 20
      dx = x + ppart(2,j,k)*dt
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
! deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp,nx&
     &,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 14 flops/particle, 8 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, cu, irc
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vi, where i = y,z
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! cu(i,j) = ith component of current density at grid point j
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of current array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ntmax, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real anx, edgelx, edgerx, dxp, amx, vy, vz, x, dx
      real scu
      dimension scu(2,MXV)
!     dimension scu(2,mx+1)
      anx = real(nx)
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,vy,vz,edgelx,edgerx,
!$OMP& scu)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
! deposit current within tile to local accumulator
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
! advance position half a time-step
      if (dt.eq.0.0) go to 30
      dx = x + ppart(2,j,k)*dt
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
! set new position
      ppart(1,j,k) = dx
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
! deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 40 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   40 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRJPPOST1L(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,mx,nxv&
     &,mx1,ipbc)
! for 1-2/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 24 flops/particle, 1 divide, 1 sqrt, 9 loads, 5 stores
! input: all, output: part, cu
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pi*gami, where i = y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = x momentum of particle n in tile m
! ppart(3,n,m) = y momentum of particle n in tile m
! ppart(4,n,m) = z momentum of particle n in tile m
! cu(i,j) = ith component of current density at grid point j
! kpic = number of particles per tile
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of current array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real edgelx, ci2, edgerx, dxp, amx, vx, vy, vz, x, dx, p2, gami
      real scu
      dimension scu(2,MXV)
!     dimension scu(2,mx+1)
      ci2 = ci*ci
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,vx,vy,vz,p2,gami,scu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
! find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! deposit current within tile to local accumulator
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
! advance position half a time-step
      if (dt.eq.0.0) go to 20
      dx = x + vx*dt
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
! deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,    &
     &idimp,nx,mx,nxv,mx1,ntmax,irc)
! for 1-2/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 24 flops/particle, 1 divide, 1 sqrt, 9 loads, 5 stores
! input: all except ncl, ihole, irc, output: part, ncl, ihole, cu, irc
! current density is approximated by values at the nearest grid points
! cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pi*gami, where i = y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = x momentum of particle n in tile m
! ppart(3,n,m) = y momentum of particle n in tile m
! ppart(4,n,m) = z momentum of particle n in tile m
! cu(i,j) = ith component of current density at grid point j
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nxv = second dimension of current array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ntmax, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real edgelx, ci2, edgerx, dxp, amx, vx, vy, vz, x, dx, p2, gami
      real anx
      real scu
      dimension scu(2,MXV)
!c     dimension scu(2,mx+1)
      ci2 = ci*ci
      anx = real(nx)
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,vx,vy,vz,p2,gami,     
!$OMP& edgelx,edgerx,scu)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
! zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
! loop over particles in tile
      do 30 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
! find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! deposit current within tile to local accumulator
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
! advance position half a time-step
      if (dt.eq.0.0) go to 30
      dx = x + vx*dt
! find particles going out of bounds
      nn = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
! set new position
      ppart(1,j,k) = dx
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
! deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 40 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   40 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GMJPPOST1L(ppart,amu,kpic,qm,nppmx,idimp,mx,nxv,mx1)
! for 1-2/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 14 flops/particle, 8 loads, 4 stores
! input: all, output: part, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
! where vj = vj(t-dt/2) and vk = vk(t-dt/2)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! amu(i,j) = ith component of momentum flux at grid point j
! kpic = number of particles per tile
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! mx = number of grids in sorting cell in x
! nxv = second dimension of flux array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, mx, nxv, mx1
      real qm
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real dxp, amx, x, vx, v1, v2
      real samu
      dimension samu(2,MXV)
!     dimension samu(2,mx+1)
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,vx,v1,v2,samu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 10 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
! deposit momentum flux within tile to local accumulator
      vx = ppart(2,j,k)
      v1 = vx*ppart(3,j,k)
      v2 = vx*ppart(4,j,k)
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
   20 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRMJPPOST1L(ppart,amu,kpic,qm,ci,nppmx,idimp,mx,nxv,mx1&
     &)
! for 1-2/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation for relativistic particles
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 23 flops/particle, 1 divide, 8 loads, 4 stores
! input: all, output: part, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
! where n = nearest grid point and dx = x-n
! and qci = qm*pj*pk*gami2, where jk = xy,xz for i = 1, 2
! where pj = pj(t-dt/2) and pk = pk(t-dt/2)
! where gami2 = 1./(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! ppart(3,n,m) = momentum py of particle n in tile m
! ppart(4,n,m) = momentum pz of particle n in tile m
! amu(i,j) = ith component of momentum flux at grid point j
! kpic = number of particles per tile
! qm = charge on particle, in units of e
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! mx = number of grids in sorting cell in x
! nxv = second dimension of flux array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, mx, nxv, mx1
      real qm, ci
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), amu(2,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real ci2, gami2, dxp, amx, x, vx, vy, vz, p2, v1, v2
      real samu
      dimension samu(2,MXV)
!     dimension samu(2,mx+1)
      ci2 = ci*ci
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,vx,vy,vz,p2,v1,v2,samu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 10 j = 1, mx+1
      samu(1,j) = 0.0
      samu(2,j) = 0.0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
c find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
! deposit momentum flux within tile to local accumulator
      v1 = vx*vy*gami2
      v2 = vx*vz*gami2
      samu(1,nn) = samu(1,nn) + v1*amx
      samu(2,nn) = samu(2,nn) + v2*amx
      samu(1,nn+1) = samu(1,nn+1) + v1*dxp
      samu(2,nn+1) = samu(2,nn+1) + v2*dxp
   20 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      amu(1,j+noff) = amu(1,j+noff) + samu(1,j)
      amu(2,j+noff) = amu(2,j+noff) + samu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      amu(1,1+noff) = amu(1,1+noff) + samu(1,1)
!$OMP ATOMIC
      amu(2,1+noff) = amu(2,1+noff) + samu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff) = amu(1,nn+noff) + samu(1,nn)
!$OMP ATOMIC
         amu(2,nn+noff) = amu(2,nn+noff) + samu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
