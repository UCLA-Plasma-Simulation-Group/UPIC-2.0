!-----------------------------------------------------------------------
! Fortran Library for depositing current density
! 2-1/2D Vector/MPI/OpenMP PIC Codes:
! VPPGJPPOST2L calculates particle current density using linear
!              interpolation and advances particle positions half a
!              time-step with various particle boundary conditions
! VPPGJPPOSTF2L calculates particle current density using linear
!               interpolation, advances particle positions half a
!               time-step with periodic boundary conditions, and
!               determines list of particles which are leaving each tile
! VPPGRJPPOST2L calculates particle current density using linear
!               interpolation for relativistic particles and advances
!               particle positions half a time-step with with various
!               particle boundary conditions
! VPPGRJPPOSTF2L calculates particle current density using linear
!                interpolation for relativistic particles, advances
!                particle positions half a time-step with periodic
!                boundary conditions, and determines list of particles
!                which are leaving each tile
! VPPGMJPPOST2L calculates particle momentum flux using linear
!               interpolation
! VPPGRMJPPOST2L calculates relativistic particle momentum flux using
!                linear interpolation
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: august 17, 2018
!-----------------------------------------------------------------------
      subroutine VPPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny&
     &,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 41 flops/particle, 17 loads, 14 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
!     dimension scu(3,MXV*MYV)
      dimension scu(3,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,5)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,mnoff,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,scu,n,s,t)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 60 m = 1, ipp
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = ppart(3,j+joff,k)
      t(j,4) = ppart(4,j+joff,k)
      t(j,5) = ppart(5,j+joff,k)
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,3)
      vy = t(j,4)
      vz = t(j,5)
!dir$ ivdep
      do 30 i = 1, lvect
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + vx*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + vy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + vz*s(j,i)
   30 continue
   40 continue
! advance position half a time-step
      if (dt.eq.0.0) go to 60
! !dir$ vector aligned
      do 50 j = 1, npblk
      dx = t(j,1) + t(j,3)*dt
      dy = t(j,2) + t(j,4)*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(4,j+joff,k) = -t(j,4)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
      endif
! set new position
      ppart(1,j+joff,k) = dx
      ppart(2,j+joff,k) = dy
   50 continue
   60 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 70 j = nps, nppp
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
! deposit current
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      dx = amx*amy
      dy = dxp*amy
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 70
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
   70 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
!dir$ ivdep
      do 80 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
   80 continue
   90 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     & + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     & + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     & + scu(3,i+lxv*(mm-1))
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine V2PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx, &
     &ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 41 flops/particle, 17 loads, 14 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
!     dimension scu(3,MXV*MYV)
      dimension scu(3,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,5)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,mnoff,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,scu,n,s,t)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 60 m = 1, ipp
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = ppart(3,j+joff,k)
      t(j,4) = ppart(4,j+joff,k)
      t(j,5) = ppart(5,j+joff,k)
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      vx = t(j,3)
      vy = t(j,4)
      vz = t(j,5)
!dir$ ivdep
      do 30 i = 1, lvect
      if (i.gt.2) nn = mm
      scu(1,i+nn) = scu(1,i+nn) + vx*s(j,i)
      scu(2,i+nn) = scu(2,i+nn) + vy*s(j,i)
      scu(3,i+nn) = scu(3,i+nn) + vz*s(j,i)
   30 continue
   40 continue
! advance position half a time-step
      if (dt.eq.0.0) go to 60
! !dir$ vector aligned
      do 50 j = 1, npblk
      dx = t(j,1) + t(j,3)*dt
      dy = t(j,2) + t(j,4)*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(4,j+joff,k) = -t(j,4)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
      endif
! set new position
      ppart(1,j+joff,k) = dx
      ppart(2,j+joff,k) = dy
   50 continue
   60 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 70 j = nps, nppp
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
! deposit current
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      dx = amx*amy
      dy = dxp*amy
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 70
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
   70 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
!dir$ ivdep
      do 80 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
   80 continue
   90 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     & + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     & + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     & + scu(3,i+lxv*(mm-1))
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,  &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 41 flops/particle, 17 loads, 14 stores
! input: all except ncl, ihole, irc,
! output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, ih, nh, nn, mm, lxv
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
!     dimension scu(3,MXV*MYV)
      dimension scu(3,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,5)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      anx = real(nx)
      any = real(ny)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,ih,nh,mnoff, &
!$OMP& x,y,vx,vy,vz,dxp,dyp,amx,amy,dx,dy,edgelx,edgely,edgerx,edgery,  &
!$OMP& scu,n,s,t)
      do 140 k = 1, mxyp1
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
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 80 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 30 j = 1, npblk
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = ppart(3,j+joff,k)
      t(j,4) = ppart(4,j+joff,k)
      t(j,5) = ppart(5,j+joff,k)
   30 continue
! deposit current within tile to local accumulator
      do 50 j = 1, npblk
      vx = t(j,3)
      vy = t(j,4)
      vz = t(j,5)
!dir$ ivdep
      do 40 i = 1, lvect
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + vx*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + vy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + vz*s(j,i)
   40 continue
   50 continue
! advance position half a time-step
      if (dt.eq.0.0) go to 80
! !dir$ vector aligned
      do 60 j = 1, npblk
      dx = t(j,1) + t(j,3)*dt
      dy = t(j,2) + t(j,4)*dt
      s(j,1) = dx
      s(j,2) = dy
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
      n(j) = mm
   60 continue
! increment counters
      do 70 j = 1, npblk
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
   70 continue
   80 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 90 j = nps, nppp
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
! deposit current
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      dx = amx*amy
      dy = dxp*amy
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 90
      dx = x + vx*dt
      dy = y + vy*dt
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
   90 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 110 j = 2, mm
!dir$ ivdep
      do 100 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
  100 continue
  110 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 120 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     & + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     & + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     & + scu(3,i+lxv*(mm-1))
      endif
  120 continue
      nn = min(mx+1,nxv-noffp)
      do 130 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  130 continue
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  140 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 150 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  150    continue
         irc = ih
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp, &
     &nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
      real scu
!     dimension scu(3,MXV*MYV)
      dimension scu(3,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,8)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      ci2 = ci*ci
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
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,mnoff,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,ux,uy,uz,p2,gami,scu,n,s,t)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 60 m = 1, ipp
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
! find inverse gamma
      ux = ppart(3,j+joff,k)
      uy = ppart(4,j+joff,k)
      uz = ppart(5,j+joff,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = ux
      t(j,4) = uy
      t(j,5) = uz
      t(j,6) = ux*gami
      t(j,7) = uy*gami
      t(j,8) = uz*gami
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,6)
      vy = t(j,7)
      vz = t(j,8)
!dir$ ivdep
      do 30 i = 1, lvect
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + vx*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + vy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + vz*s(j,i)
   30 continue
   40 continue
! advance position half a time-step
      if (dt.eq.0.0) go to 60
! !dir$ vector aligned
      do 50 j = 1, npblk
      dx = t(j,1) + t(j,6)*dt
      dy = t(j,2) + t(j,7)*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(4,j+joff,k) = -t(j,4)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(3,j+joff,k) = -t(j,3)
         endif
      endif
! set new position
      ppart(1,j+joff,k) = dx
      ppart(2,j+joff,k) = dy
   50 continue
   60 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 70 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      ux = ppart(3,j,k)
      uy = ppart(4,j,k)
      uz = ppart(5,j,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 70
      dx = x + vx*dt
      dy = y + vy*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(4,j,k) = -uy
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(3,j,k) = -ux
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   70 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
!dir$ ivdep
      do 80 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
   80 continue
   90 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     & + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     & + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     & + scu(3,i+lxv*(mm-1))
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt, &
     &ci,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! if dt /= 0, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
! input: all except ncl, ihole, irc,
! output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k) = ith component of current density at grid point (j,kk),
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
!       returns new ntmax required
! optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, ih, nh, lxv
      real ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
!     dimension scu(3,MXV*MYV)
      dimension scu(3,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,8)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,ih,nh,mnoff, &
!$OMP& x,y,vx,vy,vz,ux,uy,uz,dxp,dyp,amx,amy,dx,dy,edgelx,edgely,edgerx,&
!$OMP& edgery,p2,gami,scu,n,s,t)
      do 140 k = 1, mxyp1
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
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 80 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 30 j = 1, npblk
! find interpolation weights
      x = ppart(1,j+joff,k)
      y = ppart(2,j+joff,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      ux = ppart(3,j+joff,k)
      uy = ppart(4,j+joff,k)
      uz = ppart(5,j+joff,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = ux
      t(j,4) = uy
      t(j,5) = uz
      t(j,6) = ux*gami
      t(j,7) = uy*gami
      t(j,8) = uz*gami
   30 continue
! deposit current within tile to local accumulator
      do 50 j = 1, npblk
      vx = t(j,6)
      vy = t(j,7)
      vz = t(j,8)
!dir$ ivdep
      do 40 i = 1, lvect
      scu(1,n(j)+mn(i)) = scu(1,n(j)+mn(i)) + vx*s(j,i)
      scu(2,n(j)+mn(i)) = scu(2,n(j)+mn(i)) + vy*s(j,i)
      scu(3,n(j)+mn(i)) = scu(3,n(j)+mn(i)) + vz*s(j,i)
   40 continue
   50 continue
! advance position half a time-step
      if (dt.eq.0.0) go to 80
! !dir$ vector aligned
      do 60 j = 1, npblk
      dx = t(j,1) + t(j,6)*dt
      dy = t(j,2) + t(j,7)*dt
      s(j,1) = dx
      s(j,2) = dy
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
      n(j) = mm
   60 continue
! increment counters
      do 70 j = 1, npblk
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
   70 continue
   80 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 90 j = nps, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
! find inverse gamma
      ux = ppart(3,j,k)
      uy = ppart(4,j,k)
      uz = ppart(5,j,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 90
      dx = x + vx*dt
      dy = y + vy*dt
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
   90 continue
! deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 110 j = 2, mm
!dir$ ivdep
      do 100 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) +   &
     & scu(3,i+lxv*(j-1))
  100 continue
  110 continue
! deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 120 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))&
     & + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))&
     & + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))&
     & + scu(3,i+lxv*(mm-1))
      endif
  120 continue
      nn = min(mx+1,nxv-noffp)
      do 130 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) +   &
     &scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1))&
     & + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1))&
     & + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1))&
     & + scu(3,nn+lxv*(j-1))
      endif
  130 continue
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,k) = ih
  140 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 150 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
  150    continue
         irc = ih
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,&
     &nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 51 flops/particle, 21 loads, 16 stores
! input: all, output: amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where vj = vj(t-dt/2) and vk = vk(t-dt/2)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = x velocity of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = y velocity of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = z velocity of particle n in partition in tile m
! at t - dt/2
! amu(i,j,k) = ith component of momentum flux at grid point j,kk
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of flux array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, v1, v2, v3, v4
      real samu
!     dimension samu(4,MXV*MYV)
      dimension samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,mnoff,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,v1,v2,v3,v4,samu,n,s,t)
      do 110 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
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
      t(j,1) = ppart(3,j+joff,k)
      t(j,2) = ppart(4,j+joff,k)
      t(j,3) = ppart(5,j+joff,k)
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
!dir$ ivdep
      do 30 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
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
! deposit momentum flux
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
   60 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
   70 continue
   80 continue
! deposit momentum flux to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 90 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
      endif
   90 continue
      nn = min(mx+1,nxv-noffp)
      do 100 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
      endif
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx&
     &,my,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates particle momentum flux
! using first-order spline interpolation for relativistic particles
! vectorizable/OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored in segmented array
! 62 flops/particle, 1 divide, 21 loads, 16 stores
! input: all, output: amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where pj = pj(t-dt/2) and pk = pk(t-dt/2)
! where gami2 = 1./(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = x momentum of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = y momentum of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = z momentum of particle n in partition in tile m
! at t - dt/2
! amu(i,j,k) = ith component of momentum flux at grid point j,kk
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of flux array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, mx, my, nxv, nypmx, mx1, mxyp1
      real qm, ci
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), amu(4,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noffp, moffp, nppp, ipp, joff, nps
      integer mnoff, i, j, k, m, nn, mm, lxv
      real ci2, gami, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, v1, v2, v3, v4
      real samu
!     dimension samu(4,MXV*MYV)
      dimension samu(4,(mx+1)*(my+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
      lxv = mx + 1
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,mnoff,x,y,vx,&
!$OMP& vy,vz,dxp,dyp,amx,amy,dx,dy,v1,v2,v3,v4,p2,gami,samu,n,s,t)
      do 110 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1,(mx+1)*(my+1)
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
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
! find inverse gamma
      vx = ppart(3,j+joff,k)
      vy = ppart(4,j+joff,k)
      vz = ppart(5,j+joff,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      n(j) = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
      t(j,1) = vx*gami
      t(j,2) = vy*gami
      t(j,3) = vz*gami
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
!dir$ ivdep
      do 30 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
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
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit momentum flux
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
   60 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + &
     & samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + &
     & samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + &
     & samu(3,i+lxv*(j-1))
      amu(4,i+noffp+nxv*(j+moffp-1)) = amu(4,i+noffp+nxv*(j+moffp-1)) + &
     & samu(4,i+lxv*(j-1))
   70 continue
   80 continue
! deposit momentum flux to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 90 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp) = amu(4,i+noffp+nxv*moffp) + samu(4,i)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)) + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)) + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)) + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)) =                              &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)) + samu(4,i+lxv*(mm-1))
      endif
   90 continue
      nn = min(mx+1,nxv-noffp)
      do 100 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + &
     &samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + &
     &samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + &
     &samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)) = amu(4,1+noffp+nxv*(j+moffp-1)) + &
     &samu(4,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(1,nn+noffp+nxv*(j+moffp-1)) + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(2,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(2,nn+noffp+nxv*(j+moffp-1)) + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(3,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(3,nn+noffp+nxv*(j+moffp-1)) + samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
         amu(4,nn+noffp+nxv*(j+moffp-1)) =                              &
     &   amu(4,nn+noffp+nxv*(j+moffp-1)) + samu(4,nn+lxv*(j-1))
      endif
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
