!-----------------------------------------------------------------------
! Fortran Library for depositing current density
! 3D Vector/MPI/OpenMP PIC Codes:
! VPPGJPPOST32L calculates particle current density using linear
!               interpolation and advances particle positions half a
!               time-step with various particle boundary conditions
! VPPGJPPOSTF32L calculates particle current density using linear
!                interpolation, advances particle positions half a
!                time-step with periodic boundary conditions, and
!                determines list of particles which are leaving each
!                tile
! VPPGRJPPOST32L calculates particle current density using linear
!                interpolation for relativistic particles and advances
!                particle positions half a time-step with with various
!                particle boundary conditions
! VPPGRJPPOSTF32L calculates particle current density using linear
!                 interpolation for relativistic particles, advances
!                 particle positions half a time-step with periodic
!                 boundary conditions, and determines list of particles
!                 which are leaving each tile
! VPPGMJPPOST32L calculates particle momentum flux using linear
!                interpolation
! VPPGRMJPPOST32L calculates relativistic particle momentum flux using
!                 linear interpolation
! SET_PVZERO3 zeros out current density array.
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: may 2, 2018
!-----------------------------------------------------------------------
      subroutine VPPGJPPOST32L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx, &
     &ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation
! in addition, particle positions are advanced a half time-step
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 69 flops/particle, 30 loads, 27 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qm, dt
      real ppart, cu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, x, y, z
      real vx, vy, vz
      real scu
!     dimension scu(3,MXV*MYV*MZ)
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,6)
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
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgelz = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      edgerz = real(nz-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,nm,lm,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx, &
!$OMP& dy,dz,scu,n,s,t)
      do 180 l = 1, mxyzp1
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
      t(j,4) = ppart(4,j+joff,l)
      t(j,5) = ppart(5,j+joff,l)
      t(j,6) = ppart(6,j+joff,l)
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
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
      dx = t(j,1) + t(j,4)*dt
      dy = t(j,2) + t(j,5)*dt
      dz = t(j,3) + t(j,6)*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(4,j+joff,l) = -t(j,4)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(5,j+joff,l) = -t(j,5)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            ppart(6,j+joff,l) = -t(j,6)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(4,j+joff,l) = -t(j,4)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(5,j+joff,l) = -t(j,5)
         endif
      endif
! set new position
      ppart(1,j+joff,l) = dx
      ppart(2,j+joff,l) = dy
      ppart(3,j+joff,l) = dz
   50 continue
   60 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 70 j = nps, nppp
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
! deposit current
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = amx*amz
      dy = amy*amz
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = dyp*amz
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dx1*amz
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      scu(1,mm) = scu(1,mm) + vx*dx
      scu(2,mm) = scu(2,mm) + vy*dx
      scu(3,mm) = scu(3,mm) + vz*dx
      dx = dyp*dzp
      scu(1,mm+1) = scu(1,mm+1) + vx*dy
      scu(2,mm+1) = scu(2,mm+1) + vy*dy
      scu(3,mm+1) = scu(3,mm+1) + vz*dy
      dy = dx1*dzp
      scu(1,mm+lxv) = scu(1,mm+lxv) + vx*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + vy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + vz*dx
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + vx*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + vy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 70
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
   70 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 100 k = 2, ll
      do 90 j = 2, mm
!dir$ ivdep
      do 80 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
   80 continue
   90 continue
  100 continue
! deposit current to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 120 j = 2, mm
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  110 continue
  120 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 150 k = 1, ll
      do 130 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  130 continue
      do 140 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  140 continue
  150 continue
      if (lm > mz) then
         do 160 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  160    continue
         do 170 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  170    continue
      endif
  180 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt,&
     &nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,    &
     &ntmax,idds,irc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation.
! in addition, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 69 flops/particle, 30 loads, 27 stores
! input: all except ncl, ihole, irc, output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
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
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, ih, nh, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, x, y, z
      real vx, vy, vz
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real scu
!     dimension scu(3,MXV*MYV*MZ)
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,6)
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
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,ih,nh,nm,lm,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,  &
!$OMP& dx1,dx,dy,dz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,scu,n,s,t)
      do 200 l = 1, mxyzp1
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
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)*(mz+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 26
      ncl(j,l) = 0
   20 continue
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 80 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 30 j = 1, npblk
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
      t(j,4) = ppart(4,j+joff,l)
      t(j,5) = ppart(5,j+joff,l)
      t(j,6) = ppart(6,j+joff,l)
   30 continue
! deposit current within tile to local accumulator
      do 50 j = 1, npblk
      vx = t(j,4)
      vy = t(j,5)
      vz = t(j,6)
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
      dx = t(j,1) + t(j,4)*dt
      dy = t(j,2) + t(j,5)*dt
      dz = t(j,3) + t(j,6)*dt
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
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
      n(j) = mm
   60 continue
! increment counters
      do 70 j = 1, npblk
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
! deposit current
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = amx*amz
      dy = amy*amz
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = dyp*amz
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dx1*amz
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      scu(1,mm) = scu(1,mm) + vx*dx
      scu(2,mm) = scu(2,mm) + vy*dx
      scu(3,mm) = scu(3,mm) + vz*dx
      dx = dyp*dzp
      scu(1,mm+1) = scu(1,mm+1) + vx*dy
      scu(2,mm+1) = scu(2,mm+1) + vy*dy
      scu(3,mm+1) = scu(3,mm+1) + vz*dy
      dy = dx1*dzp
      scu(1,mm+lxv) = scu(1,mm+lxv) + vx*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + vy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + vz*dx
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + vx*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + vy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 90
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
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
   90 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 120 k = 2, ll
      do 110 j = 2, mm
!dir$ ivdep
      do 100 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
  100 continue
  110 continue
  120 continue
! deposit current to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 140 j = 2, mm
      do 130 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  130 continue
  140 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 170 k = 1, ll
      do 150 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  150 continue
      do 160 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  160 continue
  170 continue
      if (lm > mz) then
         do 180 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  180    continue
         do 190 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  190    continue
      endif
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  200 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 210 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  210    continue
         irc = ih
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRJPPOST32L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,&
     &nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! in addition, particle positions are advanced a half time-step
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
! input: all, output: ppart, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = x momentum of particle n in partition in tile m
! ppart(5,n,m) = y momentum of particle n in partition in tile m
! ppart(6,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, x, y, z
      real vx, vy, vz, ux, uy, uz, p2, gami
      real scu
!     dimension scu(3,MXV*MYV*MZ)
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
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
      ci2 = ci*ci
! set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgelz = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      edgerz = real(nz-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,nm,lm,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx, &
!$OMP& vy,vz,ux,uy,uz,p2,gami,scu,n,s,t)
      do 180 l = 1, mxyzp1
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
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      ux = ppart(4,j+joff,l)
      uy = ppart(5,j+joff,l)
      uz = ppart(6,j+joff,l)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
      t(j,4) = ux
      t(j,5) = uy
      t(j,6) = uz
      t(j,7) = ux*gami
      t(j,8) = uy*gami
      t(j,9) = uz*gami
   20 continue
! deposit current within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,7)
      vy = t(j,8)
      vz = t(j,9)
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
      dx = t(j,1) + t(j,7)*dt
      dy = t(j,2) + t(j,8)*dt
      dz = t(j,3) + t(j,9)*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(4,j+joff,l) = -t(j,4)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(5,j+joff,l) = -t(j,5)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            ppart(6,j+joff,l) = -t(j,6)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            ppart(4,j+joff,l) = -t(j,4)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            ppart(5,j+joff,l) = -t(j,5)
         endif
      endif
! set new position
      ppart(1,j+joff,l) = dx
      ppart(2,j+joff,l) = dy
      ppart(3,j+joff,l) = dz
   50 continue
   60 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 70 j = nps, nppp
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
! find inverse gamma
      ux = ppart(4,j,l)
      uy = ppart(5,j,l)
      uz = ppart(6,j,l)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = dyp*amz
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dx1*amz
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      scu(1,mm) = scu(1,mm) + vx*dx
      scu(2,mm) = scu(2,mm) + vy*dx
      scu(3,mm) = scu(3,mm) + vz*dx
      dx = dyp*dzp
      scu(1,mm+1) = scu(1,mm+1) + vx*dy
      scu(2,mm+1) = scu(2,mm+1) + vy*dy
      scu(3,mm+1) = scu(3,mm+1) + vz*dy
      dy = dx1*dzp
      scu(1,mm+lxv) = scu(1,mm+lxv) + vx*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + vy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + vz*dx
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + vx*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + vy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 70
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -uy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -uz
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -uy
         endif
      endif
! set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   70 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 100 k = 2, ll
      do 90 j = 2, mm
!dir$ ivdep
      do 80 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
   80 continue
   90 continue
  100 continue
! deposit charge to interior points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 120 j = 2, mm
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  110 continue
  120 continue
! deposit current to edge points in global array
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 150 k = 1, ll
      do 130 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  130 continue
      do 140 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  140 continue
  150 continue
      if (lm > mz) then
         do 160 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  160    continue
         do 170 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  170    continue
      endif
  180 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRJPPOSTF32L(ppart,cu,kpic,ncl,ihole,noff,nyzp,qm,dt&
     &,ci,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &ntmax,idds,irc)
! for 3d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! in addition, particle positions are advanced a half time-step
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
! input: all except ncl, ihole, irc, output: ppart, cu, ncl, ihole, irc
! current density is approximated by values at the nearest grid points
! cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = x momentum of particle n in partition in tile m
! ppart(5,n,m) = y momentum of particle n in partition in tile m
! ppart(6,n,m) = z momentum of particle n in partition in tile m
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
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
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      integer nppmx, idimp, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), cu(3,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, ih, nh, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, x, y, z
      real vx, vy, vz, ux, uy, uz, ci2, p2, gami
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real scu
!     dimension scu(3,MXV*MYV*MZ)
      dimension scu(3,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,9)
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
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,ih,nh,nm,lm,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy, &
!$OMP& dz,vx,vy,vz,ux,uy,uz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,p2&
!$OMP& ,gami,scu,n,s,t)
      do 200 l = 1, mxyzp1
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
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)*(mz+1)
      scu(1,j) = 0.0
      scu(2,j) = 0.0
      scu(3,j) = 0.0
   10 continue
! clear counters
      do 20 j = 1, 26
      ncl(j,l) = 0
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
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      ux = ppart(4,j+joff,l)
      uy = ppart(5,j+joff,l)
      uz = ppart(6,j+joff,l)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
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
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
      t(j,4) = ux
      t(j,5) = uy
      t(j,6) = uz
      t(j,7) = ux*gami
      t(j,8) = uy*gami
      t(j,9) = uz*gami
   30 continue
! deposit current within tile to local accumulator
      do 50 j = 1, npblk
      vx = t(j,7)
      vy = t(j,8)
      vz = t(j,9)
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
      dx = t(j,1) + t(j,7)*dt
      dy = t(j,2) + t(j,8)*dt
      dz = t(j,3) + t(j,9)*dt
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
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
      n(j) = mm
   60 continue
! increment counters
      do 70 j = 1, npblk
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
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      ux = ppart(4,j,l)
      uy = ppart(5,j,l)
      uz = ppart(6,j,l)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = dyp*amz
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dx1*amz
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      dx = amx*dzp
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
      mm = nn + lxyv
      dy = amy*dzp
      scu(1,mm) = scu(1,mm) + vx*dx
      scu(2,mm) = scu(2,mm) + vy*dx
      scu(3,mm) = scu(3,mm) + vz*dx
      dx = dyp*dzp
      scu(1,mm+1) = scu(1,mm+1) + vx*dy
      scu(2,mm+1) = scu(2,mm+1) + vy*dy
      scu(3,mm+1) = scu(3,mm+1) + vz*dy
      dy = dx1*dzp
      scu(1,mm+lxv) = scu(1,mm+lxv) + vx*dx
      scu(2,mm+lxv) = scu(2,mm+lxv) + vy*dx
      scu(3,mm+lxv) = scu(3,mm+lxv) + vz*dx
      scu(1,mm+1+lxv) = scu(1,mm+1+lxv) + vx*dy
      scu(2,mm+1+lxv) = scu(2,mm+1+lxv) + vy*dy
      scu(3,mm+1+lxv) = scu(3,mm+1+lxv) + vz*dy
! advance position half a time-step
      if (dt.eq.0.0) go to 90
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
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
   90 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 120 k = 2, ll
      do 110 j = 2, mm
!dir$ ivdep
      do 100 i = 2, nn
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(1,i+lxv*(j-1)+lxyv*(k-1))
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(2,i+lxv*(j-1)+lxyv*(k-1))
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                  &
     & scu(3,i+lxv*(j-1)+lxyv*(k-1))
  100 continue
  110 continue
  120 continue
! deposit current to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 140 j = 2, mm
      do 130 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(1,i+lxv*(j-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(2,i+lxv*(j-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                        &
     &cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + scu(3,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  130 continue
  140 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 170 k = 1, ll
      do 150 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                        &
     &cu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + scu(3,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  150 continue
      do 160 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                  &
     &cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                    &
     &+ scu(3,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =              &
     &   cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                &
     &   + scu(3,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  160 continue
  170 continue
      if (lm > mz) then
         do 180 i = 2, nn
!$OMP ATOMIC
         cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                    &
     &   cu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                      &
     &   + scu(3,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  180    continue
         do 190 j = 1, mm
!$OMP ATOMIC
         cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =              &
     &   cu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))                &
     &   + scu(3,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =          &
     &      cu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))            &
     &      + scu(3,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  190    continue
      endif
! set error and end of file flag
      if (nh.gt.0) irc = 1
      ihole(1,1,l) = ih
  200 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) then
         ih = 0
         do 210 l = 1, mxyzp1
         ih = max(ih,ihole(1,1,l))
  210    continue
         irc = ih
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGMJPPOST32L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my&
     &,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle momentum flux
! using first-order linear interpolation
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 121 flops/particle, 52 loads, 48 stores
! input: all, output: ppart, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where vj = vj(t-dt/2) and vk = vk(t-dt/2)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      real ppart, amu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, x, y, z
      real vx, vy, vz, v1, v2, v3, v4, v5, v6
      real samu
!     dimension samu(6,MXV*MYV*MZV)
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
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
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,nm,lm,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx, &
!$OMP& dy,v1,v2,v3,v4,v5,v6,samu,n,s,t)
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
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
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
      t(j,1) = ppart(4,j+joff,l)
      t(j,2) = ppart(5,j+joff,l)
      t(j,3) = ppart(6,j+joff,l)
   20 continue
! deposit momentum flux within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
!dir$ ivdep
      do 30 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
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
! deposit momentum flux
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = amx*amz
      dy = amy*amz
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
   60 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 90 k = 2, ll
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
   70 continue
   80 continue
   90 continue
! deposit momentum flux to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 110 j = 2, mm
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  100 continue
  110 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 140 k = 1, ll
      do 120 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  120 continue
      do 130 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  130 continue
  140 continue
      if (lm > mz) then
         do 150 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  150    continue
         do 160 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  160    continue
      endif
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPGRMJPPOST32L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp, &
     &mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates particle momentum flux
! using first-order linear interpolation for relativistic particles
! for distributed data, with 2D spatial decomposition
! vectorizable/OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 134 flops/particle, 1 divide, 52 loads, 48 stores
! input: all, output: ppart, amu
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
! amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
! amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
! amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
! amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
! amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
! amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
! amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! and qci = qm*pj*pk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
! where pj = pj(t-dt/2) and pk = pk(t-dt/2)
! where gami2 = 1./(1.+sum(pi**2)*ci*ci)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum vz of particle n in partition in tile m
! at t - dt/2
! cu(i,j,k,l) = ith component of current density at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic = number of particles per tile
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! qm = charge on particle, in units of e
! ci = reciprical of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! mx/my/mz = number of grids in sorting cell in x/y/z
! nxv = second dimension of current array, must be >= nx+1
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
      real qm, ci
      real ppart, amu
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), amu(6,nxv*nypmx*nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxyp1, noffp, moffp, loffp, nppp, ipp, joff, nps
      integer i, j, k, l, m, mnoff, lnoff, nn, mm, ll, nm, lm
      integer lxv, lxyv, nxyv
      real ci2, gami, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy
      real x, y, z, vx, vy, vz, p2, v1, v2, v3, v4, v5, v6
      real samu
!     dimension samu(6,MXV*MYV*MZV)
      dimension samu(6,(mx+1)*(my+1)*(mz+1))
! scratch arrays
      integer n, mn
      real s, t
!dir$ attributes align : 64 :: n, mn, s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
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
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,m,noffp,moffp,loffp,nppp,ipp,joff,nps,nn,mm,ll,  &
!$OMP& mnoff,lnoff,nm,lm,x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx, &
!$OMP& dy,v1,v2,v3,v4,v5,v6,p2,gami,samu,n,s,t)
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
      samu(1,j) = 0.0
      samu(2,j) = 0.0
      samu(3,j) = 0.0
      samu(4,j) = 0.0
      samu(5,j) = 0.0
      samu(6,j) = 0.0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
      x = ppart(1,j+joff,l)
      y = ppart(2,j+joff,l)
      z = ppart(3,j+joff,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      vx = ppart(4,j+joff,l)
      vy = ppart(5,j+joff,l)
      vz = ppart(6,j+joff,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
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
      t(j,1) = vx*gami
      t(j,2) = vy*gami
      t(j,3) = vz*gami
   20 continue
! deposit momentum flux within tile to local accumulator
      do 40 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      vz = t(j,3)
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
!dir$ ivdep
      do 30 i = 1, lvect
      samu(1,n(j)+mn(i)) = samu(1,n(j)+mn(i)) + v1*s(j,i)
      samu(2,n(j)+mn(i)) = samu(2,n(j)+mn(i)) + v2*s(j,i)
      samu(3,n(j)+mn(i)) = samu(3,n(j)+mn(i)) + v3*s(j,i)
      samu(4,n(j)+mn(i)) = samu(4,n(j)+mn(i)) + v4*s(j,i)
      samu(5,n(j)+mn(i)) = samu(5,n(j)+mn(i)) + v5*s(j,i)
      samu(6,n(j)+mn(i)) = samu(6,n(j)+mn(i)) + v6*s(j,i)
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
! find inverse gamma
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + lxv*(mm - mnoff) + lxyv*(ll - lnoff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit momentum flux
      dx = amx*amz
      dy = amy*amz
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx
      samu(4,nn) = samu(4,nn) + v4*dx
      samu(5,nn) = samu(5,nn) + v5*dx
      samu(6,nn) = samu(6,nn) + v6*dx
      dx = dyp*amz
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy
      samu(4,nn+1) = samu(4,nn+1) + v4*dy
      samu(5,nn+1) = samu(5,nn+1) + v5*dy
      samu(6,nn+1) = samu(6,nn+1) + v6*dy
      dy = dx1*amz
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx
      samu(4,nn+lxv) = samu(4,nn+lxv) + v4*dx
      samu(5,nn+lxv) = samu(5,nn+lxv) + v5*dx
      samu(6,nn+lxv) = samu(6,nn+lxv) + v6*dx
      dx = amx*dzp
      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy
      samu(4,nn+1+lxv) = samu(4,nn+1+lxv) + v4*dy
      samu(5,nn+1+lxv) = samu(5,nn+1+lxv) + v5*dy
      samu(6,nn+1+lxv) = samu(6,nn+1+lxv) + v6*dy
      mm = nn + lxyv
      dy = amy*dzp
      samu(1,mm) = samu(1,mm) + v1*dx
      samu(2,mm) = samu(2,mm) + v2*dx
      samu(3,mm) = samu(3,mm) + v3*dx
      samu(4,mm) = samu(4,mm) + v4*dx
      samu(5,mm) = samu(5,mm) + v5*dx
      samu(6,mm) = samu(6,mm) + v6*dx
      dx = dyp*dzp
      samu(1,mm+1) = samu(1,mm+1) + v1*dy
      samu(2,mm+1) = samu(2,mm+1) + v2*dy
      samu(3,mm+1) = samu(3,mm+1) + v3*dy
      samu(4,mm+1) = samu(4,mm+1) + v4*dy
      samu(5,mm+1) = samu(5,mm+1) + v5*dy
      samu(6,mm+1) = samu(6,mm+1) + v6*dy
      dy = dx1*dzp
      samu(1,mm+lxv) = samu(1,mm+lxv) + v1*dx
      samu(2,mm+lxv) = samu(2,mm+lxv) + v2*dx
      samu(3,mm+lxv) = samu(3,mm+lxv) + v3*dx
      samu(4,mm+lxv) = samu(4,mm+lxv) + v4*dx
      samu(5,mm+lxv) = samu(5,mm+lxv) + v5*dx
      samu(6,mm+lxv) = samu(6,mm+lxv) + v6*dx
      samu(1,mm+1+lxv) = samu(1,mm+1+lxv) + v1*dy
      samu(2,mm+1+lxv) = samu(2,mm+1+lxv) + v2*dy
      samu(3,mm+1+lxv) = samu(3,mm+1+lxv) + v3*dy
      samu(4,mm+1+lxv) = samu(4,mm+1+lxv) + v4*dy
      samu(5,mm+1+lxv) = samu(5,mm+1+lxv) + v5*dy
      samu(6,mm+1+lxv) = samu(6,mm+1+lxv) + v6*dy
   60 continue
! deposit momentum flux to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 90 k = 2, ll
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(1,i+lxv*(j-1)+lxyv*(k-1))
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(2,i+lxv*(j-1)+lxyv*(k-1))
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(3,i+lxv*(j-1)+lxyv*(k-1))
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(4,i+lxv*(j-1)+lxyv*(k-1))
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(5,i+lxv*(j-1)+lxyv*(k-1))
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) +                 &
     & samu(6,i+lxv*(j-1)+lxyv*(k-1))
   70 continue
   80 continue
   90 continue
! deposit momentum flux to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 110 j = 2, mm
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(1,i+lxv*(j-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(2,i+lxv*(j-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(3,i+lxv*(j-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(4,i+lxv*(j-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(5,i+lxv*(j-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) =                       &
     &amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*loffp) + samu(6,i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,i+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,i+lxv*(j-1)+lxyv*(lm-1))
      endif
  100 continue
  110 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 140 k = 1, ll
      do 120 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(1,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(1,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(2,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(2,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(3,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(3,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(4,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(4,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(5,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(5,i+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) =                       &
     &amu(6,i+noffp+nxv*moffp+nxyv*(k+loffp-1)) + samu(6,i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,i+lxv*(mm-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,i+lxv*(mm-1)+lxyv*(k-1))
      endif
  120 continue
      do 130 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(1,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(2,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(3,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(4,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(5,1+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
      amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =                 &
     &amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))                   &
     &+ samu(6,1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(1,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(2,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(3,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(4,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(5,nm+lxv*(j-1)+lxyv*(k-1))
!$OMP ATOMIC
         amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1)) =             &
     &   amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(k+loffp-1))               &
     &   + samu(6,nm+lxv*(j-1)+lxyv*(k-1))
      endif
  130 continue
  140 continue
      if (lm > mz) then
         do 150 i = 2, nn
!$OMP ATOMIC
         amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(1,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(1,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(2,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(2,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(3,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(3,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(4,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(4,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(5,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(5,i+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1)) =                   &
     &   amu(6,i+noffp+nxv*moffp+nxyv*(lm+loffp-1))                     &
     &   + samu(6,i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,i+lxv*(mm-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,i+noffp+nxv*(mm+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  150    continue
         do 160 j = 1, mm
!$OMP ATOMIC
         amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(1,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(1,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(2,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(2,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(3,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(3,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(4,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(4,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(5,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(5,1+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
         amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =             &
     &   amu(6,1+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))               &
     &   + samu(6,1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(1,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(1,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(2,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(2,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(3,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(3,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(4,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(4,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(5,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(5,nm+lxv*(j-1)+lxyv*(lm-1))
!$OMP ATOMIC
            amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1)) =         &
     &      amu(6,nm+noffp+nxv*(j+moffp-1)+nxyv*(lm+loffp-1))           &
     &      + samu(6,nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  160    continue
      endif
  170 continue
!$OMP END PARALLEL DO
      return
      end
