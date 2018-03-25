!-----------------------------------------------------------------------
! Fortran library for spectral field solvers
! 3D MPI/OpenMP PIC Code:
! MPPOIS332 solves 3d poisson's equation for smoothed electric field
! MPADDQEI32 adds electron and ion densities
! MPPCUPERP32 calculates the transverse current in fourier space
! MIPPBPOISP332 solves 3d poisson's equation for unsmoothed magnetic
!               field
! MPPMAXWEL32 solves 3d maxwell's equation for unsmoothed transverse
!             electric and magnetic fields
! MPPEMFIELD32 adds and smooths or copies and smooths complex vector
!              fields in fourier space
! MPADDCUEI32 adds electron and ion current densities
! MPADDAMUI32 adds electron and ion momentum flux densities
! MPPBADDEXT32 adds constant to magnetic field for 3d code
! MPPADDVRFIELD32 calculates a = b + c
! MPPBBPOISP332 solves 3d poisson's equation in fourier space for
!               smoothed magnetic field
! MPPDCUPERP32 calculates transverse part of the derivative of the
!              current density from the momentum flux
! MPPADCUPERP32 calculates transverse part of the derivative of the
!               current density from the momentum flux and acceleration
!               density
! MPPEPOISP332 solves 3d poisson's equation in fourier space for
!              smoothed or unsmoothed transverse electric field
! MPPOTP32 solves 3d poisson's equation for potential
! MPPELFIELD32 solves 3d poisson's equation for unsmoothed electric
!              field
! MPPDIVF32 calculates the divergence in fourier space
! MPPGRADF32 calculates the gradient in fourier space
! MPPCURLF32 calculates the curl in fourier space
! MPPAVPOT332 calculates 3d vector potential from magnetic field
! MCUAVE33 averages current in fourier space for 3d code
! MPPAVRPOT332 solves 3d poisson's equation for the radiative part of	
!              the vector potential
! MPPAPOTP32 solves 3d poisson's equation for vector potential
! MPPETFIELD332 solves 3d poisson's equation in fourier space for
!               unsmoothed transverse electric field
! MPPSMOOTH32 provides a 3d scalar smoothing function
! MPPSMOOTH332 provides a 3d vector smoothing function
! PPRDMODES32 extracts lowest order scalar modes from packed array
!             stores them into a location in an unpacked array
! PPWRMODES32 extracts lowest order scalar modes from a location in an
!             unpacked array and stores them into a packed array
! PPRDVMODES32 extracts lowest order vector modes from packed array
!              stores them into a location in an unpacked array
! PPWRVMODES32 extracts lowest order vector modes from a location in an
!              unpacked array and stores them into a packed array
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 15, 2018
!-----------------------------------------------------------------------
      subroutine MPPOIS332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,  &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions for distributed data,
! with 2D spatial decomposition and OpenMP
! for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz,
!                       kxyp,kyzp,nzhd
! output: ffc
! for isign =/ 0, input: q,ffc,isign,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,
!                        kyzp,nzhd
! output: fxyz,we
! approximate flop count is:
! 62*nxc*nyc*nzc + 33*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the equation used is:
! fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
! fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
! fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
! fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
! fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
! fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
! q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
! fxyz(1,l,j,k) = x component of force/charge
! fxyz(2,l,j,k) = y component of force/charge
! fxyz(3,l,j,k) = z component of force/charge
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! if isign = 0, form factor array is prepared
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! ax/ay/az = half-width of particle in x/y/z direction
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! electric field energy is also calculated, using
! we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
!    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer isign, nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      integer nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nzv,kxyp,kyzp), fxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero, zt1, zt2
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      if (isign.ne.0) go to 40
      if (kstrt.gt.(nvpy*nvpz)) return
! prepare form factor array
      do 30 k = 1, kyzps
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*real(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyps
      dkx = dnx*real(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*real(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.0) then
         ffc(l,j,k) = cmplx(affp,1.0)
      else
         ffc(l,j,k) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
      return
! calculate force/charge and sum field energy
   40 sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 160
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,at1,at2,at3,at4,zt1,zt2,wp)     &
!$OMP& REDUCTION(+:sum1)
      do 60 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 50 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,k))*aimag(ffc(l,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,k)),-real(q(l,j,k)))
            zt2 = cmplx(aimag(q(l1,j,k)),-real(q(l1,j,k)))
            fxyz(1,l,j,k) = at2*zt1
            fxyz(2,l,j,k) = at3*zt1
            fxyz(3,l,j,k) = at4*zt1
            fxyz(1,l1,j,k) = at2*zt2
            fxyz(2,l1,j,k) = at3*zt2
            fxyz(3,l1,j,k) = -at4*zt2
            wp = wp + at1*(q(l,j,k)*conjg(q(l,j,k))                     &
     &              + q(l1,j,k)*conjg(q(l1,j,k)))
   50       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,k))*aimag(ffc(1,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(aimag(q(1,j,k)),-real(q(1,j,k)))
            fxyz(1,1,j,k) = at2*zt1
            fxyz(2,1,j,k) = at3*zt1
            fxyz(3,1,j,k) = zero
            fxyz(1,l1,j,k) = zero
            fxyz(2,l1,j,k) = zero
            fxyz(3,l1,j,k) = zero
            wp = wp + at1*(q(1,j,k)*conjg(q(1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   60 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,at1,at3,at4,zt1,zt2,wp)         &
!$OMP& REDUCTION(+:sum2)
      do 90 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 70 l = 2, nzh
               l1 = nz2 - l
               at1 = real(ffc(l,1,k))*aimag(ffc(l,1,k))
               at3 = dky*at1
               at4 = dnz*real(l - 1)*at1
               zt1 = cmplx(aimag(q(l,1,k)),-real(q(l,1,k)))
               zt2 = cmplx(aimag(q(l1,1,k)),-real(q(l1,1,k)))
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = at3*zt1
               fxyz(3,l,1,k) = at4*zt1
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = at3*zt2
               fxyz(3,l1,1,k) = -at4*zt2
               wp = wp + at1*(q(l,1,k)*conjg(q(l,1,k))                  &
     &                 + q(l1,1,k)*conjg(q(l1,1,k)))
   70          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = real(ffc(1,1,k))*aimag(ffc(1,1,k))
               at3 = dky*at1
               zt1 = cmplx(aimag(q(1,1,k)),-real(q(1,1,k)))
               fxyz(1,1,1,k) = zero
               fxyz(2,1,1,k) = at3*zt1
               fxyz(3,1,1,k) = zero
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = zero
               fxyz(3,l1,1,k) = zero
               wp = wp + at1*(q(1,1,k)*conjg(q(1,1,k)))
! throw away kx = nx/2
            else
               do 80 l = 1, nz
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = zero
               fxyz(3,l,1,k) = zero
   80          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   90 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,at1,at2,at4,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum2)
         do 110 j = 1, kxyps
         dkx = dnx*real(j + joff)
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 100 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,1))*aimag(ffc(l,j,1))
            at2 = dkx*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,1)),-real(q(l,j,1)))
            zt2 = cmplx(aimag(q(l1,j,1)),-real(q(l1,j,1)))
            fxyz(1,l,j,1) = at2*zt1
            fxyz(2,l,j,1) = zero
            fxyz(3,l,j,1) = at4*zt1
            fxyz(1,l1,j,1) = at2*zt2
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = -at4*zt2
            wp = wp + at1*(q(l,j,1)*conjg(q(l,j,1))                     &
     &              + q(l1,j,1)*conjg(q(l1,j,1)))
  100       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,1))*aimag(ffc(1,j,1))
            at2 = dkx*at1
            zt1 = cmplx(aimag(q(1,j,1)),-real(q(1,j,1)))
            fxyz(1,1,j,1) = at2*zt1
            fxyz(2,1,j,1) = zero
            fxyz(3,1,j,1) = zero
            fxyz(1,l1,j,1) = zero
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = zero
            wp = wp + at1*(q(1,j,1)*conjg(q(1,j,1)))
         endif
         sum2 = sum2 + wp
  110    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,1,1))*aimag(ffc(l,1,1))
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,1,1)),-real(q(l,1,1)))
            fxyz(1,l,1,1) = zero
            fxyz(2,l,1,1) = zero
            fxyz(3,l,1,1) = at4*zt1
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
            wp = wp + at1*(q(l,1,1)*conjg(q(l,1,1)))
  120       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fxyz(1,1,1,1) = zero
            fxyz(2,1,1,1) = zero
            fxyz(3,1,1,1) = zero
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 140 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 130 l = 1, nz
            fxyz(1,l,j,k1) = zero
            fxyz(2,l,j,k1) = zero
            fxyz(3,l,j,k1) = zero
  130       continue
         endif
  140    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 150 l = 1, nz
            fxyz(1,l,1,k1) = zero
            fxyz(2,l,1,k1) = zero
            fxyz(3,l,1,k1) = zero
  150       continue
         endif
      endif
      sum1 = sum1 + sum2
  160 continue
      we = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPADDQEI32(qe,qi,nyzp,nx,nxe,nypmx,nzpmx,idds)
! adds electron and ion densities
! assumes guard cells have already been added
! qe/qi = charge density for electrons/ions
! nyzp(1) = number of primary gridpoints in y in particle partition m
! nyzp(2) = number of primary gridpoints in z in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real qe, qi
      integer nyzp
      dimension qe(nxe,nypmx,nzpmx), qi(nxe,nypmx,nzpmx)
      dimension nyzp(idds)
      integer j, k, l
!$OMP PARALLEL DO PRIVATE(j,k,l)
      do 30 l = 1, nyzp(2)
      do 20 k = 1, nyzp(1)
      do 10 j = 1, nx
      qe(j,k,l) = qe(j,k,l) + qi(j,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCUPERP32(cu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! this subroutine calculates the transverse current in fourier space
! for distributed data with 2D spatial decomposition
! input: all output: cu
! approximate flop count is:
! 100*nxc*nyc*nzc + 36*(nxc*nyc + nxc*nzc + nyc*nzc)
! and (nx/2)*nyc*nzc divides
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the transverse current is calculated using the equation:
! cux(kx,ky,kz) = cux(kx,ky,kz) - kx*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
!                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! cuy(kx,ky,kz) = cuy(kx,ky,kz) - ky*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
!                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! cuz(kx,ky,kz) = cuz(kx,ky,kz) - kz*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
!                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! cux(kx=pi) = cuy(kx=pi) = cuz(kx=pi) = 0,
! cux(ky=pi) = cuy(ky=pi) = cux(ky=pi) = 0,
! cux(kz=pi) = cuy(kz=pi) = cuz(kz=pi) = 0,
! cux(kx=0,ky=0,kz=0) = cuy(kx=0,ky=0,kz=0) = cuz(kx=0,ky=0,kz=0) = 0.
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of procs in y/z
! nzv = second dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex cu
      dimension cu(3,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, dkx2, dky2, dkxy2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate transverse part of current
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,zt1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            zt1 = at1*(dkx*cu(1,l,j,k) + dky*cu(2,l,j,k)                &
     &          + dkz*cu(3,l,j,k))
            cu(1,l,j,k) = cu(1,l,j,k) - dkx*zt1
            cu(2,l,j,k) = cu(2,l,j,k) - dky*zt1
            cu(3,l,j,k) = cu(3,l,j,k) - dkz*zt1
            zt1 = at1*(dkx*cu(1,l1,j,k) + dky*cu(2,l1,j,k)              &
     &          - dkz*cu(3,l1,j,k))
            cu(1,l1,j,k) = cu(1,l1,j,k) - dkx*zt1
            cu(2,l1,j,k) = cu(2,l1,j,k) - dky*zt1
            cu(3,l1,j,k) = cu(3,l1,j,k) + dkz*zt1
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            zt1 = at1*(dkx*cu(1,1,j,k) + dky*cu(2,1,j,k))
            cu(1,1,j,k) = cu(1,1,j,k) - dkx*zt1
            cu(2,1,j,k) = cu(2,1,j,k) - dky*zt1
            cu(1,l1,j,k) = zero
            cu(2,l1,j,k) = zero
            cu(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,dky2,at1,zt1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               zt1 = at1*(dky*cu(2,l,1,k) + dkz*cu(3,l,1,k))
               cu(2,l,1,k) = cu(2,l,1,k) - dky*zt1
               cu(3,l,1,k) = cu(3,l,1,k) - dkz*zt1
               zt1 = at1*(dky*cu(2,l1,1,k) - dkz*cu(3,l1,1,k))
               cu(2,l1,1,k) = cu(2,l1,1,k) - dky*zt1
               cu(3,l1,1,k) = cu(3,l1,1,k) + dkz*zt1
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               cu(2,1,1,k) = zero
               cu(1,l1,1,k) = zero
               cu(2,l1,1,k) = zero
               cu(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               cu(1,l,1,k) = zero
               cu(2,l,1,k) = zero
               cu(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,dkx2,at1,zt1)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            zt1 = at1*(dkx*cu(1,l,j,1) + dkz*cu(3,l,j,1))
            cu(1,l,j,1) = cu(1,l,j,1) - dkx*zt1
            cu(3,l,j,1) = cu(3,l,j,1) - dkz*zt1
            zt1 = at1*(dkx*cu(1,l1,j,1) - dkz*cu(3,l1,j,1))
            cu(1,l1,j,1) = cu(1,l1,j,1) - dkx*zt1
            cu(3,l1,j,1) = cu(3,l1,j,1) + dkz*zt1
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            cu(1,1,j,1) = zero
            cu(1,l1,j,1) = zero
            cu(2,l1,j,1) = zero
            cu(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            cu(3,l,1,1) = zero
            cu(1,l1,1,1) = zero
            cu(2,l1,1,1) = zero
            cu(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            cu(1,1,1,1) = zero
            cu(2,1,1,1) = zero
            cu(3,1,1,1) = zero
            cu(1,l1,1,1) = zero
            cu(2,l1,1,1) = zero
            cu(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            cu(1,l,j,k1) = zero
            cu(2,l,j,k1) = zero
            cu(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            cu(1,l,1,k1) = zero
            cu(2,l,1,k1) = zero
            cu(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MIPPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,   &
     &nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! magnetic field with periodic boundary conditions, for distributed data
! with 2D spatial decomposition and OpenMP
! input: cu,ffc,ci,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
! output: bxyz, wm
! approximate flop count is:
! 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the magnetic field is calculated using the equations:
! bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz)),
! by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz)),
! bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz)),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
! bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
! bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
! cu(l,j,k) = complex current density for fourier mode jj-1,kk-1,l-1
! bxyz(i,l,j,k) = i component of complex magnetic field
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! ci = reciprical of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real ci, wm
      complex cu, bxyz, ffc
      dimension cu(3,nzv,kxyp,kyzp), bxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real ci2, dnx, dny, dnz, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,at1,at2,at3,at4,zt1,zt2,zt3,wp) &
!$OMP& REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,k))
            zt1 = cmplx(-aimag(cu(3,l,j,k)),real(cu(3,l,j,k)))
            zt2 = cmplx(-aimag(cu(2,l,j,k)),real(cu(2,l,j,k)))
            zt3 = cmplx(-aimag(cu(1,l,j,k)),real(cu(1,l,j,k)))
            bxyz(1,l,j,k) = at3*zt1 - at4*zt2
            bxyz(2,l,j,k) = at4*zt3 - at2*zt1
            bxyz(3,l,j,k) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(cu(3,l1,j,k)),real(cu(3,l1,j,k)))
            zt2 = cmplx(-aimag(cu(2,l1,j,k)),real(cu(2,l1,j,k)))
            zt3 = cmplx(-aimag(cu(1,l1,j,k)),real(cu(1,l1,j,k)))
            bxyz(1,l1,j,k) = at3*zt1 + at4*zt2
            bxyz(2,l1,j,k) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,k) = at2*zt2 - at3*zt3
            wp = wp + at1*(cu(1,l,j,k)*conjg(cu(1,l,j,k))               &
     &              + cu(2,l,j,k)*conjg(cu(2,l,j,k))                    &
     &              + cu(3,l,j,k)*conjg(cu(3,l,j,k))                    &
     &              + cu(1,l1,j,k)*conjg(cu(1,l1,j,k))                  &
     &              + cu(2,l1,j,k)*conjg(cu(2,l1,j,k))                  &
     &              + cu(3,l1,j,k)*conjg(cu(3,l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = at1*aimag(ffc(1,j,k))
            zt1 = cmplx(-aimag(cu(3,1,j,k)),real(cu(3,1,j,k)))
            zt2 = cmplx(-aimag(cu(2,1,j,k)),real(cu(2,1,j,k)))
            zt3 = cmplx(-aimag(cu(1,1,j,k)),real(cu(1,1,j,k)))
            bxyz(1,1,j,k) = at3*zt1
            bxyz(2,1,j,k) = -at2*zt1
            bxyz(3,1,j,k) = at2*zt2 - at3*zt3
            bxyz(1,l1,j,k) = zero
            bxyz(2,l1,j,k) = zero
            bxyz(3,l1,j,k) = zero
            wp = wp + at1*(cu(1,1,j,k)*conjg(cu(1,1,j,k))               &
     &              + cu(2,1,j,k)*conjg(cu(2,1,j,k))                    &
     &              + cu(3,1,j,k)*conjg(cu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,at1,at3,at4,zt1,zt2,zt3,wp)     &
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = ci2*real(ffc(l,1,k))
               at3 = dky*at1
               at4 = dnz*real(l - 1)*at1
               at1 = at1*aimag(ffc(l,1,k))
               zt1 = cmplx(-aimag(cu(3,l,1,k)),real(cu(3,l,1,k)))
               zt2 = cmplx(-aimag(cu(2,l,1,k)),real(cu(2,l,1,k)))
               zt3 = cmplx(-aimag(cu(1,l,1,k)),real(cu(1,l,1,k)))
               bxyz(1,l,1,k) = at3*zt1 - at4*zt2
               bxyz(2,l,1,k) = at4*zt3
               bxyz(3,l,1,k) = -at3*zt3
               zt1 = cmplx(-aimag(cu(3,l1,1,k)),real(cu(3,l1,1,k)))
               zt2 = cmplx(-aimag(cu(2,l1,1,k)),real(cu(2,l1,1,k)))
               zt3 = cmplx(-aimag(cu(1,l1,1,k)),real(cu(1,l1,1,k)))
               bxyz(1,l1,1,k) = at3*zt1 + at4*zt2
               bxyz(2,l1,1,k) = -at4*zt3
               bxyz(3,l1,1,k) = -at3*zt3
               wp = wp + at1*(cu(1,l,1,k)*conjg(cu(1,l,1,k))            &
     &                 + cu(2,l,1,k)*conjg(cu(2,l,1,k))                 &
     &                 + cu(3,l,1,k)*conjg(cu(3,l,1,k))                 &
     &                 + cu(1,l1,1,k)*conjg(cu(1,l1,1,k))               &
     &                 + cu(2,l1,1,k)*conjg(cu(2,l1,1,k))               &
     &                 + cu(3,l1,1,k)*conjg(cu(3,l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = ci2*real(ffc(1,1,k))
               at3 = dky*at1
               at1 = at1*aimag(ffc(1,1,k))
               zt1 = cmplx(-aimag(cu(3,1,1,k)),real(cu(3,1,1,k)))
               zt3 = cmplx(-aimag(cu(1,1,1,k)),real(cu(1,1,1,k)))
               bxyz(1,1,1,k) = at3*zt1
               bxyz(2,1,1,k) = zero
               bxyz(3,1,1,k) = -at3*zt3
               bxyz(1,l1,1,k) = zero
               bxyz(2,l1,1,k) = zero
               bxyz(3,l1,1,k) = zero
               wp = wp + at1*(cu(1,1,1,k)*conjg(cu(1,1,1,k))            &
     &                 + cu(2,1,1,k)*conjg(cu(2,1,1,k))                 &
     &                 + cu(3,1,1,k)*conjg(cu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k) = zero
               bxyz(2,l,1,k) = zero
               bxyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,at1,at2,at4,zt1,zt2,zt3,wp)        &
!$OMP& REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,1))
            at2 = dkx*at1
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,1))
            zt1 = cmplx(-aimag(cu(3,l,j,1)),real(cu(3,l,j,1)))
            zt2 = cmplx(-aimag(cu(2,l,j,1)),real(cu(2,l,j,1)))
            zt3 = cmplx(-aimag(cu(1,l,j,1)),real(cu(1,l,j,1)))
            bxyz(1,l,j,1) = -at4*zt2
            bxyz(2,l,j,1) = at4*zt3 - at2*zt1
            bxyz(3,l,j,1) = at2*zt2
            zt1 = cmplx(-aimag(cu(3,l1,j,1)),real(cu(3,l1,j,1)))
            zt2 = cmplx(-aimag(cu(2,l1,j,1)),real(cu(2,l1,j,1)))
            zt3 = cmplx(-aimag(cu(1,l1,j,1)),real(cu(1,l1,j,1)))
            bxyz(1,l1,j,1) = at4*zt2
            bxyz(2,l1,j,1) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,1) = at2*zt2
            wp = wp + at1*(cu(1,l,j,1)*conjg(cu(1,l,j,1))               &
     &              + cu(2,l,j,1)*conjg(cu(2,l,j,1))                    &
     &              + cu(3,l,j,1)*conjg(cu(3,l,j,1))                    &
     &              + cu(1,l1,j,1)*conjg(cu(1,l1,j,1))                  &
     &              + cu(2,l1,j,1)*conjg(cu(2,l1,j,1))                  &
     &              + cu(3,l1,j,1)*conjg(cu(3,l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,1))
            at2 = dkx*at1
            at1 = at1*aimag(ffc(1,j,1))
            zt1 = cmplx(-aimag(cu(3,1,j,1)),real(cu(3,1,j,1)))
            zt2 = cmplx(-aimag(cu(2,1,j,1)),real(cu(2,1,j,1)))
            bxyz(1,1,j,1) = zero
            bxyz(2,1,j,1) = -at2*zt1
            bxyz(3,1,j,1) = at2*zt2
            bxyz(1,l1,j,1) = zero
            bxyz(2,l1,j,1) = zero
            bxyz(3,l1,j,1) = zero
            wp = wp + at1*(cu(1,1,j,1)*conjg(cu(1,1,j,1))               &
     &              + cu(2,1,j,1)*conjg(cu(2,1,j,1))                    &
     &              + cu(3,1,j,1)*conjg(cu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,1,1))
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,1,1))
            zt2 = cmplx(-aimag(cu(2,l,1,1)),real(cu(2,l,1,1)))
            zt3 = cmplx(-aimag(cu(1,l,1,1)),real(cu(1,l,1,1)))
            bxyz(1,l,1,1) = -at4*zt2
            bxyz(2,l,1,1) = at4*zt3
            bxyz(3,l,1,1) = zero
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
            wp = wp + at1*(cu(1,l,1,1)*conjg(cu(1,l,1,1))               &
     &              + cu(2,l,1,1)*conjg(cu(2,l,1,1))                    &
     &              + cu(3,l,1,1)*conjg(cu(3,l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1) = zero
            bxyz(2,1,1,1) = zero
            bxyz(3,1,1,1) = zero
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            bxyz(1,l,j,k1) = zero
            bxyz(2,l,j,k1) = zero
            bxyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            bxyz(1,l,1,k1) = zero
            bxyz(2,l,1,k1) = zero
            bxyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      wm = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,nz,&
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d maxwell's equation in fourier space for
! transverse electric and magnetic fields with periodic boundary
! conditions, for distributed data with 2D spatial decomposition
! with OpenMP
! input: all, output: wf, wm, exyz, bxyz
! approximate flop count is:
! 679*nxc*nyc*nzc + 149*(nxc*nyc + nxc*nzc + nyc*nzc)
! plus nxc*nyc*nzc divides
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz - 1, and
! nvpy/nvpz = number of procs in y/z
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky,kz) = bx(kx,ky,kz) - .5*dt*sqrt(-1)*
!                (ky*ez(kx,ky,kz)-kz*ey(kx,ky,kz))
! by(kx,ky,kz) = by(kx,ky,kz) - .5*dt*sqrt(-1)*
!               (kz*ex(kx,ky,kz)-kx*ez(kx,ky,kz))
! bz(kx,ky,kz) = bz(kx,ky,kz) - .5*dt*sqrt(-1)*
!               (kx*ey(kx,ky,kz)-ky*ex(kx,ky,kz))
! the electric field is then updated a whole step using the equations:
! ex(kx,ky,kz) = ex(kx,ky,kz) + c2*dt*sqrt(-1)*
!  (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz)) - affp*dt*cux(kx,ky,kz)*s(kx,ky,kz)
! ey(kx,ky,kz) = ey(kx,ky,kz) + c2*dt*sqrt(-1)*
!  (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz)) - affp*dt*cuy(kx,ky,kz)*s(kx,ky,kz)
! ez(kx,ky,kz) = ez(kx,ky,kz) + c2*dt*sqrt(-1)*
!  (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz)) - affp*dt*cuz(kx,ky,kz)*s(kx,ky,kz)
! the magnetic field is finally updated the remaining half step with
! the new electric field and the previous magnetic field equations.
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
! and s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
! j,k,l = fourier mode numbers, except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
! ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
! ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
! and similarly for bx, by, bz.
! exyz(i,l,j,k) = i component of complex transverse electric field
! bxyz(i,l,j,k) = i component of complex magnetic field
! cu(i,l,j,k) = i component of complex current density
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s,
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
! magnetic field energy is also calculated, using
! wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real affp, ci, dt, wf, wm
      complex exyz, bxyz, cu, ffc
      dimension exyz(3,nzv,kxyp,kyzp), bxyz(3,nzv,kxyp,kyzp)
      dimension cu(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dth, c2, cdt, adt, anorm, dkx, dky, dkz, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws, sum1, sum2, sum3, sum4
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum1 = 0.0d0
      sum2 = 0.0d0
! calculate the electromagnetic fields
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkz,dky,dkx,afdt,zt1,zt2,zt3,zt4,zt5,zt6&
!$OMP& ,zt7,zt8,zt9,ws,wp) REDUCTION(+:sum1,sum2)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      ws = 0.0d0
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            afdt = adt*aimag(ffc(l,j,k))
! update magnetic field half time step, ky > 0, kz > 0
            zt1 = cmplx(-aimag(exyz(3,l,j,k)),real(exyz(3,l,j,k)))
            zt2 = cmplx(-aimag(exyz(2,l,j,k)),real(exyz(2,l,j,k)))
            zt3 = cmplx(-aimag(exyz(1,l,j,k)),real(exyz(1,l,j,k)))
            zt4 = bxyz(1,l,j,k) - dth*(dky*zt1 - dkz*zt2)
            zt5 = bxyz(2,l,j,k) - dth*(dkz*zt3 - dkx*zt1)
            zt6 = bxyz(3,l,j,k) - dth*(dkx*zt2 - dky*zt3)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,j,k) + cdt*(dky*zt1 - dkz*zt2)               &
     &          - afdt*cu(1,l,j,k)
            zt8 = exyz(2,l,j,k) + cdt*(dkz*zt3 - dkx*zt1)               &
     &          - afdt*cu(2,l,j,k)
            zt9 = exyz(3,l,j,k) + cdt*(dkx*zt2 - dky*zt3)               &
     &          - afdt*cu(3,l,j,k)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,j,k) = zt7
            exyz(2,l,j,k) = zt8
            exyz(3,l,j,k) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)            &
     &                     + zt9*conjg(zt9))
            zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,l,j,k) = zt4
            bxyz(2,l,j,k) = zt5
            bxyz(3,l,j,k) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)            &
     &                     + zt6*conjg(zt6))
! update magnetic field half time step, ky > 0, kz < 0
            zt1 = cmplx(-aimag(exyz(3,l1,j,k)),real(exyz(3,l1,j,k)))
            zt2 = cmplx(-aimag(exyz(2,l1,j,k)),real(exyz(2,l1,j,k)))
            zt3 = cmplx(-aimag(exyz(1,l1,j,k)),real(exyz(1,l1,j,k)))
            zt4 = bxyz(1,l1,j,k) - dth*(dky*zt1 + dkz*zt2)
            zt5 = bxyz(2,l1,j,k) + dth*(dkz*zt3 + dkx*zt1)
            zt6 = bxyz(3,l1,j,k) - dth*(dkx*zt2 - dky*zt3)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l1,j,k) + cdt*(dky*zt1 + dkz*zt2)              &
     &          - afdt*cu(1,l1,j,k)
            zt8 = exyz(2,l1,j,k) - cdt*(dkz*zt3 + dkx*zt1)              &
     &          - afdt*cu(2,l1,j,k)
            zt9 = exyz(3,l1,j,k) + cdt*(dkx*zt2 - dky*zt3)              &
     &          - afdt*cu(3,l1,j,k)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l1,j,k) = zt7
            exyz(2,l1,j,k) = zt8
            exyz(3,l1,j,k) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)            &
     &                     + zt9*conjg(zt9))
            zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
            zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,l1,j,k) = zt4
            bxyz(2,l1,j,k) = zt5
            bxyz(3,l1,j,k) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)            &
     &                     + zt6*conjg(zt6))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            afdt = adt*aimag(ffc(1,j,k))
! update magnetic field half time step, ky > 0
            zt1 = cmplx(-aimag(exyz(3,1,j,k)),real(exyz(3,1,j,k)))
            zt2 = cmplx(-aimag(exyz(2,1,j,k)),real(exyz(2,1,j,k)))
            zt3 = cmplx(-aimag(exyz(1,1,j,k)),real(exyz(1,1,j,k)))
            zt4 = bxyz(1,1,j,k) - dth*(dky*zt1)
            zt5 = bxyz(2,1,j,k) + dth*(dkx*zt1)
            zt6 = bxyz(3,1,j,k) - dth*(dkx*zt2 - dky*zt3)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,1,j,k) + cdt*(dky*zt1) - afdt*cu(1,1,j,k)
            zt8 = exyz(2,1,j,k) - cdt*(dkx*zt1) - afdt*cu(2,1,j,k)
            zt9 = exyz(3,1,j,k) + cdt*(dkx*zt2 - dky*zt3)               &
     &          - afdt*cu(3,1,j,k)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,1,j,k) = zt7
            exyz(2,1,j,k) = zt8
            exyz(3,1,j,k) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)            &
     &                     + zt9*conjg(zt9))
            zt4 = zt4 - dth*(dky*zt1)
            zt5 = zt5 + dth*(dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,1,j,k) = zt4
            bxyz(2,1,j,k) = zt5
            bxyz(3,1,j,k) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)            &
     &                     + zt6*conjg(zt6))
            bxyz(1,l1,j,k) = zero
            bxyz(2,l1,j,k) = zero
            bxyz(3,l1,j,k) = zero
            exyz(1,l1,j,k) = zero
            exyz(2,l1,j,k) = zero
            exyz(3,l1,j,k) = zero
         endif
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum3 = 0.0d0
      sum4 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dkz,dky,afdt,zt1,zt2,zt3,zt4,zt5,zt6&
!$OMP& ,zt7,zt8,zt9,ws,wp) REDUCTION(+:sum3,sum4)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         ws = 0.0d0
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               afdt = adt*aimag(ffc(l,1,k))
! update magnetic field half time step, kz > 0
               zt1 = cmplx(-aimag(exyz(3,l,1,k)),real(exyz(3,l,1,k)))
               zt2 = cmplx(-aimag(exyz(2,l,1,k)),real(exyz(2,l,1,k)))
               zt3 = cmplx(-aimag(exyz(1,l,1,k)),real(exyz(1,l,1,k)))
               zt4 = bxyz(1,l,1,k) - dth*(dky*zt1 - dkz*zt2)
               zt5 = bxyz(2,l,1,k) - dth*(dkz*zt3)
               zt6 = bxyz(3,l,1,k) + dth*(dky*zt3)
! update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt2 = cmplx(-aimag(zt5),real(zt5))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,l,1,k) + cdt*(dky*zt1 - dkz*zt2)            &
     &             - afdt*cu(1,l,1,k)
               zt8 = exyz(2,l,1,k) + cdt*(dkz*zt3) - afdt*cu(2,l,1,k)
               zt9 = exyz(3,l,1,k) - cdt*(dky*zt3) - afdt*cu(3,l,1,k)
! update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt2 = cmplx(-aimag(zt8),real(zt8))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,l,1,k) = zt7
               exyz(2,l,1,k) = zt8
               exyz(3,l,1,k) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)         &
     &                        + zt9*conjg(zt9))
               zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
               zt5 = zt5 - dth*(dkz*zt3)
               zt6 = zt6 + dth*(dky*zt3) 
               bxyz(1,l,1,k) = zt4
               bxyz(2,l,1,k) = zt5
               bxyz(3,l,1,k) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)         &
     &                        + zt6*conjg(zt6))
! update magnetic field half time step, kz < 0
               zt1 = cmplx(-aimag(exyz(3,l1,1,k)),real(exyz(3,l1,1,k)))
               zt2 = cmplx(-aimag(exyz(2,l1,1,k)),real(exyz(2,l1,1,k)))
               zt3 = cmplx(-aimag(exyz(1,l1,1,k)),real(exyz(1,l1,1,k)))
               zt4 = bxyz(1,l1,1,k) - dth*(dky*zt1 + dkz*zt2)
               zt5 = bxyz(2,l1,1,k) + dth*(dkz*zt3)
               zt6 = bxyz(3,l1,1,k) + dth*(dky*zt3)
! update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt2 = cmplx(-aimag(zt5),real(zt5))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,l1,1,k) + cdt*(dky*zt1 + dkz*zt2)           &
     &             - afdt*cu(1,l1,1,k)
               zt8 = exyz(2,l1,1,k) - cdt*(dkz*zt3) - afdt*cu(2,l1,1,k)
               zt9 = exyz(3,l1,1,k) - cdt*(dky*zt3) - afdt*cu(3,l1,1,k)
! update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt2 = cmplx(-aimag(zt8),real(zt8))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,l1,1,k) = zt7
               exyz(2,l1,1,k) = zt8
               exyz(3,l1,1,k) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)         &
     &                        + zt9*conjg(zt9))
               zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
               zt5 = zt5 + dth*(dkz*zt3)
               zt6 = zt6 + dth*(dky*zt3)
               bxyz(1,l1,1,k) = zt4
               bxyz(2,l1,1,k) = zt5
               bxyz(3,l1,1,k) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)         &
     &                        + zt6*conjg(zt6))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               afdt = adt*aimag(ffc(1,1,k))
! update magnetic field half time step
               zt1 = cmplx(-aimag(exyz(3,1,1,k)),real(exyz(3,1,1,k)))
               zt3 = cmplx(-aimag(exyz(1,1,1,k)),real(exyz(1,1,1,k)))
               zt4 = bxyz(1,1,1,k) - dth*(dky*zt1)
               zt6 = bxyz(3,1,1,k) + dth*(dky*zt3)
! update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,1,k)
               zt9 = exyz(3,1,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,1,k)
! update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,1,1,k) = zt7
               exyz(2,1,1,k) = zero
               exyz(3,1,1,k) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
               zt4 = zt4 - dth*(dky*zt1)
               zt6 = zt6 + dth*(dky*zt3)
               bxyz(1,1,1,k) = zt4
               bxyz(2,1,1,k) = zero
               bxyz(3,1,1,k) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
               bxyz(1,l1,1,k) = zero
               bxyz(2,l1,1,k) = zero
               bxyz(3,l1,1,k) = zero
               exyz(1,l1,1,k) = zero
               exyz(2,l1,1,k) = zero
               exyz(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k) = zero
               bxyz(2,l,1,k) = zero
               bxyz(3,l,1,k) = zero
               exyz(1,l,1,k) = zero
               exyz(2,l,1,k) = zero
               exyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum3 = sum3 + ws
         sum4 = sum4 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum3
      sum2 = sum2 + sum4
! mode numbers ky = 0, ny/2
      sum3 = 0.0d0
      sum4 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,k1,l1,dkz,dkx,afdt,zt1,zt2,zt3,zt4,zt5,zt6&
!$OMP& ,zt7,zt8,zt9,ws,wp) REDUCTION(+:sum3,sum4)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         ws = 0.0d0
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            afdt = adt*aimag(ffc(l,j,1))
            zt1 = cmplx(-aimag(exyz(3,l,j,1)),real(exyz(3,l,j,1)))
            zt2 = cmplx(-aimag(exyz(2,l,j,1)),real(exyz(2,l,j,1)))
            zt3 = cmplx(-aimag(exyz(1,l,j,1)),real(exyz(1,l,j,1)))
            zt4 = bxyz(1,l,j,1) + dth*(dkz*zt2)
            zt5 = bxyz(2,l,j,1) - dth*(dkz*zt3 - dkx*zt1)
            zt6 = bxyz(3,l,j,1) - dth*(dkx*zt2)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,j,1) - cdt*(dkz*zt2) - afdt*cu(1,l,j,1)
            zt8 = exyz(2,l,j,1) + cdt*(dkz*zt3 - dkx*zt1)               &
     &          - afdt*cu(2,l,j,1)
            zt9 = exyz(3,l,j,1) + cdt*(dkx*zt2) - afdt*cu(3,l,j,1)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,j,1) = zt7
            exyz(2,l,j,1) = zt8
            exyz(3,l,j,1) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)            &
     &                     + zt9*conjg(zt9))
            zt4 = zt4 + dth*(dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,l,j,1) = zt4
            bxyz(2,l,j,1) = zt5
            bxyz(3,l,j,1) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)            &
     &                     + zt6*conjg(zt6))
! update magnetic field half time step, kz < 0
            zt1 = cmplx(-aimag(exyz(3,l1,j,1)),real(exyz(3,l1,j,1)))
            zt2 = cmplx(-aimag(exyz(2,l1,j,1)),real(exyz(2,l1,j,1)))
            zt3 = cmplx(-aimag(exyz(1,l1,j,1)),real(exyz(1,l1,j,1)))
            zt4 = bxyz(1,l1,j,1) - dth*(dkz*zt2)
            zt5 = bxyz(2,l1,j,1) + dth*(dkz*zt3 + dkx*zt1)
            zt6 = bxyz(3,l1,j,1) - dth*(dkx*zt2)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l1,j,1) + cdt*(dkz*zt2) - afdt*cu(1,l1,j,1)
            zt8 = exyz(2,l1,j,1) - cdt*(dkz*zt3 + dkx*zt1)              &
     &          - afdt*cu(2,l1,j,1)
            zt9 = exyz(3,l1,j,1) + cdt*(dkx*zt2) - afdt*cu(3,l1,j,1)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l1,j,1) = zt7
            exyz(2,l1,j,1) = zt8
            exyz(3,l1,j,1) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)            &
     &                     + zt9*conjg(zt9))
            zt4 = zt4 - dth*(dkz*zt2)
            zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,l1,j,1) = zt4
            bxyz(2,l1,j,1) = zt5
            bxyz(3,l1,j,1) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)            &
     &                     + zt6*conjg(zt6))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            afdt = adt*aimag(ffc(1,j,1))
! update magnetic field half time step
            zt1 = cmplx(-aimag(exyz(3,1,j,1)),real(exyz(3,1,j,1)))
            zt2 = cmplx(-aimag(exyz(2,1,j,1)),real(exyz(2,1,j,1)))
            zt5 = bxyz(2,1,j,1) + dth*(dkx*zt1)
            zt6 = bxyz(3,1,j,1) - dth*(dkx*zt2)
! update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt8 = exyz(2,1,j,1) - cdt*(dkx*zt1) - afdt*cu(2,1,j,1)
            zt9 = exyz(3,1,j,1) + cdt*(dkx*zt2) - afdt*cu(3,1,j,1)
! update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            exyz(1,1,j,1) = zero
            exyz(2,1,j,1) = zt8
            exyz(3,1,j,1) = zt9
            ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
            zt5 = zt5 + dth*(dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,1,j,1) = zero
            bxyz(2,1,j,1) = zt5
            bxyz(3,1,j,1) = zt6
            wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
            bxyz(1,l1,j,1) = zero
            bxyz(2,l1,j,1) = zero
            bxyz(3,l1,j,1) = zero
            exyz(1,l1,j,1) = zero
            exyz(2,l1,j,1) = zero
            exyz(3,l1,j,1) = zero
         endif
         sum3 = sum3 + ws
         sum4 = sum4 + wp
   70    continue
!$OMP END PARALLEL DO
         ws = 0.0d0
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            afdt = adt*aimag(ffc(l,1,1))
! update magnetic field half time step
            zt2 = cmplx(-aimag(exyz(2,l,1,1)),real(exyz(2,l,1,1)))
            zt3 = cmplx(-aimag(exyz(1,l,1,1)),real(exyz(1,l,1,1)))
            zt4 = bxyz(1,l,1,1) + dth*(dkz*zt2)
            zt5 = bxyz(2,l,1,1) - dth*(dkz*zt3)
! update electric field whole time step
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,1,1) - cdt*(dkz*zt2) - afdt*cu(1,l,1,1)
            zt8 = exyz(2,l,1,1) + cdt*(dkz*zt3) - afdt*cu(2,l,1,1)
! update magnetic field half time step and store electric field
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,1,1) = zt7
            exyz(2,l,1,1) = zt8
            exyz(3,l,1,1) = zero
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8))
            zt4 = zt4 + dth*(dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3)
            bxyz(1,l,1,1) = zt4
            bxyz(2,l,1,1) = zt5
            bxyz(3,l,1,1) = zero
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5))
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1) = zero
            bxyz(2,1,1,1) = zero
            bxyz(3,1,1,1) = zero
            exyz(1,1,1,1) = zero
            exyz(2,1,1,1) = zero
            exyz(3,1,1,1) = zero
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
         endif
         sum3 = sum3 + ws
         sum4 = sum4 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            bxyz(1,l,j,k1) = zero
            bxyz(2,l,j,k1) = zero
            bxyz(3,l,j,k1) = zero
            exyz(1,l,j,k1) = zero
            exyz(2,l,j,k1) = zero
            exyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            bxyz(1,l,1,k1) = zero
            bxyz(2,l,1,k1) = zero
            bxyz(3,l,1,k1) = zero
            exyz(1,l,1,k1) = zero
            exyz(2,l,1,k1) = zero
            exyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum3
      sum2 = sum2 + sum4
  120 continue      
      wf = real(nx)*real(ny)*real(nz)*sum1
      wm = real(nx)*real(ny)*real(nz)*c2*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nvpy,  &
     &nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine either adds complex vector fields if isign > 0
! or copies complex vector fields if isign < 0
! includes additional smoothing
! uses OpenMP
      implicit none
      integer isign, nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      integer nzhd
      complex fxyz, exyz, ffc
      dimension fxyz(3,nzv,kxyp,kyzp), exyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer i, j, k, l, nxh, nzh, nz2, js, ks, kxyps, kyzps, l1, kk
      real at1
      nxh = nx/2
      nzh = max(1,nz/2)
      nz2 = nz + 2
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nxh-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
! add the fields
      if (isign.gt.0) then
! mode numbers 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,l,kk,l1,at1)
         do 30 kk = 1, kyzps*kxyps
         k = (kk - 1)/kxyps
         j = kk - kxyps*k
         k = k + 1
         do 20 l = 2, nzh
         l1 = nz2 - l
         at1 = aimag(ffc(l,j,k))
         do 10 i = 1, 3
         fxyz(i,l,j,k) = fxyz(i,l,j,k) + exyz(i,l,j,k)*at1
         fxyz(i,l1,j,k) = fxyz(i,l1,j,k) + exyz(i,l1,j,k)*at1
   10    continue
   20    continue
   30    continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kz = 0, nz/2
!$OMP PARALLEL DO PRIVATE(i,j,k,l1,at1)
         do 60 k = 1, kyzps
         do 50 j = 1, kxyps
         l1 = nzh + 1
         at1 = aimag(ffc(1,j,k))
         do 40 i = 1, 3
         fxyz(i,1,j,k) = fxyz(i,1,j,k) + exyz(i,1,j,k)*at1
         fxyz(i,l1,j,k) = fxyz(i,l1,j,k) + exyz(i,l1,j,k)*at1
   40    continue
   50    continue
   60    continue
!$OMP END PARALLEL DO
! copy the fields
      else if (isign.lt.0) then
! mode numbers 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,l,kk,l1,at1)
         do 90 kk = 1, kyzps*kxyps
         k = (kk - 1)/kxyps
         j = kk - kxyps*k
         k = k + 1
         do 80 l = 2, nzh
         l1 = nz2 - l
         at1 = aimag(ffc(l,j,k))
         do 70 i = 1, 3
         fxyz(i,l,j,k) = exyz(i,l,j,k)*at1
         fxyz(i,l1,j,k) = exyz(i,l1,j,k)*at1
   70    continue
   80    continue
   90    continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kz = 0, nz/2
!$OMP PARALLEL DO PRIVATE(i,j,k,l1,at1)
         do 120 k = 1, kyzps
         do 110 j = 1, kxyps
         l1 = nzh + 1
         at1 = aimag(ffc(1,j,k))
         do 100 i = 1, 3
         fxyz(i,1,j,k) = exyz(i,1,j,k)*at1
         fxyz(i,l1,j,k) = exyz(i,l1,j,k)*at1
  100    continue
  110    continue
  120    continue
!$OMP END PARALLEL DO
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPADDCUEI32(cue,cui,nyzp,nx,nxe,nypmx,nzpmx,idds)
! adds electron and ion current densities
! assumes guard cells have already been added
! cue/cui = current density for electrons/ions
! nyzp(1) = number of primary gridpoints in y in particle partition m
! nyzp(2) = number of primary gridpoints in z in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real cue, cui
      integer nyzp
      dimension cue(3,nxe,nypmx,nzpmx), cui(3,nxe,nypmx,nzpmx)
      dimension nyzp(idds)
      integer j, k, l
!$OMP PARALLEL DO PRIVATE(j,k,l)
      do 30 l = 1, nyzp(2)
      do 20 k = 1, nyzp(1)
      do 10 j = 1, nx
      cue(1,j,k,l) = cue(1,j,k,l) + cui(1,j,k,l)
      cue(2,j,k,l) = cue(2,j,k,l) + cui(2,j,k,l)
      cue(3,j,k,l) = cue(3,j,k,l) + cui(3,j,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPADDAMUI32(amu,amui,nyzp,nx,nxe,nypmx,nzpmx,idds)
! adds electron and ion momentum flux densities
! assumes guard cells have already been added
! amu/amui = momentum flux density for electrons/ions
! nyzp(1) = number of primary gridpoints in y in particle partition m
! nyzp(2) = number of primary gridpoints in z in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real amu, amui
      integer nyzp
      dimension amu(6,nxe,nypmx,nzpmx), amui(6,nxe,nypmx,nzpmx)
      dimension nyzp(idds)
      integer j, k, l
!$OMP PARALLEL DO PRIVATE(j,k,l)
      do 30 l = 1, nyzp(2)
      do 20 k = 1, nyzp(1)
      do 10 j = 1, nx
      amu(1,j,k,l) = amu(1,j,k,l) + amui(1,j,k,l)
      amu(2,j,k,l) = amu(2,j,k,l) + amui(2,j,k,l)
      amu(3,j,k,l) = amu(3,j,k,l) + amui(3,j,k,l)
      amu(4,j,k,l) = amu(4,j,k,l) + amui(4,j,k,l)
      amu(5,j,k,l) = amu(5,j,k,l) + amui(5,j,k,l)
      amu(6,j,k,l) = amu(6,j,k,l) + amui(6,j,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPBADDEXT32(bxyz,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx, &
     &idds)
! adds constant to magnetic field for 3d code
! for distributed data with 2D decomposition and OpenMP
! bxyz = magnetic field
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      real omx, omy, omz
      integer nyzp
      real bxyz
      dimension bxyz(3,nxe,nypmx,nzpmx)
      dimension nyzp(idds)
! local data
      integer j, k, l, ll
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
      do 20 ll = 1, nyzp(1)*nyzp(2)
      l = (ll - 1)/nyzp(1)
      k = ll - nyzp(1)*l
      l = l + 1
      do 10 j = 1, nx
      bxyz(1,j,k,l) = bxyz(1,j,k,l) + omx
      bxyz(2,j,k,l) = bxyz(2,j,k,l) + omy
      bxyz(3,j,k,l) = bxyz(3,j,k,l) + omz
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADDVRFIELD32(a,b,c,ndim,nxe,nypmx,nzpmx)
! this subroutine calculates a = b + c for distributed real vector field
! with 2D decomposition and OpenMP
      implicit none
      integer ndim, nxe, nypmx, nzpmx
      real a, b, c
      dimension a(ndim,nxe,nypmx,nzpmx)
      dimension b(ndim,nxe,nypmx,nzpmx), c(ndim,nxe,nypmx,nzpmx)
! local data
      integer i, j, k, l, ll
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
      do 30 ll = 1, nypmx*nzpmx
      l = (ll - 1)/nypmx
      k = ll - nypmx*l
      l = l + 1
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k,l) = b(i,j,k,l) + c(i,j,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPBBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,   &
     &nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! magnetic field (or convolution of magnetic field over particle shape)
! with periodic boundary conditions, for distributed data,
! with 2D spatial decomposition and OpenMP
! input: cu,ffc,ci,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
! output: bxyz,wm
! approximate flop count is:
! 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the magnetic field is calculated using the equations:
! bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz))*s(kx,ky,kz),
! by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz))*s(kx,ky,kz),
! bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
!                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz))*s(kx,ky,kz),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
! bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
! bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
! cu(l,j,k) = complex current density for fourier mode jj-1,kk-1,l-1
! bxyz(1,l,j,k) = x component of complex magnetic field
! bxyz(2,l,j,k) = y component of complex magnetic field
! bxyz(3,l,j,k) = z component of complex magnetic field
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! ci = reciprical of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real ci, wm
      complex cu, bxyz, ffc
      dimension cu(3,nzv,kxyp,kyzp), bxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real ci2, dnx, dny, dnz, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,at1,at2,at3,at4,zt1,zt2,zt3,wp) &
!$OMP& REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,k))*aimag(ffc(l,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(-aimag(cu(3,l,j,k)),real(cu(3,l,j,k)))
            zt2 = cmplx(-aimag(cu(2,l,j,k)),real(cu(2,l,j,k)))
            zt3 = cmplx(-aimag(cu(1,l,j,k)),real(cu(1,l,j,k)))
            bxyz(1,l,j,k) = at3*zt1 - at4*zt2
            bxyz(2,l,j,k) = at4*zt3 - at2*zt1
            bxyz(3,l,j,k) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(cu(3,l1,j,k)),real(cu(3,l1,j,k)))
            zt2 = cmplx(-aimag(cu(2,l1,j,k)),real(cu(2,l1,j,k)))
            zt3 = cmplx(-aimag(cu(1,l1,j,k)),real(cu(1,l1,j,k)))
            bxyz(1,l1,j,k) = at3*zt1 + at4*zt2
            bxyz(2,l1,j,k) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,k) = at2*zt2 - at3*zt3
            wp = wp + at1*(cu(1,l,j,k)*conjg(cu(1,l,j,k))               &
     &              + cu(2,l,j,k)*conjg(cu(2,l,j,k))                    &
     &              + cu(3,l,j,k)*conjg(cu(3,l,j,k))                    &
     &              + cu(1,l1,j,k)*conjg(cu(1,l1,j,k))                  &
     &              + cu(2,l1,j,k)*conjg(cu(2,l1,j,k))                  &
     &              + cu(3,l1,j,k)*conjg(cu(3,l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,k))*aimag(ffc(1,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(-aimag(cu(3,1,j,k)),real(cu(3,1,j,k)))
            zt2 = cmplx(-aimag(cu(2,1,j,k)),real(cu(2,1,j,k)))
            zt3 = cmplx(-aimag(cu(1,1,j,k)),real(cu(1,1,j,k)))
            bxyz(1,1,j,k) = at3*zt1
            bxyz(2,1,j,k) = -at2*zt1
            bxyz(3,1,j,k) = at2*zt2 - at3*zt3
            bxyz(1,l1,j,k) = zero
            bxyz(2,l1,j,k) = zero
            bxyz(3,l1,j,k) = zero
            wp = wp + at1*(cu(1,1,j,k)*conjg(cu(1,1,j,k))               &
     &              + cu(2,1,j,k)*conjg(cu(2,1,j,k))                    &
     &              + cu(3,1,j,k)*conjg(cu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,at1,at3,at4,zt1,zt2,zt3,wp)     &
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = ci2*real(ffc(l,1,k))*aimag(ffc(l,1,k))
               at3 = dky*at1
               at4 = dnz*real(l - 1)*at1
               zt1 = cmplx(-aimag(cu(3,l,1,k)),real(cu(3,l,1,k)))
               zt2 = cmplx(-aimag(cu(2,l,1,k)),real(cu(2,l,1,k)))
               zt3 = cmplx(-aimag(cu(1,l,1,k)),real(cu(1,l,1,k)))
               bxyz(1,l,1,k) = at3*zt1 - at4*zt2
               bxyz(2,l,1,k) = at4*zt3
               bxyz(3,l,1,k) = -at3*zt3
               zt1 = cmplx(-aimag(cu(3,l1,1,k)),real(cu(3,l1,1,k)))
               zt2 = cmplx(-aimag(cu(2,l1,1,k)),real(cu(2,l1,1,k)))
               zt3 = cmplx(-aimag(cu(1,l1,1,k)),real(cu(1,l1,1,k)))
               bxyz(1,l1,1,k) = at3*zt1 + at4*zt2
               bxyz(2,l1,1,k) = -at4*zt3
               bxyz(3,l1,1,k) = -at3*zt3
               wp = wp + at1*(cu(1,l,1,k)*conjg(cu(1,l,1,k))            &
     &                 + cu(2,l,1,k)*conjg(cu(2,l,1,k))                 &
     &                 + cu(3,l,1,k)*conjg(cu(3,l,1,k))                 &
     &                 + cu(1,l1,1,k)*conjg(cu(1,l1,1,k))               &
     &                 + cu(2,l1,1,k)*conjg(cu(2,l1,1,k))               &
     &                 + cu(3,l1,1,k)*conjg(cu(3,l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = ci2*real(ffc(1,1,k))*aimag(ffc(1,1,k))
               at3 = dky*at1
               zt1 = cmplx(-aimag(cu(3,1,1,k)),real(cu(3,1,1,k)))
               zt3 = cmplx(-aimag(cu(1,1,1,k)),real(cu(1,1,1,k)))
               bxyz(1,1,1,k) = at3*zt1
               bxyz(2,1,1,k) = zero
               bxyz(3,1,1,k) = -at3*zt3
               bxyz(1,l1,1,k) = zero
               bxyz(2,l1,1,k) = zero
               bxyz(3,l1,1,k) = zero
               wp = wp + at1*(cu(1,1,1,k)*conjg(cu(1,1,1,k))            &
     &                 + cu(2,1,1,k)*conjg(cu(2,1,1,k))                 &
     &                 + cu(3,1,1,k)*conjg(cu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k) = zero
               bxyz(2,l,1,k) = zero
               bxyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,at1,at2,at4,zt1,zt2,zt3,wp)        &
!$OMP& REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,1))*aimag(ffc(l,j,1))
            at2 = dkx*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(-aimag(cu(3,l,j,1)),real(cu(3,l,j,1)))
            zt2 = cmplx(-aimag(cu(2,l,j,1)),real(cu(2,l,j,1)))
            zt3 = cmplx(-aimag(cu(1,l,j,1)),real(cu(1,l,j,1)))
            bxyz(1,l,j,1) = -at4*zt2
            bxyz(2,l,j,1) = at4*zt3 - at2*zt1
            bxyz(3,l,j,1) = at2*zt2
            zt1 = cmplx(-aimag(cu(3,l1,j,1)),real(cu(3,l1,j,1)))
            zt2 = cmplx(-aimag(cu(2,l1,j,1)),real(cu(2,l1,j,1)))
            zt3 = cmplx(-aimag(cu(1,l1,j,1)),real(cu(1,l1,j,1)))
            bxyz(1,l1,j,1) = at4*zt2
            bxyz(2,l1,j,1) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,1) = at2*zt2
            wp = wp + at1*(cu(1,l,j,1)*conjg(cu(1,l,j,1))               &
     &              + cu(2,l,j,1)*conjg(cu(2,l,j,1))                    &
     &              + cu(3,l,j,1)*conjg(cu(3,l,j,1))                    &
     &              + cu(1,l1,j,1)*conjg(cu(1,l1,j,1))                  &
     &              + cu(2,l1,j,1)*conjg(cu(2,l1,j,1))                  &
     &              + cu(3,l1,j,1)*conjg(cu(3,l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,1))*aimag(ffc(1,j,1))
            at2 = dkx*at1
            zt1 = cmplx(-aimag(cu(3,1,j,1)),real(cu(3,1,j,1)))
            zt2 = cmplx(-aimag(cu(2,1,j,1)),real(cu(2,1,j,1)))
            bxyz(1,1,j,1) = zero
            bxyz(2,1,j,1) = -at2*zt1
            bxyz(3,1,j,1) = at2*zt2
            bxyz(1,l1,j,1) = zero
            bxyz(2,l1,j,1) = zero
            bxyz(3,l1,j,1) = zero
            wp = wp + at1*(cu(1,1,j,1)*conjg(cu(1,1,j,1))               &
     &              + cu(2,1,j,1)*conjg(cu(2,1,j,1))                    &
     &              + cu(3,1,j,1)*conjg(cu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,1,1))*aimag(ffc(l,1,1))
            at4 = dnz*real(l - 1)*at1
            zt2 = cmplx(-aimag(cu(2,l,1,1)),real(cu(2,l,1,1)))
            zt3 = cmplx(-aimag(cu(1,l,1,1)),real(cu(1,l,1,1)))
            bxyz(1,l,1,1) = -at4*zt2
            bxyz(2,l,1,1) = at4*zt3
            bxyz(3,l,1,1) = zero
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
            wp = wp + at1*(cu(1,l,1,1)*conjg(cu(1,l,1,1))               &
     &              + cu(2,l,1,1)*conjg(cu(2,l,1,1))                    &
     &              + cu(3,l,1,1)*conjg(cu(3,l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1) = zero
            bxyz(2,1,1,1) = zero
            bxyz(3,1,1,1) = zero
            bxyz(1,l1,1,1) = zero
            bxyz(2,l1,1,1) = zero
            bxyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            bxyz(1,l,j,k1) = zero
            bxyz(2,l,j,k1) = zero
            bxyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            bxyz(1,l,1,k1) = zero
            bxyz(2,l,1,k1) = zero
            bxyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      wm = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,&
     &kyzp)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux
! in 3d with periodic boundary conditions.
! for distributed data, with 2D spatial decomposition and OpenMP
! input: all except dcu, output: dcu
! approximate flop count is:
! 220*nxc*nyc*nzc + 70*(nxc*nyc + nxc*nzc + nyc*nzc)
! and (nx/2)*nyc*nzc divides
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky,kz) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
! dcu(2,kx,ky,kz) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
! dcu(3,kx,ky,kz) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
! dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
! dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
! dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
! the transverse part is calculated using the equation:
! dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
!                 - kx*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
!                 - ky*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
!                 - kz*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! on output:
! dcu(i,l,j,k) = transverse part of complex derivative of current for
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! amu(1,l,j,k) = xx component of complex momentum flux
! amu(2,l,j,k) = xy component of complex momentum flux
! amu(3,l,j,k) = xz component of complex momentum flux
! amu(4,l,j,k) = yy component of complex momentum flux
! amu(5,l,j,k) = yz component of complex momentum flux
! amu(6,l,j,k) = zz component of complex momentum flux
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of procs in y/z
! nzv = second dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex dcu, amu
      dimension dcu(3,nzv,kxyp,kyzp), amu(6,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, dkx2, dky2, dkxy2, at1
      complex zero, zt1, zt2, zt3, zt4, zt5
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate transverse part of current
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,zt1,zt2,zt3, &
!$OMP& zt4,zt5)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            zt1 = cmplx(aimag(amu(1,l,j,k)),-real(amu(1,l,j,k)))
            zt2 = cmplx(aimag(amu(2,l,j,k)),-real(amu(2,l,j,k)))
            zt3 = cmplx(aimag(amu(3,l,j,k)),-real(amu(3,l,j,k)))
            zt1 = dkx*zt1 + dky*zt2 + dkz*zt3
            zt4 = cmplx(aimag(amu(4,l,j,k)),-real(amu(4,l,j,k)))
            zt5 = cmplx(aimag(amu(5,l,j,k)),-real(amu(5,l,j,k)))
            zt2 = dkx*zt2 + dky*zt4 + dkz*zt5
            zt4 = cmplx(aimag(amu(6,l,j,k)),-real(amu(6,l,j,k)))
            zt3 = dkx*zt3 + dky*zt5 + dkz*zt4
            zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3)
            dcu(1,l,j,k) = zt1 - dkx*zt4
            dcu(2,l,j,k) = zt2 - dky*zt4
            dcu(3,l,j,k) = zt3 - dkz*zt4
            zt1 = cmplx(aimag(amu(1,l1,j,k)),-real(amu(1,l1,j,k)))
            zt2 = cmplx(aimag(amu(2,l1,j,k)),-real(amu(2,l1,j,k)))
            zt3 = cmplx(aimag(amu(3,l1,j,k)),-real(amu(3,l1,j,k)))
            zt1 = dkx*zt1 + dky*zt2 - dkz*zt3
            zt4 = cmplx(aimag(amu(4,l1,j,k)),-real(amu(4,l1,j,k)))
            zt5 = cmplx(aimag(amu(5,l1,j,k)),-real(amu(5,l1,j,k)))
            zt2 = dkx*zt2 + dky*zt4 - dkz*zt5
            zt4 = cmplx(aimag(amu(6,l1,j,k)),-real(amu(6,l1,j,k)))
            zt3 = dkx*zt3 + dky*zt5 - dkz*zt4
            zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3)
            dcu(1,l1,j,k) = zt1 - dkx*zt4
            dcu(2,l1,j,k) = zt2 - dky*zt4
            dcu(3,l1,j,k) = zt3 + dkz*zt4
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            zt1 = cmplx(aimag(amu(1,1,j,k)),-real(amu(1,1,j,k)))
            zt2 = cmplx(aimag(amu(2,1,j,k)),-real(amu(2,1,j,k)))
            zt3 = cmplx(aimag(amu(3,1,j,k)),-real(amu(3,1,j,k)))
            zt1 = dkx*zt1 + dky*zt2
            zt4 = cmplx(aimag(amu(4,1,j,k)),-real(amu(4,1,j,k)))
            zt5 = cmplx(aimag(amu(5,1,j,k)),-real(amu(5,1,j,k)))
            zt2 = dkx*zt2 + dky*zt4
            zt3 = dkx*zt3 + dky*zt5
            zt4 = at1*(dkx*zt1 + dky*zt2)
            dcu(1,1,j,k) = zt1 - dkx*zt4
            dcu(2,1,j,k) = zt2 - dky*zt4
            dcu(3,1,j,k) = zt3
            dcu(1,l1,j,k) = zero
            dcu(2,l1,j,k) = zero
            dcu(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,dky2,at1,zt1,zt2,zt3,zt4,zt5&
!$OMP& )
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               zt2 = cmplx(aimag(amu(2,l,1,k)),-real(amu(2,l,1,k)))
               zt3 = cmplx(aimag(amu(3,l,1,k)),-real(amu(3,l,1,k)))
               zt1 = dky*zt2 + dkz*zt3
               zt4 = cmplx(aimag(amu(4,l,1,k)),-real(amu(4,l,1,k)))
               zt5 = cmplx(aimag(amu(5,l,1,k)),-real(amu(5,l,1,k)))
               zt2 = dky*zt4 + dkz*zt5
               zt4 = cmplx(aimag(amu(6,l,1,k)),-real(amu(6,l,1,k)))
               zt3 = dky*zt5 + dkz*zt4
               zt4 = at1*(dky*zt2 + dkz*zt3)
               dcu(1,l,1,k) = zt1
               dcu(2,l,1,k) = zt2 - dky*zt4
               dcu(3,l,1,k) = zt3 - dkz*zt4
               zt2 = cmplx(aimag(amu(2,l1,1,k)),-real(amu(2,l1,1,k)))
               zt3 = cmplx(aimag(amu(3,l1,1,k)),-real(amu(3,l1,1,k)))
               zt1 = dky*zt2 - dkz*zt3
               zt4 = cmplx(aimag(amu(4,l1,1,k)),-real(amu(4,l1,1,k)))
               zt5 = cmplx(aimag(amu(5,l1,1,k)),-real(amu(5,l1,1,k)))
               zt2 = dky*zt4 - dkz*zt5
               zt4 = cmplx(aimag(amu(6,l1,1,k)),-real(amu(6,l1,1,k)))
               zt3 = dky*zt5 - dkz*zt4
               zt4 = at1*(dky*zt2 - dkz*zt3)
               dcu(1,l1,1,k) = zt1
               dcu(2,l1,1,k) = zt2 - dky*zt4
               dcu(3,l1,1,k) = zt3 + dkz*zt4
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt2 = cmplx(aimag(amu(2,1,1,k)),-real(amu(2,1,1,k)))
               zt1 = dky*zt2
               zt5 = cmplx(aimag(amu(5,1,1,k)),-real(amu(5,1,1,k)))
               zt3 = dky*zt5
               dcu(1,1,1,k) = zt1
               dcu(2,1,1,k) = zero
               dcu(3,1,1,k) = zt3
               dcu(1,l1,1,k) = zero
               dcu(2,l1,1,k) = zero
               dcu(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               dcu(1,l,1,k) = zero
               dcu(2,l,1,k) = zero
               dcu(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,dkx2,at1,zt1,zt2,zt3,zt4,zt5)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            zt1 = cmplx(aimag(amu(1,l,j,1)),-real(amu(1,l,j,1)))
            zt2 = cmplx(aimag(amu(2,l,j,1)),-real(amu(2,l,j,1)))
            zt3 = cmplx(aimag(amu(3,l,j,1)),-real(amu(3,l,j,1)))
            zt1 = dkx*zt1 + dkz*zt3
            zt5 = cmplx(aimag(amu(5,l,j,1)),-real(amu(5,l,j,1)))
            zt2 = dkx*zt2 + dkz*zt5
            zt4 = cmplx(aimag(amu(6,l,j,1)),-real(amu(6,l,j,1)))
            zt3 = dkx*zt3 + dkz*zt4
            zt4 = at1*(dkx*zt1 + dkz*zt3)
            dcu(1,l,j,1) = zt1 - dkx*zt4
            dcu(2,l,j,1) = zt2
            dcu(3,l,j,1) = zt3 - dkz*zt4
            zt1 = cmplx(aimag(amu(1,l1,j,1)),-real(amu(1,l1,j,1)))
            zt2 = cmplx(aimag(amu(2,l1,j,1)),-real(amu(2,l1,j,1)))
            zt3 = cmplx(aimag(amu(3,l1,j,1)),-real(amu(3,l1,j,1)))
            zt1 = dkx*zt1 - dkz*zt3
            zt5 = cmplx(aimag(amu(5,l1,j,1)),-real(amu(5,l1,j,1)))
            zt2 = dkx*zt2 - dkz*zt5
            zt4 = cmplx(aimag(amu(6,l1,j,1)),-real(amu(6,l1,j,1)))
            zt3 = dkx*zt3 - dkz*zt4
            zt4 = at1*(dkx*zt1 - dkz*zt3)
            dcu(1,l1,j,1) = zt1 - dkx*zt4
            dcu(2,l1,j,1) = zt2
            dcu(3,l1,j,1) = zt3 + dkz*zt4
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt2 = cmplx(aimag(amu(2,1,j,1)),-real(amu(2,1,j,1)))
            zt3 = cmplx(aimag(amu(3,1,j,1)),-real(amu(3,1,j,1)))
            zt2 = dkx*zt2
            zt3 = dkx*zt3
            dcu(1,1,j,1) = zero
            dcu(2,1,j,1) = zt2
            dcu(3,1,j,1) = zt3
            dcu(1,l1,j,1) = zero
            dcu(2,l1,j,1) = zero
            dcu(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt3 = cmplx(aimag(amu(3,l,1,1)),-real(amu(3,l,1,1)))
            zt1 = dkz*zt3
            zt5 = cmplx(aimag(amu(5,l,1,1)),-real(amu(5,l,1,1)))
            zt2 = dkz*zt5
            dcu(1,l,1,1) = zt1
            dcu(2,l,1,1) = zt2
            dcu(3,l,1,1) = zero
            dcu(1,l1,1,1) = zero
            dcu(2,l1,1,1) = zero
            dcu(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            dcu(1,1,1,1) = zero
            dcu(2,1,1,1) = zero
            dcu(3,1,1,1) = zero
            dcu(1,l1,1,1) = zero
            dcu(2,l1,1,1) = zero
            dcu(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            dcu(1,l,j,k1) = zero
            dcu(2,l,j,k1) = zero
            dcu(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            dcu(1,l,1,k1) = zero
            dcu(2,l,1,k1) = zero
            dcu(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADCUPERP32(dcu,amu,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp&
     &,kyzp)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux and acceleration density
! in 3d with periodic boundary conditions.
! for distributed data, with 2D spatial decomposition and OpenMP
! input: all output: dcu
! approximate flop count is:
! 244*nxc*nyc*nzc + 82*(nxc*nyc + nxc*nzc + nyc*nzc)
! and (nx/2)*nyc*nzc divides
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
!                -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
! dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
!                -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
! dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
!                 -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
! dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
! dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
! dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
! the transverse part is calculated using the equation:
! dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
!                 - kx*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
!                 - ky*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
!                 - kz*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
!                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! on input:
! dcu(i,l,j,k) = complex acceleration density for fourier mode (j-1,k-1)
! on output:
! dcu(i,l,j,k) = transverse part of complex derivative of current for
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! amu(1,l,j,k) = xx component of complex momentum flux
! amu(2,l,j,k) = xy component of complex momentum flux
! amu(3,l,j,k) = xz component of complex momentum flux
! amu(4,l,j,k) = yy component of complex momentum flux
! amu(5,l,j,k) = yz component of complex momentum flux
! amu(6,l,j,k) = zz component of complex momentum flux
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of procs in y/z
! nzv = second dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex dcu, amu
      dimension dcu(3,nzv,kxyp,kyzp), amu(6,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, dkx2, dky2, dkxy2, at1
      complex zero, zt1, zt2, zt3, zt4, zt5
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate transverse part of current
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,zt1,zt2,zt3, &
!$OMP& zt4,zt5)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            zt1 = cmplx(aimag(amu(1,l,j,k)),-real(amu(1,l,j,k)))
            zt2 = cmplx(aimag(amu(2,l,j,k)),-real(amu(2,l,j,k)))
            zt3 = cmplx(aimag(amu(3,l,j,k)),-real(amu(3,l,j,k)))
            zt1 = dcu(1,l,j,k) + dkx*zt1 + dky*zt2 + dkz*zt3
            zt4 = cmplx(aimag(amu(4,l,j,k)),-real(amu(4,l,j,k)))
            zt5 = cmplx(aimag(amu(5,l,j,k)),-real(amu(5,l,j,k)))
            zt2 = dcu(2,l,j,k) + dkx*zt2 + dky*zt4 + dkz*zt5
            zt4 = cmplx(aimag(amu(6,l,j,k)),-real(amu(6,l,j,k)))
            zt3 = dcu(3,l,j,k) + dkx*zt3 + dky*zt5 + dkz*zt4
            zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3)
            dcu(1,l,j,k) = zt1 - dkx*zt4
            dcu(2,l,j,k) = zt2 - dky*zt4
            dcu(3,l,j,k) = zt3 - dkz*zt4
            zt1 = cmplx(aimag(amu(1,l1,j,k)),-real(amu(1,l1,j,k)))
            zt2 = cmplx(aimag(amu(2,l1,j,k)),-real(amu(2,l1,j,k)))
            zt3 = cmplx(aimag(amu(3,l1,j,k)),-real(amu(3,l1,j,k)))
            zt1 = dcu(1,l1,j,k) + dkx*zt1 + dky*zt2 - dkz*zt3
            zt4 = cmplx(aimag(amu(4,l1,j,k)),-real(amu(4,l1,j,k)))
            zt5 = cmplx(aimag(amu(5,l1,j,k)),-real(amu(5,l1,j,k)))
            zt2 = dcu(2,l1,j,k) + dkx*zt2 + dky*zt4 - dkz*zt5
            zt4 = cmplx(aimag(amu(6,l1,j,k)),-real(amu(6,l1,j,k)))
            zt3 = dcu(3,l1,j,k) + dkx*zt3 + dky*zt5 - dkz*zt4
            zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3)
            dcu(1,l1,j,k) = zt1 - dkx*zt4
            dcu(2,l1,j,k) = zt2 - dky*zt4
            dcu(3,l1,j,k) = zt3 + dkz*zt4
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            zt1 = cmplx(aimag(amu(1,1,j,k)),-real(amu(1,1,j,k)))
            zt2 = cmplx(aimag(amu(2,1,j,k)),-real(amu(2,1,j,k)))
            zt3 = cmplx(aimag(amu(3,1,j,k)),-real(amu(3,1,j,k)))
            zt1 = dcu(1,1,j,k) + dkx*zt1 + dky*zt2
            zt4 = cmplx(aimag(amu(4,1,j,k)),-real(amu(4,1,j,k)))
            zt5 = cmplx(aimag(amu(5,1,j,k)),-real(amu(5,1,j,k)))
            zt2 = dcu(2,1,j,k) + dkx*zt2 + dky*zt4
            zt3 = dcu(3,1,j,k) + dkx*zt3 + dky*zt5
            zt4 = at1*(dkx*zt1 + dky*zt2)
            dcu(1,1,j,k) = zt1 - dkx*zt4
            dcu(2,1,j,k) = zt2 - dky*zt4
            dcu(3,1,j,k) = zt3
            dcu(1,l1,j,k) = zero
            dcu(2,l1,j,k) = zero
            dcu(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,dky2,at1,zt1,zt2,zt3,zt4,zt5&
!$OMP& )
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               zt2 = cmplx(aimag(amu(2,l,1,k)),-real(amu(2,l,1,k)))
               zt3 = cmplx(aimag(amu(3,l,1,k)),-real(amu(3,l,1,k)))
               zt1 = dcu(1,l,1,k) + dky*zt2 + dkz*zt3
               zt4 = cmplx(aimag(amu(4,l,1,k)),-real(amu(4,l,1,k)))
               zt5 = cmplx(aimag(amu(5,l,1,k)),-real(amu(5,l,1,k)))
               zt2 = dcu(2,l,1,k) + dky*zt4 + dkz*zt5
               zt4 = cmplx(aimag(amu(6,l,1,k)),-real(amu(6,l,1,k)))
               zt3 = dcu(3,l,1,k) + dky*zt5 + dkz*zt4
               zt4 = at1*(dky*zt2 + dkz*zt3)
               dcu(1,l,1,k) = zt1
               dcu(2,l,1,k) = zt2 - dky*zt4
               dcu(3,l,1,k) = zt3 - dkz*zt4
               zt2 = cmplx(aimag(amu(2,l1,1,k)),-real(amu(2,l1,1,k)))
               zt3 = cmplx(aimag(amu(3,l1,1,k)),-real(amu(3,l1,1,k)))
               zt1 = dcu(1,l1,1,k) + dky*zt2 - dkz*zt3
               zt4 = cmplx(aimag(amu(4,l1,1,k)),-real(amu(4,l1,1,k)))
               zt5 = cmplx(aimag(amu(5,l1,1,k)),-real(amu(5,l1,1,k)))
               zt2 = dcu(2,l1,1,k) + dky*zt4 - dkz*zt5
               zt4 = cmplx(aimag(amu(6,l1,1,k)),-real(amu(6,l1,1,k)))
               zt3 = dcu(3,l1,1,k) + dky*zt5 - dkz*zt4
               zt4 = at1*(dky*zt2 - dkz*zt3)
               dcu(1,l1,1,k) = zt1
               dcu(2,l1,1,k) = zt2 - dky*zt4
               dcu(3,l1,1,k) = zt3 + dkz*zt4
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt2 = cmplx(aimag(amu(2,1,1,k)),-real(amu(2,1,1,k)))
               zt1 = dcu(1,1,1,k) + dky*zt2
               zt5 = cmplx(aimag(amu(5,1,1,k)),-real(amu(5,1,1,k)))
               zt3 = dcu(3,1,1,k) + dky*zt5
               dcu(1,1,1,k) = zt1
               dcu(2,1,1,k) = zero
               dcu(3,1,1,k) = zt3
               dcu(1,l1,1,k) = zero
               dcu(2,l1,1,k) = zero
               dcu(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               dcu(1,l,1,k) = zero
               dcu(2,l,1,k) = zero
               dcu(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,dkx2,at1,zt1,zt2,zt3,zt4,zt5)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            zt1 = cmplx(aimag(amu(1,l,j,1)),-real(amu(1,l,j,1)))
            zt2 = cmplx(aimag(amu(2,l,j,1)),-real(amu(2,l,j,1)))
            zt3 = cmplx(aimag(amu(3,l,j,1)),-real(amu(3,l,j,1)))
            zt1 = dcu(1,l,j,1) + dkx*zt1 + dkz*zt3
            zt5 = cmplx(aimag(amu(5,l,j,1)),-real(amu(5,l,j,1)))
            zt2 = dcu(2,l,j,1) + dkx*zt2 + dkz*zt5
            zt4 = cmplx(aimag(amu(6,l,j,1)),-real(amu(6,l,j,1)))
            zt3 = dcu(3,l,j,1) + dkx*zt3 + dkz*zt4
            zt4 = at1*(dkx*zt1 + dkz*zt3)
            dcu(1,l,j,1) = zt1 - dkx*zt4
            dcu(2,l,j,1) = zt2
            dcu(3,l,j,1) = zt3 - dkz*zt4
            zt1 = cmplx(aimag(amu(1,l1,j,1)),-real(amu(1,l1,j,1)))
            zt2 = cmplx(aimag(amu(2,l1,j,1)),-real(amu(2,l1,j,1)))
            zt3 = cmplx(aimag(amu(3,l1,j,1)),-real(amu(3,l1,j,1)))
            zt1 = dcu(1,l1,j,1) + dkx*zt1 - dkz*zt3
            zt5 = cmplx(aimag(amu(5,l1,j,1)),-real(amu(5,l1,j,1)))
            zt2 = dcu(2,l1,j,1) + dkx*zt2 - dkz*zt5
            zt4 = cmplx(aimag(amu(6,l1,j,1)),-real(amu(6,l1,j,1)))
            zt3 = dcu(3,l1,j,1) + dkx*zt3 - dkz*zt4
            zt4 = at1*(dkx*zt1 - dkz*zt3)
            dcu(1,l1,j,1) = zt1 - dkx*zt4
            dcu(2,l1,j,1) = zt2
            dcu(3,l1,j,1) = zt3 + dkz*zt4
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt2 = cmplx(aimag(amu(2,1,j,1)),-real(amu(2,1,j,1)))
            zt3 = cmplx(aimag(amu(3,1,j,1)),-real(amu(3,1,j,1)))
            zt2 = dcu(2,1,j,1) + dkx*zt2
            zt3 = dcu(3,1,j,1) + dkx*zt3
            dcu(1,1,j,1) = zero
            dcu(2,1,j,1) = zt2
            dcu(3,1,j,1) = zt3
            dcu(1,l1,j,1) = zero
            dcu(2,l1,j,1) = zero
            dcu(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt3 = cmplx(aimag(amu(3,l,1,1)),-real(amu(3,l,1,1)))
            zt1 = dcu(1,l,1,1) + dkz*zt3
            zt5 = cmplx(aimag(amu(5,l,1,1)),-real(amu(5,l,1,1)))
            zt2 = dcu(2,l,1,1) + dkz*zt5
            dcu(1,l,1,1) = zt1
            dcu(2,l,1,1) = zt2
            dcu(3,l,1,1) = zero
            dcu(1,l1,1,1) = zero
            dcu(2,l1,1,1) = zero
            dcu(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            dcu(1,1,1,1) = zero
            dcu(2,1,1,1) = zero
            dcu(3,1,1,1) = zero
            dcu(1,l1,1,1) = zero
            dcu(2,l1,1,1) = zero
            dcu(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            dcu(1,l,j,k1) = zero
            dcu(2,l,j,k1) = zero
            dcu(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            dcu(1,l,1,k1) = zero
            dcu(2,l,1,k1) = zero
            dcu(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEPOISP332(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf&
     &,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! transverse electric field (or convolution of transverse electric field
! over particle shape), with periodic boundary conditions.
! for distributed data, with 2D spatial decomposition and OpenMP
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! for isign = 0, output: ffe
! input: isign,ax,ay,az,affp,wp0,nx,ny,nz,kstrt,nvpy,nvpz,kxyp,kyzp,nzhd
! for isign /= 0, output: exyz,wf
! input: dcu,ffe,isign,affp,ci,nx,ny,nz,kstrt,nzv,kxyp,kyzp,nzhd
! approximate flop count is:
! 128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! if isign = 0, form factor array is prepared
! if isign = -1, smoothed transverse electric field is calculated
! using the equation:
! ex(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcux(kx,ky,kz)*s(kx,ky,kz)
! ey(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuy(kx,ky,kz)*s(kx,ky,kz)
! ez(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuz(kx,ky,kz)*s(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
! ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
! ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
! if isign = 1, unsmoothed transverse electric field is calculated
! using the equation:
! ex(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcux(kx,ky,kz)
! ey(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuy(kx,ky,kz)
! ez(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuz(kx,ky,kz)
! dcu(l,j,k) = transverse part of complex derivative of current for
! fourier mode jj-1,kk-1,l-1
! exyz(1,l,j,k) = x component of complex transverse electric field
! exyz(2,l,j,k) = y component of complex transverse electric field
! exyz(3,l,j,k) = z component of complex transverse electric field
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffe(l,j,k)) = finite-size particle shape factor s
! real(ffe(l,j,k)) = potential green's function g
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! ax/ay/az = half-width of particle in x/y/z direction
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! wp0 = normalized total plasma frequency squared
! ci = reciprical of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*nz*sum((affp/((kx**2+ky**2+kz**2)*ci*ci)**2)
!    |dcu(kx,ky,kz)*s(kx,ky,kz)|**2)
! this expression is valid only if the derivative of current is
! divergence-free
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer isign, nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      integer nzhd
      real ax, ay, az, affp, wp0, ci, wf
      complex dcu, exyz, ffe
      dimension dcu(3,nzv,kxyp,kyzp), exyz(3,nzv,kxyp,kyzp)
      dimension ffe(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, ci2, wpc
      real at1, at2, at3, at4, at5, at6
      complex zero
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.(nvpy*nvpz)) return
      wpc = wp0*ci2
! prepare form factor array
      do 30 k = 1, kyzps
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*real(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyps
      dkx = dnx*real(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*real(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.0) then
         ffe(l,j,k) = cmplx(affp,1.0)
      else
         ffe(l,j,k) = cmplx(affp*at6/(at5 + wpc*at6*at6),at6)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 170
! calculate smoothed transverse electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 160
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 60 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 50 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,k))
            at1 = at2*aimag(ffe(l,j,k))
            at2 = at2*at2
            exyz(1,l,j,k) = at1*dcu(1,l,j,k)
            exyz(2,l,j,k) = at1*dcu(2,l,j,k)
            exyz(3,l,j,k) = at1*dcu(3,l,j,k)
            exyz(1,l1,j,k) = at1*dcu(1,l1,j,k)
            exyz(2,l1,j,k) = at1*dcu(2,l1,j,k)
            exyz(3,l1,j,k) = at1*dcu(3,l1,j,k)
            wp = wp + at2*(dcu(1,l,j,k)*conjg(dcu(1,l,j,k))             &
     &              + dcu(2,l,j,k)*conjg(dcu(2,l,j,k))                  &
     &              + dcu(3,l,j,k)*conjg(dcu(3,l,j,k))                  &
     &              + dcu(1,l1,j,k)*conjg(dcu(1,l1,j,k))                &
     &              + dcu(2,l1,j,k)*conjg(dcu(2,l1,j,k))                &
     &              + dcu(3,l1,j,k)*conjg(dcu(3,l1,j,k)))
   50       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,k))
            at1 = at2*aimag(ffe(1,j,k))
            at2 = at2*at2
            exyz(1,1,j,k) = at1*dcu(1,1,j,k)
            exyz(2,1,j,k) = at1*dcu(2,1,j,k)
            exyz(3,1,j,k) = at1*dcu(3,1,j,k)
            exyz(1,l1,j,k) = zero
            exyz(2,l1,j,k) = zero
            exyz(3,l1,j,k) = zero
            wp = wp + at2*(dcu(1,1,j,k)*conjg(dcu(1,1,j,k))             &
     &              + dcu(2,1,j,k)*conjg(dcu(2,1,j,k))                  &
     &              + dcu(3,1,j,k)*conjg(dcu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   60 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 90 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 70 l = 2, nzh
               l1 = nz2 - l
               at2 = -ci2*real(ffe(l,1,k))
               at1 = at2*aimag(ffe(l,1,k))
               at2 = at2*at2
               exyz(1,l,1,k) = at1*dcu(1,l,1,k)
               exyz(2,l,1,k) = at1*dcu(2,l,1,k)
               exyz(3,l,1,k) = at1*dcu(3,l,1,k)
               exyz(1,l1,1,k) = at1*dcu(1,l1,1,k)
               exyz(2,l1,1,k) = at1*dcu(2,l1,1,k)
               exyz(3,l1,1,k) = at1*dcu(3,l1,1,k)
               wp = wp + at2*(dcu(1,l,1,k)*conjg(dcu(1,l,1,k))          &
     &                 + dcu(2,l,1,k)*conjg(dcu(2,l,1,k))               &
     &                 + dcu(3,l,1,k)*conjg(dcu(3,l,1,k))               &
     &                 + dcu(1,l1,1,k)*conjg(dcu(1,l1,1,k))             &
     &                 + dcu(2,l1,1,k)*conjg(dcu(2,l1,1,k))             &
     &                 + dcu(3,l1,1,k)*conjg(dcu(3,l1,1,k)))
   70          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = -ci2*real(ffe(1,1,k))
               at1 = at2*aimag(ffe(1,1,k))
               at2 = at2*at2
               exyz(1,1,1,k) = at1*dcu(1,1,1,k)
               exyz(2,1,1,k) = at1*dcu(2,1,1,k)
               exyz(3,1,1,k) = at1*dcu(3,1,1,k)
               exyz(1,l1,1,k) = zero
               exyz(2,l1,1,k) = zero
               exyz(3,l1,1,k) = zero
               wp = wp + at2*(dcu(1,1,1,k)*conjg(dcu(1,1,1,k))          &
     &                 + dcu(2,1,1,k)*conjg(dcu(2,1,1,k))               &
     &                 + dcu(3,1,1,k)*conjg(dcu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 80 l = 1, nz
               exyz(1,l,1,k) = zero
               exyz(2,l,1,k) = zero
               exyz(3,l,1,k) = zero
   80          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   90 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 110 j = 1, kxyps
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 100 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,1))
            at1 = at2*aimag(ffe(l,j,1))
            at2 = at2*at2
            exyz(1,l,j,1) = at1*dcu(1,l,j,1)
            exyz(2,l,j,1) = at1*dcu(2,l,j,1)
            exyz(3,l,j,1) = at1*dcu(3,l,j,1)
            exyz(1,l1,j,1) = at1*dcu(1,l1,j,1)
            exyz(2,l1,j,1) = at1*dcu(2,l1,j,1)
            exyz(3,l1,j,1) = at1*dcu(3,l1,j,1)
            wp = wp + at2*(dcu(1,l,j,1)*conjg(dcu(1,l,j,1))             &
     &              + dcu(2,l,j,1)*conjg(dcu(2,l,j,1))                  &
     &              + dcu(3,l,j,1)*conjg(dcu(3,l,j,1))                  &
     &              + dcu(1,l1,j,1)*conjg(dcu(1,l1,j,1))                &
     &              + dcu(2,l1,j,1)*conjg(dcu(2,l1,j,1))                &
     &              + dcu(3,l1,j,1)*conjg(dcu(3,l1,j,1)))
  100       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,1))
            at1 = at2*aimag(ffe(1,j,1))
            at2 = at2*at2
            exyz(1,1,j,1) = at1*dcu(1,1,j,1)
            exyz(2,1,j,1) = at1*dcu(2,1,j,1)
            exyz(3,1,j,1) = at1*dcu(3,1,j,1)
            exyz(1,l1,j,1) = zero
            exyz(2,l1,j,1) = zero
            exyz(3,l1,j,1) = zero
            wp = wp + at2*(dcu(1,1,j,1)*conjg(dcu(1,1,j,1))             &
     &              + dcu(2,1,j,1)*conjg(dcu(2,1,j,1))                  &
     &              + dcu(3,1,j,1)*conjg(dcu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
  110    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,1,1))
            at1 = at2*aimag(ffe(l,1,1))
            at2 = at2*at2
            exyz(1,l,1,1) = at1*dcu(1,l,1,1)
            exyz(2,l,1,1) = at1*dcu(2,l,1,1)
            exyz(3,l,1,1) = at1*dcu(3,l,1,1)
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
            wp = wp + at2*(dcu(1,l,1,1)*conjg(dcu(1,l,1,1))             &
     &              + dcu(2,l,1,1)*conjg(dcu(2,l,1,1))                  &
     &              + dcu(3,l,1,1)*conjg(dcu(3,l,1,1)))
  120       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            exyz(1,1,1,1) = zero
            exyz(2,1,1,1) = zero
            exyz(3,1,1,1) = zero
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 140 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 130 l = 1, nz
            exyz(1,l,j,k1) = zero
            exyz(2,l,j,k1) = zero
            exyz(3,l,j,k1) = zero
  130       continue
         endif
  140    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 150 l = 1, nz
            exyz(1,l,1,k1) = zero
            exyz(2,l,1,k1) = zero
            exyz(3,l,1,k1) = zero
  150       continue
         endif
      endif
      sum1 = sum1 + sum2
  160 continue
      wf = real(nx)*real(ny)*real(nz)*sum1/affp
      return
! calculate unsmoothed transverse electric field and sum field energy
  170 sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 330
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 190 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 180 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,k))
            at1 = at2*at2
            exyz(1,l,j,k) = at2*dcu(1,l,j,k)
            exyz(2,l,j,k) = at2*dcu(2,l,j,k)
            exyz(3,l,j,k) = at2*dcu(3,l,j,k)
            exyz(1,l1,j,k) = at2*dcu(1,l1,j,k)
            exyz(2,l1,j,k) = at2*dcu(2,l1,j,k)
            exyz(3,l1,j,k) = at2*dcu(3,l1,j,k)
            wp = wp + at1*(dcu(1,l,j,k)*conjg(dcu(1,l,j,k))             &
     &              + dcu(2,l,j,k)*conjg(dcu(2,l,j,k))                  &
     &              + dcu(3,l,j,k)*conjg(dcu(3,l,j,k))                  &
     &              + dcu(1,l1,j,k)*conjg(dcu(1,l1,j,k))                &
     &              + dcu(2,l1,j,k)*conjg(dcu(2,l1,j,k))                &
     &              + dcu(3,l1,j,k)*conjg(dcu(3,l1,j,k)))
  180       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,k))
            at1 = at2*at2
            exyz(1,1,j,k) = at2*dcu(1,1,j,k)
            exyz(2,1,j,k) = at2*dcu(2,1,j,k)
            exyz(3,1,j,k) = at2*dcu(3,1,j,k)
            exyz(1,l1,j,k) = zero
            exyz(2,l1,j,k) = zero
            exyz(3,l1,j,k) = zero
            wp = wp + at1*(dcu(1,1,j,k)*conjg(dcu(1,1,j,k))             &
     &              + dcu(2,1,j,k)*conjg(dcu(2,1,j,k))                  &
     &              + dcu(3,1,j,k)*conjg(dcu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
  190 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 220 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 200 l = 2, nzh
               l1 = nz2 - l
               at2 = -ci2*real(ffe(l,1,k))
               at1 = at2*at2
               exyz(1,l,1,k) = at2*dcu(1,l,1,k)
               exyz(2,l,1,k) = at2*dcu(2,l,1,k)
               exyz(3,l,1,k) = at2*dcu(3,l,1,k)
               exyz(1,l1,1,k) = at2*dcu(1,l1,1,k)
               exyz(2,l1,1,k) = at2*dcu(2,l1,1,k)
               exyz(3,l1,1,k) = at2*dcu(3,l1,1,k)
               wp = wp + at1*(dcu(1,l,1,k)*conjg(dcu(1,l,1,k))          &
     &                + dcu(2,l,1,k)*conjg(dcu(2,l,1,k))                &
     &                + dcu(3,l,1,k)*conjg(dcu(3,l,1,k))                &
     &                + dcu(1,l1,1,k)*conjg(dcu(1,l1,1,k))              &
     &                + dcu(2,l1,1,k)*conjg(dcu(2,l1,1,k))              &
     &                + dcu(3,l1,1,k)*conjg(dcu(3,l1,1,k)))
  200          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = -ci2*real(ffe(1,1,k))
               at1 = at2*at2
               exyz(1,1,1,k) = at2*dcu(1,1,1,k)
               exyz(2,1,1,k) = at2*dcu(2,1,1,k)
               exyz(3,1,1,k) = at2*dcu(3,1,1,k)
               exyz(1,l1,1,k) = zero
               exyz(2,l1,1,k) = zero
               exyz(3,l1,1,k) = zero
               wp = wp + at1*(dcu(1,1,1,k)*conjg(dcu(1,1,1,k))          &
     &                 + dcu(2,1,1,k)*conjg(dcu(2,1,1,k))               &
     &                 + dcu(3,1,1,k)*conjg(dcu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 210 l = 1, nz
               exyz(1,l,1,k) = zero
               exyz(2,l,1,k) = zero
               exyz(3,l,1,k) = zero
  210          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
  220 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 240 j = 1, kxyps
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 230 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,1))
            at1 = at2*at2
            exyz(1,l,j,1) = at2*dcu(1,l,j,1)
            exyz(2,l,j,1) = at2*dcu(2,l,j,1)
            exyz(3,l,j,1) = at2*dcu(3,l,j,1)
            exyz(1,l1,j,1) = at2*dcu(1,l1,j,1)
            exyz(2,l1,j,1) = at2*dcu(2,l1,j,1)
            exyz(3,l1,j,1) = at2*dcu(3,l1,j,1)
            wp = wp + at1*(dcu(1,l,j,1)*conjg(dcu(1,l,j,1))             &
     &              + dcu(2,l,j,1)*conjg(dcu(2,l,j,1))                  &
     &              + dcu(3,l,j,1)*conjg(dcu(3,l,j,1))                  &
     &              + dcu(1,l1,j,1)*conjg(dcu(1,l1,j,1))                &
     &              + dcu(2,l1,j,1)*conjg(dcu(2,l1,j,1))                &
     &              + dcu(3,l1,j,1)*conjg(dcu(3,l1,j,1)))
  230       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,1))
            at1 = at2*at2
            exyz(1,1,j,1) = at2*dcu(1,1,j,1)
            exyz(2,1,j,1) = at2*dcu(2,1,j,1)
            exyz(3,1,j,1) = at2*dcu(3,1,j,1)
            exyz(1,l1,j,1) = zero
            exyz(2,l1,j,1) = zero
            exyz(3,l1,j,1) = zero
            wp = wp + at1*(dcu(1,1,j,1)*conjg(dcu(1,1,j,1))             &
     &              + dcu(2,1,j,1)*conjg(dcu(2,1,j,1))                  &
     &              + dcu(3,1,j,1)*conjg(dcu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
  240    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 250 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,1,1))
            at1 = at2*at2
            exyz(1,l,1,1) = at2*dcu(1,l,1,1)
            exyz(2,l,1,1) = at2*dcu(2,l,1,1)
            exyz(3,l,1,1) = at2*dcu(3,l,1,1)
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
            wp = wp + at1*(dcu(1,l,1,1)*conjg(dcu(1,l,1,1))             &
     &              + dcu(2,l,1,1)*conjg(dcu(2,l,1,1))                  &
     &              + dcu(3,l,1,1)*conjg(dcu(3,l,1,1)))
  250       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            exyz(1,1,1,1) = zero
            exyz(2,1,1,1) = zero
            exyz(3,1,1,1) = zero
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 270 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 260 l = 1, nz
            exyz(1,l,j,k1) = zero
            exyz(2,l,j,k1) = zero
            exyz(3,l,j,k1) = zero
  260       continue
         endif
  270    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 320 l = 1, nz
            exyz(1,l,1,k1) = zero
            exyz(2,l,1,k1) = zero
            exyz(3,l,1,k1) = zero
  320       continue
         endif
      endif
      sum1 = sum1 + sum2
  330 continue
      wf = real(nx)*real(ny)*real(nz)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp&
     &,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! potential, with periodic boundary conditions for distributed data,
! with 2D spatial decomposition and OpenMP
! input: q,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
! output: pot,we
! approximate flop count is:
! 41*nxc*nyc*nzc + 21*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! potential is calculated using the equation:
! pot(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! pot(kx=pi) = 0, pot(ky=pi) = 0, pot(kz=pi) = 0, and
! pot(kx=0,ky=0,kz=0) = 0.
! q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
! pot(l,j,k) = complex potential
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! electric field energy is also calculated, using
! we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
!    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
! where affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real we
      complex q, pot, ffc
      dimension q(nzv,kxyp,kyzp), pot(nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real at1, at2
      complex zero
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,k))
            at1 = at2*aimag(ffc(l,j,k))
            pot(l,j,k) = at2*q(l,j,k)
            pot(l1,j,k) = at2*q(l1,j,k)
            wp = wp + at1*(q(l,j,k)*conjg(q(l,j,k))                     &
     &              + q(l1,j,k)*conjg(q(l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,k))
            at1 = at2*aimag(ffc(1,j,k))
            pot(1,j,k) = at2*q(1,j,k)
            pot(l1,j,k) = zero
            wp = wp + at1*(q(1,j,k)*conjg(q(1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at2 = real(ffc(l,1,k))
               at1 = at2*aimag(ffc(l,1,k))
               pot(l,1,k) = at2*q(l,1,k)
               pot(l1,1,k) = at2*q(l1,1,k)
               wp = wp + at1*(q(l,1,k)*conjg(q(l,1,k))                  &
     &                 + q(l1,1,k)*conjg(q(l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = real(ffc(1,1,k))
               at1 = at2*aimag(ffc(1,1,k))
               pot(1,1,k) = at2*q(1,1,k)
               pot(l1,1,k) = zero
               wp = wp + at1*(q(1,1,k)*conjg(q(1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               pot(l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,1))
            at1 = at2*aimag(ffc(l,j,1))
            pot(l,j,1) = at2*q(l,j,1)
            pot(l1,j,1) = at2*q(l1,j,1)
            wp = wp + at1*(q(l,j,1)*conjg(q(l,j,1))                     &
     &              + q(l1,j,1)*conjg(q(l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,1))
            at1 = at2*aimag(ffc(1,j,1))
            pot(1,j,1) = at2*q(1,j,1)
            pot(l1,j,1) = zero
            wp = wp + at1*(q(1,j,1)*conjg(q(1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,1,1))
            at1 = at2*aimag(ffc(l,1,1))
            pot(l,1,1) = at2*q(l,1,1)
            pot(l1,1,1) = zero
            wp = wp + at1*(q(l,1,1)*conjg(q(l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1) = zero
            pot(l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            pot(l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            pot(l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      we = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPELFIELD32(q,fxyz,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! unsmoothed longitudinal electric field, with periodic boundary
! conditions for distributed data
! input: q,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
! output: fxyz,we
! approximate flop count is:
! 62*nxc*nyc*nzc + 33*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the equation used is:
! fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz),
! fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz),
! fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
! fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
! fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
! fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
! q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
! fxyz(1,l,j,k) = x component of complex electric field,
! fxyz(2,l,j,k) = y component of complex electric field,
! fxyz(3,l,j,k) = z component of complex electric field,
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! electric field energy is also calculated, using
! we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
!    |q(kx,ky,kz)*s(kx,ky,kz)|**2),
! where affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real we
      complex q, fxyz, ffc
      dimension q(nzv,kxyp,kyzp), fxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate unsmoothed longitudinal electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,at1,at2,at3,at4,zt1,zt2,wp)     &
!$OMP& REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,k))
            zt1 = cmplx(aimag(q(l,j,k)),-real(q(l,j,k)))
            zt2 = cmplx(aimag(q(l1,j,k)),-real(q(l1,j,k)))
            fxyz(1,l,j,k) = at2*zt1
            fxyz(2,l,j,k) = at3*zt1
            fxyz(3,l,j,k) = at4*zt1
            fxyz(1,l1,j,k) = at2*zt2
            fxyz(2,l1,j,k) = at3*zt2
            fxyz(3,l1,j,k) = -at4*zt2
            wp = wp + at1*(q(l,j,k)*conjg(q(l,j,k))                     &
     &              + q(l1,j,k)*conjg(q(l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = at1*aimag(ffc(1,j,k))
            zt1 = cmplx(aimag(q(1,j,k)),-real(q(1,j,k)))
            fxyz(1,1,j,k) = at2*zt1
            fxyz(2,1,j,k) = at3*zt1
            fxyz(3,1,j,k) = zero
            fxyz(1,l1,j,k) = zero
            fxyz(2,l1,j,k) = zero
            fxyz(3,l1,j,k) = zero
            wp = wp + at1*(q(1,j,k)*conjg(q(1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,at1,at3,at4,zt1,zt2,wp)         &
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = real(ffc(l,1,k))
               at3 = dky*at1
               at4 = dnz*real(l - 1)*at1
               at1 = at1*aimag(ffc(l,1,k))
               zt1 = cmplx(aimag(q(l,1,k)),-real(q(l,1,k)))
               zt2 = cmplx(aimag(q(l1,1,k)),-real(q(l1,1,k)))
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = at3*zt1
               fxyz(3,l,1,k) = at4*zt1
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = at3*zt2
               fxyz(3,l1,1,k) = -at4*zt2
               wp = wp + at1*(q(l,1,k)*conjg(q(l,1,k))                  &
     &                 + q(l1,1,k)*conjg(q(l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = real(ffc(1,1,k))
               at3 = dky*at1
               at1 = at1*aimag(ffc(1,1,k))
               zt1 = cmplx(aimag(q(1,1,k)),-real(q(1,1,k)))
               fxyz(1,1,1,k) = zero
               fxyz(2,1,1,k) = at3*zt1
               fxyz(3,1,1,k) = zero
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = zero
               fxyz(3,l1,1,k) = zero
               wp = wp + at1*(q(1,1,k)*conjg(q(1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = zero
               fxyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,at1,at2,at4,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,1))
            at2 = dkx*at1
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,1))
            zt1 = cmplx(aimag(q(l,j,1)),-real(q(l,j,1)))
            zt2 = cmplx(aimag(q(l1,j,1)),-real(q(l1,j,1)))
            fxyz(1,l,j,1) = at2*zt1
            fxyz(2,l,j,1) = zero
            fxyz(3,l,j,1) = at4*zt1
            fxyz(1,l1,j,1) = at2*zt2
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = -at4*zt2
            wp = wp + at1*(q(l,j,1)*conjg(q(l,j,1))                     &
     &              + q(l1,j,1)*conjg(q(l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,1))
            at2 = dkx*at1
            at1 = at1*aimag(ffc(1,j,1))
            zt1 = cmplx(aimag(q(1,j,1)),-real(q(1,j,1)))
            fxyz(1,1,j,1) = at2*zt1
            fxyz(2,1,j,1) = zero
            fxyz(3,1,j,1) = zero
            fxyz(1,l1,j,1) = zero
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = zero
            wp = wp + at1*(q(1,j,1)*conjg(q(1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,1,1))
            at4 = dnz*real(l - 1)*at1
            at1 = at1*aimag(ffc(l,1,1))
            zt1 = cmplx(aimag(q(l,1,1)),-real(q(l,1,1)))
            fxyz(1,l,1,1) = zero
            fxyz(2,l,1,1) = zero
            fxyz(3,l,1,1) = at4*zt1
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
            wp = wp + at1*(q(l,1,1)*conjg(q(l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fxyz(1,1,1,1) = zero
            fxyz(2,1,1,1) = zero
            fxyz(3,1,1,1) = zero
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            fxyz(1,l,j,k1) = zero
            fxyz(2,l,j,k1) = zero
            fxyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            fxyz(1,l,1,k1) = zero
            fxyz(2,l,1,k1) = zero
            fxyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      we = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! this subroutine calculates the divergence in fourier space
! for distributed data with 2D spatial decomposition and OpenMP
! input: all except df, output: df
! approximate flop count is:
! 35*nxc*nyc*nzc + 16*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the divergence is calculated using the equations:
! df(kx,ky,kz) = sqrt(-1)*(kx*fx(kx,ky,kz)+ky*fy(kx,ky,kz)
!                       +kz*fz(kx,ky,kz))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! df(kx=pi) = 0, df(ky=pi) = 0, df(kz=pi) = 0
! and df(kx=0,ky=0,kz=0) = 0.
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex f, df
      dimension f(3,nzv,kxyp,kyzp), df(nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate the divergence
      if (kstrt.gt.(nvpy*nvpz)) return
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1)
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkx*f(1,l,j,k) + dky*f(2,l,j,k) + dkz*f(3,l,j,k)
            df(l,j,k) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,k) + dky*f(2,l1,j,k) - dkz*f(3,l1,j,k)
            df(l1,j,k) = cmplx(-aimag(zt1),real(zt1))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,k) + dky*f(2,1,j,k)
            df(1,j,k) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = dky*f(2,l,1,k) + dkz*f(3,l,1,k)
               df(l,1,k) = cmplx(-aimag(zt1),real(zt1))
               zt1 = dky*f(2,l1,1,k) - dkz*f(3,l1,1,k)
               df(l1,1,k) = cmplx(-aimag(zt1),real(zt1))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = dky*f(2,1,1,k)
               df(1,1,k) = cmplx(-aimag(zt1),real(zt1))
               df(l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               df(l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkx*f(1,l,j,1) + dkz*f(3,l,j,1)
            df(l,j,1) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,1) - dkz*f(3,l1,j,1)
            df(l1,j,1) = cmplx(-aimag(zt1),real(zt1))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,1)
            df(1,j,1) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkz*f(3,l,1,1)
            df(l,1,1) = cmplx(-aimag(zt1),real(zt1))
            df(l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            df(1,1,1) = zero
            df(l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            df(l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            df(l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! this subroutine calculates the gradient in fourier space
! for distributed data with 2D spatial decomposition and OpenMP
! input: all except f, output: f
! approximate flop count is:
! 30*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the gradient is calculated using the equations:
! fx(kx,ky,kz) = sqrt(-1)*kx*df(kx,ky,kz)
! fy(kx,ky,kz) = sqrt(-1)*ky*df(kx,ky,kz)
! fz(kx,ky,kz) = sqrt(-1)*kz*df(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
! fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
! fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
! fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex df, f
      dimension df(nzv,kxyp,kyzp), f(3,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate the gradient
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,j,k)),real(df(l,j,k)))
            f(1,l,j,k) = dkx*zt1
            f(2,l,j,k) = dky*zt1
            f(3,l,j,k) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,k)),real(df(l1,j,k)))
            f(1,l1,j,k) = dkx*zt1
            f(2,l1,j,k) = dky*zt1
            f(3,l1,j,k) = -dkz*zt1
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,k)),real(df(1,j,k)))
            f(1,1,j,k) = dkx*zt1
            f(2,1,j,k) = dky*zt1
            f(3,1,j,k) = zero
            f(1,l1,j,k) = zero
            f(2,l1,j,k) = zero
            f(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = cmplx(-aimag(df(l,1,k)),real(df(l,1,k)))
               f(1,l,1,k) = zero
               f(2,l,1,k) = dky*zt1
               f(3,l,1,k) = dkz*zt1
               zt1 = cmplx(-aimag(df(l1,1,k)),real(df(l1,1,k)))
               f(1,l1,1,k) = zero
               f(2,l1,1,k) = dky*zt1
               f(3,l1,1,k) = -dkz*zt1
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(df(1,1,k)),real(df(1,1,k)))
               f(1,1,1,k) = zero
               f(2,1,1,k) = dky*zt1
               f(3,1,1,k) = zero
               f(1,l1,1,k) = zero
               f(2,l1,1,k) = zero
               f(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               f(1,l,1,k) = zero
               f(2,l,1,k) = zero
               f(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,j,1)),real(df(l,j,1)))
            f(1,l,j,1) = dkx*zt1
            f(2,l,j,1) = zero
            f(3,l,j,1) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,1)),real(df(l1,j,1)))
            f(1,l1,j,1) = dkx*zt1
            f(2,l1,j,1) = zero
            f(3,l1,j,1) = -dkz*zt1
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,1)),real(df(1,j,1)))
            f(1,1,j,1) = dkx*zt1
            f(2,1,j,1) = zero
            f(3,1,j,1) = zero
            f(1,l1,j,1) = zero
            f(2,l1,j,1) = zero
            f(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,1,1)),real(df(l,1,1)))
            f(1,l,1,1) = zero
            f(2,l,1,1) = zero
            f(3,l,1,1) = dkz*zt1
            f(1,l1,1,1) = zero
            f(2,l1,1,1) = zero
            f(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            f(1,1,1,1) = zero
            f(2,1,1,1) = zero
            f(3,1,1,1) = zero
            f(1,l1,1,1) = zero
            f(2,l1,1,1) = zero
            f(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            f(1,l,j,k1) = zero
            f(2,l,j,k1) = zero
            f(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            f(1,l,1,k1) = zero
            f(2,l,1,k1) = zero
            f(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
! this subroutine calculates the curl in fourier space
! for distributed data with 2D spatial decomposition and OpenMP
! input: all except g, output: g
! approximate flop count is:
! 86*nxc*nyc*nzc + 32*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the curl is calculated using the equations:
! gx(kx,ky,kz) = sqrt(-1)*(ky*fz(kx,ky,kz)-kz*fy(kx,ky,kz))
! gy(kx,ky,kz) = sqrt(-1)*(kz*fx(kx,ky,kz)-kx*fz(kx,ky,kz))
! gz(kx,ky,kz) = sqrt(-1)*(kx*fy(kx,ky,kz)-ky*fx(kx,ky,kz))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! gx(kx=pi) = gy(kx=pi) = gz(kx=pi) = 0,
! gx(ky=pi) = gy(ky=pi) = gx(ky=pi) = 0,
! gx(kz=pi) = gy(kz=pi) = gz(kz=pi) = 0,
! gx(kx=0,ky=0,kz=0) = gy(kx=0,ky=0,kz=0) = gz(kx=0,ky=0,kz=0) = 0.
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex f, g
      dimension f(3,nzv,kxyp,kyzp), g(3,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate the curl
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1,zt2,zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,k)),real(f(3,l,j,k)))
            zt2 = cmplx(-aimag(f(2,l,j,k)),real(f(2,l,j,k)))
            zt3 = cmplx(-aimag(f(1,l,j,k)),real(f(1,l,j,k)))
            g(1,l,j,k) = dky*zt1 - dkz*zt2
            g(2,l,j,k) = dkz*zt3 - dkx*zt1
            g(3,l,j,k) = dkx*zt2 - dky*zt3
            zt1 = cmplx(-aimag(f(3,l1,j,k)),real(f(3,l1,j,k)))
            zt2 = cmplx(-aimag(f(2,l1,j,k)),real(f(2,l1,j,k)))
            zt3 = cmplx(-aimag(f(1,l1,j,k)),real(f(1,l1,j,k)))
            g(1,l1,j,k) = dky*zt1 + dkz*zt2
            g(2,l1,j,k) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,k) = dkx*zt2 - dky*zt3
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,k)),real(f(3,1,j,k)))
            zt2 = cmplx(-aimag(f(2,1,j,k)),real(f(2,1,j,k)))
            zt3 = cmplx(-aimag(f(1,1,j,k)),real(f(1,1,j,k)))
            g(1,1,j,k) = dky*zt1
            g(2,1,j,k) = -dkx*zt1
            g(3,1,j,k) = dkx*zt2 - dky*zt3
            g(1,l1,j,k) = zero
            g(2,l1,j,k) = zero
            g(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1,zt2,zt3)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = cmplx(-aimag(f(3,l,1,k)),real(f(3,l,1,k)))
               zt2 = cmplx(-aimag(f(2,l,1,k)),real(f(2,l,1,k)))
               zt3 = cmplx(-aimag(f(1,l,1,k)),real(f(1,l,1,k)))
               g(1,l,1,k) = dky*zt1 - dkz*zt2
               g(2,l,1,k) = dkz*zt3
               g(3,l,1,k) = -dky*zt3
               zt1 = cmplx(-aimag(f(3,l1,1,k)),real(f(3,l1,1,k)))
               zt2 = cmplx(-aimag(f(2,l1,1,k)),real(f(2,l1,1,k)))
               zt3 = cmplx(-aimag(f(1,l1,1,k)),real(f(1,l1,1,k)))
               g(1,l1,1,k) = dky*zt1 + dkz*zt2
               g(2,l1,1,k) = -dkz*zt3
               g(3,l1,1,k) = -dky*zt3
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(f(3,1,1,k)),real(f(3,1,1,k)))
               zt3 = cmplx(-aimag(f(1,1,1,k)),real(f(1,1,1,k)))
               g(1,1,1,k) = dky*zt1
               g(2,1,1,k) = zero
               g(3,1,1,k) = -dky*zt3
               g(1,l1,1,k) = zero
               g(2,l1,1,k) = zero
               g(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               g(1,l,1,k) = zero
               g(2,l,1,k) = zero
               g(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,1)),real(f(3,l,j,1)))
            zt2 = cmplx(-aimag(f(2,l,j,1)),real(f(2,l,j,1)))
            zt3 = cmplx(-aimag(f(1,l,j,1)),real(f(1,l,j,1)))
            g(1,l,j,1) = -dkz*zt2
            g(2,l,j,1) = dkz*zt3 - dkx*zt1
            g(3,l,j,1) = dkx*zt2
            zt1 = cmplx(-aimag(f(3,l1,j,1)),real(f(3,l1,j,1)))
            zt2 = cmplx(-aimag(f(2,l1,j,1)),real(f(2,l1,j,1)))
            zt3 = cmplx(-aimag(f(1,l1,j,1)),real(f(1,l1,j,1)))
            g(1,l1,j,1) = dkz*zt2
            g(2,l1,j,1) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,1) = dkx*zt2
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,1)),real(f(3,1,j,1)))
            zt2 = cmplx(-aimag(f(2,1,j,1)),real(f(2,1,j,1)))
            g(1,1,j,1) = zero
            g(2,1,j,1) = -dkx*zt1
            g(3,1,j,1) = dkx*zt2
            g(1,l1,j,1) = zero
            g(2,l1,j,1) = zero
            g(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt2 = cmplx(-aimag(f(2,l,1,1)),real(f(2,l,1,1)))
            zt3 = cmplx(-aimag(f(1,l,1,1)),real(f(1,l,1,1)))
            g(1,l,1,1) = -dkz*zt2
            g(2,l,1,1) = dkz*zt3
            g(3,l,1,1) = zero
            g(1,l1,1,1) = zero
            g(2,l1,1,1) = zero
            g(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            g(1,1,1,1) = zero
            g(2,1,1,1) = zero
            g(3,1,1,1) = zero
            g(1,l1,1,1) = zero
            g(2,l1,1,1) = zero
            g(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            g(1,l,j,k1) = zero
            g(2,l,j,k1) = zero
            g(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            g(1,l,1,k1) = zero
            g(2,l,1,k1) = zero
            g(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp&
     &,kyzp)
! this subroutine calculates 3d vector potential from magnetic field
! in fourier space with periodic boundary conditions,
! for distributed data, with 2D spatial decomposition and OpenMP
! input: bxyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp
! output: axyz
! approximate flop count is:
! 99*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the vector potential is calculated using the equations:
! ax(kx,ky,kz) = sqrt(-1)*
!                (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! ay(kx,ky,kz) = sqrt(-1)*
!                (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! az(kx,ky,kz) = sqrt(-1)*
!                (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
! ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
! ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
! bxyz(i,l,j,k) = i component of complex magnetic field
! axyz(i,l,j,k) = i component of complex vector potential
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex bxyz, axyz
      dimension bxyz(3,nzv,kxyp,kyzp), axyz(3,nzv,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, dkx2, dky2, dkxy2
      real at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate vector potential
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,at2,at3,at4, &
!$OMP& zt1,zt2,zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,k)),real(bxyz(3,l,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,k)),real(bxyz(2,l,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,k)),real(bxyz(1,l,j,k)))
            axyz(1,l,j,k) = at3*zt1 - at4*zt2
            axyz(2,l,j,k) = at4*zt3 - at2*zt1
            axyz(3,l,j,k) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(bxyz(3,l1,j,k)),real(bxyz(3,l1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,k)),real(bxyz(2,l1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,k)),real(bxyz(1,l1,j,k)))
            axyz(1,l1,j,k) = at3*zt1 + at4*zt2
            axyz(2,l1,j,k) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,k) = at2*zt2 - at3*zt3
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(-aimag(bxyz(3,1,j,k)),real(bxyz(3,1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,k)),real(bxyz(2,1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,1,j,k)),real(bxyz(1,1,j,k)))
            axyz(1,1,j,k) = at3*zt1
            axyz(2,1,j,k) = -at2*zt1
            axyz(3,1,j,k) = at2*zt2 - at3*zt3
            axyz(1,l1,j,k) = zero
            axyz(2,l1,j,k) = zero
            axyz(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dky2,dkz,at1,at3,at4,zt1,zt2,zt3&
!$OMP& )
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               at3 = dky*at1
               at4 = dkz*at1
               zt1 = cmplx(-aimag(bxyz(3,l,1,k)),real(bxyz(3,l,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l,1,k)),real(bxyz(2,l,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l,1,k)),real(bxyz(1,l,1,k)))
               axyz(1,l,1,k) = at3*zt1 - at4*zt2
               axyz(2,l,1,k) = at4*zt3
               axyz(3,l,1,k) = -at3*zt3
               zt1 = cmplx(-aimag(bxyz(3,l1,1,k)),real(bxyz(3,l1,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l1,1,k)),real(bxyz(2,l1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l1,1,k)),real(bxyz(1,l1,1,k)))
               axyz(1,l1,1,k) = at3*zt1 + at4*zt2
               axyz(2,l1,1,k) = -at4*zt3
               axyz(3,l1,1,k) = -at3*zt3
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at3 = 1.0/dky
               zt1 = cmplx(-aimag(bxyz(3,1,1,k)),real(bxyz(3,1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,1,1,k)),real(bxyz(1,1,1,k)))
               axyz(1,1,1,k) = at3*zt1
               axyz(2,1,1,k) = zero
               axyz(3,1,1,k) = -at3*zt3
               axyz(1,l1,1,k) = zero
               axyz(2,l1,1,k) = zero
               axyz(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k) = zero
               axyz(2,l,1,k) = zero
               axyz(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkx2,dkz,at1,at2,at4,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            at2 = dkx*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,1)),real(bxyz(3,l,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,1)),real(bxyz(2,l,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,1)),real(bxyz(1,l,j,1)))
            axyz(1,l,j,1) = -at4*zt2
            axyz(2,l,j,1) = at4*zt3 - at2*zt1
            axyz(3,l,j,1) = at2*zt2
            zt1 = cmplx(-aimag(bxyz(3,l1,j,1)),real(bxyz(3,l1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,1)),real(bxyz(2,l1,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,1)),real(bxyz(1,l1,j,1)))
            axyz(1,l1,j,1) = at4*zt2
            axyz(2,l1,j,1) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,1) = at2*zt2
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = 1.0/dkx
            zt1 = cmplx(-aimag(bxyz(3,1,j,1)),real(bxyz(3,1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,1)),real(bxyz(2,1,j,1)))
            axyz(1,1,j,1) = zero
            axyz(2,1,j,1) = -at2*zt1
            axyz(3,1,j,1) = at2*zt2
            axyz(1,l1,j,1) = zero
            axyz(2,l1,j,1) = zero
            axyz(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at4 = 1.0/dkz
            zt2 = cmplx(-aimag(bxyz(2,l,1,1)),real(bxyz(2,l,1,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,1,1)),real(bxyz(1,l,1,1)))
            axyz(1,l,1,1) = -at4*zt2
            axyz(2,l,1,1) = at4*zt3
            axyz(3,l,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1) = zero
            axyz(2,1,1,1) = zero
            axyz(3,1,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1) = zero
            axyz(2,l,j,k1) = zero
            axyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1) = zero
            axyz(2,l,1,k1) = zero
            axyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MCUAVE33(cuave,cunew,cuold,nz,kxyp,kyzp,nzv)
! this subroutine averages current in fourier space for 2-1/2d code
! input: all except cuave, output: cuave
! cunew(i,l,j,k),cuold(i,l,j,k) = complex current densities to be
! averaged
! cuave(i,l,j,k) = average complex current density
! for component i, all for fourier mode (j-1,k-1,l-1)
! nz = system length in z direction
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzv = second dimension of field arrays, must be >= nz
      implicit none
      integer nz, kxyp, kyzp, nzv
      complex cuave, cunew, cuold
      dimension cuave(3,nzv,kxyp,kyzp), cuold(3,nzv,kxyp,kyzp)
      dimension cunew(3,nzv,kxyp,kyzp)
! local data
      integer j, k, l
!$OMP PARALLEL DO PRIVATE(j,k,l)
      do 30 k = 1, kyzp
      do 20 j = 1, kxyp
      do 10 l = 1, nz
      cuave(1,l,j,k) = 0.5*(cunew(1,l,j,k) + cuold(1,l,j,k))
      cuave(2,l,j,k) = 0.5*(cunew(2,l,j,k) + cuold(2,l,j,k))
      cuave(3,l,j,k) = 0.5*(cunew(3,l,j,k) + cuold(3,l,j,k))
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAVRPOT332(axyz,bxyz,ffc,affp,ci,nx,ny,nz,kstrt,nvpy,&
     &nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for the
! radiative part of the vector potential
! with periodic boundary conditions,
! for distributed data, with 2D spatial decomposition and OpenMP
! input: axyz,bxyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp
! output: axyz
! approximate flop count is:
! 105*nxc*nyc*nzc + 42*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! the vector potential is calculated using the equations:
! ax(kx,ky,kz) = sqrt(-1)*
! (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz) - affp*ci2*cux(kx,ky,kz)*s(kx,ky,kz))
! /(kx*kx+ky*ky+kz*kz)
! ay(kx,ky,kz) = sqrt(-1)*
! (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz) + affp*ci2*cuy(kx,ky,kz)*s(kx,ky,kz))
! /(kx*kx+ky*ky+kz*kz)
! az(kx,ky,kz) = sqrt(-1)*
! (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz) - affp*ci2*cuz(kx,ky,kz)*s(kx,ky,kz))
! /(kx*kx+ky*ky+kz*kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers, except for
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
! ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
! ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
! axyz(i,l,j,k) = on entry, complex current density cu
! axyz(i,l,j,k) = on exit, complex current radiative vector potential
! bxyz(i,l,j,k) = complex magnetic field
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! ci = reciprical of velocity of light
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real affp, ci
      complex axyz, bxyz, ffc
      dimension axyz(3,nzv,kxyp,kyzp), bxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, afc2, dkx, dky, dkz, dkx2, dky2, dkxy2
      real at1, at2
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      afc2 = affp*ci*ci
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate the radiative vector potential
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL                                                          &
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,at2,zt1,zt2, &
!$OMP& zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            at2 = afc2*aimag(ffc(l,j,k))
            zt1 = cmplx(-aimag(bxyz(3,l,j,k)),real(bxyz(3,l,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,k)),real(bxyz(2,l,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,k)),real(bxyz(1,l,j,k)))
            axyz(1,l,j,k) = at1*(dky*zt1 - dkz*zt2 - at2*axyz(1,l,j,k))
            axyz(2,l,j,k) = at1*(dkz*zt3 - dkx*zt1 - at2*axyz(2,l,j,k))
            axyz(3,l,j,k) = at1*(dkx*zt2 - dky*zt3 - at2*axyz(3,l,j,k))
            zt1 = cmplx(-aimag(bxyz(3,l1,j,k)),real(bxyz(3,l1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,k)),real(bxyz(2,l1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,k)),real(bxyz(1,l1,j,k)))
            axyz(1,l1,j,k) = at1*(dky*zt1 + dkz*zt2                     &
     &                     - at2*axyz(1,l1,j,k))
            axyz(2,l1,j,k) = -at1*(dkz*zt3 + dkx*zt1                    &
     &                     + at2*axyz(2,l1,j,k))
            axyz(3,l1,j,k) = at1*(dkx*zt2 - dky*zt3                     &
     &                     - at2*axyz(3,l1,j,k))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            at2 = afc2*aimag(ffc(1,j,k))
            zt1 = cmplx(-aimag(bxyz(3,1,j,k)),real(bxyz(3,1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,k)),real(bxyz(2,1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,1,j,k)),real(bxyz(1,1,j,k)))
            axyz(1,1,j,k) = at1*(dky*zt1 - at2*axyz(1,1,j,k))
            axyz(2,1,j,k) = -at1*(dkx*zt1 + at2*axyz(2,1,j,k))
            axyz(3,1,j,k) = at1*(dkx*zt2 - dky*zt3 - at2*axyz(3,1,j,k))
            axyz(1,l1,j,k) = zero
            axyz(2,l1,j,k) = zero
            axyz(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,dky2,at1,at2,zt1,zt2,zt3)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               at2 = afc2*aimag(ffc(l,1,k))
               zt1 = cmplx(-aimag(bxyz(3,l,1,k)),real(bxyz(3,l,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l,1,k)),real(bxyz(2,l,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l,1,k)),real(bxyz(1,l,1,k)))
               axyz(1,l,1,k) = at1*(dky*zt1 - dkz*zt2                   &
     &                       - at2*axyz(1,l,1,k))
               axyz(2,l,1,k) = at1*(dkz*zt3 - at2*axyz(2,l,1,k))
               axyz(3,l,1,k) = -at1*(dky*zt3 + at2*axyz(3,l,1,k))
               zt1 = cmplx(-aimag(bxyz(3,l1,1,k)),real(bxyz(3,l1,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l1,1,k)),real(bxyz(2,l1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l1,1,k)),real(bxyz(1,l1,1,k)))
               axyz(1,l1,1,k) = at1*(dky*zt1 + dkz*zt2                  &
     &                        - at2*axyz(1,l1,1,k))
               axyz(2,l1,1,k) = -at1*(dkz*zt3 + at2*axyz(2,l1,1,k))
               axyz(3,l1,1,k) = -at1*(dky*zt3 + at2*axyz(3,l1,1,k))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = 1.0/(dky*dky)
               at2 = afc2*aimag(ffc(1,1,k))
               zt1 = cmplx(-aimag(bxyz(3,1,1,k)),real(bxyz(3,1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,1,1,k)),real(bxyz(1,1,1,k)))
               axyz(1,1,1,k) = at1*(dky*zt1 - at2*axyz(1,1,1,k))
               axyz(2,1,1,k) = zero
               axyz(3,1,1,k) = -at1*(dky*zt3 + at2*axyz(3,1,1,k))
               axyz(1,l1,1,k) = zero
               axyz(2,l1,1,k) = zero
               axyz(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k) = zero
               axyz(2,l,1,k) = zero
               axyz(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,dkx2,at1,at2,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            at2 = afc2*aimag(ffc(l,j,1))
            zt1 = cmplx(-aimag(bxyz(3,l,j,1)),real(bxyz(3,l,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,1)),real(bxyz(2,l,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,1)),real(bxyz(1,l,j,1)))
            axyz(1,l,j,1) = -at1*(dkz*zt2 + at2*axyz(1,l,j,1))
            axyz(2,l,j,1) = at1*(dkz*zt3 - dkx*zt1 - at2*axyz(2,l,j,1))
            axyz(3,l,j,1) = at1*(dkx*zt2 - at2*axyz(3,l,j,1))
            zt1 = cmplx(-aimag(bxyz(3,l1,j,1)),real(bxyz(3,l1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,1)),real(bxyz(2,l1,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,1)),real(bxyz(1,l1,j,1)))
            axyz(1,l1,j,1) = at1*(dkz*zt2 - at2*axyz(1,l1,j,1))
            axyz(2,l1,j,1) = -at1*(dkz*zt3 + dkx*zt1                    &
     &                     + at2*axyz(2,l1,j,1))
            axyz(3,l1,j,1) = at1*(dkx*zt2 - at2*axyz(3,l1,j,1))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/(dkx*dkx)
            at2 = afc2*aimag(ffc(1,j,1))
            zt1 = cmplx(-aimag(bxyz(3,1,j,1)),real(bxyz(3,1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,1)),real(bxyz(2,1,j,1)))
            axyz(1,1,j,1) = zero
            axyz(2,1,j,1) = -at1*(dkx*zt1 + at2*axyz(2,1,j,1))
            axyz(3,1,j,1) = at1*(dkx*zt2 - at2*axyz(3,1,j,1))
            axyz(1,l1,j,1) = zero
            axyz(2,l1,j,1) = zero
            axyz(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz)
            at2 = afc2*aimag(ffc(l,1,1))
            zt2 = cmplx(-aimag(bxyz(2,l,1,1)),real(bxyz(2,l,1,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,1,1)),real(bxyz(1,l,1,1)))
            axyz(1,l,1,1) = -at1*(dkz*zt2 + at2*axyz(1,l,1,1))
            axyz(2,l,1,1) = at1*(dkz*zt3 - at2*axyz(2,l,1,1))
            axyz(3,l,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1) = zero
            axyz(2,1,1,1) = zero
            axyz(3,1,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1) = zero
            axyz(2,l,j,k1) = zero
            axyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1) = zero
            axyz(2,l,1,k1) = zero
            axyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAPOTP32(cu,axyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz, &
     &nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! for vector potential with periodic boundary conditions
! for distributed data, with 2D spatial decomposition and OpenMP
! input: cu,ffc,ci,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
! output: axyz,wm
! approximate flop count is:
! 128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! vector potential is calculated using the equation:
! ax(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cux(kx,ky,kz)
! ay(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuy(kx,ky,kz)
! az(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuz(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
! ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
! ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
! cu(l,j,k) = complex current density for fourier mode jj-1,kk-1,l-1
! axyz(1,l,j,k) = x component of complex vector potential
! axyz(2,l,j,k) = y component of complex vector potential
! axyz(3,l,j,k) = z component of complex vector potential
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! ci = reciprical of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real ci, wm
      complex cu, axyz, ffc
      dimension cu(3,nzv,kxyp,kyzp), axyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real ci2, at1, at2
      complex zero
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      ci2 = ci*ci
! calculate vector potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,j,k))
            at1 = at2*aimag(ffc(l,j,k))
            axyz(1,l,j,k) = at2*cu(1,l,j,k)
            axyz(2,l,j,k) = at2*cu(2,l,j,k)
            axyz(3,l,j,k) = at2*cu(3,l,j,k)
            axyz(1,l1,j,k) = at2*cu(1,l1,j,k)
            axyz(2,l1,j,k) = at2*cu(2,l1,j,k)
            axyz(3,l1,j,k) = at2*cu(3,l1,j,k)
            wp = wp + at1*(cu(1,l,j,k)*conjg(cu(1,l,j,k))               &
     &              + cu(2,l,j,k)*conjg(cu(2,l,j,k))                    &
     &              + cu(3,l,j,k)*conjg(cu(3,l,j,k))                    &
     &              + cu(1,l1,j,k)*conjg(cu(1,l1,j,k))                  &
     &              + cu(2,l1,j,k)*conjg(cu(2,l1,j,k))                  &
     &              + cu(3,l1,j,k)*conjg(cu(3,l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = ci2*real(ffc(1,j,k))
            at1 = at2*aimag(ffc(1,j,k))
            axyz(1,1,j,k) = at2*cu(1,1,j,k)
            axyz(2,1,j,k) = at2*cu(2,1,j,k)
            axyz(3,1,j,k) = at2*cu(3,1,j,k)
            axyz(1,l1,j,k) = zero
            axyz(2,l1,j,k) = zero
            axyz(3,l1,j,k) = zero
            wp = wp + at1*(cu(1,1,j,k)*conjg(cu(1,1,j,k))               &
     &              + cu(2,1,j,k)*conjg(cu(2,1,j,k))                    &
     &              + cu(3,1,j,k)*conjg(cu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at2 = ci2*real(ffc(l,1,k))
               at1 = at2*aimag(ffc(l,1,k))
               axyz(1,l,1,k) = at2*cu(1,l,1,k)
               axyz(2,l,1,k) = at2*cu(2,l,1,k)
               axyz(3,l,1,k) = at2*cu(3,l,1,k)
               axyz(1,l1,1,k) = at2*cu(1,l1,1,k)
               axyz(2,l1,1,k) = at2*cu(2,l1,1,k)
               axyz(3,l1,1,k) = at2*cu(3,l1,1,k)
               wp = wp + at1*(cu(1,l,1,k)*conjg(cu(1,l,1,k))            &
     &                 + cu(2,l,1,k)*conjg(cu(2,l,1,k))                 &
     &                 + cu(3,l,1,k)*conjg(cu(3,l,1,k))                 &
     &                 + cu(1,l1,1,k)*conjg(cu(1,l1,1,k))               &
     &                 + cu(2,l1,1,k)*conjg(cu(2,l1,1,k))               &
     &                 + cu(3,l1,1,k)*conjg(cu(3,l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = ci2*real(ffc(1,1,k))
               at1 = at2*aimag(ffc(1,1,k))
               axyz(1,1,1,k) = at2*cu(1,1,1,k)
               axyz(2,1,1,k) = at2*cu(2,1,1,k)
               axyz(3,1,1,k) = at2*cu(3,1,1,k)
               axyz(1,l1,1,k) = zero
               axyz(2,l1,1,k) = zero
               axyz(3,l1,1,k) = zero
               wp = wp + at1*(cu(1,1,1,k)*conjg(cu(1,1,1,k))            &
     &                 + cu(2,1,1,k)*conjg(cu(2,1,1,k))                 &
     &                 + cu(3,1,1,k)*conjg(cu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k) = zero
               axyz(2,l,1,k) = zero
               axyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,j,1))
            at1 = at2*aimag(ffc(l,j,1))
            axyz(1,l,j,1) = at2*cu(1,l,j,1)
            axyz(2,l,j,1) = at2*cu(2,l,j,1)
            axyz(3,l,j,1) = at2*cu(3,l,j,1)
            axyz(1,l1,j,1) = at2*cu(1,l1,j,1)
            axyz(2,l1,j,1) = at2*cu(2,l1,j,1)
            axyz(3,l1,j,1) = at2*cu(3,l1,j,1)
            wp = wp + at1*(cu(1,l,j,1)*conjg(cu(1,l,j,1))               &
     &              + cu(2,l,j,1)*conjg(cu(2,l,j,1))                    &
     &              + cu(3,l,j,1)*conjg(cu(3,l,j,1))                    &
     &              + cu(1,l1,j,1)*conjg(cu(1,l1,j,1))                  &
     &              + cu(2,l1,j,1)*conjg(cu(2,l1,j,1))                  &
     &              + cu(3,l1,j,1)*conjg(cu(3,l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = ci2*real(ffc(1,j,1))
            at1 = at2*aimag(ffc(1,j,1))
            axyz(1,1,j,1) = at2*cu(1,1,j,1)
            axyz(2,1,j,1) = at2*cu(2,1,j,1)
            axyz(3,1,j,1) = at2*cu(3,1,j,1)
            axyz(1,l1,j,1) = zero
            axyz(2,l1,j,1) = zero
            axyz(3,l1,j,1) = zero
            wp = wp + at1*(cu(1,1,j,1)*conjg(cu(1,1,j,1))               &
     &              + cu(2,1,j,1)*conjg(cu(2,1,j,1))                    &
     &              + cu(3,1,j,1)*conjg(cu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,1,1))
            at1 = at2*aimag(ffc(l,1,1))
            axyz(1,l,1,1) = at2*cu(1,l,1,1)
            axyz(2,l,1,1) = at2*cu(2,l,1,1)
            axyz(3,l,1,1) = at2*cu(3,l,1,1)
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
            wp = wp + at1*(cu(1,l,1,1)*conjg(cu(1,l,1,1))               &
     &              + cu(2,l,1,1)*conjg(cu(2,l,1,1))                    &
     &              + cu(3,l,1,1)*conjg(cu(3,l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1) = zero
            axyz(2,1,1,1) = zero
            axyz(3,1,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1) = zero
            axyz(2,l,j,k1) = zero
            axyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1) = zero
            axyz(2,l,1,k1) = zero
            axyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      wm = real(nx)*real(ny)*real(nz)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPETFIELD332(dcu,exyz,ffe,affp,ci,wf,nx,ny,nz,kstrt,  &
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
! this subroutine solves 3d poisson's equation in fourier space for
! unsmoothed transverse electric field, with periodic boundary
! conditions, for distributed data, with 2D spatial decomposition
! and with OpenMP
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! input: dcu,ffe,affp,ci,nx,ny,nz,kstrt,nzv,kxyp,kyzp,nzhd
! output: exyz,wf
! approximate flop count is:
! 128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! unsmoothed transverse electric field is calculated using the equation:
! ex(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcux(kx,ky,kz)
! ey(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuy(kx,ky,kz)
! ez(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuz(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
! ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
! ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
! dcu(l,j,k) = transverse part of complex derivative of current for
! fourier mode jj-1,kk-1,l-1
! exyz(1,l,j,k) = x component of complex transverse electric field
! exyz(2,l,j,k) = y component of complex transverse electric field
! exyz(3,l,j,k) = z component of complex transverse electric field
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffe(l,j,k)) = finite-size particle shape factor s
! real(ffe(l,j,k)) = potential green's function g
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! affp = normalization constant = nx*ny*nz/np,
! where np=number of particles
! ci = reciprical of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*nz*sum((affp/((kx**2+ky**2+kz**2)*ci*ci)**2)
!    |dcu(kx,ky,kz)*s(kx,ky,kz)|**2)
! this expression is valid only if the derivative of current is
! divergence-free
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real affp, ci, wf
      complex dcu, exyz, ffe
      dimension dcu(3,nzv,kxyp,kyzp), exyz(3,nzv,kxyp,kyzp)
      dimension ffe(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real ci2, at1, at2
      complex zero
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      ci2 = ci*ci
! calculate unsmoothed transverse electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,k))
            at1 = at2*at2
            exyz(1,l,j,k) = at2*dcu(1,l,j,k)
            exyz(2,l,j,k) = at2*dcu(2,l,j,k)
            exyz(3,l,j,k) = at2*dcu(3,l,j,k)
            exyz(1,l1,j,k) = at2*dcu(1,l1,j,k)
            exyz(2,l1,j,k) = at2*dcu(2,l1,j,k)
            exyz(3,l1,j,k) = at2*dcu(3,l1,j,k)
            wp = wp + at1*(dcu(1,l,j,k)*conjg(dcu(1,l,j,k))             &
     &              + dcu(2,l,j,k)*conjg(dcu(2,l,j,k))                  &
     &              + dcu(3,l,j,k)*conjg(dcu(3,l,j,k))                  &
     &              + dcu(1,l1,j,k)*conjg(dcu(1,l1,j,k))                &
     &              + dcu(2,l1,j,k)*conjg(dcu(2,l1,j,k))                &
     &              + dcu(3,l1,j,k)*conjg(dcu(3,l1,j,k)))
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,k))
            at1 = at2*at2
            exyz(1,1,j,k) = at2*dcu(1,1,j,k)
            exyz(2,1,j,k) = at2*dcu(2,1,j,k)
            exyz(3,1,j,k) = at2*dcu(3,1,j,k)
            exyz(1,l1,j,k) = zero
            exyz(2,l1,j,k) = zero
            exyz(3,l1,j,k) = zero
            wp = wp + at1*(dcu(1,1,j,k)*conjg(dcu(1,1,j,k))             &
     &              + dcu(2,1,j,k)*conjg(dcu(2,1,j,k))                  &
     &              + dcu(3,1,j,k)*conjg(dcu(3,1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at2 = -ci2*real(ffe(l,1,k))
               at1 = at2*at2
               exyz(1,l,1,k) = at2*dcu(1,l,1,k)
               exyz(2,l,1,k) = at2*dcu(2,l,1,k)
               exyz(3,l,1,k) = at2*dcu(3,l,1,k)
               exyz(1,l1,1,k) = at2*dcu(1,l1,1,k)
               exyz(2,l1,1,k) = at2*dcu(2,l1,1,k)
               exyz(3,l1,1,k) = at2*dcu(3,l1,1,k)
               wp = wp + at1*(dcu(1,l,1,k)*conjg(dcu(1,l,1,k))          &
     &                + dcu(2,l,1,k)*conjg(dcu(2,l,1,k))                &
     &                + dcu(3,l,1,k)*conjg(dcu(3,l,1,k))                &
     &                + dcu(1,l1,1,k)*conjg(dcu(1,l1,1,k))              &
     &                + dcu(2,l1,1,k)*conjg(dcu(2,l1,1,k))              &
     &                + dcu(3,l1,1,k)*conjg(dcu(3,l1,1,k)))
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = -ci2*real(ffe(1,1,k))
               at1 = at2*at2
               exyz(1,1,1,k) = at2*dcu(1,1,1,k)
               exyz(2,1,1,k) = at2*dcu(2,1,1,k)
               exyz(3,1,1,k) = at2*dcu(3,1,1,k)
               exyz(1,l1,1,k) = zero
               exyz(2,l1,1,k) = zero
               exyz(3,l1,1,k) = zero
               wp = wp + at1*(dcu(1,1,1,k)*conjg(dcu(1,1,1,k))          &
     &                 + dcu(2,1,1,k)*conjg(dcu(2,1,1,k))               &
     &                 + dcu(3,1,1,k)*conjg(dcu(3,1,1,k)))
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               exyz(1,l,1,k) = zero
               exyz(2,l,1,k) = zero
               exyz(3,l,1,k) = zero
   40          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
! mode numbers ky = 0, ny/2
      sum2 = 0.0d0
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,j,1))
            at1 = at2*at2
            exyz(1,l,j,1) = at2*dcu(1,l,j,1)
            exyz(2,l,j,1) = at2*dcu(2,l,j,1)
            exyz(3,l,j,1) = at2*dcu(3,l,j,1)
            exyz(1,l1,j,1) = at2*dcu(1,l1,j,1)
            exyz(2,l1,j,1) = at2*dcu(2,l1,j,1)
            exyz(3,l1,j,1) = at2*dcu(3,l1,j,1)
            wp = wp + at1*(dcu(1,l,j,1)*conjg(dcu(1,l,j,1))             &
     &              + dcu(2,l,j,1)*conjg(dcu(2,l,j,1))                  &
     &              + dcu(3,l,j,1)*conjg(dcu(3,l,j,1))                  &
     &              + dcu(1,l1,j,1)*conjg(dcu(1,l1,j,1))                &
     &              + dcu(2,l1,j,1)*conjg(dcu(2,l1,j,1))                &
     &              + dcu(3,l1,j,1)*conjg(dcu(3,l1,j,1)))
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = -ci2*real(ffe(1,j,1))
            at1 = at2*at2
            exyz(1,1,j,1) = at2*dcu(1,1,j,1)
            exyz(2,1,j,1) = at2*dcu(2,1,j,1)
            exyz(3,1,j,1) = at2*dcu(3,1,j,1)
            exyz(1,l1,j,1) = zero
            exyz(2,l1,j,1) = zero
            exyz(3,l1,j,1) = zero
            wp = wp + at1*(dcu(1,1,j,1)*conjg(dcu(1,1,j,1))             &
     &              + dcu(2,1,j,1)*conjg(dcu(2,1,j,1))                  &
     &              + dcu(3,1,j,1)*conjg(dcu(3,1,j,1)))
         endif
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at2 = -ci2*real(ffe(l,1,1))
            at1 = at2*at2
            exyz(1,l,1,1) = at2*dcu(1,l,1,1)
            exyz(2,l,1,1) = at2*dcu(2,l,1,1)
            exyz(3,l,1,1) = at2*dcu(3,l,1,1)
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
            wp = wp + at1*(dcu(1,l,1,1)*conjg(dcu(1,l,1,1))             &
     &              + dcu(2,l,1,1)*conjg(dcu(2,l,1,1))                  &
     &              + dcu(3,l,1,1)*conjg(dcu(3,l,1,1)))
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            exyz(1,1,1,1) = zero
            exyz(2,1,1,1) = zero
            exyz(3,1,1,1) = zero
            exyz(1,l1,1,1) = zero
            exyz(2,l1,1,1) = zero
            exyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            exyz(1,l,j,k1) = zero
            exyz(2,l,j,k1) = zero
            exyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            exyz(1,l,1,k1) = zero
            exyz(2,l,1,k1) = zero
            exyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      sum1 = sum1 + sum2
  120 continue
      wf = real(nx)*real(ny)*real(nz)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,&
     &kyzp,nzhd)
! this subroutine provides a 3d smoothing function,
! in fourier space, with periodic boundary conditions
! for distributed data, with 2D spatial decomposition and OpenMP
! input: q,ffc,nx,ny,nz,kstrt,nzv,kxyp,kyzp,nzhd, output: qs
! approximate flop count is:
! 8*nxc*nyc*nzc + 4*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! smoothing is calculated using the equation:
! qs(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! qs(kx=pi) = 0, qs(ky=pi) = 0, qs(kz=pi) = 0, and
! qs(kx=0,ky=0,kz=0) = 0.
! q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
! qs(l,j,k) = complex smoothed charge density,
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      complex q, qs, ffc
      dimension q(nzv,kxyp,kyzp), qs(nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate smoothing
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k))
            qs(l,j,k) = at1*q(l,j,k)
            qs(l1,j,k) = at1*q(l1,j,k)
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k))
            qs(1,j,k) = at1*q(1,j,k)
            qs(l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k))
               qs(l,1,k) = at1*q(l,1,k)
               qs(l1,1,k) = at1*q(l1,1,k)
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k))
               qs(1,1,k) = at1*q(1,1,k)
               qs(l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               qs(l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1)
         do 70 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1))
            qs(l,j,1) = at1*q(l,j,1)
            qs(l1,j,1) = at1*q(l1,j,1)
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1))
            qs(1,j,1) = at1*q(1,j,1)
            qs(l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1))
            qs(l,1,1) = at1*q(l,1,1)
            qs(l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1))
            qs(1,1,1) = cmplx(at1*real(q(1,1,1)),0.0)
            qs(l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            qs(l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            qs(l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  &
     &kxyp,kyzp,nzhd)
! this subroutine provides a 3d vector smoothing function
! in fourier space, with periodic boundary conditions
! or distributed data, with 2D spatial decomposition and OpenMP
! input: cu,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
! output: cus
! approximate flop count is:
! 24*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
! where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
! nvpy/nvpz = number of procs in y/z
! smoothing is calculated using the equation:
! cusx(kx,ky,kz) = cux(kx,ky,kz)*s(kx,ky,kz)
! cusy(kx,ky,kz) = cuy(kx,ky,kz)*s(kx,ky,kz)
! cusz(kx,ky,kz) = cuz(kx,ky,kz)*s(kx,ky,kz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
! j,k,l = fourier mode numbers,
! g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
! s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
! cusx(kx=pi) = cusy(kx=pi) = cusz(kx=pi) = 0,
! cusx(ky=pi) = cusy(ky=pi) = cusx(ky=pi) = 0,
! cusx(kz=pi) = cusy(kz=pi) = cusz(kz=pi) = 0,
! cusx(kx=0,ky=0,kz=0) = cusy(kx=0,ky=0,kz=0) = cusz(kx=0,ky=0,kz=0) = 0
! cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
! cus(1,l,j,k,m) = x component of complex smoothed current density
! cus(2,l,j,k,m) = y component of complex smoothed current density
! cus(3,l,j,k,m) = z component of complex smoothed current density
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! aimag(ffc(l,j,k)) = finite-size particle shape factor s
! real(ffc(l,j,k)) = potential green's function g
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      complex cu, cus, ffc
      dimension cu(3,nzv,kxyp,kyzp), cus(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
! calculate smoothing
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k))
            cus(1,l,j,k) = at1*cu(1,l,j,k)
            cus(2,l,j,k) = at1*cu(2,l,j,k)
            cus(3,l,j,k) = at1*cu(3,l,j,k)
            cus(1,l1,j,k) = at1*cu(1,l1,j,k)
            cus(2,l1,j,k) = at1*cu(2,l1,j,k)
            cus(3,l1,j,k) = at1*cu(3,l1,j,k)
   10       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k))
            cus(1,1,j,k) = at1*cu(1,1,j,k)
            cus(2,1,j,k) = at1*cu(2,1,j,k)
            cus(3,1,j,k) = at1*cu(3,1,j,k)
            cus(1,l1,j,k) = zero
            cus(2,l1,j,k) = zero
            cus(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
! mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if (js.eq.0) then
! keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k))
               cus(1,l,1,k) = at1*cu(1,l,1,k)
               cus(2,l,1,k) = at1*cu(2,l,1,k)
               cus(3,l,1,k) = at1*cu(3,l,1,k)
               cus(1,l1,1,k) = at1*cu(1,l1,1,k)
               cus(2,l1,1,k) = at1*cu(2,l1,1,k)
               cus(3,l1,1,k) = at1*cu(3,l1,1,k)
   30          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k))
               cus(1,1,1,k) = at1*cu(1,1,1,k)
               cus(2,1,1,k) = at1*cu(2,1,1,k)
               cus(3,1,1,k) = at1*cu(3,1,1,k)
               cus(1,l1,1,k) = zero
               cus(2,l1,1,k) = zero
               cus(3,l1,1,k) = zero
! throw away kx = nx/2
            else
               do 40 l = 1, nz
               cus(1,l,1,k) = zero
               cus(2,l,1,k) = zero
               cus(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1)
         do 70 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1))
            cus(1,l,j,1) = at1*cu(1,l,j,1)
            cus(2,l,j,1) = at1*cu(2,l,j,1)
            cus(3,l,j,1) = at1*cu(3,l,j,1)
            cus(1,l1,j,1) = at1*cu(1,l1,j,1)
            cus(2,l1,j,1) = at1*cu(2,l1,j,1)
            cus(3,l1,j,1) = at1*cu(3,l1,j,1)
   60       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1))
            cus(1,1,j,1) = at1*cu(1,1,j,1)
            cus(2,1,j,1) = at1*cu(2,1,j,1)
            cus(3,1,j,1) = at1*cu(3,1,j,1)
            cus(1,l1,j,1) = zero
            cus(2,l1,j,1) = zero
            cus(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1))
            cus(1,l,1,1) = at1*cu(1,l,1,1)
            cus(2,l,1,1) = at1*cu(2,l,1,1)
            cus(3,l,1,1) = at1*cu(3,l,1,1)
            cus(1,l1,1,1) = zero
            cus(2,l1,1,1) = zero
            cus(3,l1,1,1) = zero
   80       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1))
            cus(1,1,1,1) = cmplx(at1*real(cu(1,1,1,1)),0.)
            cus(2,1,1,1) = cmplx(at1*real(cu(2,1,1,1)),0.)
            cus(3,1,1,1) = cmplx(at1*real(cu(3,1,1,1)),0.)
            cus(1,l1,1,1) = zero
            cus(2,l1,1,1) = zero
            cus(3,l1,1,1) = zero
         endif
      endif
! throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            cus(1,l,j,k1) = zero
            cus(2,l,j,k1) = zero
            cus(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            cus(1,l,1,k1) = zero
            cus(2,l,1,k1) = zero
            cus(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRDMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,    &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! this subroutine extracts lowest order modes from packed complex array
! pot and stores them into a location in an unpacked complex array pott
! modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
! and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
! and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
! where js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
! and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! except kx = NX/2 is stored at location kxyp+1 when js=0,
! and ky = NY/2 is stored at location 1 when ks=nvp/2.
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! modeszd = first dimension of array pott,
! where modeszd  = min(2*modesz-1,nz)
! modesypd = third dimension of array pott, modesypd >= kyzp
! modesxpd = second dimension of array pott,
! modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
! in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, kstrt, nvpy, nvpz, nzv
      integer kxyp, kyzp, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp), pott(modeszd,modesxpd,modesypd)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))&
     & then
         if (k1.gt.nyh) k1 = k1 - ny
         do 20 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pott(2*l-2,j,k) = pot(l,j,k)
            pott(2*l-1,j,k) = pot(l1,j,k)
   10       continue
! mode numbers kz = 0, nz/2
            pott(1,j,k) = pot(1,j,k)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,j,k) = pot(l1,j,k)
            endif
         endif
   20    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            if (k1.gt.0) then
               do 30 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,1,k) = pot(l,1,k)
               pott(2*l-1,1,k) = pot(l1,1,k)
   30          continue
! mode numbers kz = 0, nz/2
               pott(1,1,k) = pot(1,1,k)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,1,k) = pot(l1,1,k)
               endif
! kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 40 l = 2, lmax
                  l1 = nz2 - l
                  pott(2*l-2,1,k) = pot(l,1,k)
                  pott(2*l-1,1,k) = pot(l1,1,k)
   40             continue
! mode numbers kz = 0, nz/2
                  pott(1,1,k) = pot(1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(nz,1,k) = pot(l1,1,k)
                  endif
               endif
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! ky = 0
      if (ks.eq.0) then
         do 70 j = 1, jmax
         if ((j+joff).gt.1) then
            do 60 l = 2, lmax
            l1 = nz2 - l
            pott(2*l-2,j,1) = pot(l,j,1)
            pott(2*l-1,j,1) = pot(l1,j,1)
   60       continue
! mode numbers kz = 0, nz/2
            pott(1,j,1) = pot(1,j,1)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,j,1) = pot(l1,j,1)
            endif
         endif
   70    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            do 80 l = 2, lmax
            pott(2*l-2,1,1) = pot(l,1,1)
            pott(2*l-1,1,1) = conjg(pot(l,1,1))
   80       continue
! mode numbers kz = 0, nz/2
            pott(1,1,1) = cmplx(real(pot(1,1,1)),0.0)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,1,1) = cmplx(real(pot(l1,1,1)),0.0)
            endif
! kx = nx/2
            if (modesx.gt.nxh) then
               do 90 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,j1,1) = conjg(pot(l1,1,1))
               pott(2*l-1,j1,1) = pot(l1,1,1)
   90          continue
! mode numbers kz = 0, nz/2
               pott(1,j1,1) = cmplx(aimag(pot(1,1,1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,j1,1) = cmplx(aimag(pot(l1,1,1)),0.0)
               endif
            endif
         endif
      endif
! ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         if (modesy.gt.nyh) then
            do 110 j = 1, jmax
            if ((j+joff).gt.1) then
               do 100 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,j,k1) = pot(l,j,k1)
               pott(2*l-1,j,k1) = pot(l1,j,k1)
  100          continue
! mode numbers kz = 0, nz/2
               pott(1,j,k1) = pot(1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,j,k1) = pot(l1,j,k1)
               endif
            endif
  110       continue
! mode numbers kx = 0, nx/2
            if (js.eq.0) then
! kx = 0
               do 120 l = 2, lmax
               pott(2*l-2,1,k1) = pot(l,1,k1)
               pott(2*l-1,1,k1) = conjg(pot(l,1,k1))
  120          continue
! mode numbers kz = 0, nz/2
               pott(1,1,k1) = cmplx(real(pot(1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,1,k1) = cmplx(real(pot(l1,1,k1)),0.0)
               endif
! kx  = nx/2
               if (modesx.gt.nxh) then
                  do 130 l = 2, lmax
                  l1 = nz2 - l
                  pott(2*l-2,j1,k1) = conjg(pot(l1,1,k1))
                  pott(2*l-1,j1,k1) = pot(l1,1,k1)
  130             continue
! mode numbers kz = 0, nz/2
                  pott(1,j1,k1) = cmplx(aimag(pot(1,1,k1)),0.0)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(nz,j1,k1) = cmplx(aimag(pot(l1,1,k1)),0.0)
                  endif
               endif
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPWRMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,    &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex array pott and stores them into a packed complex
! array pot
! modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
! and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
! and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
! where js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
! and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! except kx = NX/2 is stored at location kxyp+1 when js=0,
! and ky = NY/2 is stored at location 1 when ks=nvp/2.
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! modeszd = first dimension of array pott,
! where modeszd  = min(2*modesz-1,nz)
! modesypd = third dimension of array pott, modesypd >= kyzp
! modesxpd = second dimension of array pott,
! modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
! in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, kstrt, nvpy, nvpz, nzv
      integer kxyp, kyzp, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp), pott(modeszd,modesxpd,modesypd)
! local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 100 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))&
     & then
         if (k1.gt.nyh) k1 = k1 - ny
         do 30 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,k) = pott(2*l-2,j,k)
            pot(l1,j,k) = pott(2*l-1,j,k)
   10       continue
            do 20 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,k) = zero
            pot(l1,j,k) = zero
   20       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            pot(1,j,k) = pott(1,j,k)
            pot(l1,j,k) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,k) = pott(nz,j,k)
            endif
         endif
   30    continue
         do 50 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 40 l = 1, nz
            pot(l,j,k) = zero
   40       continue
         endif
   50    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               pot(l,1,k) = pott(2*l-2,1,k)
               pot(l1,1,k) = pott(2*l-1,1,k)
   60          continue
               do 70 l = lmax+1, nzh
               l1 = nz2 - l
               pot(l,1,k) = zero
               pot(l1,1,k) = zero
   70          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               pot(1,1,k) = pott(1,1,k)
               pot(l1,1,k) = zero
               if (modesz.gt.nzh) then
                  pot(l1,1,k) = pott(nz,1,k)
               endif
! kx = nx/2
            else
               do 80 l = 1, nz
               pot(l,1,k) = zero
   80          continue
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  pot(l,1,k) = pott(2*l-2,1,k)
                  pot(l1,1,k) = pott(2*l-1,1,k)
   90             continue
! mode numbers kz = 0, nz/2
                  pot(1,1,k) = pott(1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k) = pott(nz,1,k)
                  endif
               endif
            endif
         endif
      endif
  100 continue
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 140 k = 1, kyzps
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 120 j = 1, kxyps
         if ((j+joff).gt.1) then
            do 110 l = 1, nz
            pot(l,j,k) = zero
  110       continue
         endif
  120    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 130 l = 1, nz
            pot(l,1,k) = zero
  130       continue
         endif
      endif
  140 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! ky = 0
      if (ks.eq.0) then
         do 170 j = 1, jmax
         if ((j+joff).gt.1) then
            do 150 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,1) = pott(2*l-2,j,1)
            pot(l1,j,1) = pott(2*l-1,j,1)
  150       continue
            do 160 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,1) = zero
            pot(l1,j,1) = zero
  160       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,j,1) = pott(1,j,1)
            pot(l1,j,1) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,1) = pott(nz,j,1)
            endif
         endif
  170    continue
         do 190 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 180 l = 1, nz
            pot(l,j,1) = zero
  180       continue
         endif
  190    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            do 200 l = 2, lmax
            l1 = nz2 - l
            pot(l,1,1) = pott(2*l-2,1,1)
            pot(l1,1,1) = zero
  200       continue
            do 210 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,1,1) = zero
            pot(l1,1,1) = zero
  210       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1) = cmplx(real(pott(1,1,1)),0.0)
            pot(l1,1,1) = zero
            if (modesz.gt.nzh) then
               pot(l1,1,1) = cmplx(real(pott(nz,1,1)),0.0)
            endif
! kx = nx/2
            if (modesx.gt.nxh) then
               do 220 l = 2, lmax
               l1 = nz2 - l
               pot(l1,1,1) = conjg(pott(2*l-2,j1,1))
  220          continue
! mode numbers kz = 0, nz/2
               pot(1,1,1) = cmplx(real(pot(1,1,1)),real(pott(1,j1,1)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,1) = cmplx(real(pot(l1,1,1)),                &
     &                                real(pott(nz,j1,1)))
               endif
            endif
         endif
      endif
! ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 240 j = 1, jmax
         if ((j+joff).gt.1) then
            do 230 l = 1, nz
            pot(l,j,k1) = zero
  230       continue
         endif
  240    continue
         do 260 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 250 l = 1, nz
            pot(l,j,k1) = zero
  250       continue
         endif
  260    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 270 l = 1, nz
            pot(l,1,k1) = zero
  270       continue
         endif
         if (modesy.gt.nyh) then
            do 290 j = 1, jmax
            if ((j+joff).gt.1) then
               do 280 l = 2, lmax
               l1 = nz2 - l
               pot(l,j,k1) = pott(2*l-2,j,k1)
               pot(l1,j,k1) = pott(2*l-1,j,k1)
  280          continue
! mode numbers kz = 0, nz/2
               pot(1,j,k1) = pott(1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,j,k1) = pott(nz,j,k1)
               endif
            endif
  290       continue
! mode numbers kx = 0, nx/2
            if (js.eq.0) then
! kx = 0
               do 300 l = 2, lmax
               pot(l,1,k1) = pott(2*l-2,1,k1)
  300          continue
! mode numbers kz = 0, nz/2
               pot(1,1,k1) = cmplx(real(pott(1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,k1) = cmplx(real(pott(nz,1,k1)),0.0)
               endif
! kx  = nx/2
               if (modesx.gt.nxh) then
                  do 310 l = 2, lmax
                  l1 = nz2 - l
                  pot(l1,1,k1) = conjg(pott(2*l-2,j1,k1))
  310             continue
! mode numbers kz = 0, nz/2
                  pot(1,1,k1) = cmplx(real(pot(1,1,k1)),                &
     &                                real(pott(1,j1,k1)))
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k1) = cmplx(real(pot(l1,1,k1)),           &
     &                                    real(pott(nz,j1,k1)))
                  endif
               endif
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRDVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! this subroutine extracts lowest order modes from packed complex vector
! array vpot and stores them into a location in an unpacked complex
! vector array vpott
! modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
! and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
! and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
! where js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
! and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! except kx = NX/2 is stored at location kxyp+1 when js=0,
! and ky = NY/2 is stored at location 1 when ks=nvp/2.
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! ndim = number of field arrays, must be >= 1
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! modeszd = second dimension of array vpott,
! where modeszd  = min(2*modesz-1,nz)
! modesypd = fourth dimension of array vpott, modesypd >= kyzp
! modesxpd = thid dimension of array vpott,
! modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
! in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, ndim, kstrt
      integer nvpy, nvpz, nzv, kxyp, kyzp, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp)
      dimension vpott(ndim,modeszd,modesxpd,modesypd)
! local data
      integer i, j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 110 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))&
     & then
         if (k1.gt.nyh) k1 = k1 - ny
         do 40 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpott(i,2*l-2,j,k) = vpot(i,l,j,k)
            vpott(i,2*l-1,j,k) = vpot(i,l1,j,k)
   10       continue
   20       continue
! mode numbers kz = 0, nz/2
            do 30 i = 1, ndim
            vpott(i,1,j,k) = vpot(i,1,j,k)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,j,k) = vpot(i,l1,j,k)
            endif
   30       continue
         endif
   40    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               do 50 i = 1, ndim
               vpott(i,2*l-2,1,k) = vpot(i,l,1,k)
               vpott(i,2*l-1,1,k) = vpot(i,l1,1,k)
   50          continue
   60          continue
! mode numbers kz = 0, nz/2
               do 70 i = 1, ndim
               vpott(i,1,1,k) = vpot(i,1,1,k)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,1,k) = vpot(i,l1,1,k)
               endif
   70          continue
! kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  do 80 i = 1, ndim
                  vpott(i,2*l-2,1,k) = vpot(i,l,1,k)
                  vpott(i,2*l-1,1,k) = vpot(i,l1,1,k)
   80             continue
   90             continue
! mode numbers kz = 0, nz/2
                  do 100 i = 1, ndim
                  vpott(i,1,1,k) = vpot(i,1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(i,nz,1,k) = vpot(i,l1,1,k)
                  endif
  100             continue
               endif
            endif
         endif
      endif
  110 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! ky = 0
      if (ks.eq.0) then
         do 150 j = 1, jmax
         if ((j+joff).gt.1) then
            do 130 l = 2, lmax
            l1 = nz2 - l
            do 120 i = 1, ndim
            vpott(i,2*l-2,j,1) = vpot(i,l,j,1)
            vpott(i,2*l-1,j,1) = vpot(i,l1,j,1)
  120       continue
  130       continue
! mode numbers kz = 0, nz/2
            do 140 i = 1, ndim
            vpott(i,1,j,1) = vpot(i,1,j,1)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,j,1) = vpot(i,l1,j,1)
            endif
  140       continue
         endif
  150    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            do 170 l = 2, lmax
            do 160 i = 1, ndim
            vpott(i,2*l-2,1,1) = vpot(i,l,1,1)
            vpott(i,2*l-1,1,1) = conjg(vpot(i,l,1,1))
  160       continue
  170       continue
! mode numbers kz = 0, nz/2
            do 180 i = 1, ndim
            vpott(i,1,1,1) = cmplx(real(vpot(i,1,1,1)),0.0)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,1,1) = cmplx(real(vpot(i,l1,1,1)),0.0)
            endif
  180       continue
! kx = nx/2
            if (modesx.gt.nxh) then
               do 200 l = 2, lmax
               l1 = nz2 - l
               do 190 i = 1, ndim
               vpott(i,2*l-2,j1,1) = conjg(vpot(i,l1,1,1))
               vpott(i,2*l-1,j1,1) = vpot(i,l1,1,1)
  190          continue
  200          continue
! mode numbers kz = 0, nz/2
               do 210 i = 1, ndim
               vpott(i,1,j1,1) = cmplx(aimag(vpot(i,1,1,1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,j1,1) = cmplx(aimag(vpot(i,l1,1,1)),0.0)
               endif
  210          continue
            endif
         endif
      endif
! ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         if (modesy.gt.nyh) then
            do 250 j = 1, jmax
            if ((j+joff).gt.1) then
               do 230 l = 2, lmax
               l1 = nz2 - l
               do 220 i = 1, ndim
               vpott(i,2*l-2,j,k1) = vpot(i,l,j,k1)
               vpott(i,2*l-1,j,k1) = vpot(i,l1,j,k1)
  220          continue
  230          continue
! mode numbers kz = 0, nz/2
               do 240 i = 1, ndim
               vpott(i,1,j,k1) = vpot(i,1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,j,k1) = vpot(i,l1,j,k1)
               endif
  240          continue
            endif
  250       continue
! mode numbers kx = 0, nx/2
            if (js.eq.0) then
! kx = 0
               do 270 l = 2, lmax
               do 260 i = 1, ndim
               vpott(i,2*l-2,1,k1) = vpot(i,l,1,k1)
               vpott(i,2*l-1,1,k1) = conjg(vpot(i,l,1,k1))
  260          continue
  270          continue
! mode numbers kz = 0, nz/2
               do 280 i = 1, ndim
               vpott(i,1,1,k1) = cmplx(real(vpot(i,1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,1,k1) = cmplx(real(vpot(i,l1,1,k1)),0.0)
               endif
  280          continue
! kx  = nx/2
               if (modesx.gt.nxh) then
                  do 300 l = 2, lmax
                  l1 = nz2 - l
                  do 290 i = 1, ndim
                  vpott(i,2*l-2,j1,k1) = conjg(vpot(i,l1,1,k1))
                  vpott(i,2*l-1,j1,k1) = vpot(i,l1,1,k1)
  290             continue
  300             continue
! mode numbers kz = 0, nz/2
                  do 310 i = 1, ndim
                  vpott(i,1,j1,k1) = cmplx(aimag(vpot(i,1,1,k1)),0.0)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(i,nz,j1,k1) = cmplx(aimag(vpot(i,l1,1,k1)),  &
     &                                         0.0)
                  endif
  310             continue
               endif
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPWRVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex vector array vpott and stores them into packed
! complex vector array vpot
! modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
! and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
! and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
! where js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
! and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! except kx = NX/2 is stored at location kxyp+1 when js=0,
! and ky = NY/2 is stored at location 1 when ks=nvp/2.
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! ndim = number of field arrays, must be >= 1
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nzv = first dimension of field arrays, must be >= nz
! kxyp/kyzp = number of complex grids in each field partition in
! x/y direction
! modeszd = second dimension of array vpott,
! where modeszd  = min(2*modesz-1,nz)
! modesypd = fourth dimension of array vpott, modesypd >= kyzp
! modesxpd = thid dimension of array vpott,
! modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
! in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, ndim, kstrt
      integer nvpy, nvpz, nzv, kxyp, kyzp, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp)
      dimension vpott(ndim,modeszd,modesxpd,modesypd)
! local data
      integer i, j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
! find processor id and offsets in y/z
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.(nvpy*nvpz)) return
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 200 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))&
     & then
         if (k1.gt.nyh) k1 = k1 - ny
         do 60 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpot(i,l,j,k) = vpott(i,2*l-2,j,k)
            vpot(i,l1,j,k) = vpott(i,2*l-1,j,k)
   10       continue
   20       continue
            do 40 l = lmax+1, nzh
            l1 = nz2 - l
            do 30 i = 1, ndim
            vpot(i,l,j,k) = zero
            vpot(i,l1,j,k) = zero
   30       continue
   40       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            do 50 i = 1, ndim
            vpot(i,1,j,k) = vpott(i,1,j,k)
            vpot(i,l1,j,k) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,k) = vpott(i,nz,j,k)
            endif
   50       continue
         endif
   60    continue
         do 90 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 80 l = 1, nz
            do 70 i = 1, ndim
            vpot(i,l,j,k) = zero
   70       continue
   80       continue
         endif
   90    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            if (k1.gt.0) then
               do 110 l = 2, lmax
               l1 = nz2 - l
               do 100 i = 1, ndim
               vpot(i,l,1,k) = vpott(i,2*l-2,1,k)
               vpot(i,l1,1,k) = vpott(i,2*l-1,1,k)
  100          continue
  110          continue
               do 130 l = lmax+1, nzh
               l1 = nz2 - l
               do 120 i = 1, ndim
               vpot(i,l,1,k) = zero
               vpot(i,l1,1,k) = zero
  120          continue
  130          continue
! mode numbers kz = 0, nz/2
               l1 = nzh + 1
               do 140 i = 1, ndim
               vpot(i,1,1,k) = vpott(i,1,1,k)
               vpot(i,l1,1,k) = zero
               if (modesz.gt.nzh) then
                  vpot(i,l1,1,k) = vpott(i,nz,1,k)
               endif
  140          continue
! kx = nx/2
            else
               do 160 l = 1, nz
               do 150 i = 1, ndim
               vpot(i,l,1,k) = zero
  150          continue
  160          continue
               if (modesx.gt.nxh) then
                  do 180 l = 2, lmax
                  l1 = nz2 - l
                  do 170 i = 1, ndim
                  vpot(i,l,1,k) = vpott(i,2*l-2,1,k)
                  vpot(i,l1,1,k) = vpott(i,2*l-1,1,k)
  170             continue
  180             continue
! mode numbers kz = 0, nz/2
                  do 190 i = 1, ndim
                  vpot(i,1,1,k) = vpott(i,1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k) = vpott(i,nz,1,k)
                  endif
  190             continue
               endif
            endif
         endif
      endif
  200 continue
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 260 k = 1, kyzps
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 230 j = 1, kxyps
         if ((j+joff).gt.1) then
            do 220 l = 1, nz
            do 210 i = 1, ndim
            vpot(i,l,j,k) = zero
  210       continue
  220       continue
         endif
  230    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 250 l = 1, nz
            do 240 i = 1, ndim
            vpot(i,l,1,k) = zero
  240       continue
  250       continue
         endif
      endif
  260 continue
!$OMP END PARALLEL DO
! mode numbers ky = 0, ny/2
! ky = 0
      if (ks.eq.0) then
         do 320 j = 1, jmax
         if ((j+joff).gt.1) then
            do 280 l = 2, lmax
            l1 = nz2 - l
            do 270 i = 1, ndim
            vpot(i,l,j,1) = vpott(i,2*l-2,j,1)
            vpot(i,l1,j,1) = vpott(i,2*l-1,j,1)
  270       continue
  280       continue
            do 300 l = lmax+1, nzh
            l1 = nz2 - l
            do 290 i = 1, ndim
            vpot(i,l,j,1) = zero
            vpot(i,l1,j,1) = zero
  290       continue
  300       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 310 i = 1, ndim
            vpot(i,1,j,1) = vpott(i,1,j,1)
            vpot(i,l1,j,1) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,1) = vpott(i,nz,j,1)
            endif
  310       continue
         endif
  320    continue
         do 350 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 340 l = 1, nz
            do 330 i = 1, ndim
            vpot(i,l,j,1) = zero
  330       continue
  340       continue
         endif
  350    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
! kx = 0
            do 370 l = 2, lmax
            l1 = nz2 - l
            do 360 i = 1, ndim
            vpot(i,l,1,1) = vpott(i,2*l-2,1,1)
            vpot(i,l1,1,1) = zero
  360       continue
  370       continue
            do 390 l = lmax+1, nzh
            l1 = nz2 - l
            do 380 i = 1, ndim
            vpot(i,l,1,1) = zero
            vpot(i,l1,1,1) = zero
  380       continue
  390       continue
! mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 400 i = 1, ndim
            vpot(i,1,1,1) = cmplx(real(vpott(i,1,1,1)),0.0)
            vpot(i,l1,1,1) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,1,1) = cmplx(real(vpott(i,nz,1,1)),0.0)
            endif
  400       continue
! kx = nx/2
            if (modesx.gt.nxh) then
               do 420 l = 2, lmax
               l1 = nz2 - l
               do 410 i = 1, ndim
               vpot(i,l1,1,1) = conjg(vpott(i,2*l-2,j1,1))
  410          continue
  420          continue
! mode numbers kz = 0, nz/2
               do 430 i = 1, ndim
               vpot(i,1,1,1) = cmplx(real(vpot(i,1,1,1)),               &
     &                               real(vpott(i,1,j1,1)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,1) = cmplx(real(vpot(i,l1,1,1)),          &
     &                                   real(vpott(i,nz,j1,1)))
               endif
  430          continue
            endif
         endif
      endif
! ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 460 j = 1, jmax
         if ((j+joff).gt.1) then
            do 450 l = 1, nz
            do 440 i = 1, ndim
            vpot(i,l,j,k1) = zero
  440       continue
  450       continue
         endif
  460    continue
         do 490 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 480 l = 1, nz
            do 470 i = 1, ndim
            vpot(i,l,j,k1) = zero
  470       continue
  480       continue
         endif
  490    continue
! mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 510 l = 1, nz
            do 500 i = 1, ndim
            vpot(i,l,1,k1) = zero
  500       continue
  510       continue
         endif
         if (modesy.gt.nyh) then
            do 550 j = 1, jmax
            if ((j+joff).gt.1) then
               do 530 l = 2, lmax
               l1 = nz2 - l
               do 520 i = 1, ndim
               vpot(i,l,j,k1) = vpott(i,2*l-2,j,k1)
               vpot(i,l1,j,k1) = vpott(i,2*l-1,j,k1)
  520          continue
  530          continue
! mode numbers kz = 0, nz/2
               do 540 i = 1, ndim
               vpot(i,1,j,k1) = vpott(i,1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,j,k1) = vpott(i,nz,j,k1)
               endif
  540          continue
            endif
  550       continue
! mode numbers kx = 0, nx/2
            if (js.eq.0) then
! kx = 0
               do 570 l = 2, lmax
               do 560 i = 1, ndim
               vpot(i,l,1,k1) = vpott(i,2*l-2,1,k1)
  560          continue
  570          continue
! mode numbers kz = 0, nz/2
               do 580 i = 1, ndim
               vpot(i,1,1,k1) = cmplx(real(vpott(i,1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,k1) = cmplx(real(vpott(i,nz,1,k1)),0.0)
               endif
  580          continue
! kx  = nx/2
               if (modesx.gt.nxh) then
                  do 600 l = 2, lmax
                  l1 = nz2 - l
                  do 590 i = 1, ndim
                  vpot(i,l1,1,k1) = conjg(vpott(i,2*l-2,j1,k1))
  590             continue
  600             continue
! mode numbers kz = 0, nz/2
                  do 610 i = 1, ndim
                  vpot(i,1,1,k1) = cmplx(real(vpot(i,1,1,k1)),          &
     &                                   real(vpott(i,1,j1,k1)))
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k1) = cmplx(real(vpot(i,l1,1,k1)),     &
     &                                       real(vpott(i,nz,j1,k1)))
                  endif
  610             continue
               endif
            endif
         endif
      endif
      return
      end
