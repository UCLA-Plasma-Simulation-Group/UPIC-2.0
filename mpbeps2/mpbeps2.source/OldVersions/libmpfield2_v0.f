!-----------------------------------------------------------------------
! Fortran library for spectral field solvers
! 2D MPI/OpenMP PIC Code:
! MPPOIS22 solves 2d poisson's equation for smoothed electric field
! MPPOIS23 solves 2-1/2d poisson's equation for smoothed electric field
! MPPADDQEI2 adds electron and ion densities
! MPPCUPERP2 calculates the transverse current in fourier space
! MIPPBPOISP23 solves 2-1/2d poisson's equation for unsmoothed magnetic
!              field
! MPPMAXWEL2 solves 2-1/2d maxwell's equation for unsmoothed transverse
!            electric and magnetic fields
! MPPEMFIELD2 adds and smooths or copies and smooths complex vector
!             fields in fourier space
! MPPADDCUEI2 adds electron and ion current densities
! MPPADDAMUI2 adds electron and ion momentum flux densities
! PPBADDEXT2 adds constant to magnetic field for 2-1/2d code
! PPADDVRFIELD2 calculates a = b + c
! MPPBBPOISP23 solves 2-1/2d poisson's equation in fourier space for
!              smoothed magnetic field
! MPPDCUPERP23 calculates transverse part of the derivative of the
!              current density from the momentum flux
! MPPADCUPERP23 Calculates transverse part of the derivative of the
!               current density from the momentum flux and acceleration
!               density
! MPPEPOISP23 solves 2-1/2d poisson's equation in fourier space for
!             smoothed or unsmoothed transverse electric field
! MPPOTP2 solves 2d poisson's equation for potential
! MPPELFIELD22 solves 2d poisson's equation for unsmoothed electric
!              field
! MPPELFIELD23 solves 2-1/2d poisson's equation for unsmoothed electric
!              field
! MPPDIVF2 calculates the divergence in fourier space
! MPPGRADF2 calculates the gradient in fourier space
! MPPCURLF2 calculates the curl in fourier space
! MPPAVPOT23 calculates 2-1/2d vector potential from magnetic field
! MCUAVE23 averages current in fourier space for 2-1/2d code
! MPPAVRPOT23 solves 2-1/2d poisson's equation for the radiative part of
!             the vector potential
! MPPAPOTP23 solves 2-1/2d poisson's equation for vector potential
! MPPETFIELD23 solves 2-1/2d poisson's equation in fourier space for
!              unsmoothed transverse electric field
! MPPSMOOTH2 provides a 2d scalar smoothing function
! MPPSMOOTH23 provides a 2-1/2d vector smoothing function
! PPRDMODES2 extracts lowest order scalar modes from packed array
!            stores them into a location in an unpacked array
! PPWRMODES2 extracts lowest order scalar modes from a location in an
!            unpacked array and stores them into a packed array
! PPRDVMODES2 extracts lowest order vector modes from packed array
!             stores them into a location in an unpacked array
! PPWRVMODES2 extracts lowest order vector modes from a location in an
!             unpacked array and stores them into a packed array
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 26, 2018
!-----------------------------------------------------------------------
      subroutine MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,&
     &kxp,nyhd)
! this subroutine solves 2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions, for distributed data.
! for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
! output: ffc
! for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd,
! output: fxy,we
! approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the equation used is:
! fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
! fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! q(k,j) = complex charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of complex force/charge,
! fxy(2,k,j) = y component of complex force/charge,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! if isign = 0, form factor array is prepared
! if isign is not equal to 0, force/charge is calculated.
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(2,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
! prepare form factor array
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffc(k,j) = cmplx(affp,1.0)
      else
         ffc(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 70
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,at1,at2,at3,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))*aimag(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   40    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))*aimag(ffc(1,j))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))*aimag(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   60    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
      endif
      sum1 = sum1 + wp
   70 continue
      we = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,&
     &kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions.  Zeros out z component.
! for distributed data.
! for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
! output: ffc
! for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd,
! output: fxy,we
! approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the equation used is:
! fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
! fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! q(k,j) = complex charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of complex force/charge,
! fxy(2,k,j) = y component of complex force/charge,
! fxy(3,k,j) = zero,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! if isign = 0, form factor array is prepared
! if isign is not equal to 0, force/charge is calculated.
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
! prepare form factor array
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffc(k,j) = cmplx(affp,1.0)
      else
         ffc(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 70
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2-
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,at1,at2,at3,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))*aimag(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(3,k,j) = zero
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   40    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))*aimag(ffc(1,j))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(3,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))*aimag(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(3,k,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   60    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(3,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   70 continue
      we = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADDQEI2(qe,qi,nyp,nx,nxe,nypmx)
! adds electron and ion densities
! assumes guard cells have already been added
! qe/qi = charge density for electrons/ions
! nyp = number of primary gridpoints in y in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
      implicit none
      integer nx, nxe, nypmx
      integer nyp
      real qe, qi
      dimension qe(nxe,nypmx), qi(nxe,nypmx)
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      qe(j,k) = qe(j,k) + qi(j,k)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
! this subroutine calculates the transverse current in fourier space
! input: all, output: cu
! approximate flop count is: 36*nxc*nyc
! and nxc*nyc divides
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the transverse current is calculated using the equation:
! cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
! and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
! cu(i,k,j) = i-th component of complex current density and
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex cu
      dimension cu(3,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate transverse part of current
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dkx2,dky,at1,zt1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j) + dky*cu(2,k,j))
         cu(1,k,j) = cu(1,k,j) - dkx*zt1
         cu(2,k,j) = cu(2,k,j) - dky*zt1
         zt1 = at1*(dkx*cu(1,k1,j) - dky*cu(2,k1,j))
         cu(1,k1,j) = cu(1,k1,j) - dkx*zt1
         cu(2,k1,j) = cu(2,k1,j) + dky*zt1
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j) = zero
         cu(1,k1,j) = zero
         cu(2,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1) = zero
         cu(1,k1,1) = zero
         cu(2,k1,1) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1) = zero
         cu(2,1,1) = zero
         cu(1,k1,1) = zero
         cu(2,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! magnetic field with periodic boundary conditions for distributed data.
! input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
! approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
! by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
! bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
! bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
! cu(i,k,j) = i-th component of complex current density and
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! ci = reciprical of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
! affp = normalization constant = nx*ny/np, where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nyv,kxp), bxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real ci2, dnx, dny, dkx, dky, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,at1,at2,at3,zt1,zt2,zt3,wp)    &
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffc(k,j))
         zt1 = cmplx(-aimag(cu(3,k,j)),real(cu(3,k,j)))
         zt2 = cmplx(-aimag(cu(2,k,j)),real(cu(2,k,j)))
         zt3 = cmplx(-aimag(cu(1,k,j)),real(cu(1,k,j)))
         bxy(1,k,j) = at2*zt1
         bxy(2,k,j) = -at3*zt1
         bxy(3,k,j) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(cu(3,k1,j)),real(cu(3,k1,j)))
         zt2 = cmplx(-aimag(cu(2,k1,j)),real(cu(2,k1,j)))
         zt3 = cmplx(-aimag(cu(1,k1,j)),real(cu(1,k1,j)))
         bxy(1,k1,j) = -at2*zt1
         bxy(2,k1,j) = -at3*zt1
         bxy(3,k1,j) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j)*conjg(cu(1,k,j))                      &
     &   + cu(2,k,j)*conjg(cu(2,k,j)) + cu(3,k,j)*conjg(cu(3,k,j))      &
     &   + cu(1,k1,j)*conjg(cu(1,k1,j)) + cu(2,k1,j)*conjg(cu(2,k1,j))  &
     &   + cu(3,k1,j)*conjg(cu(3,k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j))
         at2 = dkx*at1
         at1 = at1*aimag(ffc(1,j))
         zt1 = cmplx(-aimag(cu(3,1,j)),real(cu(3,1,j)))
         zt2 = cmplx(-aimag(cu(2,1,j)),real(cu(2,1,j)))
         bxy(1,1,j) = zero
         bxy(2,1,j) = -at2*zt1
         bxy(3,1,j) = at2*zt2
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         wp = wp + at1*(cu(1,1,j)*conjg(cu(1,1,j))                      &
     &   + cu(2,1,j)*conjg(cu(2,1,j)) + cu(3,1,j)*conjg(cu(3,1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,1))
         at2 = dky*at1
         at1 = at1*aimag(ffc(k,1))
         zt1 = cmplx(-aimag(cu(3,k,1)),real(cu(3,k,1)))
         zt2 = cmplx(-aimag(cu(1,k,1)),real(cu(1,k,1)))
         bxy(1,k,1) = at2*zt1
         bxy(2,k,1) = zero
         bxy(3,k,1) = -at2*zt2
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         wp = wp + at1*(cu(1,k,1)*conjg(cu(1,k,1))                      &
     &   + cu(2,k,1)*conjg(cu(2,k,1)) + cu(3,k,1)*conjg(cu(3,k,1)))
   30    continue
         k1 = nyh + 1
         bxy(1,1,1) = zero
         bxy(2,1,1) = zero
         bxy(3,1,1) = zero
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,&
     &nyv,kxp,nyhd)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse electric and magnetic fields with periodic boundary
! conditions.
! input: all, output: wf, wm, exy, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
! by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
! bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
! the electric field is then updated a whole step using the equations:
! ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
!                       - affp*dt*cux(kx,ky)*s(kx,ky)
! ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
!                       - affp*dt*cuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*dt*cuz(kx,ky)*s(kx,ky)
! the magnetic field is finally updated the remaining half step with
! the new electric field and the previous magnetic field equations.
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! and similarly for bx, by, bz.
! cu(i,k,j) = i-th component of complex current density and
! exy(i,k,j) = i-th component of complex electric field,
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffc(k,j)) = finite-size particle shape factor s
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
! magnetic field energy is also calculated, using
! wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(3,nyv,kxp), bxy(3,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dth, c2, cdt, adt, anorm, dkx, dky, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws, sum1, sum2
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum1 = 0.0d0
      sum2 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,afdt,zt1,zt2,zt3,zt4,zt5,zt6,  &
!$OMP& zt7,zt8,zt9,ws,wp) REDUCTION(+:sum1,sum2)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      ws = 0.0d0
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffc(k,j))
! update magnetic field half time step, ky > 0
         zt1 = cmplx(-aimag(exy(3,k,j)),real(exy(3,k,j)))
         zt2 = cmplx(-aimag(exy(2,k,j)),real(exy(2,k,j)))
         zt3 = cmplx(-aimag(exy(1,k,j)),real(exy(1,k,j)))
         zt4 = bxy(1,k,j) - dth*(dky*zt1)
         zt5 = bxy(2,k,j) + dth*(dkx*zt1)
         zt6 = bxy(3,k,j) - dth*(dkx*zt2 - dky*zt3)
! update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,j) + cdt*(dky*zt1) - afdt*cu(1,k,j)
         zt8 = exy(2,k,j) - cdt*(dkx*zt1) - afdt*cu(2,k,j)
         zt9 = exy(3,k,j) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,k,j)
! update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,j) = zt7
         exy(2,k,j) = zt8
         exy(3,k,j) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + zt9*conjg(zt9))
         zt4 = zt4 - dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
         bxy(1,k,j) = zt4
         bxy(2,k,j) = zt5
         bxy(3,k,j) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + zt6*conjg(zt6))
! update magnetic field half time step, ky < 0
         zt1 = cmplx(-aimag(exy(3,k1,j)),real(exy(3,k1,j)))
         zt2 = cmplx(-aimag(exy(2,k1,j)),real(exy(2,k1,j)))
         zt3 = cmplx(-aimag(exy(1,k1,j)),real(exy(1,k1,j)))
         zt4 = bxy(1,k1,j) + dth*(dky*zt1)
         zt5 = bxy(2,k1,j) + dth*(dkx*zt1)
         zt6 = bxy(3,k1,j) - dth*(dkx*zt2 + dky*zt3)
! update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k1,j) - cdt*(dky*zt1) - afdt*cu(1,k1,j)
         zt8 = exy(2,k1,j) - cdt*(dkx*zt1) - afdt*cu(2,k1,j)
         zt9 = exy(3,k1,j) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,k1,j)
! update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k1,j) = zt7
         exy(2,k1,j) = zt8
         exy(3,k1,j) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + zt9*conjg(zt9))
         zt4 = zt4 + dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
         bxy(1,k1,j) = zt4
         bxy(2,k1,j) = zt5
         bxy(3,k1,j) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + zt6*conjg(zt6))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         afdt = adt*aimag(ffc(1,j))
! update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,1,j)),real(exy(3,1,j)))
         zt2 = cmplx(-aimag(exy(2,1,j)),real(exy(2,1,j)))
         zt5 = bxy(2,1,j) + dth*(dkx*zt1)
         zt6 = bxy(3,1,j) - dth*(dkx*zt2)
! update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt8 = exy(2,1,j) - cdt*(dkx*zt1) - afdt*cu(2,1,j)
         zt9 = exy(3,1,j) + cdt*(dkx*zt2) - afdt*cu(3,1,j)
! update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         exy(1,1,j) = zero
         exy(2,1,j) = zt8
         exy(3,1,j) = zt9
         ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2)
         bxy(1,1,j) = zero
         bxy(2,1,j) = zt5
         bxy(3,1,j) = zt6
         wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffc(k,1))
! update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,k,1)),real(exy(3,k,1)))
         zt3 = cmplx(-aimag(exy(1,k,1)),real(exy(1,k,1)))
         zt4 = bxy(1,k,1) - dth*(dky*zt1)
         zt6 = bxy(3,k,1) + dth*(dky*zt3)
! update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,1) + cdt*(dky*zt1) - afdt*cu(1,k,1)
         zt9 = exy(3,k,1) - cdt*(dky*zt3) - afdt*cu(3,k,1)
! update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,1) = zt7
         exy(2,k,1) = zero
         exy(3,k,1) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
         zt4 = zt4 - dth*(dky*zt1)
         zt6 = zt6 + dth*(dky*zt3)
         bxy(1,k,1) = zt4
         bxy(2,k,1) = zero
         bxy(3,k,1) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         bxy(1,1,1) = zero
         bxy(2,1,1) = zero
         bxy(3,1,1) = zero
         exy(1,1,1) = zero
         exy(2,1,1) = zero
         exy(3,1,1) = zero
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   40 continue
      wf = real(nx)*real(ny)*sum1
      wm = real(nx)*real(ny)*c2*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine either adds complex vector fields if isign > 0
! or copies complex vector fields if isign < 0
! includes additional smoothing
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nyv,kxp), exy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer i, nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      if (kstrt.gt.nxh) return
! add the fields
      if (isign.gt.0) then
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1,at1)
         do 40 j = 1, kxps
         do 20 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         do 10 i = 1, 3
         fxy(i,k,j) = fxy(i,k,j) + exy(i,k,j)*at1
         fxy(i,k1,j) = fxy(i,k1,j) + exy(i,k1,j)*at1
   10    continue
   20    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         do 30 i = 1, 3
         fxy(i,1,j) = fxy(i,1,j) + exy(i,1,j)*at1
         fxy(i,k1,j) = fxy(i,k1,j) + exy(i,k1,j)*at1
   30    continue
   40    continue
!$OMP END PARALLEL DO
! copy the fields
      else if (isign.lt.0) then
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1,at1)
         do 80 j = 1, kxps
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         do 50 i = 1, 3
         fxy(i,k,j) = exy(i,k,j)*at1
         fxy(i,k1,j) = exy(i,k1,j)*at1
   50    continue
   60    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         do 70 i = 1, 3
         fxy(i,1,j) = exy(i,1,j)*at1
         fxy(i,k1,j) = exy(i,k1,j)*at1
   70    continue
   80    continue
!$OMP END PARALLEL DO
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADDCUEI23(cue,cui,nyp,nx,nxe,nypmx)
! adds electron and ion current densities
! assumes guard cells have already been added
! cue/cui = current density for electrons/ions
! nyp = number of primary gridpoints in y in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
      implicit none
      integer nx, nxe, nypmx
      integer nyp
      real cue, cui
      dimension cue(3,nxe,nypmx), cui(3,nxe,nypmx)
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      cue(1,j,k) = cue(1,j,k) + cui(1,j,k)
      cue(2,j,k) = cue(2,j,k) + cui(2,j,k)
      cue(3,j,k) = cue(3,j,k) + cui(3,j,k)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADDAMUI23(amu,amui,nyp,nx,nxe,nypmx)
! adds electron and ion momentum flux densities
! assumes guard cells have already been added
! amu/amui = momentum flux density for electrons/ions
! nyp = number of primary gridpoints in y in particle partition m
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
      implicit none
      integer nx, nxe, nypmx
      integer nyp
      real amu, amui
      dimension amu(4,nxe,nypmx), amui(4,nxe,nypmx)
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      amu(1,j,k) = amu(1,j,k) + amui(1,j,k)
      amu(2,j,k) = amu(2,j,k) + amui(2,j,k)
      amu(3,j,k) = amu(3,j,k) + amui(3,j,k)
      amu(4,j,k) = amu(4,j,k) + amui(4,j,k)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx)
! adds constant to magnetic field for 2-1/2d code
! bxy = magnetic field
! nyp = number of primary (complete) gridpoints in particle partition
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z
! nx = system length in x direction
! nxe = first dimension of field array, must be >= nx
! nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer nyp, nx, nxe, nypmx
      real omx, omy, omz
      real bxy
      dimension bxy(3,nxe,nypmx)
! local data
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 k = 1, nyp
      do 10 j = 1, nx
      bxy(1,j,k) = bxy(1,j,k) + omx
      bxy(2,j,k) = bxy(2,j,k) + omy
      bxy(3,j,k) = bxy(3,j,k) + omz
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPADDVRFIELD2(a,b,c,ndim,nxe,nypmx)
! this subroutine calculates a = b + c for distributed real vector field
      implicit none
      integer ndim, nxe, nypmx
      real a, b, c
      dimension a(ndim,nxe,nypmx)
      dimension b(ndim,nxe,nypmx), c(ndim,nxe,nypmx)
! local data
      integer i, j, k
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 30 k = 1, nypmx
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k) = b(i,j,k) + c(i,j,k)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPBBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! magnetic field (or convolution of magnetic field over particle shape)
! with periodic boundary conditions for distributed data.
! input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,nyhd, output: bxy,wm
! approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
! by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
! bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
!             s(kx,ky),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
! bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
! cu(i,k,j) = i-th component of complex current density and
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*|cu(kx,ky)*s(kx,ky|**2)
! affp = normalization constant = nx*ny/np, where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nyv,kxp), bxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real ci2, dnx, dny, dkx, dky, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,at1,at2,at3,zt1,zt2,zt3,wp)    &
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,j))*aimag(ffc(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,k,j)),real(cu(3,k,j)))
         zt2 = cmplx(-aimag(cu(2,k,j)),real(cu(2,k,j)))
         zt3 = cmplx(-aimag(cu(1,k,j)),real(cu(1,k,j)))
         bxy(1,k,j) = at2*zt1
         bxy(2,k,j) = -at3*zt1
         bxy(3,k,j) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(cu(3,k1,j)),real(cu(3,k1,j)))
         zt2 = cmplx(-aimag(cu(2,k1,j)),real(cu(2,k1,j)))
         zt3 = cmplx(-aimag(cu(1,k1,j)),real(cu(1,k1,j)))
         bxy(1,k1,j) = -at2*zt1
         bxy(2,k1,j) = -at3*zt1
         bxy(3,k1,j) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j)*conjg(cu(1,k,j))                      &
     &   + cu(2,k,j)*conjg(cu(2,k,j)) + cu(3,k,j)*conjg(cu(3,k,j))      &
     &   + cu(1,k1,j)*conjg(cu(1,k1,j)) + cu(2,k1,j)*conjg(cu(2,k1,j))  &
     &   + cu(3,k1,j)*conjg(cu(3,k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j))*aimag(ffc(1,j))
         at2 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,1,j)),real(cu(3,1,j)))
         zt2 = cmplx(-aimag(cu(2,1,j)),real(cu(2,1,j)))
         bxy(1,1,j) = zero
         bxy(2,1,j) = -at2*zt1
         bxy(3,1,j) = at2*zt2
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         wp = wp + at1*(cu(1,1,j)*conjg(cu(1,1,j))                      &
     &   + cu(2,1,j)*conjg(cu(2,1,j)) + cu(3,1,j)*conjg(cu(3,1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,1))*aimag(ffc(k,1))
         at2 = dky*at1
         zt1 = cmplx(-aimag(cu(3,k,1)),real(cu(3,k,1)))
         zt2 = cmplx(-aimag(cu(1,k,1)),real(cu(1,k,1)))
         bxy(1,k,1) = at2*zt1
         bxy(2,k,1) = zero
         bxy(3,k,1) = -at2*zt2
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         wp = wp + at1*(cu(1,k,1)*conjg(cu(1,k,1))                      &
     &   + cu(2,k,1)*conjg(cu(2,k,1)) + cu(3,k,1)*conjg(cu(3,k,1)))
   30    continue
         k1 = nyh + 1
         bxy(1,1,1) = zero
         bxy(2,1,1) = zero
         bxy(3,1,1) = zero
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux
! in 2-1/2d with periodic boundary conditions.
! input: all, output: dcu
! approximate flop count is: 45*nxc*nyc
! and nxc*nyc divides
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
! dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
! dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
! the transverse part is calculated using the equation:
! dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!               (kx*kx+ky*ky)
! dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!               (kx*kx+ky*ky)
! on output:
! dcu(i,k,j) = i-th component of transverse part of complex derivative
! of current for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! amu(1,k,j) = xx-yy component of complex momentum flux
! amu(2,k,j) = xy component of complex momentum flux
! amu(3,k,j) = zx component of complex momentum flux
! amu(4,k,j) = zy component of complex momentum flux
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex dcu, amu
      dimension dcu(3,nyv,kxp), amu(4,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate transverse part of current
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,k1,dkx,dkx2,dky,dky2,dkxy,dkxy2,at1,zt1,zt2,zt3)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j)),-real(amu(1,k,j)))
         zt2 = cmplx(aimag(amu(2,k,j)),-real(amu(2,k,j)))
         zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
         dcu(1,k,j) = dky*zt3
         dcu(2,k,j) = -dkx*zt3
         zt1 = cmplx(aimag(amu(3,k,j)),-real(amu(3,k,j)))
         zt2 = cmplx(aimag(amu(4,k,j)),-real(amu(4,k,j)))
         dcu(3,k,j) = dkx*zt1 + dky*zt2
         zt1 = cmplx(aimag(amu(1,k1,j)),-real(amu(1,k1,j)))
         zt2 = cmplx(aimag(amu(2,k1,j)),-real(amu(2,k1,j)))
         zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
         dcu(1,k1,j) = dky*zt3
         dcu(2,k1,j) = dkx*zt3
         zt1 = cmplx(aimag(amu(3,k1,j)),-real(amu(3,k1,j)))
         zt2 = cmplx(aimag(amu(4,k1,j)),-real(amu(4,k1,j)))
         dcu(3,k1,j) = dkx*zt1 - dky*zt2
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j)),-real(amu(2,1,j)))
         dcu(1,1,j) = zero
         dcu(2,1,j) = dkx*zt2
         zt1 = cmplx(aimag(amu(3,1,j)),-real(amu(3,1,j)))
         dcu(3,1,j) = dkx*zt1
         dcu(1,k1,j) = zero
         dcu(2,k1,j) = zero
         dcu(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1)),-real(amu(2,k,1)))
         dcu(1,k,1) = dky*zt2
         dcu(2,k,1) = zero
         zt2 = cmplx(aimag(amu(4,k,1)),-real(amu(4,k,1)))
         dcu(3,k,1) = dky*zt2
         dcu(1,k1,1) = zero
         dcu(2,k1,1) = zero
         dcu(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1) = zero
         dcu(2,1,1) = zero
         dcu(3,1,1) = zero
         dcu(1,k1,1) = zero
         dcu(2,k1,1) = zero
         dcu(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux and acceleration density
! in 2-1/2d with periodic boundary conditions.
! input: all, output: dcu
! approximate flop count is: 65*nxc*nyc
! and nxc*nyc divides
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
! dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
! dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
! the transverse part is calculated using the equation:
! dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!               (kx*kx+ky*ky)
! dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!               (kx*kx+ky*ky)
! on input:
! dcu(i,j,k) = complex acceleration density for fourier mode (jj-1,k-1)
! on output:
! dcu(i,k,j) = i-th component of transverse part of complex derivative
! of current for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! amu(1,k,j) = xx-yy component of complex momentum flux
! amu(2,k,j) = xy component of complex momentum flux
! amu(3,k,j) = zx component of complex momentum flux
! amu(4,k,j) = zy component of complex momentum flux
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex dcu, amu
      dimension dcu(3,nyv,kxp), amu(4,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate transverse part of current
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,k1,dkx,dkx2,dky,dky2,dkxy,dkxy2,at1,zt1,zt2,zt3)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j)),-real(amu(1,k,j)))
         zt2 = cmplx(aimag(amu(2,k,j)),-real(amu(2,k,j)))
         zt3 = at1*(dky*dcu(1,k,j) - dkx*dcu(2,k,j) + dkxy*zt1          &
     &       + dkxy2*zt2)
         dcu(1,k,j) = dky*zt3
         dcu(2,k,j) = -dkx*zt3
         zt1 = cmplx(aimag(amu(3,k,j)),-real(amu(3,k,j)))
         zt2 = cmplx(aimag(amu(4,k,j)),-real(amu(4,k,j)))
         dcu(3,k,j) = dcu(3,k,j) + dkx*zt1 + dky*zt2
         zt1 = cmplx(aimag(amu(1,k1,j)),-real(amu(1,k1,j)))
         zt2 = cmplx(aimag(amu(2,k1,j)),-real(amu(2,k1,j)))
         zt3 = at1*(dky*dcu(1,k1,j) + dkx*dcu(2,k1,j)                   &
     &       + dkxy*zt1 - dkxy2*zt2)
         dcu(1,k1,j) = dky*zt3
         dcu(2,k1,j) = dkx*zt3
         zt1 = cmplx(aimag(amu(3,k1,j)),-real(amu(3,k1,j)))
         zt2 = cmplx(aimag(amu(4,k1,j)),-real(amu(4,k1,j)))
         dcu(3,k1,j) = dcu(3,k1,j) + dkx*zt1 - dky*zt2
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j)),-real(amu(2,1,j)))
         dcu(1,1,j) = zero
         dcu(2,1,j) = dcu(2,1,j) + dkx*zt2
         zt1 = cmplx(aimag(amu(3,1,j)),-real(amu(3,1,j)))
         dcu(3,1,j) = dcu(3,1,j) + dkx*zt1
         dcu(1,k1,j) = zero
         dcu(2,k1,j) = zero
         dcu(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1)),-real(amu(2,k,1)))
         dcu(1,k,1) = dcu(1,k,1) + dky*zt2
         dcu(2,k,1) = zero
         zt2 = cmplx(aimag(amu(4,k,1)),-real(amu(4,k,1)))
         dcu(3,k,1) = dcu(3,k,1) + dky*zt2
         dcu(1,k1,1) = zero
         dcu(2,k1,1) = zero
         dcu(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1) = zero
         dcu(2,1,1) = zero
         dcu(1,k1,1) = zero
         dcu(2,k1,1) = zero
         dcu(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx, &
     &ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! transverse electric field (or convolution of transverse electric field
! over particle shape), with periodic boundary conditions.
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,kstrt,nyv,kxp,nyhd,
! output: ffe
! for isign /= 0, input: dcu,ffe,isign,affp,ci,nx,ny,kstrt,nyv,kxp,nyhd,
! output: exy,wf
! approximate flop count is: 59*nxc*nyc + 32*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! if isign = 0, form factor array is prepared
! if isign = -1, smoothed transverse electric field is calculated
! using the equations:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! if isign = 1, unsmoothed transverse electric field is calculated
! using the equations:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
! dcu(i,k,j) = i-th component of transverse part of complex derivative
! of current,
! exy(i,k,j) = i-th component of complex transverse electric field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffe(k,j)) = finite-size particle shape factor s
! real(ffe(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! wp0 = normalized total plasma frequency squared
! ci = reciprical of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
!    |dcu(kx,ky)*s(kx,ky)|**2)
! this expression is valid only if the derivative of current is
! divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(3,nyv,kxp), exy(3,nyv,kxp)
      dimension ffe(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
      wpc = wp0*ci2
! prepare form factor array
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffe(k,j) = cmplx(affp,1.0)
      else
         ffe(k,j) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
! calculate smoothed transverse electric field and sum field energy
      sum1 = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 70
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 50 j = 1, kxps
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*aimag(ffe(k,j))
         at2 = at2*at2
         exy(1,k,j) = at1*dcu(1,k,j)
         exy(2,k,j) = at1*dcu(2,k,j)
         exy(3,k,j) = at1*dcu(3,k,j)
         exy(1,k1,j) = at1*dcu(1,k1,j)
         exy(2,k1,j) = at1*dcu(2,k1,j)
         exy(3,k1,j) = at1*dcu(3,k1,j)
         wp = wp + at2*(dcu(1,k,j)*conjg(dcu(1,k,j))                    &
     &   + dcu(2,k,j)*conjg(dcu(2,k,j)) + dcu(3,k,j)*conjg(dcu(3,k,j))  &
     &   + dcu(1,k1,j)*conjg(dcu(1,k1,j))                               &
     &   + dcu(2,k1,j)*conjg(dcu(2,k1,j))                               &
     &   + dcu(3,k1,j)*conjg(dcu(3,k1,j)))
   40    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*aimag(ffe(1,j))
         at2 = at2*at2
         exy(1,1,j) = at1*dcu(1,1,j)
         exy(2,1,j) = at1*dcu(2,1,j)
         exy(3,1,j) = at1*dcu(3,1,j)
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
         wp = wp + at2*(dcu(1,1,j)*conjg(dcu(1,1,j))                    &
     &   + dcu(2,1,j)*conjg(dcu(2,1,j)) + dcu(3,1,j)*conjg(dcu(3,1,j)))
      endif
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*aimag(ffe(k,1))
         at2 = at2*at2
         exy(1,k,1) = at1*dcu(1,k,1)
         exy(2,k,1) = at1*dcu(2,k,1)
         exy(3,k,1) = at1*dcu(3,k,1)
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
         wp = wp + at2*(dcu(1,k,1)*conjg(dcu(1,k,1))                    &
     &   + dcu(2,k,1)*conjg(dcu(2,k,1)) + dcu(3,k,1)*conjg(dcu(3,k,1)))
   60    continue
         k1 = nyh + 1
         exy(1,1,1) = zero
         exy(2,1,1) = zero
         exy(3,1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   70 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
! calculate unsmoothed transverse electric field and sum field energy
   80 sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 120
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 100 j = 1, kxps
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 90 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*at2
         exy(1,k,j) = at2*dcu(1,k,j)
         exy(2,k,j) = at2*dcu(2,k,j)
         exy(3,k,j) = at2*dcu(3,k,j)
         exy(1,k1,j) = at2*dcu(1,k1,j)
         exy(2,k1,j) = at2*dcu(2,k1,j)
         exy(3,k1,j) = at2*dcu(3,k1,j)
         wp = wp + at1*(dcu(1,k,j)*conjg(dcu(1,k,j))                    &
     &   + dcu(2,k,j)*conjg(dcu(2,k,j)) + dcu(3,k,j)*conjg(dcu(3,k,j))  &
     &   + dcu(1,k1,j)*conjg(dcu(1,k1,j))                               &
     &   + dcu(2,k1,j)*conjg(dcu(2,k1,j))                               &
     &   + dcu(3,k1,j)*conjg(dcu(3,k1,j)))
   90    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*at2
         exy(1,1,j) = at2*dcu(1,1,j)
         exy(2,1,j) = at2*dcu(2,1,j)
         exy(3,1,j) = at2*dcu(3,1,j)
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
         wp = wp + at1*(dcu(1,1,j)*conjg(dcu(1,1,j))                    &
     &   + dcu(2,1,j)*conjg(dcu(2,1,j)) + dcu(3,1,j)*conjg(dcu(3,1,j)))
      endif
      sum1 = sum1 + wp
  100 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*at2
         exy(1,k,1) = at2*dcu(1,k,1)
         exy(2,k,1) = at2*dcu(2,k,1)
         exy(3,k,1) = at2*dcu(3,k,1)
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
         wp = wp + at1*(dcu(1,k,1)*conjg(dcu(1,k,1))                    &
     &   + dcu(2,k,1)*conjg(dcu(2,k,1)) + dcu(3,k,1)*conjg(dcu(3,k,1)))
  110    continue
         k1 = nyh + 1
         exy(1,1,1) = zero
         exy(2,1,1) = zero
         exy(3,1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
  120 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2d poisson's equation in fourier space for
! potential, with periodic boundary conditions for distributed data.
! input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: pot,we
! approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! potential is calculated using the equation:
! pot(kx,ky) = g(kx,ky)*q(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! pot(kx=pi) = pot(ky=pi) = 0, and pot(kx=0,ky=0) = 0.
! q(k,j) = complex charge density for fourier mode (jj-1,k-1)
! pot(k,j) = complex potential
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real we
      complex q, pot, ffc
      dimension q(nyv,kxp), pot(nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1, at2
      complex zero
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
! calculate potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 j = 1, kxps
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,j))
         at1 = at2*aimag(ffc(k,j))
         pot(k,j) = at2*q(k,j)
         pot(k1,j) = at2*q(k1,j)
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = real(ffc(1,j))
         at1 = at2*aimag(ffc(1,j))
         pot(1,j) = at2*q(1,j)
         pot(k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,1))
         at1 = at2*aimag(ffc(k,1))
         pot(k,1) = at2*q(k,1)
         pot(k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   30    continue
         k1 = nyh + 1
         pot(1,1) = zero
         pot(k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      we = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPELFIELD22(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2d poisson's equation in fourier space for
! unsmoothed longitudinal electric field, with periodic boundary
! conditions for distributed data
! input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: fxy,we
! approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the equation used is:
! fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
! fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
! fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! q(k,j) = complex charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of complex electric field,
! fxy(2,k,j) = y component of complex electric field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2),
! where affp = normalization constant = nx*ny/np, where np=number of
! particles
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(2,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, at1, at2, at3
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate unsmoothed longitudinal electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,at1,at2,at3,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         at1 = at1*aimag(ffc(k,j))
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))
         at3 = dkx*at1
         at1 = at1*aimag(ffc(1,j))
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         at1 = at1*aimag(ffc(k,1))
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   30    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      we = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPELFIELD23(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! unsmoothed longitudinal electric field, with periodic boundary
! conditions for distributed data
! input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: fxy,we
! approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the equation used is:
! fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
! fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
! fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! q(k,j) = complex charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of complex electric field,
! fxy(2,k,j) = y component of complex electric field,
! fxy(3,k,j) = zero,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np, where np=number of
! particles
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, at1, at2, at3
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate unsmoothed longitudinal electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,at1,at2,at3,zt1,zt2,wp)            &
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         at1 = at1*aimag(ffc(k,j))
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(3,k,j) = zero
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))
         at3 = dkx*at1
         at1 = at1*aimag(ffc(1,j))
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(3,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         at1 = at1*aimag(ffc(k,1))
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(3,k,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   30    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(3,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      we = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
! this subroutine calculates the divergence in fourier space
! input: all except df, output: df
! approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the divergence is calculated using the equation:
! df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp
      complex f, df
      dimension f(ndim,nyv,kxp), df(nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate the divergence
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = dkx*f(1,k,j) + dky*f(2,k,j)
         df(k,j) = cmplx(-aimag(zt1),real(zt1))
         zt1 = dkx*f(1,k1,j) - dky*f(2,k1,j)
         df(k1,j) = cmplx(-aimag(zt1),real(zt1))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         df(1,j) = dkx*cmplx(-aimag(f(1,1,j)),real(f(1,1,j)))
         df(k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         df(k,1) = dky*cmplx(-aimag(f(2,k,1)),real(f(2,k,1)))
         df(k1,1) = zero
   30    continue
         k1 = nyh + 1
         df(1,1) = zero
         df(k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
! this subroutine calculates the gradient in fourier space
! input: all except f, output: f
! approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the gradient is calculated using the equations:
! fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
! fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
! and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp
      complex df, f
      dimension df(nyv,kxp), f(ndim,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate the gradient
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(df(k,j)),real(df(k,j)))
         f(1,k,j) = dkx*zt1
         f(2,k,j) = dky*zt1
         zt1 = cmplx(-aimag(df(k1,j)),real(df(k1,j)))
         f(1,k1,j) = dkx*zt1
         f(2,k1,j) = -dky*zt1
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         f(1,1,j) = dkx*cmplx(-aimag(df(1,j)),real(df(1,j)))
         f(2,1,j) = zero
         f(1,k1,j) = zero
         f(2,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         f(1,k,1) = zero
         f(2,k,1) = dky*cmplx(-aimag(df(k,1)),real(df(k,1)))
         f(1,k1,1) = zero
         f(2,k1,1) = zero
   30    continue
         k1 = nyh + 1
         f(1,1,1) = zero
         f(2,1,1) = zero
         f(1,k1,1) = zero
         f(2,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
! this subroutine calculates the curl in fourier space
! input: all except g, output: g
! approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the curl is calculated using the equations:
! gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
! gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
! gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
! and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex f, g
      dimension f(3,nyv,kxp), g(3,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate the curl
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1,zt2,zt3)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(f(3,k,j)),real(f(3,k,j)))
         zt2 = cmplx(-aimag(f(2,k,j)),real(f(2,k,j)))
         zt3 = cmplx(-aimag(f(1,k,j)),real(f(1,k,j)))
         g(1,k,j) = dky*zt1
         g(2,k,j) = -dkx*zt1
         g(3,k,j) = dkx*zt2 - dky*zt3
         zt1 = cmplx(-aimag(f(3,k1,j)),real(f(3,k1,j)))
         zt2 = cmplx(-aimag(f(2,k1,j)),real(f(2,k1,j)))
         zt3 = cmplx(-aimag(f(1,k1,j)),real(f(1,k1,j)))
         g(1,k1,j) = -dky*zt1
         g(2,k1,j) = -dkx*zt1
         g(3,k1,j) = dkx*zt2 + dky*zt3
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt1 = cmplx(-aimag(f(3,1,j)),real(f(3,1,j)))
         zt2 = cmplx(-aimag(f(2,1,j)),real(f(2,1,j)))
         g(1,1,j) = zero
         g(2,1,j) = -dkx*zt1
         g(3,1,j) = dkx*zt2
         g(1,k1,j) = zero
         g(2,k1,j) = zero
         g(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(f(3,k,1)),real(f(3,k,1)))
         zt2 = cmplx(-aimag(f(1,k,1)),real(f(1,k,1)))
         g(1,k,1) = dky*zt1
         g(2,k,1) = zero
         g(3,k,1) = -dky*zt2
         g(1,k1,1) = zero
         g(2,k1,1) = zero
         g(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         g(1,1,1) = zero
         g(2,1,1) = zero
         g(3,1,1) = zero
         g(1,k1,1) = zero
         g(2,k1,1) = zero
         g(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
! this subroutine calculates 2-1/2d vector potential from magnetic field
! in fourier space with periodic boundary conditions
! for distributed data.
! input: bxy,nx,ny,kstrt,nyv,kxp, output: axy
! approximate flop count is: 38*nxc*nyc + 10*(nxc + nyc)
! and nxc*nyc divides
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the vector potential is calculated using the equations:
! ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
! ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
! az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
! ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
! bxy(i,k,j) = i-th component of complex magnetic field,
! axy(i,k,j) = i-th component of complex vector potential,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex bxy, axy
      dimension bxy(3,nyv,kxp), axy(3,nyv,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, at1, at2, at3
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate vector potential
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dkx2,dky,at1,at2,at3,zt1,zt2,zt3)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(bxy(3,k,j)),real(bxy(3,k,j)))
         zt2 = cmplx(-aimag(bxy(2,k,j)),real(bxy(2,k,j)))
         zt3 = cmplx(-aimag(bxy(1,k,j)),real(bxy(1,k,j)))
         axy(1,k,j) = at2*zt1
         axy(2,k,j) = -at3*zt1
         axy(3,k,j) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(bxy(3,k1,j)),real(bxy(3,k1,j)))
         zt2 = cmplx(-aimag(bxy(2,k1,j)),real(bxy(2,k1,j)))
         zt3 = cmplx(-aimag(bxy(1,k1,j)),real(bxy(1,k1,j)))
         axy(1,k1,j) = -at2*zt1
         axy(2,k1,j) = -at3*zt1
         axy(3,k1,j) = at3*zt2 + at2*zt3
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = 1.0/dkx
         zt1 = cmplx(-aimag(bxy(3,1,j)),real(bxy(3,1,j)))
         zt2 = cmplx(-aimag(bxy(2,1,j)),real(bxy(2,1,j)))
         axy(1,1,j) = zero
         axy(2,1,j) = -at2*zt1
         axy(3,1,j) = at2*zt2
         axy(1,k1,j) = zero
         axy(2,k1,j) = zero
         axy(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at2 = 1.0/dky
         zt1 = cmplx(-aimag(bxy(3,k,1)),real(bxy(3,k,1)))
         zt2 = cmplx(-aimag(bxy(1,k,1)),real(bxy(1,k,1)))
         axy(1,k,1) = at2*zt1
         axy(2,k,1) = zero
         axy(3,k,1) = -at2*zt2
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1) = zero
         axy(2,1,1) = zero
         axy(3,1,1) = zero
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MCUAVE23(cuave,cunew,cuold,ny,kxp,nyv)
! this subroutine averages current in fourier space for 2-1/2d code
! input: all except cuave, output: cuave
! cunew(i,k,j),cuold(i,k,j) = complex current densities to be averaged
! cuave(i,k,j) = average complex current density
! for component i, all for fourier mode (j-1,k-1)
! ny = system length in y direction
! kxp = number of data values per block
! nyv = first dimension of field arrays, must be >= ny
      implicit none
      integer ny, kxp, nyv
      complex cuave, cunew, cuold
      dimension cuave(3,nyv,kxp), cuold(3,nyv,kxp), cunew(3,nyv,kxp)
! local data
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 j = 1, kxp
      do 10 k = 1, ny
      cuave(1,k,j) = 0.5*(cunew(1,k,j) + cuold(1,k,j))
      cuave(2,k,j) = 0.5*(cunew(2,k,j) + cuold(2,k,j))
      cuave(3,k,j) = 0.5*(cunew(3,k,j) + cuold(3,k,j))
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,   &
     &nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! the radiative part of the vector potential
! with periodic boundary conditions, for distributed data.
! input: all, output: axy
! approximate flop count is: 68*nxc*nyc + 20*(nxc + nyc)
! and nxc*nyc divides
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the radiative vector potential is updated using the equations:
! ax(kx,ky) = (sqrt(-1)*ky*bz(kx,ky)
!                       - affp*ci2*cux(kx,ky)*s(kx,ky)/(kx*kx+ky*ky)
! ay(kx,ky) = -(sqrt(-1)*kx*bz(kx,ky)
!                       + affp*ci2*cuy(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
! az(kx,ky) = (sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*ci2*cuz(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, ci2 = ci*ci
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
! ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
! axy(i,k,j) = on entry, i-th component of complex current density cu,
! axy(i,k,j) = on exit, i-th component of complex radiative vector
! potential,
! bxy(i,k,j) = i-th component of complex magnetic field,
! aimag(ffc(k,j)) = finite-size particle shape factor s,
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! all for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp = number of data values per block
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci
      complex axy, bxy, ffc
      dimension axy(3,nyv,kxp), bxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, afc2, dkx, dkx2, dky, at1, at2
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      afc2 = affp*ci*ci
      zero = cmplx(0.,0.)
! calculate the radiative vector potential
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dkx2,dky,at1,at2,zt1,zt2,zt3)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         at2 = afc2*aimag(ffc(k,j))
! update radiative vector potential, ky > 0
         zt1 = cmplx(-aimag(bxy(3,k,j)),real(bxy(3,k,j)))
         zt2 = cmplx(-aimag(bxy(2,k,j)),real(bxy(2,k,j)))
         zt3 = cmplx(-aimag(bxy(1,k,j)),real(bxy(1,k,j)))
         axy(1,k,j) = at1*(dky*zt1 - at2*axy(1,k,j))
         axy(2,k,j) = -at1*(dkx*zt1 + at2*axy(2,k,j))
         axy(3,k,j) = at1*((dkx*zt2 - dky*zt3) - at2*axy(3,k,j))
! update radiative vector potential, ky < 0
         zt1 = cmplx(-aimag(bxy(3,k1,j)),real(bxy(3,k1,j)))
         zt2 = cmplx(-aimag(bxy(2,k1,j)),real(bxy(2,k1,j)))
         zt3 = cmplx(-aimag(bxy(1,k1,j)),real(bxy(1,k1,j)))
         axy(1,k1,j) = -at1*(dky*zt1 + at2*axy(1,k1,j))
         axy(2,k1,j) = -at1*(dkx*zt1 + at2*axy(2,k1,j))
         axy(3,k1,j) = at1*((dkx*zt2 + dky*zt3) - at2*axy(3,k1,j))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = 1.0/dkx2
         at2 = afc2*aimag(ffc(1,j))
! update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,1,j)),real(bxy(3,1,j)))
         zt2 = cmplx(-aimag(bxy(2,1,j)),real(bxy(2,1,j)))
         axy(1,1,j) = zero
         axy(2,1,j) = -at1*(dkx*zt1 + at2*axy(2,1,j))
         axy(3,1,j) = at1*(dkx*zt2 - at2*axy(3,1,j))
         axy(1,k1,j) = zero
         axy(2,k1,j) = zero
         axy(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky)
         at2 = afc2*aimag(ffc(k,1))
! update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,k,1)),real(bxy(3,k,1)))
         zt3 = cmplx(-aimag(bxy(1,k,1)),real(bxy(1,k,1)))
         axy(1,k,1) = at1*(dky*zt1 - at2*axy(1,k,1))
         axy(2,k,1) = zero
         axy(3,k,1) = -at1*(dky*zt3 + at2*axy(3,k,1))
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1) = zero
         axy(2,1,1) = zero
         axy(3,1,1) = zero
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAPOTP23(cu,axy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! vector potential with periodic boundary conditions
! for distributed data.
! input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,nyhd, output: axy,wm
! approximate flop count is: 63*nxc*nyc + 33*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! vector potential is calculated using the equation:
! ax(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
! ay(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
! az(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
! ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
! ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
! cu(i,k,j) = i-th component of complex current density and
! axy(i,k,j) = i-th component of complex vector potential,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstart - 1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real ci, wm
      complex cu, axy, ffc
      dimension cu(3,nyv,kxp), axy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real ci2, at1, at2
      complex zero
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate vector potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 j = 1, kxps
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,j))
         at1 = at2*aimag(ffc(k,j))
         axy(1,k,j) = at2*cu(1,k,j)
         axy(2,k,j) = at2*cu(2,k,j)
         axy(3,k,j) = at2*cu(3,k,j)
         axy(1,k1,j) = at2*cu(1,k1,j)
         axy(2,k1,j) = at2*cu(2,k1,j)
         axy(3,k1,j) = at2*cu(3,k1,j)
         wp = wp + at1*(cu(1,k,j)*conjg(cu(1,k,j))                      &
     &   + cu(2,k,j)*conjg(cu(2,k,j)) + cu(3,k,j)*conjg(cu(3,k,j))      &
     &   + cu(1,k1,j)*conjg(cu(1,k1,j)) + cu(2,k1,j)*conjg(cu(2,k1,j))  &
     &   + cu(3,k1,j)*conjg(cu(3,k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = ci2*real(ffc(1,j))
         at1 = at2*aimag(ffc(1,j))
         axy(1,1,j) = at2*cu(1,1,j)
         axy(2,1,j) = at2*cu(2,1,j)
         axy(3,1,j) = at2*cu(3,1,j)
         axy(1,k1,j) = zero
         axy(2,k1,j) = zero
         axy(3,k1,j) = zero
         wp = wp + at1*(cu(1,1,j)*conjg(cu(1,1,j))                      &
     &   + cu(2,1,j)*conjg(cu(2,1,j)) + cu(3,1,j)*conjg(cu(3,1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,1))
         at1 = at2*aimag(ffc(k,1))
         axy(1,k,1) = at2*cu(1,k,1)
         axy(2,k,1) = at2*cu(2,k,1)
         axy(3,k,1) = at2*cu(3,k,1)
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
         wp = wp + at1*(cu(1,k,1)*conjg(cu(1,k,1))                      &
     &   + cu(2,k,1)*conjg(cu(2,k,1)) + cu(3,k,1)*conjg(cu(3,k,1)))
   30    continue
         k1 = nyh + 1
         axy(1,1,1) = zero
         axy(2,1,1) = zero
         axy(3,1,1) = zero
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPETFIELD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,kxp&
     &,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! unsmoothed transverse electric field, with periodic boundary
! conditions for distributed data
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! input: dcu,ffe,affp,ci,nx,ny,kstrt,nyv,kxp,nyhd, output: exy,wf
! approximate flop count is: 59*nxc*nyc + 32*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! unsmoothed transverse electric field is calculated using the equation:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! dcu(i,k,j) = i-th component of transverse part of complex derivative
! of current,
! exy(i,k,j) = i-th component of complex transverse electric field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffe(k,j)) = finite-size particle shape factor s
! real(ffe(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
!    |dcu(kx,ky)*s(kx,ky)|**2)
! this expression is valid only if the derivative of current is
! divergence-free
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp = number of data values per block
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, wf
      complex dcu, exy, ffe
      dimension dcu(3,nyv,kxp), exy(3,nyv,kxp)
      dimension ffe(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real ci2, at1, at2
      complex zero
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate unsmoothed transverse electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 j = 1, kxps
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*at2
         exy(1,k,j) = at2*dcu(1,k,j)
         exy(2,k,j) = at2*dcu(2,k,j)
         exy(3,k,j) = at2*dcu(3,k,j)
         exy(1,k1,j) = at2*dcu(1,k1,j)
         exy(2,k1,j) = at2*dcu(2,k1,j)
         exy(3,k1,j) = at2*dcu(3,k1,j)
         wp = wp + at1*(dcu(1,k,j)*conjg(dcu(1,k,j))                    &
     &   + dcu(2,k,j)*conjg(dcu(2,k,j)) + dcu(3,k,j)*conjg(dcu(3,k,j))  &
     &   + dcu(1,k1,j)*conjg(dcu(1,k1,j))                               &
     &   + dcu(2,k1,j)*conjg(dcu(2,k1,j))                               &
     &   + dcu(3,k1,j)*conjg(dcu(3,k1,j)))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*at2
         exy(1,1,j) = at2*dcu(1,1,j)
         exy(2,1,j) = at2*dcu(2,1,j)
         exy(3,1,j) = at2*dcu(3,1,j)
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
         wp = wp + at1*(dcu(1,1,j)*conjg(dcu(1,1,j))                    &
     &   + dcu(2,1,j)*conjg(dcu(2,1,j)) + dcu(3,1,j)*conjg(dcu(3,1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*at2
         exy(1,k,1) = at2*dcu(1,k,1)
         exy(2,k,1) = at2*dcu(2,k,1)
         exy(3,k,1) = at2*dcu(3,k,1)
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
         wp = wp + at1*(dcu(1,k,1)*conjg(dcu(1,k,1))                    &
     &   + dcu(2,k,1)*conjg(dcu(2,k,1)) + dcu(3,k,1)*conjg(dcu(3,k,1)))
  30    continue
         k1 = nyh + 1
         exy(1,1,1) = zero
         exy(2,1,1) = zero
         exy(3,1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine provides a 2d scalar smoothing function
! in fourier space, with periodic boundary conditions
! for distributed data.
! input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: qs
! approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! smoothing is calculated using the equation:
! qs(kx,ky) = q(kx,ky)*s(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! qs(kx=pi) = qs(ky=pi) = 0, and qs(kx=0,ky=0) = 0.
! q(k,j) = complex charge density
! qs(k,j) = complex smoothed charge density,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      complex q, qs, ffc
      dimension q(nyv,kxp), qs(nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
! calculate smoothing
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
      do 20 j = 1, kxps
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         qs(k,j) = at1*q(k,j)
         qs(k1,j) = at1*q(k1,j)
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         qs(1,j) = at1*q(1,j)
         qs(k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1))
         qs(k,1) = at1*q(k,1)
         qs(k1,1) = zero
   30    continue
         k1 = nyh + 1
         qs(1,1) = cmplx(aimag(ffc(1,1))*real(q(1,1)),0.)
         qs(k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
! this subroutine provides a 2d vector smoothing function
! in fourier space, with periodic boundary conditions.
! for distributed data.
! input: cu,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: cus
! approximate flop count is: 12*nxc*nyc + 6*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! smoothing is calculated using the equation:
! cusx(kx,ky) = cux(kx,ky)*s(kx,ky)
! cusy(kx,ky) = cuy(kx,ky)*s(kx,ky)
! cusz(kx,ky) = cuz(kx,ky)*s(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! cusx(kx=pi) = cusy(kx=pi) = cusz(kx=pi) = 0,
! cusx(ky=pi) = cusy(ky=pi) = cusz(ky=pi) = 0,
! cusx(kx=0,ky=0) = cusy(kx=0,ky=0) = cusz(kx=0,ky=0) = 0.
! cu(i,k,j) = i-th component of complex current density and
! cus(i,k,j) = i-th component of complex smoothed current density
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! kxp = number of data values per block
! kstrt = starting data block number
! aimag(ffc(k,j)) = finite-size particle shape factor s
! real(ffc(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      complex cu, cus, ffc
      dimension cu(3,nyv,kxp), cus(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
! calculate smoothing
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
      do 20 j = 1, kxps
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         cus(1,k,j) = at1*cu(1,k,j)
         cus(2,k,j) = at1*cu(2,k,j)
         cus(3,k,j) = at1*cu(3,k,j)
         cus(1,k1,j) = at1*cu(1,k1,j)
         cus(2,k1,j) = at1*cu(2,k1,j)
         cus(3,k1,j) = at1*cu(3,k1,j)
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         cus(1,1,j) = at1*cu(1,1,j)
         cus(2,1,j) = at1*cu(2,1,j)
         cus(3,1,j) = at1*cu(3,1,j)
         cus(1,k1,j) = zero
         cus(2,k1,j) = zero
         cus(3,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1))
         cus(1,k,1) = at1*cu(1,k,1)
         cus(2,k,1) = at1*cu(2,k,1)
         cus(3,k,1) = at1*cu(3,k,1)
         cus(1,k1,1) = zero
         cus(2,k1,1) = zero
         cus(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         at1 = aimag(ffc(1,1))
         cus(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
         cus(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
         cus(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
         cus(1,k1,1) = zero
         cus(2,k1,1) = zero
         cus(3,k1,1) = zero
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp, &
     &modesxpd,modesyd)
! this subroutine extracts lowest order modes from packed complex array
! pot and stores them into a location in an unpacked complex array pott
! modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
! and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
! except kx = NX/2 is stored at location kxp+1 when idproc=0.
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! kstrt = starting data block number
! nyv = first dimension of input array pot, nyv >= ny
! kxp = number of data values per block
! modesyd = second dimension of array pott,
! where modesyd >= min(2*modesy-1,ny)
! modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
! unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp), pott(modesyd,modesxpd)
! local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, j1, k1, ks, joff
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(j,k,k1)
      do 20 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pott(2*k-2,j) = pot(k,j)
         pott(2*k-1,j) = pot(k1,j)
   10    continue
! mode numbers ky = 0, ny/2
         pott(1,j) = pot(1,j)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(ny,j) = pot(k1,j)
         endif
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, kmax
         k1 = ny2 - k
         pott(2*k-2,1) = pot(k,1)
         pott(2*k-1,1) = conjg(pot(k,1))
         if (modesx.gt.nxh) then
            pott(2*k-2,j1) = conjg(pot(k1,1))
            pott(2*k-1,j1) = pot(k1,1)
         endif
   30    continue
         pott(1,1) = cmplx(real(pot(1,1)),0.0)
         if (modesx.gt.nxh) then
            pott(1,j1) = cmplx(aimag(pot(1,1)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(ny,1) = cmplx(real(pot(k1,1)),0.0)
            if (modesx.gt.nxh) then
               pott(ny,j1) = cmplx(aimag(pot(k1,1)),0.0)
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp, &
     &modesxpd,modesyd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex array pott and stores them into a packed complex
! array pot
! modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
! and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
! except kx = NX/2 is stored at location kxp+1 when idproc=0.
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! kstrt = starting data block number
! nyv = first dimension of input array pot, nyv >= ny
! kxp = number of data values per block
! modesyd = second dimension of array pott,
! where modesyd  >= min(2*modesy-1,ny)
! modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
! unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp), pott(modesyd,modesxpd)
! local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, j1, k1, ks, joff
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(j,k,k1)
      do 30 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pot(k,j) = pott(2*k-2,j)
         pot(k1,j) = pott(2*k-1,j)
   10    continue
         do 20 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,j) = zero
         pot(k1,j) = zero
   20    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         pot(1,j) = pott(1,j)
         pot(k1,j) = zero
         if (modesy.gt.nyh) then
            pot(k1,j) = pott(ny,j)
         endif
      endif
   30 continue
!$OMP END PARALLEL DO
      do 50 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         pot(k,j) = zero
         pot(k1,j) = zero
   40    continue
         k1 = nyh + 1
! mode numbers ky = 0, ny/2
         pot(1,j) = zero
         pot(k1,j) = zero
      endif
   50 continue
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         pot(k,1) = pott(2*k-2,1)
         pot(k1,1) = zero
         if (modesx.gt.nxh) then
            pot(k1,1) = conjg(pott(2*k-2,j1))
         endif
   60    continue
         do 70 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,1) = zero
         pot(k1,1) = zero
   70    continue
         k1 = nyh + 1
         pot(1,1) = cmplx(real(pott(1,1)),0.0)
         pot(k1,1) = zero
         if (modesx.gt.nxh) then
            pot(1,1) = cmplx(real(pot(1,1)),real(pott(1,j1)))
         endif
         if (modesy.gt.nyh) then
            pot(k1,1) = cmplx(real(pott(ny,1)),0.0)
            if (modesx.gt.nxh) then
               pot(k1,1) = cmplx(real(pot(k1,1)),real(pott(ny,j1)))
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt, &
     &nyv,kxp,modesxpd,modesyd)
! this subroutine extracts lowest order modes from packed complex vector
! array vpot and stores them into a location in an unpacked complex
! vector array vpott
! modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
! and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
! except kx = NX/2 is stored at location kxp+1 when idproc=0.
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! ndim = number of field arrays, must be >= 1
! kstrt = starting data block number
! nyv = second dimension of input array vpot, nyv >= ny
! kxp = number of data values per block
! modesyd = third dimension of array vpott,
! where modesyd >= min(2*modesy-1,ny)
! modesxpd = fourth dimension of array vpott,
! modesxpd >= min(modesx,kxp),  unless modesx = nx/2+1,
! in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, ndim, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp)
      dimension vpott(ndim,modesyd,modesxpd)
! local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, j1, k1, ks, joff
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 40 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         k1 = ny2 - k
         do 10 i = 1, ndim
         vpott(i,2*k-2,j) = vpot(i,k,j)
         vpott(i,2*k-1,j) = vpot(i,k1,j)
   10    continue
   20    continue
! mode numbers ky = 0, ny/2
         do 30 i = 1, ndim
         vpott(i,1,j) = vpot(i,1,j)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(i,ny,j) = vpot(i,k1,j)
         endif
   30    continue
      endif
   40 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         do 50 i = 1, ndim
         vpott(i,2*k-2,1) = vpot(i,k,1)
         vpott(i,2*k-1,1) = conjg(vpot(i,k,1))
         if (modesx.gt.nxh) then
            vpott(i,2*k-2,j1) = conjg(vpot(i,k1,1))
            vpott(i,2*k-1,j1) = vpot(i,k1,1)
         endif
   50    continue
   60    continue
         do 70 i = 1, ndim
         vpott(i,1,1) = cmplx(real(vpot(i,1,1)),0.0)
         if (modesx.gt.nxh) then
            vpott(i,1,j1) = cmplx(aimag(vpot(i,1,1)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(i,ny,1) = cmplx(real(vpot(i,k1,1)),0.0)
            if (modesx.gt.nxh) then
               vpott(i,ny,j1) = cmplx(aimag(vpot(i,k1,1)),0.0)
            endif
         endif
   70    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt, &
     &nyv,kxp,modesxpd,modesyd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex vector array vpott and stores them into a packed
! complex vector array vpot
! modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
! and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
! except kx = NX/2 is stored at location kxp+1 when idproc=0.
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! ndim = number of field arrays, must be >= 1
! kstrt = starting data block number
! nyv = second dimension of input array vpot, nyv >= ny
! kxp = number of data values per block
! modesyd = third dimension of array vpott,
! where modesyd  >= min(2*modesy-1,ny)
! modesxpd = fourth dimension of array vpott,
! modesxpd >= min(modesx,kxp) unless modesx = nx/2+1,
! in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, ndim, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp)
      dimension vpott(ndim,modesyd,modesxpd)
! local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, j1, k1, ks, joff
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.nxh) return
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 60 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         do 10 i = 1, ndim
         k1 = ny2 - k
         vpot(i,k,j) = vpott(i,2*k-2,j)
         vpot(i,k1,j) = vpott(i,2*k-1,j)
   10    continue
   20    continue
         do 40 k = kmax+1, nyh
         k1 = ny2 - k
         do 30 i = 1, ndim
         vpot(i,k,j) = zero
         vpot(i,k1,j) = zero
   30    continue
   40    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 50 i = 1, ndim
         vpot(i,1,j) = vpott(i,1,j)
         vpot(i,k1,j) = zero
         if (modesy.gt.nyh) then
            vpot(i,k1,j) = vpott(i,ny,j)
         endif
   50    continue
      endif
   60 continue
!$OMP END PARALLEL DO
      do 100 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 i = 1, ndim
         vpot(i,k,j) = zero
         vpot(i,k1,j) = zero
   70    continue
   80    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 90 i = 1, ndim
         vpot(i,1,j) = zero
         vpot(i,k1,j) = zero
   90    continue
      endif
  100 continue
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 120 k = 2, kmax
         k1 = ny2 - k
         do 110 i = 1, ndim
         vpot(i,k,1) = vpott(i,2*k-2,1)
         vpot(i,k1,1) = zero
         if (modesx.gt.nxh) then
            vpot(i,k1,1) = conjg(vpott(i,2*k-2,j1))
         endif
  110    continue
  120    continue
         do 140 k = kmax+1, nyh
         k1 = ny2 - k
         do 130 i = 1, ndim
         vpot(i,k,1) = zero
         vpot(i,k1,1) = zero
  130    continue
  140    continue
         k1 = nyh + 1
         do 150 i = 1, ndim
         vpot(i,1,1) = cmplx(real(vpott(i,1,1)),0.0)
         vpot(i,k1,1) = zero
         if (modesx.gt.nxh) then
            vpot(i,1,1) = cmplx(real(vpot(i,1,1)),real(vpott(i,1,j1)))
         endif
         if (modesy.gt.nyh) then
            vpot(i,k1,1) = cmplx(real(vpott(i,ny,1)),0.0)
            if (modesx.gt.nxh) then
               vpot(i,k1,1) = cmplx(real(vpot(i,k1,1)),                 &
     &                              real(vpott(i,ny,j1)))
            endif
         endif
  150    continue
      endif
      return
      end
