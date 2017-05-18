!-----------------------------------------------------------------------
! Fortran library for spectral field solvers with dirichlet boundary
! conditions
! 2D MPI/OpenMP PIC Code:
! MPPOISD22 solves 2d poisson's equation for smoothed electric field
! MPPOISD23 solves 2-1/2d poisson's equation for smoothed electric field
! MPPCUPERPD2 calculates the transverse current in fourier space
! MIPPBPOISD23 solves 2-1/2d poisson's equation for unsmoothed magnetic
!              field
! MPPMAXWELD2 solves 2-1/2d maxwell's equation for unsmoothed transverse
!             electric and magnetic fields
! MPPEMFIELDR2 adds and smooths or copies and smooths real vector fields
!              in fourier space
! MPPBBPOISP23 solves 2-1/2d poisson's equation in fourier space for
!              smoothed magnetic field
! MPPDCUPERP23 calculates transverse part of the derivative of the
!              current density from the momentum flux
! MPPADCUPERPD23 calculates transverse part of the derivative of the
!                current density from the momentum flux and acceleration
!                density
! MPPEPOISP23 solves 2-1/2d poisson's equation in fourier space for
!             smoothed or unsmoothed transverse electric field
! MPPOTPD2 solves 2d poisson's equation for potential
! MPPELFIELDD22 solves 2d poisson's equation for unsmoothed electric
!               field
! MPPELFIELDD23 solves 2-1/2d poisson's equation for unsmoothed electric
!               field
! MPPDIVFD2 calculates the divergence in fourier space
! MPPGRADFD2 calculates the gradient in fourier space
! MPPCURLFD2 calculates the curl in fourier space
! MPPAVPOTD23 calculates 2-1/2d vector potential from magnetic field
! MPPAPOTD23 solves 2-1/2d poisson's equation for vector potential
! MPPETFIELDD23 solves 2-1/2d poisson's equation in fourier space for
!               unsmoothed transverse electric field
! MPPSMOOTHD2 provides a 2d scalar smoothing function
! MPPSMOOTHD23 provides a 2-1/2d vector smoothing function
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: april 27, 2017
!-----------------------------------------------------------------------
      subroutine MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,nyd)
! this subroutine solves 2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,nyd
!               output: ffd
! for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,nyd
!                 output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed force/charge,
! fxy(2,k,j) = y component of transformed force/charge,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! if isign = 0, form factor array is prepared
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real ax, ay, affp, we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(2,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      if (isign.ne.0) go to 30
! prepare form factor array
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-0.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffd(k,j) = cmplx(affp,1.0)
      else
         ffd(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 80
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, ny
         at1 = real(ffd(k,j))*aimag(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         wp = wp + at1*q(k,j)**2
   40    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 60 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
   60    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 70 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
   70 continue
   80 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! Zeros out z component
! for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,nyd
!               output: ffd
! for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,nyd
!                 output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed force/charge,
! fxy(2,k,j) = y component of transformed force/charge,
! fxy(3,k,j) = zero,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! if isign = 0, form factor array is prepared
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real ax, ay, affp, we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      if (isign.ne.0) go to 30
! prepare form factor array
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-0.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffd(k,j) = cmplx(affp,1.0)
      else
         ffd(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 80
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, ny
         at1 = real(ffd(k,j))*aimag(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         fxy(3,k,j) = 0.0
         wp = wp + at1*q(k,j)**2
   40    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(3,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      fxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 60 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
         fxy(3,k,1) = 0.0
   60    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 70 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
      fxy(3,k,kxp2s+1) = 0.0
   70 continue
   80 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates the transverse current in fourier space
! with dirichlet boundary conditions (zero potential).
! input: all, output: cu
! approximate flop count is: 10*nx*ny and nx*ny divides
! the transverse current is calculated using the equation:
! cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! cu(i,k,j) = i-th component of transformed current density and
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real cu
      dimension cu(3,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, dkx2, at1, at2, at3, at4
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate transverse part of current
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dkx2,dky,at1,at2,at3,at4)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         at2 = at1*(dkx*cu(1,k,j) + dky*cu(2,k,j))
         at3 = cu(1,k,j) - dkx*at2
         at4 = cu(2,k,j) - dky*at2
         cu(1,k,j) = at3
         cu(2,k,j) = at4
   10    continue
! mode numbers ky = 0, ny
         cu(1,1,j) = 0.0
         cu(1,ny+1,j) = 0.0
         cu(2,ny+1,j) = 0.0
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         cu(2,k,1) = 0.0
   30    continue
         cu(1,1,1) = 0.0
         cu(2,1,1) = 0.0
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      cu(1,k,kxp2s+1) = 0.0
      cu(2,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MIPPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! magnetic field with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: cu,ffd,ci,nx,ny,kstrt,nyv,kxp2,nyd, output: bxy,wm
! approximate flop count is: 20*nx*ny
! magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci*g(kx,ky)*ky*cuz(kx,ky),
! by(kx,ky) = -ci*ci*g(kx,ky)*kx*cuz(kx,ky),
! bz(kx,ky) = ci*ci*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! cu(i,k,j) = i-th component of transformed current density and
! bxy(i,k,j) = i-th component of transformed magnetic field,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)) = finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real ci, wm
      real cu, bxy
      dimension cu(3,nyv,kxp2+1), bxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, ci2, at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffd(k,j))
         bxy(1,k,j) = at2*cu(3,k,j)
         bxy(2,k,j) = -at3*cu(3,k,j)
         bxy(3,k,j) = at3*cu(2,k,j) - at2*cu(1,k,j)
         wp = wp + 2.0*at1*(cu(1,k,j)**2 + cu(2,k,j)**2 + cu(3,k,j)**2)
   10    continue
! mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j))
         at2 = dkx*at1
         at1 = at1*aimag(ffd(1,j))
         bxy(1,1,j) = 0.0
         bxy(2,1,j) = 0.0
         bxy(3,1,j) = at2*cu(2,1,j)
         wp = wp + at1*(cu(2,1,j)**2 + cu(3,1,j)**2)
      endif
      bxy(1,ny+1,j) = 0.0
      bxy(2,ny+1,j) = 0.0
      bxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,1))
         at2 = dky*at1
         at1 = at1*aimag(ffd(k,1))
         bxy(1,k,1) = 0.0
         bxy(2,k,1) = 0.0
         bxy(3,k,1) = -at2*cu(1,k,1)
         wp = wp + at1*(cu(1,k,1)**2 + cu(3,k,1)**2)
   30    continue
         bxy(1,1,1) = 0.0
         bxy(2,1,1) = 0.0
         bxy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      bxy(1,k,kxp2s+1) = 0.0
      bxy(2,k,kxp2s+1) = 0.0
      bxy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt&
     &,nyv,kxp2,nyd)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse electric and magnetic fields with dirichlet boundary
! conditions (zero potential).
! using fast sine/cosine transforms for distributed data.
! input: all, output: wf, wm, exy, bxy
! approximate flop count is: 61*nx*ny
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky) = bx(kx,ky) - .5*dt*ky*ez(kx,ky)
! by(kx,ky) = by(kx,ky) + .5*dt*kx*ez(kx,ky)
! bz(kx,ky) = bz(kx,ky) - .5*dt*(kx*ey(kx,ky)-ky*ex(kx,ky))
! the electric field is then updated a whole step using the equations:
! ex(kx,ky) = ex(kx,ky) - c2*dt*ky*bz(kx,ky)
!                       - affp*dt*cux(kx,ky)*s(kx,ky)
! ey(kx,ky) = ey(kx,ky) + c2*dt*kx*bz(kx,ky)
!                       - affp*dt*cuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = ez(kx,ky) - c2*dt*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*dt*cuz(kx,ky)*s(kx,ky)
! the magnetic field is finally updated the remaining half step with
! the new electric field and the previous magnetic field equations.
! where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! modes nx and ny are zeroed out
! cu(i,k,j) = i-th component of transformed current density and
! exy(i,k,j) = i-th component of transformed electric field,
! bxy(i,k,j) = i-th component of transformed magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! aimag(ffd(k,j)) = finite-size particle shape factor s
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprocal of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
! magnetic field energy is also calculated, using
! wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
! nx/ny = system length in x/y direction
! kxp2 = number of data values per block
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real affp, ci, dt, wf, wm
      real cu, exy, bxy
      dimension exy(3,nyv,kxp2+1), bxy(3,nyv,kxp2+2), cu(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dth, c2, cdt, adt, anorm, dkx, dky, afdt
      real at4, at5, at6, at7, at8, at9
      double precision wp, ws, sum1, sum2
      if (ci.le.0.0) return
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum1 = 0.0d0
      sum2 = 0.0d0
      if (kstrt.gt.nx) go to 50
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky,afdt,at4,at5,at6,at7,at8,at9,ws,  
!$OMP& wp)
!$OMP& REDUCTION(+:sum1,sum2)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      ws = 0.0d0
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffd(k,j))
         at7 = exy(1,k,j)
         at8 = exy(2,k,j)
         at9 = exy(3,k,j)
! update magnetic field half time step, ky > 0
         at4 = bxy(1,k,j) - dth*(dky*at9)
         at5 = bxy(2,k,j) + dth*(dkx*at9)
         at6 = bxy(3,k,j) + dth*(dky*at7 - dkx*at8)
! update electric field whole time step
         at7 = at7 - cdt*(dky*at6) - afdt*cu(1,k,j)
         at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,k,j)
         at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,k,j)
! update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dky*at7 - dkx*at8)
         ws = ws + 2.0*anorm*(at7*at7 + at8*at8 + at9*at9)
         wp = wp + 2.0*anorm*(at4*at4 + at5*at5 + at6*at6)
         exy(1,k,j) = at7
         exy(2,k,j) = at8
         exy(3,k,j) = at9
         bxy(1,k,j) = at4
         bxy(2,k,j) = at5
         bxy(3,k,j) = at6
   10    continue
! mode numbers ky = 0, ny
         afdt = adt*aimag(ffd(1,j))
         at8 = exy(2,1,j)
         at9 = exy(3,1,j)
! update magnetic field half time step, ky > 0
         at5 = bxy(2,1,j) + dth*(dkx*at9)
         at6 = bxy(3,1,j) - dth*(dkx*at8)
! update electric field whole time step
         at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,1,j)
         at9 = at9 - cdt*(dkx*at5) - afdt*cu(3,1,j)
! update magnetic field half time step and store electric field
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 - dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j) = 0.0
         exy(2,1,j) = at8
         exy(3,1,j) = at9
         bxy(1,1,j) = 0.0
         bxy(2,1,j) = at5
         bxy(3,1,j) = at6
      endif
      exy(1,ny+1,j) = 0.0
      exy(2,ny+1,j) = 0.0
      exy(3,ny+1,j) = 0.0
      bxy(1,ny+1,j) = 0.0
      bxy(2,ny+1,j) = 0.0
      bxy(3,ny+1,j) = 0.0
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffd(k,1))
         at7 = exy(1,k,1)
         at9 = exy(3,k,1)
! update magnetic field half time step, ky > 0
         at4 = bxy(1,k,1) - dth*(dky*at9)
         at6 = bxy(3,k,1) + dth*(dky*at7)
! update electric field whole time step
         at7 = at7 - cdt*(dky*at6) - afdt*cu(1,k,1)
         at9 = at9 + cdt*(dky*at4) - afdt*cu(3,k,1)
! update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at6 = at6 + dth*(dky*at7)
         ws = ws + anorm*(at7*at7 + at9*at9)
         wp = wp + anorm*(at4*at4 + at6*at6)
         exy(1,k,1) = at7
         exy(2,k,1) = 0.0
         exy(3,k,1) = at9
         bxy(1,k,1) = at4
         bxy(2,k,1) = 0.0
         bxy(3,k,1) = at6
   30    continue
! zero out kx = ky = 0 mode
         exy(1,1,1) = 0.0
         exy(2,1,1) = 0.0
         exy(3,1,1) = 0.0
         bxy(1,1,1) = 0.0
         bxy(2,1,1) = 0.0
         bxy(3,1,1) = 0.0
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      exy(1,k,kxp2s+1) = 0.0
      exy(2,k,kxp2s+1) = 0.0
      exy(3,k,kxp2s+1) = 0.0
      bxy(1,k,kxp2s+1) = 0.0
      bxy(2,k,kxp2s+1) = 0.0
      bxy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wf = real(nx)*real(ny)*sum1
      wm = real(nx)*real(ny)*c2*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,nyd&
     &)
! this subroutine either adds complex vector fields if isign > 0
! or copies complex vector fields if isign < 0
! includes additional smoothing
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real fxy, exy
      complex ffd
      dimension fxy(3,nyv,kxp2+1), exy(3,nyv,kxp2+1)
      dimension ffd(nyd,kxp2)
! local data
      integer i, j, k, ks, ny1, kxp2s
      real at1
      ks = kstrt - 1
      ny1 = ny + 1
      kxp2s = min(kxp2,max(0,nx-kxp2*ks))
      if (kstrt.gt.nx) return
! smooth the em field and add
      if (isign.gt.0) then
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(i,j,k,at1)
         do 40 j = 1, kxp2s
         do 20 k = 1, ny
         at1 = aimag(ffd(k,j))
         do 10 i = 1, 3
         fxy(i,k,j) = fxy(i,k,j) + exy(i,k,j)*at1
   10    continue
   20    continue
! mode numbers ky = 0, ny
         do 30 i = 1, 3
         fxy(i,ny+1,j) = fxy(i,ny+1,j) + exy(i,ny+1,j)
   30    continue
   40    continue
!$OMP END PARALLEL DO
         do 60 k = 1, ny1
         do 50 i = 1, 3
         fxy(i,k,kxp2s+1) = fxy(i,k,kxp2s+1) + exy(i,k,kxp2s+1)
   50    continue
   60    continue
! copy and smooth the magnetic fields
      else if (isign.lt.0) then
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(i,j,k,at1)
         do 100 j = 1, kxp2s
         do 80 k = 1, ny
         at1 = aimag(ffd(k,j))
         do 70 i = 1, 3
         fxy(i,k,j) = exy(i,k,j)*at1
   70    continue
   80    continue
! mode numbers ky = 0, ny
         do 90 i = 1, 3
         fxy(i,ny+1,j) = exy(i,ny+1,j)
   90    continue
  100    continue
!$OMP END PARALLEL DO
         do 120 k = 1, ny1
         do 110 i = 1, 3
         fxy(i,k,kxp2s+1) = exy(i,k,kxp2s+1)
  110    continue
  120    continue
! copy the electric fields
      else if (isign.eq.0) then
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do 160 j = 1, kxp2s
         do 140 k = 1, ny
         do 130 i = 1, 3
         fxy(i,k,j) = exy(i,k,j)
  130    continue
  140    continue
! mode numbers ky = 0, ny
         do 150 i = 1, 3
         fxy(i,ny+1,j) = exy(i,ny+1,j)
  150    continue
  160    continue
!$OMP END PARALLEL DO
         do 180 k = 1, ny1
         do 170 i = 1, 3
         fxy(i,k,kxp2s+1) = exy(i,k,kxp2s+1)
  170    continue
  180    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPBBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! smoothed magnetic field (or convolution of magnetic field over 
! particle shape) with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,nyd, output: bxy,wm
! approximate flop count is: 20*nx*ny
! magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci**g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
! by(kx,ky) = -ci*ci*sg(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
! bz(kx,ky) = ci*ci*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*s(kx,ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! cu(i,k,j) = i-th component of transformed current density and
! bxy(i,k,j) = i-th component of transformed magnetic field,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! aimag(ffd(k,j)) = finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number=
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real ci, wm
      real cu, bxy
      dimension cu(3,nyv,kxp2+1), bxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, ci2, at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      ci2 = ci*ci
! calculate smoothed magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,j))*aimag(ffd(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j) = at2*cu(3,k,j)
         bxy(2,k,j) = -at3*cu(3,k,j)
         bxy(3,k,j) = at3*cu(2,k,j) - at2*cu(1,k,j)
         wp = wp + 2.0*at1*(cu(1,k,j)**2 + cu(2,k,j)**2 + cu(3,k,j)**2)
   10    continue
! mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j))*aimag(ffd(1,j))
         at2 = dkx*at1
         bxy(1,1,j) = 0.0
         bxy(2,1,j) = 0.0
         bxy(3,1,j) = at2*cu(2,1,j)
         wp = wp + at1*(cu(2,1,j)**2 + cu(3,1,j)**2)
      endif
      bxy(1,ny+1,j) = 0.0
      bxy(2,ny+1,j) = 0.0
      bxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,1))*aimag(ffd(k,1))
         at2 = dky*at1
         bxy(1,k,1) = 0.0
         bxy(2,k,1) = 0.0
         bxy(3,k,1) = -at2*cu(1,k,1)
         wp = wp + at1*(cu(1,k,1)**2 + cu(3,k,1)**2)
   30    continue
         bxy(1,1,1) = 0.0
         bxy(2,1,1) = 0.0
         bxy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      bxy(1,k,kxp2s+1) = 0.0
      bxy(2,k,kxp2s+1) = 0.0
      bxy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux in 2-1/2d
! with dirichlet boundary conditions (zero potential).
! input: all except dcu, output: dcu
! approximate flop count is: 16*nx*ny and nx*ny divides
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky) = -(kx*vx*vx-ky*vx*vy)
! dcu(2,kx,ky) = (kx*vx*vy-ky*vy*vy)
! dcu(3,kx,ky) = (kx*vx*vz+ky*vy*vz)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! the transverse part is calculated using the equation:
! dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!                (kx*kx+ky*ky)
! dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!                (kx*kx+ky*ky)
! on output:
! dcu(i,k,j) = i-th component of transverse part of transformed
! derivative of current
! amu(1,k,j) = xx-yy component of transformed momentum flux
! amu(2,k,j) = xy component of transformed momentum flux
! amu(3,k,j) = zx component of transformed momentum flux
! amu(4,k,j) = zy component of transformed momentum flux
! all for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real dcu, amu
      dimension dcu(3,nyv,kxp2+1), amu(4,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1, at3
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate transverse part of current
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,dkx,dkx2,dky,dky2,dkxy,dkxy2,at1,at3)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         at3 = at1*(-dkxy*amu(1,k,j) + dkxy2*amu(2,k,j))
         dcu(1,k,j) = dky*at3
         dcu(2,k,j) = -dkx*at3
         dcu(3,k,j) = dkx*amu(3,k,j) + dky*amu(4,k,j)
   10    continue
! mode numbers ky = 0, ny
         dcu(1,1,j) = 0.0
         dcu(2,1,j) = dkx*amu(2,1,j)
         dcu(3,1,j) = dkx*amu(3,1,j)
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         dcu(1,k,1) = dky*amu(2,k,1)
         dcu(2,k,1) = 0.0
         dcu(3,k,1) = dky*amu(4,k,1)
   30    continue
         dcu(1,1,1) = 0.0
         dcu(2,1,1) = 0.0
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      dcu(1,k,kxp2s+1) = 0.0
      dcu(2,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPADCUPERPD23(dcu,amu,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux and acceleration density
! in 2-1/2d with dirichlet boundary conditions (zero potential).
! input: all, output: dcu
! approximate flop count is: 21*nx*ny and nx*ny divides
! the derivative of the current is calculated using the equations:
! dcu(1,kx,ky) = dcu(1,kx,ky)-(kx*vx*vx-ky*vx*vy)
! dcu(2,kx,ky) = dcu(2,kx,ky)+(kx*vx*vy-ky*vy*vy)
! dcu(3,kx,ky) = dcu(3,kx,ky)+(kx*vx*vz+ky*vy*vz)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! the transverse part is calculated using the equation:
! dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!                (kx*kx+ky*ky)
! dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
!                (kx*kx+ky*ky)
! on input:
! dcu(i,j,k) = i-th component of transformed acceleration density
! on output:
! dcu(i,k,j) = i-th component of transverse part of transformed
! derivative of current
! amu(1,k,j) = xx-yy component of transformed momentum flux
! amu(2,k,j) = xy component of transformed momentum flux
! amu(3,k,j) = zx component of transformed momentum flux
! amu(4,k,j) = zy component of transformed momentum flux
! all for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real dcu, amu
      dimension dcu(3,nyv,kxp2+1), amu(4,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1, at3
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate transverse part of current
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,dkx,dkx2,dky,dky2,dkxy,dkxy2,at1,at3)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         at3 = at1*(dky*dcu(1,k,j) - dkx*dcu(2,k,j) - dkxy*amu(1,k,j)   &
     &       + dkxy2*amu(2,k,j))
         dcu(1,k,j) = dky*at3
         dcu(2,k,j) = -dkx*at3
         dcu(3,k,j) = dcu(3,k,j) + dkx*amu(3,k,j) + dky*amu(4,k,j)
   10    continue
! mode numbers ky = 0, ny
         dcu(1,1,j) = 0.0
         dcu(2,1,j) = dcu(2,1,j) + dkx*amu(2,1,j)
         dcu(3,1,j) = dcu(3,1,j) + dkx*amu(3,1,j)
      endif
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         dcu(1,k,1) = dcu(1,k,1) + dky*amu(2,k,1)
         dcu(2,k,1) = 0.0
         dcu(3,k,1) = dcu(3,k,1) + dky*amu(4,k,1)
   30    continue
         dcu(1,1,1) = 0.0
         dcu(2,1,1) = 0.0
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      dcu(1,k,kxp2s+1) = 0.0
      dcu(2,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPEPOISD23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx, &
     &ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! transverse electric field (or convolution of transverse electric field
! over particle shape),
! with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,kstrt,nyv,kxp2,nyd,
! output: ffe
! for isign /= 0, input: dcu,ffe,isign,affp,ci,nx,ny,kstrt,nyv,kxp2,nyh,
! output: exy,wf
! approximate flop count is: 15*nx*ny
! if isign = 0, form factor array is prepared
! if isign = -1, smoothed transverse electric field is calculated
! using the equations:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! if isign = 1, unsmoothed transverse electric field is calculated
! using the equations:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
! dcu(i,k,j) = i-th component of transverse part of transformed
! derivative of current,
! exy(i,k,j) = i-th component of transformed transverse electric field,
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffe(k,j)) = finite-size particle shape factor s
! real(ffe(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprocal of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
!    |dcu(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real ax, ay, affp, wp0, ci, wf
      real dcu, exy
      dimension dcu(3,nyv,kxp2+1), exy(3,nyv,kxp2+1)
      complex ffe
      dimension ffe(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
! prepare form factor array
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-0.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffe(k,j) = cmplx(affp,1.0)
      else
         ffe(k,j) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
! calculate smoothed transverse electric field and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 80
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, ny
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*aimag(ffe(k,j))
         at2 = at2*at2
         exy(1,k,j) = at1*dcu(1,k,j)
         exy(2,k,j) = at1*dcu(2,k,j)
         exy(3,k,j) = at1*dcu(3,k,j)
         wp = wp + 2.0*at2*(dcu(1,k,j)**2 + dcu(2,k,j)**2               &
     &   + dcu(3,k,j)**2)
   40    continue
! mode numbers ky = 0, ny
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*aimag(ffe(1,j))
         at2 = at2*at2
         exy(1,1,j) = 0.0
         exy(2,1,j) = at1*dcu(2,1,j)
         exy(3,1,j) = at1*dcu(3,1,j)
         wp = wp + at2*(dcu(2,1,j)**2 + dcu(3,1,j)**2)
      endif
      exy(1,ny+1,j) = 0.0
      exy(2,ny+1,j) = 0.0
      exy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 60 k = 2, ny
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*aimag(ffe(k,1))
         at2 = at2*at2
         exy(1,k,1) = at1*dcu(1,k,1)
         exy(2,k,1) = 0.0
         exy(3,k,1) = at1*dcu(3,k,1)
         wp = wp + at2*(dcu(1,k,1)**2 + dcu(3,k,1)**2)
   60    continue
         exy(1,1,1) = 0.0
         exy(2,1,1) = 0.0
         exy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 70 k = 1, ny1
      exy(1,k,kxp2s+1) = 0.0
      exy(2,k,kxp2s+1) = 0.0
      exy(3,k,kxp2s+1) = 0.0
   70 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
! calculate unsmoothed transverse electric field and sum field energy
   80 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 130
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 100 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 90 k = 2, ny
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*at2
         exy(1,k,j) = at2*dcu(1,k,j)
         exy(2,k,j) = at2*dcu(2,k,j)
         exy(3,k,j) = at2*dcu(3,k,j)
         wp = wp + 2.0*at1*(dcu(1,k,j)**2 + dcu(2,k,j)**2               &
     &   + dcu(3,k,j)**2)
   90    continue
! mode numbers ky = 0, ny
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*at2
         exy(1,1,j) = 0.0
         exy(2,1,j) = at2*dcu(2,1,j)
         exy(3,1,j) = at1*dcu(3,1,j)
         wp = wp + at1*(dcu(2,1,j)**2 + dcu(3,1,j)**2)
      endif
      exy(1,ny+1,j) = 0.0
      exy(2,ny+1,j) = 0.0
      exy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
  100 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 110 k = 2, ny
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*at2
         exy(1,k,1) = at2*dcu(1,k,1)
         exy(2,k,1) = 0.0
         exy(3,k,1) = at1*dcu(3,k,1)
         wp = wp + at1*(dcu(1,k,1)**2 + dcu(3,k,1)**2)
  110    continue
         exy(1,1,1) = 0.0
         exy(2,1,1) = 0.0
         exy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 120 k = 1, ny1
      exy(1,k,kxp2s+1) = 0.0
      exy(2,k,kxp2s+1) = 0.0
      exy(3,k,kxp2s+1) = 0.0
  120 continue
  130 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOTPD2(q,pot,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2d poisson's equation in fourier space for
! potential, with dirichlet boundary conditions (zero potential),
! using fast sine transforms for distributed data.
! input: q,ffd,nx,ny,kstrt,nyv,kxp2,nyd, output: pot,we
! approximate flop count is: 5*nx*ny
! potential is calculated using the equation:
! fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! pot(k,j) = transformed potential,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real we
      real q, pot
      dimension q(nyv,kxp2+1), pot(nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
! calculate potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at2 = real(ffd(k,j))
         at1 = at2*aimag(ffd(k,j))
         at3 = at2*q(k,j)
         pot(k,j) = at3
         wp = wp + at1*q(k,j)**2
  10     continue
      endif
! mode numbers ky = 0, ny
      pot(1,j) = 0.0
      pot(ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         pot(k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      pot(k,kxp2s+1) = 0.0
   40 continue
   50 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPELFIELDD22(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2d poisson's equation in fourier space for
! unsmoothed longitudinal electric field, with dirichlet boundary
! conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: q,ffd,nx,ny,kstrt,ny2d,kxp2,nyd, output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*q(kx,ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed electric field,
! fxy(2,k,j) = y component of transformed electric field,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(2,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate unsmoothed longitudinal electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at1 = real(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         at1 = at1*aimag(ffd(k,j))
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         wp = wp + at1*q(k,j)**2
   10    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPELFIELDD23(q,fxy,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! unsmoothed longitudinal electric field, with dirichlet boundary
! conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! Zeros out z component
! input: q,ffd,nx,ny,kstrt,ny2d,kxp2,nyd, output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed electric field,
! fxy(2,k,j) = y component of transformed electric field,
! fxy(3,k,j) = zero,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate unsmoothed longitudinal electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at1 = real(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         at1 = at1*aimag(ffd(k,j))
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         fxy(3,k,j) = 0.0
         wp = wp + at1*q(k,j)**2
   10    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(3,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      fxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
         fxy(3,k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
      fxy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2)
! this subroutine calculates the divergence in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the charge density from the electric field
! input: all except df, output: df
! approximate flop count is: 6*nx*ny
! the divergence is calculated using the equation:
! df(kx,ky) = -(kx*fx(kx,ky)+ky*fy(kx,ky))
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! modes nx and ny are zeroed out
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv ,kxp2
      real f, df
      dimension f(ndim,nyv,kxp2+1), df(nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      if (ndim.lt.2) return
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the divergence
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         df(k,j) = -(dkx*f(1,k,j) + dky*f(2,k,j))
   10    continue
      endif
! mode numbers ky = 0, ny
      df(1,j) = 0.0
      df(ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         df(k,1) = 0.0
   30    continue
      endif
      do 40 k = 1, ny1
      df(k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2)
! this subroutine calculates the gradient in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the electric field from the potential
! input: all except f, output: f
! approximate flop count is: 4*nx*ny
! the gradient is calculated using the equations:
! fx(kx,ky) = kx*df(kx,ky)
! fy(kx,ky) = ky*df(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! modes nx and ny are zeroed out
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp2
      real df, f
      dimension df(nyv,kxp2+1), f(ndim,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the gradient
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         f(1,k,j) = dkx*df(k,j)
         f(2,k,j) = dky*df(k,j)
   10    continue
      endif
! mode numbers ky = 0, ny
      f(1,1,j) = 0.0
      f(2,1,j) = 0.0
      f(1,ny+1,j) = 0.0
      f(2,ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         f(1,k,1) = 0.0
         f(2,k,1) = 0.0
   30    continue
      endif
      do 40 k = 1, ny1
      f(1,k,kxp2s+1) = 0.0
      f(2,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates the curl in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the magnetic field from the vector potential
! input: all except g, output: g
! approximate flop count is: 8*nx*ny
! the curl is calculated using the equations:
! gx(kx,ky) = ky*fz(kx,ky)
! gy(kx,ky) = -kx*fz(kx,ky)
! gz(kx,ky) = (kx*fy(kx,ky)-ky*fx(kx,ky))
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real f, g
      dimension f(3,nyv,kxp2+1), g(3,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the curl
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         g(1,k,j) = dky*f(3,k,j)
         g(2,k,j) = -dkx*f(3,k,j)
         g(3,k,j) = dkx*f(2,k,j) - dky*f(1,k,j)
   10    continue
! mode numbers ky = 0, ny
         g(1,1,j) = 0.0
         g(2,1,j) = 0.0
         g(3,1,j) = dkx*f(2,1,j)
      endif
      g(1,ny+1,j) = 0.0
      g(2,ny+1,j) = 0.0
      g(3,ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         g(1,k,1) = 0.0
         g(2,k,1) = 0.0
         g(3,k,1) = -dky*f(1,k,1)
   30    continue
         g(1,1,1) = 0.0
         g(2,1,1) = 0.0
         g(3,1,1) = 0.0
      endif
      do 40 k = 1, ny1
      g(1,k,kxp2s+1) = 0.0
      g(2,k,kxp2s+1) = 0.0
      g(3,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates 2-1/2d vector potential from magnetic field
! in fourier space with dirichlet boundary conditions (zero potential).
! using fast sine/cosine transforms for distributed data.
! input: bxy,nx,ny,kstrt,ny2d,kxp2, output: axy
! approximate flop count is: 12*nx*ny and nx*ny divides!
! the vector potential is calculated using the equations:
! ax(kx,ky) = -(ky*bz(kx,ky))/(kx*kx+ky*ky)
! ay(kx,ky) = (kx*bz(kx,ky))/(kx*kx+ky*ky)
! az(kx,ky) = -(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers
! modes nx and ny are zeroed out
! bxy(i,k,j) = i-th component of transformed magnetic field,
! axy(i,k,j) = i-th component of transformed vector potential,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! nx/ny = system length in x/y direction
! kxp2 = number of data values per block
! kstrt = starting data block number
! nyv = second dimension of field arrays, must be >= ny+1
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real bxy, axy
      dimension bxy(3,nyv,kxp2+1), axy(3,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, dkx2, at1, at2, at3, at4, at5, at6
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate vector potential
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky,dkx2,at1,at2,at3,at4,at5,at6)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         at4 = bxy(3,k,j)
         at5 = bxy(2,k,j)
         at6 = bxy(1,k,j)
         axy(1,k,j) = -at2*at4
         axy(2,k,j) = at3*at4
         axy(3,k,j) = at2*at6 - at3*at5
   10    continue
! mode numbers ky = 0, ny
         at2 = 1.0/dkx
         at4 = bxy(3,1,j)
         at5 = bxy(2,1,j)
         axy(1,1,j) = 0.0
         axy(2,1,j) = at2*at4
         axy(3,1,j) = -at2*at5
      endif
      axy(1,ny+1,j) = 0.0
      axy(2,ny+1,j) = 0.0
      axy(3,ny+1,j) = 0.0
   20 continue
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         at2 = 1.0/dky
         at4 = bxy(3,k,1)
         at6 = bxy(1,k,1)
         axy(1,k,1) = -at2*at4
         axy(2,k,1) = 0.0
         axy(3,k,1) = at2*at6
   30    continue
         axy(1,1,1) = 0.0
         axy(2,1,1) = 0.0
         axy(3,1,1) = 0.0
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      axy(1,k,kxp2s+1) = 0.0
      axy(2,k,kxp2s+1) = 0.0
      axy(3,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPAPOTD23(cu,axy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! vector potential with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,nyd, output: axy,wm
! approximate flop count is: 14*nx*ny
! vector potential is calculated using the equation:
! bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
! by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
! bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! cu(i,k,j) = i-th component of transformed current density and
! axy(i,k,j) = i-th component of transformed vector potential,
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)) = finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real ci, wm
      real cu, axy
      dimension cu(3,nyv,kxp2+1), axy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real ci2, at1, at2
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      ci2 = ci*ci
! calculate vector potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at2 = ci2*real(ffd(k,j))
         at1 = at2*aimag(ffd(k,j))
         axy(1,k,j) = at2*cu(1,k,j)
         axy(2,k,j) = at2*cu(2,k,j)
         axy(3,k,j) = at2*cu(3,k,j)
         wp = wp + 2.0*at1*(cu(1,k,j)**2 + cu(2,k,j)**2 + cu(3,k,j)**2)
   10    continue
! mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j))
         at1 = at2*aimag(ffd(1,j))
         axy(1,1,j) = 0.0
         axy(2,1,j) = at2*cu(2,1,j)
         axy(3,1,j) = at2*cu(3,1,j)
         wp = wp + at1*(cu(2,1,j)**2 + cu(3,1,j)**2)
      endif
      axy(1,ny+1,j) = 0.0
      axy(2,ny+1,j) = 0.0
      axy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         at2 = ci2*real(ffd(k,1))
         at1 = at2*aimag(ffd(k,1))
         axy(1,k,1) = at2*cu(1,k,1)
         axy(2,k,1) = 0.0
         axy(3,k,1) = at2*cu(3,k,1)
         wp = wp + at1*(cu(1,k,1)**2 + cu(3,k,1)**2)
   30    continue
         axy(1,1,1) = 0.0
         axy(2,1,1) = 0.0
         axy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      axy(1,k,kxp2s+1) = 0.0
      axy(2,k,kxp2s+1) = 0.0
      axy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPETFIELDD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,  &
     &kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! unsmoothed transverse electric field with dirichlet boundary
! conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: dcu,ffe,ci,nx,ny,kstrt,ny2d,kxp2,nyd, output: exy,wf
! approximate flop count is: 15*nx*ny
! unsmoothed transverse electric field is calculated using the equation:
! ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
! ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
! ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! dcu(i,k,j) = i-th component of transverse part of transformed
! derivative of current,
! exy(i,k,j) = i-th component of transformed transverse electric field,
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffe(k,j)) = finite-size particle shape factor s
! real(ffe(k,j)) = potential green's function g
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprocal of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
!    |dcu(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real affp, ci, wf
      real dcu, exy
      dimension dcu(3,nyv,kxp2+1), exy(3,nyv,kxp2+1)
      complex ffe
      dimension ffe(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real ci2, at1, at2
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      ci2 = ci*ci
! calculate unsmoothed transverse electric field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at2 = -ci2*real(ffe(k,j))
         at1 = at2*at2
         exy(1,k,j) = at2*dcu(1,k,j)
         exy(2,k,j) = at2*dcu(2,k,j)
         exy(3,k,j) = at2*dcu(3,k,j)
         wp = wp + 2.0*at1*(dcu(1,k,j)**2 + dcu(2,k,j)**2               &
     &   + dcu(3,k,j)**2)
   10    continue
! mode numbers ky = 0, ny
         at2 = -ci2*real(ffe(1,j))
         at1 = at2*at2
         exy(1,1,j) = 0.0
         exy(2,1,j) = at2*dcu(2,1,j)
         exy(3,1,j) = at2*dcu(3,1,j)
         wp = wp + at1*(dcu(2,1,j)**2 + dcu(3,1,j)**2)
      endif
      exy(1,ny+1,j) = 0.0
      exy(2,ny+1,j) = 0.0
      exy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         at2 = -ci2*real(ffe(k,1))
         at1 = at2*at2
         exy(1,k,1) = at2*dcu(1,k,1)
         exy(2,k,1) = 0.0
         exy(3,k,1) = at2*dcu(3,k,1)
         wp = wp + at1*(dcu(1,k,1)**2 + dcu(3,k,1)**2)
   30    continue
         exy(1,1,1) = 0.0
         exy(2,1,1) = 0.0
         exy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      exy(1,k,kxp2s+1) = 0.0
      exy(2,k,kxp2s+1) = 0.0
      exy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wf = real(nx)*real(ny)*sum1/affp
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTHD2(q,qs,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine provides a 2d scalar smoothing function
! in fourier space, with dirichlet boundary conditions (zero potential),
! using fast sine transforms for distributed data.
! input: q,ffd,nx,ny,kstrt,nyv,kxp2,nyd, output: qs
! approximate flop count is: 1*nx*ny
! smoothing is calculated using the equation:
! qs(kx,ky) = q(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! qs(k,j) = transformed smoothed charge density,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)) = finite-size particle shape factor s
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real q, qs
      dimension q(nyv,kxp2+1), qs(nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real at1, at2
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
! calculate smoothing
      if (kstrt.gt.nx) return
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2)
      do 20 j = 1, kxp2s
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at1 = aimag(ffd(k,j))
         at2 = at1*q(k,j)
         qs(k,j) = at2
  10     continue
      endif
! mode numbers ky = 0, ny
      qs(1,j) = 0.0
      qs(ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         qs(k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      qs(k,kxp2s+1) = 0.0
   40 continue
      return
      end

!-----------------------------------------------------------------------
      subroutine MPPSMOOTHD23(cu,cus,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine provides a 2d vector smoothing function
! using fast sine/cosine transforms for distributed data.
! input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,nyd, output: cus
! approximate flop count is: 3*nx*ny
! smoothing is calculated using the equation:
! cusx(kx,ky) = cux(kx,ky)*s(kx,ky)
! cusy(kx,ky) = cuy(kx,ky)*s(kx,ky)
! cusz(kx,ky) = cuz(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! cu(i,k,j) = i-th component of transformed current density and
! cus(i,k,j) = i-th component of transformed smoothed current density
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)) = finite-size particle shape factor s
! for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real cu, cus
      dimension cu(3,nyv,kxp2+1), cus(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real at1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
! calculate smoothing
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1)
      do 20 j = 1, kxp2s
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at1 = aimag(ffd(k,j))
         cus(1,k,j) = at1*cu(1,k,j)
         cus(2,k,j) = at1*cu(2,k,j)
         cus(3,k,j) = at1*cu(3,k,j)
   10    continue
! mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j))
         cus(1,1,j) = 0.0
         cus(2,1,j) = at1*cu(2,1,j)
         cus(3,1,j) = 0.0
      endif
      cus(1,ny+1,j) = 0.0
      cus(2,ny+1,j) = 0.0
      cus(3,ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         at1 = aimag(ffd(k,1))
         cus(1,k,1) = at1*cu(1,k,1)
         cus(2,k,1) = 0.0
         cus(3,k,1) = 0.0
   30    continue
         cus(1,1,1) = 0.0
         cus(2,1,1) = 0.0
         cus(3,1,1) = 0.0
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      cus(1,k,kxp2s+1) = 0.0
      cus(2,k,kxp2s+1) = 0.0
      cus(3,k,kxp2s+1) = 0.0
   40 continue
      return
      end
