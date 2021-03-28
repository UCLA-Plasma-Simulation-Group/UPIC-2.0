!-----------------------------------------------------------------------
! Fortran library for spectral field solvers
! 1D OpenMP PIC Code:
! POIS1 solves 1d poisson's equation for smoothed electric field
! ADDQEI1 adds electron and ion densities
! IBPOIS13 solves 1-2/2d poisson's equation for unsmoothed magnetic
!          field
! MAXWEL1 solves 1-2/2d maxwell's equation for unsmoothed transverse
!         electric and magnetic fields using verlet algorithm
! AMAXWEL1 solves 1-2/2d maxwell's equation for unsmoothed transverse
!         electric and magnetic fields using analytic algorithm due to
!         Irving Haber
! EMFIELD1 merges complex vector fields in fourier space
! BMFIELD1 copies complex vector fields in fourier space
! ADDCUEI13 adds electron and ion current densities
! EADDEXT1 add external traveling wave field to electric field for
!          1d code
! EADDEXT13 add external traveling and external circularly polarized
!           wave fields to electric field for 1-2/2d code
! BADDEXT1 adds constant to magnetic field for 1-2/2d code
! ADDVRFIELD13 calculates a = b + c
! BBPOIS13 solves 1-2/2d poisson's equation in fourier space for
!          smoothed magnetic field
! DCUPERP13 calculates transverse part of the derivative of the
!           current density from the momentum flux
! ADCUPERP13 calculates transverse part of the derivative of the
!            current density from the momentum flux and acceleration
!            density
! EPOIS13 solves 1-2/2d poisson's equation in fourier space for
!         smoothed or unsmoothed transverse electric field
! POTP1 solves 1d poisson's equation for potential
! ELFIELD1 solves 1d poisson's equation for unsmoothed electric field
! DIVF1 calculates the divergence in fourier space
! GRADF1 calculates the gradient in fourier space
! CURLF1 calculates the curl in fourier space
! AVPOT13 calculates 1-2/2d vector potential from magnetic field
! CUAVE13 averages current in fourier space for 1-2/2d code
! AVRPOT13 solves 1-2/2d poisson's equation for the radiative part of
!          the vector potential
! APOTP13 solves 1-2/2d poisson's equation for vector potential
! ETFIELD13 solves 1-2/2d poisson's equation in fourier space for
!           unsmoothed transverse electric field
! SMOOTH1 provides a 1d scalar smoothing function
! SMOOTH13 provides a 1-2/2d vector smoothing function
! RDMODES1 extracts lowest order scalar modes from packed array
!          stores them into a location in an unpacked array
! WRMODES1 extracts lowest order scalar modes from a location in an
!          unpacked array and stores them into a packed array
! RDVMODES1 extracts lowest order vector modes from packed array
!           stores them into a location in an unpacked array
! WRVMODES1 extracts lowest order vector modes from a location in an
!           unpacked array and stores them into a packed array
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 3, 2020
!-----------------------------------------------------------------------
      subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
! this subroutine solves 1d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions.
! for isign = 0, input: isign,ax,affp,nx, output: ffc
! for isign  /= 0, input: q,ffc,isign,nx, output: fx,we
! approximate flop count is: 6*nx
! equation used is:
! fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
! g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
! fx(k=0) = fx(k=pi) = 0.
! cmplx(q(2*j-1),q(2*j)) = complex charge density for fourier mode j-1
! cmplx(fx(2*j-1),fx(2*j)) = complex force/charge for fourier mode j-1
! if isign = 0, form factor array is prepared
! if isign is not equal to 0, force/charge is calculated
! ffc(2*j) = finite-size particle shape factor s for fourier mode j-1
! ffc(2*j-1) = potential green's function g for fourier mode j-1
! ax = half-width of particle in x direction
! affp = normalization constant = nx/np, where np = number of particles
! electric field energy is also calculated, using
! we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
! nx = system length in x direction
      implicit none
      integer isign, nx
      real ax, affp, we
      real q, fx, ffc
      dimension q(nx), fx(nx), ffc(nx)
! local data
      integer j, nxh
      real dnx, dkx, at1, at2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      if (isign.ne.0) go to 20
! prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      ffc(2*j) = exp(-.5*(dkx*ax)**2)
      ffc(2*j-1) = affp*ffc(2*j)/(dkx*dkx)
   10 continue
      ffc(1) = affp
      ffc(2) = 1.0
      return
! calculate force/charge and sum field energy
   20 wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at1 = ffc(2*j-1)*ffc(2*j)
      at2 = dnx*real(j - 1)*at1
      fx(2*j-1) = at2*q(2*j)
      fx(2*j) = -at2*q(2*j-1)
      wp = wp + at1*(q(2*j-1)**2 + q(2*j)**2)
   30 continue
! mode number kx = 0
      fx(1) = 0.
      fx(2) = 0.
      we = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine ADDQEI1(qe,qi,nx,nxe)
! adds electron and ion densities
! assumes guard cells have already been added
! qe/qi = charge density for electrons/ions
! nx = system length in x/y direction
! nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      integer nx, nxe
      real qe, qi
      dimension qe(nxe), qi(nxe)
! local data
      integer j
      do 10 j = 1, nx
      qe(j) = qe(j) + qi(j)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! magnetic field, with periodic boundary conditions.
! input: cu,ffc,ci,nx,nxv, output: byz,wm
! approximate flop count is: 29*nxc
! where nxc = nx/2 - 1
! the magnetic field is calculated using the equations:
! by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
! bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
! where kx = 2pi*j/nx, and j = fourier mode number,
! g(kx) = (affp/kx**2)*s(kx),
! s(kx) = exp(-((kx*ax)**2)/2), except for
! by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
! cu(i,j) = complex current density for fourier mode (j-1)
! byz(i,j) = i component of complex magnetic field
! all for fourier mode (j-1)
! aimag(ffc(j)) = finite-size particle shape factor s
! for fourier mode (j-1)
! real(ffc(j)) = potential green's function g
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
! affp = normalization constant = nx/np, where np=number of particles
! this expression is valid only if the current is divergence-free
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real dnx, ci2, at1, at2
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
! this subroutine solves 1d maxwell's equation in fourier space for
! transverse electric and magnetic fields with periodic boundary
! conditions.
! input: all, output: wf, wm, eyz, byz
! approximate flop count is: 87*nxc
! where nxc = nx/2 - 1
! the magnetic field is first updated half a step using the equations:
! by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
! bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
! the electric field is then updated a whole step using the equations:
! ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
! ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
! the magnetic field is finally updated the remaining half step with
! the new electric field and the previous magnetic field equations.
! where kx = 2pi*j/nx, c2 = 1./(ci*ci)
! and s(kx) = exp(-((kx*ax)**2)
! j = fourier mode numbers, except for
! ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0) = 0.
! and similarly for by, bz.
! cu(i,j) = complex current density
! eyz(i,j) = complex transverse electric field
! byz(i,j) = complex magnetic field
! for component i, all for fourier mode (j-1)
! real(ffc(1)) = affp = normalization constant = nx/np,
! where np=number of particles
! aimag(ffc(j)) = finite-size particle shape factor s,
! s(kx) = exp(-((kx*ax)**2)/2)
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*sum((1/affp)*|eyz(kx)|**2)
! magnetic field energy is also calculated, using
! wm = nx*sum((c2/affp)*|byz(kx)|**2)
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
! nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      complex zero, zt1, zt2, zt5, zt6, zt8, zt9
      double precision wp, ws
      if (ci.le.0.0) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j))
! update magnetic field half time step
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
      zt5 = byz(1,j) + dth*(dkx*zt1)
      zt6 = byz(2,j) - dth*(dkx*zt2)
! update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = eyz(1,j) - cdt*(dkx*zt1) - afdt*cu(1,j)
      zt9 = eyz(2,j) + cdt*(dkx*zt2) - afdt*cu(2,j)
! update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      eyz(1,j) = zt8
      eyz(2,j) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      byz(1,j) = zt5
      byz(2,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
! mode numbers kx = 0, nx/2
!     eyz(1,1) = zero
!     eyz(2,1) = zero
!     byz(1,1) = zero
!     byz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*c2*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine AMAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
! this subroutine solves 1d maxwell's equation in fourier space for
! transverse electric and magnetic fields with periodic boundary
! conditions, using an analytic scheme due to Irving Haber
! input: all, output: wf, wm, eyz, byz
! approximate flop count is: 91*nxc plus 1 div, 1 sin, 1 cos, 1 tan
! per particle, where nxc = nx/2 - 1
! equations being solved are:
! (c*Bn)/dt = -ick X En and
! En/dt = ick X (c*Bn) - JTn+1/2
! Note that the normalization of the E and B fields differ:
! B is normalized to the dimensionless cyclotron frequency.
! solutions are given by:
! En+1 = C*En + iS*(k X (c*Bn)) - S*JTn+1/2/c
! c*Bn+1 = C*(c*Bn) - iS*(k X En)/c + iS*T*(k X JTn+1/2)
! where En = input eyz, En+1 = output eyz
! Bn = input byz, Bn+1 = output byz, JTn+1/2 = affp*cu*s(kx)
! C = cos(k*c*dt),  S = sin(k*c*dt)/k, T = tan(k*c*dt/2)/kc
! k = kx = 2pi*j/nx, c = 1.0/ci
! and s(kx) = exp(-(kx*ax)**2/2)
! j = fourier mode numbers, except for
! ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0) = 0.
! and similarly for by, bz.
! cu(i,j) = complex current density
! eyz(i,j) = complex transverse electric field
! byz(i,j) = complex magnetic field
! for component i, all for fourier mode (j-1)
! real(ffc(1)) = affp = normalization constant = nx/np,
! where np=number of particles
! aimag(ffc(j)) = finite-size particle shape factor s
! ci = reciprocal of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*sum((1/affp)*|eyz(kx)|**2)
! magnetic field energy is also calculated, using
! wm = nx*sum((c*c/affp)*|byz(kx)|**2)
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
! nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh)
      dimension ffc(nxhd)
! local data
      integer j, nxh
      real dnx, dth, cc, cdth, affp, anorm, dkx
      real t2, t, c, s, sc, afs, aft
      complex zero, zt1, zt2, zt5, zt6, zt7, zt8
      real at1
      double precision wp, ws
      if (ci.le.0.0) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      dth = 0.5*dt
      cc = 1.0/ci
      cdth = cc*dth
      zero = cmplx(0.0,0.0)
      affp = real(ffc(1))
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
! calculate the electromagnetic fields
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      aft = affp*aimag(ffc(j))*ci
      t = dkx*cdth
      t2 = 1.0/dkx
      s = t + t
      c = cos(s)
      s = sin(s)*t2
      t = tan(t)*t2
      afs = s*aft
      sc = s*cc
      aft = aft*t
      s = s*ci
! calculate iB
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
! update electric field
      zt5 = c*eyz(1,j) - sc*(dkx*zt1) - afs*cu(1,j)
      zt6 = c*eyz(2,j) + sc*(dkx*zt2) - afs*cu(2,j)
! calculate iE
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
! store electric field and calculate energy
      eyz(1,j) = zt5
      eyz(2,j) = zt6
      at1 = anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      ws = ws + dble(at1)
! calculate ijperp
      zt7 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt8 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
! update magnetic field 
      zt5 = c*byz(1,j) + s*dkx*(zt1 - aft*zt7)
      zt6 = c*byz(2,j) - s*dkx*(zt2 - aft*zt8)
! store magnetic field and calculate energy
      byz(1,j) = zt5
      byz(2,j) = zt6
      at1 = anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
   10 continue
! mode number kx = 0
!     afs = affp*aimag(ffc(1))*dt
!     eyz(1,1) = eyz(1,1) - afs*cu(1,1)
!     eyz(2,1) = eyz(2,1) - afs*cu(2,1)
      byz(1,1) = zero
      byz(2,1) = zero
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*(cc*cc)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
! this subroutine merges complex vector fields
! includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fxyz, fx, eyz, ffc
      dimension fxyz(3,nxvh), fx(nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
! local data
      integer j, nxh
      real at1
      nxh = nx/2
! add the fields
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = eyz(1,j)*at1
      fxyz(3,j) = eyz(2,j)*at1
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
! this subroutine copies complex vector fields
! includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fyz, eyz, ffc
      dimension fyz(2,nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
! local data
      integer j, nxh
      real at1
      nxh = nx/2
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fyz(1,j) = eyz(1,j)*at1
      fyz(2,j) = eyz(2,j)*at1
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ADDCUEI13(cue,cui,nx,nxe)
! adds electron and ion current densities
! assumes guard cells have already been added
! cue/cui = current density for electrons/ions
! nx = system length in x/y direction
! nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      integer nx, nxe
      real cue, cui
      dimension cue(2,nxe), cui(2,nxe)
! local data
      integer j
      do 10 j = 1, nx
      cue(1,j) = cue(1,j) + cui(1,j)
      cue(2,j) = cue(2,j) + cui(2,j)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine EADDEXT1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx,nxe)
! add external traveling wave field to electric field
! fxe(x) = fxe(x) + e1*sin(k0*x + freq*time) + e2*cos(k0*x - freq*time)
! where
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = el0*((toff+trmp-time)/trmp), e2 = er0*((toff+trmp-time)/trmp)
!                                            if toff < time < toff+trmp
!     e1 = 0,               e2 = 0,               if time > toff+trmp
! if toff < 0 => toff = + infinity
      implicit none
      integer nx, nxe
      real amodex, freq, time, trmp, toff, el0, er0
      real fxe
      dimension fxe(nxe)
! local data
      integer j
      real tr, at, ft, dkx, xk
      if ((el0==0.0).and.(er0==0.0)) return
      tr = toff + trmp
! ramp up
      if (time < trmp) then
         at = time/trmp
      else if ((toff >= 0.0).and.(time > toff)) then
! ramp down
         if (time < tr) then
            at = (tr - time)/trmp
! shutdown
         else
            at = 0.0
         endif
! constant amplitude
      else
         at = 1.0
      endif
      ft = freq*time
      dkx = (6.28318530717959*amodex)/real(nx)
      if (at > 0.0) then
         do 10 j = 1, nx
         xk = dkx*real(j - 1)
         fxe(j) = fxe(j) + at*(er0*cos(xk - ft) + el0*sin(xk + ft))
   10    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine EADDEXT13(fxyze,amodex,freq,time,trmp,toff,el0,er0,ey0,&
     &ez0,nx,nxe)
! for 1-2/2d code
! add external traveling wave field to electrostatic electric field
! fxyze(1,x) = fxyze(1,x) + e1*sin(k0*x + freq*time)
!                         + e2*cos(k0*x - freq*time)
! where
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = el0*((toff+trmp-time)/trmp), e2 = er0*((toff+trmp-time)/trmp)
!                                            if toff < time < toff+trmp
!     e1 = 0,               e2 = 0,               if time > toff+trmp
! if toff < 0 => toff = + infinity
! add external circularly polarized wave field to electromagnetic
! electric field:
! fxyze(2,x) = fxyze(2,x) + e3*cos(k0*x + freq*time)
! fxyze(3,x) = fxyze(3,x) + e4*sin(k0*x - freq*time)
! where
!     e3 = ey0*(time/trmp), e4 = ez0*(time/trmp), if time < trmp
!     e3 = ey0,             e4 = ez0,             if trmp < time < toff

!     e3 = ey0*((toff+trmp-time)/trmp), e4 = ez0*((toff+trmp-time)/trmp)
!                                            if toff < time < toff+trmp
!     e3 = 0,               e4 = 0,               if time > toff+trmp
! if toff < 0 => toff = + infinity
      implicit none
      integer nx, nxe
      real amodex, freq, time, trmp, toff, el0, er0, ey0, ez0
      real fxyze
      dimension fxyze(3,nxe)
! local data
      integer j
      real tr, at, ft, dkx, xk
      if ((el0==0.0).and.(er0==0.0).and.(ey0==0.0).and.(ez0==0.0))      &
     &return
      tr = toff + trmp
! ramp up
      if (time < trmp) then
         at = time/trmp
      else if ((toff >= 0.0).and.(time > toff)) then
! ramp down
         if (time < tr) then
            at = (tr - time)/trmp
! shutdown
         else
            at = 0.0
         endif
! constant amplitude
      else
         at = 1.0
      endif
      ft = freq*time
      dkx = (6.28318530717959*amodex)/real(nx)
      if (at > 0.0) then
         do 10 j = 1, nx
         xk = dkx*real(j - 1)
         fxyze(1,j) = fxyze(1,j) + at*(er0*cos(xk - ft)                 &
     &                               + el0*sin(xk + ft))
         fxyze(2,j) = fxyze(2,j) + at*(ey0*cos(xk - ft))
         fxyze(3,j) = fxyze(3,j) + at*(ez0*sin(xk - ft))
   10    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine BADDEXT1(byz,omy,omz,nx,nxe)
! adds constant to magnetic field for 1-2/2d code
! byz = magnetic field
! omy/omz = magnetic field electron cyclotron frequency in y/z 
! nx = system length in x direction
! nxe = second dimension of magnetic field array, nxe must be >= nx
      implicit none
      real byz, omy, omz
      integer nx, nxe
      dimension byz(2,nxe)
! local data
      integer j
      do 10 j = 1, nx
      byz(1,j) = byz(1,j) + omy
      byz(2,j) = byz(2,j) + omz
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ADDVRFIELD13(fxyze,eyze,fxe,nxe)
! this subroutine merges longitudinal and transverse electric fields
! fxyze(1,:) = fxe
! fxyze(2:3,:) = eyze
      implicit none
      integer nxe
      real fxyze, fxe, eyze
      dimension fxyze(3,nxe), eyze(2,nxe), fxe(nxe)
! local data
      integer j
      do 10 j = 1, nxe
      fxyze(1,j) = fxe(j)
      fxyze(2,j) = eyze(1,j)
      fxyze(3,j) = eyze(2,j)
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! magnetic field (or convolution of magnetic field over particle shape)
! with periodic boundary conditions.
! input: cu,ffc,ci,nx,nxvh,nxhd output: byz,wm
! approximate flop count is: 30*nxc
! where nxc = nx/2 - 1
! the magnetic field is calculated using the equations:
! by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx)*s(kx),
! bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx)*s(kx),
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! g(kx) = (affp/kx**2)*s(kx),
! s(kx) = exp(-((kx*ax)**2+)/2), except for
! by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
! cu(i,j) = complex current density for fourier mode (j-1)
! byz(1,j) = y component of complex magnetic field
! byz(2,j) = z component of complex magnetic field
! all for fourier mode (j-1)
! aimag(ffc(j)) = finite-size particle shape factor s
! for fourier mode (j-1)
! real(ffc(j)) = potential green's function g
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
! affp = normalization constant = nx/np, where np=number of particles
! this expression is valid only if the current is divergence-free
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real dnx, ci2, at1, at2
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j))*aimag(ffc(j))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine DCUPERP13(dcu,amu,nx,nxvh)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux
! in 1-2/2d with periodic boundary conditions.
! the transverse part of the derivative of the current is calculated
! using the equations:
! dcu(1,kx) = -sqrt(-1)*kx*vx*vy
! dcu(2,kx) = -sqrt(-1)*kx*vx*vz
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
! amu(1,j) = xy component of complex momentum flux
! amu(2,j) = xz component of complex momentum flux
! all for fourier mode (j-1)
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
! local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dkx*zt1
   10 continue
      dcu(1,1) = zero
      dcu(2,1) = zero
      return
      end
!-----------------------------------------------------------------------
      subroutine ADCUPERP13(dcu,amu,nx,nxvh)
! this subroutine calculates transverse part of the derivative of
! the current density from the momentum flux and acceleration density
! in 1-2/2d with periodic boundary conditions.
! the derivative of the current is calculated using the equations:
! dcu(1,kx) = dcu(1,kx)-sqrt(-1)*kx*vx*vy
! dcu(2,kx) = dcu(2,kx)-sqrt(-1)*kx*vx*vz
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
! on input:
! dcu(i,j) = complex acceleration density for fourier mode (j-1)
! on output:
! dcu(i,j) = transverse part of complex derivative of current for
! fourier mode (j-1)
! amu(1,j) = xy component of complex momentum flux
! amu(2,j) = xz component of complex momentum flux
! all for fourier mode (j-1)
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
! local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dcu(1,j) + dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dcu(2,j) + dkx*zt1
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,   &
     &nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! transverse electric field (or convolution of transverse electric field
! over particle shape), with periodic boundary conditions.
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! for isign = 0, input: isign,ax,affp,wp0,nx,nxvh, output:ffe
! for isign =/ 0, input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
! approximate flop count is: 25*nxc
! where nxc = nx/2 - 1
! if isign = 0, form factor array is prepared
! if isign = -1, smoothed transverse electric field is calculated
! using the equation:
! ey(kx) = -ci*ci*g(kx)*dcuy(kx)*s(kx)
! ez(kx) = -ci*ci*g(kx)*dcuz(kx)*s(kx)
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
! s(kx) = exp(-((kx*ax)**2+)/2), except for
! ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0,) = 0.
! if isign = 1, unsmoothed transverse electric field is calculated
! using the equation:
! ey(kx) = -ci*ci*g(kx)*dcuy(kx)
! ez(kx) = -ci*ci*g(kx)*dcuz(kx)
! dcu(i,j) = transverse part of complex derivative of current for
! fourier mode (j-1)
! eyz(1,j) = y component of complex transverse electric field
! eyz(2,j) = z component of complex transverse electric field
! all for fourier mode (j-1)
! aimag(ffe(j)) = finite-size particle shape factor s
! for fourier mode (j-1)
! real(ffe(j)) = potential green's function g
! for fourier mode (j-1)
! ax = half-width of particle in x direction
! affp = normalization constant = nx/np, where np=number of particles
! wp0 = normalized total plasma frequency squared
! ci = reciprocal of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
! this expression is valid only if the derivative of current is
! divergence-free
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer isign, nx, nxvh, nxhd
      real ax, affp, wp0, ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxhd)
! local data
      integer nxh, j
      real dnx, ci2, wpc, dkx, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
      wpc = wp0*ci2
! prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffe(j) = cmplx(affp*at2/(dkx*dkx+ wpc*at2*at2),at2)
   10 continue
      ffe(1) = cmplx(affp,1.0)
      return
! calculate smoothed transverse electric field and sum field energy
   20 if (isign.gt.0) go to 40
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*aimag(ffe(j))
      at2 = at2*at2
      eyz(1,j) = at1*dcu(1,j)
      eyz(2,j) = at1*dcu(2,j)
      wp = wp + at2*(dcu(1,j)*conjg(dcu(1,j))                           &
     &             + dcu(2,j)*conjg(dcu(2,j)))
   30 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
! calculate unsmoothed transverse electric field and sum field energy
   40 wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 50 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j))                           &
     &             + dcu(2,j)*conjg(dcu(2,j)))
   50 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
      end
!-----------------------------------------------------------------------
      subroutine POTP1(q,pot,ffc,we,nx,nxvh,nxhd)
! this subroutine solves 1d poisson's equation in fourier space for
! potential with periodic boundary conditions.
! input: q,ffc,nx,nxvh,nxhd, output: pot,we
! approximate flop count is: 3*nx
! potential is calculated using the equation:
! pot(k) = g(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
! g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
! pot(k=0) = pot(k=pi) = 0.
! q(j) = complex charge density for fourier mode j-1
! pot(j) = complex potential for fourier mode j-1
! real(ffc(j)) = potential green's function g for fourier mode j-1
! imag(ffc(j)) = finite-size particle shape factor s for fourier mode j-1
! electric field energy is also calculated, using
! we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
! where affp = normalization constant = nx/np,
! where np = number of particles
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
! nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real we
      complex q, pot, ffc
      dimension q(nxvh), pot(nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real at1, at2
      double precision wp
      nxh = nx/2
! calculate potential and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = real(ffc(j))
      at1 = at2*aimag(ffc(j))
      pot(j) = at2*q(j)
      wp = wp + at1*(q(j)*conjg(q(j)))
   10 continue
      pot(1) = cmplx(0.0,0.0)
      we = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine ELFIELD1(q,fx,ffc,we,nx,nxvh,nxhd)
! this subroutine solves 1d poisson's equation in fourier space for
! unsmoothed longitudinal electric field with periodic boundary
! conditions.
! input: q,ffc,nx,,nxvh,nxhd, output: fx,we
! approximate flop count is: 6*nx
! equation used is:
! fx(k) = -sqrt(-1)*k*g(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
! g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
! fx(k=0) = fx(k=pi) = 0.
! q(j) = complex charge density for fourier mode j-1
! fx(j) = complex electric field for fourier mode j-1
! real(ffc(j))= finite-size particle shape factor s for fourier mode j-1
! aimag(ffc(j)) = potential green's function g for fourier mode j-1
! electric field energy is also calculated, using
! we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
! affp = normalization constant = nx/np, where np = number of particles
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real we
      complex q, fx, ffc
      dimension q(nxvh), fx(nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real dnx, at1, at2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
! calculate force/charge and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      fx(j) = at2*cmplx(aimag(q(j)),-real(q(j)))
      wp = wp + at1*(q(j)*conjg(q(j)))
   10 continue
! mode number kx = 0
      fx(1) = cmplx(0.0,0.0)
      we = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine DIVF1(f,df,nx,ndim,nxvh)
! this subroutine calculates the divergence in fourier space
! input: all except df, output: df
! approximate flop count is: 15*nxc
! where nxc = nx/2 - 1
! the divergence is calculated using the equation:
! df(kx) = sqrt(-1)*kx*fx(kx)
! where kx = 2pi*j/nx, and j = fourier mode number,
! except for df(kx=pi) = 0.
! nx = system length in x direction
! ndim = number of field arrays
! nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, ndim, nxvh
      complex f, df
      dimension f(ndim,nxvh), df(nxvh)
! local data
      integer nxh, j
      real dnx, dkx
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
! calculate the divergence
! mode numbers 0 < kx < nx/2
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      df(j) = dkx*cmplx(-aimag(f(1,j)),real(f(1,j)))
   40 continue
      df(1) = cmplx(0.0,0.0)
      return
      end
!-----------------------------------------------------------------------
      subroutine GRADF1(df,f,nx,ndim,nxvh)
! this subroutine calculates the gradient in fourier space
! input: all except f, output: f
! approximate flop count is: 5*nxc
! where nxc = nx/2 - 1
! the gradient is calculated using the equations:
! fx(kx) = sqrt(-1)*kx*df(kx)
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! except for fx(kx=pi) = 0,
! nx = system length in x direction
! ndim = number of field arrays
! nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, ndim, nxvh
      complex df, f
      dimension df(nxvh), f(ndim,nxvh)
! local data
      integer nxh, i, j
      real dnx, dkx
      complex zero
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
! calculate the gradient
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      f(1,j) = dkx*cmplx(-aimag(df(j)),real(df(j)))
   10 continue
      f(1,1) = zero
      if (ndim.eq.1) return
! handle case of ndim > 1
      do 30 j = 1, nxh
      do 20 i = 2, ndim
      f(i,j) = zero
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine CURLF1(f,g,nx,nxvh)
! this subroutine calculates the curl in fourier space
! input: all except g, output: g
! approximate flop count is: 10*nxc
! where nxc = nx/2 - 1
! the curl is calculated using the equations:
! gy(kx) = -sqrt(-1)*kx*fz(kx)
! gz(kx) = sqrt(-1)*kx*fy(kx)
! where kx = 2pi*j/nx, and j = fourier mode number,
! except for gy(kx=pi) = gz(kx=pi) = 0,
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex f, g
      dimension f(2,nxvh), g(2,nxvh)
! local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
! calculate the curl
! mode numbers 0 < kx < nx/2
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(f(2,j)),real(f(2,j)))
      zt2 = cmplx(-aimag(f(1,j)),real(f(1,j)))
      g(1,j) = -dkx*zt1
      g(2,j) = dkx*zt2
   40 continue
      g(1,1) = zero
      g(2,1) = zero
      return
      end
!-----------------------------------------------------------------------
      subroutine AVPOT13(byz,ayz,nx,nxvh)
! this subroutine calculates 1-2/2d vector potential from magnetic field
! in fourier space with periodic boundary conditions.
! input: byz, nx, nxvh, output: ayz
! approximate flop count is: 10*nxc and nxc divides,
! where nxc = nx/2 - 1
! the vector potential is calculated using the equations:
! ay(kx) = -sqrt(-1)*bz(kx))/kx
! az(kx) = sqrt(-1)*by(kx)/kx
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
! byz(i,j) = i component of complex magnetic field
! ayz(i,j) = i component of complex vector potential
! all for fourier mode (j-1)
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex byz, ayz
      dimension byz(2,nxvh), ayz(2,nxvh)
! local data
      integer j, nxh
      real dnx, dkx, at2
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
! calculate vector potential
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = 1.0/dkx
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
      ayz(1,j) = -at2*zt1
      ayz(2,j) = at2*zt2
   10 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
!-----------------------------------------------------------------------
      subroutine CUAVE13(cuave,cunew,cuold,nx,nxvh)
! this subroutine averages current in fourier space for 1-2/2d code
! input: all except cuave, output: cuave
! cunew(i,j),cuold(i,j) = complex current densities to be averaged
! cuave(i,j) = average complex current density
! for component i, all for fourier mode (j-1)
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex cuave, cunew, cuold
      dimension cuave(2,nxvh), cuold(2,nxvh), cunew(2,nxvh)
! local data
      integer nxh, j
      nxh = nx/2
      do 10 j = 1, nxh
      cuave(1,j) = 0.5*(cunew(1,j) + cuold(1,j))
      cuave(2,j) = 0.5*(cunew(2,j) + cuold(2,j))
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine AVRPOT13(ayz,byz,ffc,ci,nx,nxvh,nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! the radiative part of the vector potential
! with periodic boundary conditions.
! input: all, output: ayz
! approximate flop count is: 24*nxc
! where nxc = nx/2 - 1
! the radiative vector potential is updated using the equations:
! ay(kx) = -(sqrt(-1)*kx*bz(kx) + affp*ci2*cuy(kx)*s(kx))/(kx*kx)
! az(kx) = (sqrt(-1)*kx*by(kx) - affp*ci2*cuz(kx)*s(kx))/(kx*kx)
! where kx = 2pi*j/nx, ci2 = ci*ci
! and s(kx) = exp(-((kx*ax)**2)
! j = fourier mode numbers, except for
! ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
! ayz(i,j) = on entry, complex current density cu
! ayz(i,j) = on exit, complex current radiative vector potential
! byz(i,j) = complex magnetic field
! for component i, all for fourier mode (j-1)
! real(ffc(1)) = affp = normalization constant = nx/np,
! where np=number of particles
! aimag(ffc()) = finite-size particle shape factor s,
! s(kx) = exp(-((kx*ax)**2)/2)
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! nx = system length in x direction
! nxvh = first dimension of field arrays, must be >= nxh
! nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci
      complex ayz, byz, ffc
      dimension ayz(2,nxvh), byz(2,nxvh)
      dimension ffc(nxhd)
! local data
      integer nxh, j
      real dnx, afc2, dkx, at1, at2
      complex zero, zt1, zt2
      if (ci.le.0.) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      afc2 = real(ffc(1))*ci*ci
      zero = cmplx(0.0,0.0)
! calculate the radiative vector potential
! mode numbers 0 < kx < nx/2
      do 20 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffc(j))
! update radiative vector potential
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
      ayz(1,j) = -at1*(dkx*zt1 + at2*ayz(1,j))
      ayz(2,j) = at1*(dkx*zt2 - at2*ayz(2,j))
   20 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
!-----------------------------------------------------------------------
      subroutine APOTP13(cu,ayz,ffc,ci,wm,nx,nxvh,nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! vector potential with periodic boundary conditions.
! input: cu,ffc,ci,nx,nxvh,nxhd, output: ayz,wm
! approximate flop count is: 23*nxc
! where nxc = nx/2 - 1
! vector potential is calculated using the equation:
! by(kx) = ci*ci*g(kx)*cuy(kx)
! bz(kx) = ci*ci*g(kx)*cuz(kx)
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! g(kx) = (affp/kx**2)*s(kx),
! s(kx) = exp(-((kx*ax)**2+)/2), except for
! by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
! cu(i,j) = complex current density for fourier mode (j-1)
! ayz(1,j) = y component of complex vector potential
! ayz(2,j) = z component of complex vector potential
! all for fourier mode (j-1)
! aimag(ffc(j)) = finite-size particle shape factor s
! real(ffc(j)) = potential green's function g
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
! where affp = normalization constant = nx/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, ayz, ffc
      dimension cu(2,nxvh), ayz(2,nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real ci2, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate vector potential and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = ci2*real(ffc(j))
      at1 = at2*aimag(ffc(j))
      ayz(1,j) = at2*cu(1,j)
      ayz(2,j) = at2*cu(2,j)
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      wm = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine ETFIELD13(dcu,eyz,ffe,ci,wf,nx,nxvh,nxhd)
! this subroutine solves 1-2/2d poisson's equation in fourier space for
! unsmoothed transverse electric field with periodic boundary conditions
! using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
! A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
! input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
! approximate flop count is: 25*nxc
! where nxc = nx/2 - 1
! unsmoothed transverse electric field is calculated using the equation:
! ey(kx) = -ci*ci*g(kx)*dcuy(kx)
! ez(kx) = -ci*ci*g(kx)*dcuz(kx)
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
! s(kx) = exp(-((kx*ax)**2+)/2), except for
! ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0,) = 0.
! dcu(i,j) = transverse part of complex derivative of current for
! fourier mode (j-1)
! eyz(1,j) = y component of complex transverse electric field
! eyz(2,j) = z component of complex transverse electric field
! all for fourier mode (j-1)
! aimag(ffe(j)) = finite-size particle shape factor s
! for fourier mode (j-1)
! real(ffe(j)) = potential green's function g
! for fourier mode (j-1)
! ci = reciprocal of velocity of light
! transverse electric field energy is also calculated, using
! wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
! affp = normalization constant = nx/np, where np = number of particles
! this expression is valid only if the derivative of current is
! divergence-free
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxhd)
! local data
      integer nxh, j
      real ci2, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate unsmoothed transverse electric field and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j))                           &
     &             + dcu(2,j)*conjg(dcu(2,j)))
   10 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
      end
!-----------------------------------------------------------------------
      subroutine SMOOTH1(q,qs,ffc,nx,nxvh,nxhd)
! this subroutine provides a 1d scalar smoothing function
! in fourier space, with periodic boundary conditions.
! input: q,ffc,nx,nxvh,nxhd, output: qs, flop count = nx
! smoothing is calculated using the equation:
! qs(k) = q(k)*s(k), where k = 2pi*j/nx, j=fourier mode,
! g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2),
! except for qs(k=pi) = 0.
! q(j) = complex charge density for fourier mode j-1
! qs(j) = complex smoothed charge density for fourier mode j-1
! real(ffc(j,k)) = potential green's function g
! aimag(ffc(j,k)) = finite-size particle shape factor s
! for fourier mode j-1
! nx = system length in x direction
! nxvh = first dimension of scalar field arrays, must be >= nx/2
! nxhd = first dimension of form factor array, must be >= nx/2
      implicit none
      integer nx, nxvh, nxhd
      complex q, qs, ffc
      dimension q(nxvh), qs(nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real at1
      nxh = nx/2
! calculate smoothing
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = aimag(ffc(j))
      qs(j) = at1*q(j)
   10 continue
      qs(1) = cmplx(aimag(ffc(1))*real(q(1)),0.0)
      return
      end
!-----------------------------------------------------------------------
      subroutine SMOOTH13(cu,cus,ffc,nx,nxvh,nxhd)
! this subroutine provides a 1d vector smoothing function
! in fourier space, with periodic boundary conditions.
! input: cu,ffc,nx,nxvh, output: cus
! approximate flop count is: 4*nxc, where nxc = nx/2 - 1
! smoothing is calculated using the equation:
! cusy(kx) = cuy(kx)*s(kx)
! cusz(kx) = cuz(kx)*s(kx)
! where kx = 2pi*j/nx, and j = fourier mode numbers,
! g(kx) = (affp/kx**2)*s(kx),
! where affp = normalization constant = nx/np,
! and np=number of particles
! s(kx) = exp(-((kx*ax)**2+)/2), except for
! cusy(kx=pi) = cusz(kx=pi) = 0
! cu(i,j) = ith-component of complex current density
! cus(i,j) = ith-component of complex smoothed charge density
! all for fourier mode (j-1)
! aimag(ffc(j)) = finite-size particle shape factor s
! real(ffc(j)) = potential green's function g
! for fourier mode (j-1)
! nx = system length in x direction
! nxvh = second dimension of field arrays, must be >= nxh
! nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      complex cu, cus, ffc
      dimension cu(2,nxvh), cus(2,nxvh), ffc(nxhd)
! local data
      integer j, nxh
      real at1
      nxh = nx/2
! calculate smoothing
! mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = aimag(ffc(j))
      cus(1,j) = at1*cu(1,j)
      cus(2,j) = at1*cu(2,j)
   10 continue
      at1 = aimag(ffc(1))
      cus(1,1) = cmplx(at1*real(cu(1,1)),0.0)
      cus(2,1) = cmplx(at1*real(cu(2,1)),0.0)
      return
      end
!-----------------------------------------------------------------------
      subroutine RDMODES1(pot,pott,nx,modesx,nxvh,modesxd)
! this subroutine extracts lowest order modes from packed complex array
! pot and stores them into a location in an unpacked complex array pott
! vector array pott
! modes stored: kx=(0,1,...,NX/2)
! nx = system length in x direction
! modesx = number of modes to store in x direction,
! where modesx <= nx/2+1
! nxvh = dimension of input array pot, nxvh >= nx/2
! modesxd = dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx, nxvh, modesxd
      complex pot, pott
      dimension pot(nxvh), pott(modesxd)
! local data
      integer nxh, jmax, j, j1
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
! mode numbers 0 < kx < nx/2
      do 10 j = 2, jmax
      pott(j) = pot(j)
   10 continue
! mode numbers kx = 0, nx/2
      pott(1) = cmplx(real(pot(1)),0.0)
      if (modesx.gt.nxh) then
         pott(j1) = cmplx(aimag(pot(1)),0.0)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex array pott and stores them into a packed complex
! array pot
! modes stored: kx=(0,1,...,NX/2)
! nx = system length in x direction
! modesx = number of modes to store in x direction,
! where modesx <= nx/2+1
! nxvh = dimension of input array pot, nxvh >= nx/2
! modesxd = dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx, nxvh, modesxd
      complex pot, pott
      dimension pot(nxvh), pott(modesxd)
! local data
      integer nxh, jmax, j, j1
      complex zero
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2
      do 10 j = 2, jmax
      pot(j) = pott(j)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j) = zero
   20 continue
! mode numbers kx = 0, nx/2
      pot(1) = cmplx(real(pott(1)),0.0)
      if (modesx.gt.nxh) then
         pot(1) = cmplx(real(pot(1)),real(pott(j1)))
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine RDVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
! this subroutine extracts lowest order modes from packed complex vector
! array vpot and stores them into a location in an unpacked complex
! vector array vpott
! modes stored: kx=(0,1,...,NX/2)
! nx = system length in x direction
! modesx = number of modes to store in x direction,
! where modesx <= nx/2+1
! ndim = number of field arrays, must be >= 1
! nxvh = second dimension of input array vpot, nxvh >= nx/2
! modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, ndim, nxvh, modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(ndim,modesxd)
! local data
      integer nxh, jmax, i, j, j1
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
! mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpott(i,j) = vpot(i,j)
   10 continue
   20 continue
! mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      vpott(i,1) = cmplx(real(vpot(i,1)),0.0)
      if (modesx.gt.nxh) then
         vpott(i,j1) = cmplx(aimag(vpot(i,1)),0.0)
      endif
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex vector array vpott and stores them into a packed
! complex vector array vpot
! modes stored: kx=(0,1,...,NX/2)
! nx = system length in x direction
! modesx = number of modes to store in x direction,
! where modesx <= nx/2+1
! ndim = number of field arrays, must be >= 1
! nxvh = second dimension of input array vpot, nxvh >= nx/2
! modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, ndim, nxvh, modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(ndim,modesxd)
! local data
      integer nxh, jmax, i, j, j1
      complex zero
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j) = vpott(i,j)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j) = zero
   30 continue
   40 continue
! mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1) = cmplx(real(vpott(i,1)),0.0)
      if (modesx.gt.nxh) then
         vpot(i,1) = cmplx(real(vpot(i,1)),real(vpott(i,j1)))
      endif
   50 continue
      return
      end
