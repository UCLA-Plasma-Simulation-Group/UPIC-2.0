!-----------------------------------------------------------------------
! Library for processing complex periodic 1d scalar data with OpenMP
! WRMODES1 extracts lowest order scalar modes from a location in an
!          unpacked array and stores them into a packed array
! WRVMODES1 extracts lowest order vector modes from a location in an
!           unpacked array and stores them into a packed array
! CSPECT1 performs frequency analysis of complex time series
! ICSPECT1 performs incremental frequency analysis of complex time
!          series for one time step
! IVCSPECT1 performs incremental frequency analysis of complex vector
!           time series for one time step
! WFFT1RINIT calculates tables needed by 1d FFTs
! FFT1RXX performs performs scalar 1d real/complex FFT
! FFT1R2X performs 2 component vector 1d real/complex FFT
! FFT1R3X performs 3 component vector 1d real/complex FFT
! DIVF1 calculates the divergence in fourier space
! written by Viktor K. Decyk, UCLA
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
!-----------------------------------------------------------------------
      subroutine CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
! this subroutine performs frequency analysis of complex time series,
! pkw(w,k,1:2) = |(1/nt)*sum {fc(t,k)*exp(sqrt(-1)*w*(t-t0))}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed, real and imaginary
! parts stored adjacent, and wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! t0 = starting time value
! dt = time step
! nt = number of input data points, 
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = dimension of input array fc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, ntd, iwd, modesxd
      real t0, dt
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(ntd,modesxd)
! local data
      integer i, j, k
      real anl, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4, cwdt, swdt
      anl = 1.0/real(nt)
! loop over frequencies
      do 30 j = 1, iw
      at3 = wm(j)*dt
      cwdt = dcos(at3)
      swdt = dsin(at3)
! loop over modes
      do 20 k = 1, modesx
      at3 = wm(j)*t0
      at1 = dcos(at3)
      at2 = -dsin(at3)
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      sum4 = 0.0d0
! loop over time
      do 10 i = 1, nt
      fr = real(fc(i,k))
      fi = aimag(fc(i,k))
      sum1 = sum1 + fr*at1
      sum2 = sum2 + fi*at1
      sum3 = sum3 + fi*at2
      sum4 = sum4 + fr*at2
      at3 = at1*cwdt - at2*swdt
      at2 = at2*cwdt + at1*swdt
      at1 = at3
   10 continue
      at1 = anl*(sum1 - sum3)
      at2 = anl*(sum2 + sum4)
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd&
     &,modesxd)
! this subroutine performs incremental frequency analysis of complex
! time series for one time step
! pkw(w,k,1:2) = |(anorm/nt)*sum {fc(k)*exp(sqrt(-1)*w*(time-t0)}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent, and
! wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! pks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! nx = system length in x direction
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
! norm = 1 for potential or norm = -1 for density gives spectrum as
!        electric field energy
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = dimension of input array fc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, nx, norm, iwd, modesxd
      real time, t0
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(modesxd)
      double precision pks
      dimension pks(4,modesxd,iwd)
! local data
      integer j, k
      real anl, dnx, anorm, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4
      anl = 1.0/real(nt)
      dnx = 6.28318530717959/real(nx)
      anorm = anl
! loop over frequencies
      do 20 j = 1, iw
! loop over modes
      do 10 k = 1, modesx
      at3 = wm(j)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
! add contribution for current time
      fr = real(fc(k))
      fi = aimag(fc(k))
      sum1 = pks(1,k,j) + fr*at1
      sum2 = pks(2,k,j) + fi*at1
      sum3 = pks(3,k,j) + fi*at2
      sum4 = pks(4,k,j) + fr*at2
! save accumulation for next time
      pks(1,k,j) = sum1
      pks(2,k,j) = sum2
      pks(3,k,j) = sum3
      pks(4,k,j) = sum4
! normalize
      if (k.gt.1) then
         if (norm==1) then
            anorm = anl*(dnx*real(k - 1))
         else if (norm==(-1)) then
            anorm = anl/(dnx*real(k - 1))
         endif
      endif
! calculate spectrum for accumulated data
      at1 = anorm*(sum1 - sum3)
      at2 = anorm*(sum2 + sum4)
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anorm*(sum1 + sum3)
      at2 = anorm*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   10 continue
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx,norm&
     &,iwd,modesxd)
! this subroutine performs incremental frequency analysis of complex
! vector time series for one time step
! vpkw(1:2,w,k,1:2) = |(anorm/nt)*sum {fvc(1:2,k)*
!                                      exp(sqrt(-1)*w*(time-t0)}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent, and
! wm(w) contains the (positive) frequencies.
! on output, vpkw(1:2,:,:,1) contains result for positive frequencies,
! vpkw(1:2,:,:,2) the negative frequencies.
! vpks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! nx = system length in x direction
! norm = (-1,0,1) = normalize with (inverse curl,null,curl) op
! norm = 1 for vector potential or norm = -1 for current density gives
!        spectrum as magnetic field energy
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = second dimension of input array fvc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, nx, norm, iwd, modesxd
      real time, t0
      real wm, vpkw
      dimension wm(iwd), vpkw(2,modesxd,iwd,2)
      complex fvc
      dimension fvc(2,modesxd)
      double precision vpks
      dimension vpks(2,4,modesxd,iwd)
! local data
      integer i, j, k
      real anl, dnx, anorm, fr, fi
      double precision at1, at2, at3, at4, sum1, sum2, sum3, sum4
      anl = 1.0/real(nt)
      dnx = 6.28318530717959/real(nx)
      anorm = anl
! loop over frequencies
      do 30 j = 1, iw
! loop over modes
      do 20 k = 1, modesx
      at3 = wm(j)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
      do 10 i = 1, 2
! add contribution for current time
      fr = real(fvc(i,k))
      fi = aimag(fvc(i,k))
      sum1 = vpks(i,1,k,j) + fr*at1
      sum2 = vpks(i,2,k,j) + fi*at1
      sum3 = vpks(i,3,k,j) + fi*at2
      sum4 = vpks(i,4,k,j) + fr*at2
! save accumulation for next time
      vpks(i,1,k,j) = sum1
      vpks(i,2,k,j) = sum2
      vpks(i,3,k,j) = sum3
      vpks(i,4,k,j) = sum4
! normalize
      if (k.gt.1) then
         if (norm==1) then
            anorm = anl*(dnx*real(k - 1))
         else if (norm==(-1)) then
            anorm = anl/(dnx*real(k - 1))
         endif
      endif
! calculate spectrum for accumulated data
      at3 = anorm*(sum1 - sum3)
      at4 = anorm*(sum2 + sum4)
      vpkw(i,k,j,1) = at3*at3 + at4*at4
      at3 = anorm*(sum1 + sum3)
      at4 = anorm*(sum2 - sum4)
      vpkw(i,k,j,2) = at3*at3 + at4*at4
   10 continue
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
! this subroutine calculates tables needed by a one dimensional
! real to complex fast fourier transform and its inverse.
! input: indx, nxhd
! output: mixup, sct
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! indx = exponent which determines length in x direction,
! where nx=2**indx
! nxhd = nx/2
! written by viktor k. decyk, ucla
      implicit none
      integer indx, nxhd
      integer mixup
      complex sct
      dimension mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh
      integer j, k, lb, ll, jb, it
      real dnx, arg
      indx1 = indx - 1
      nx = 2**indx
      nxh = nx/2
! bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
! sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 30 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
! this subroutine performs a one dimensional real to complex fast
! fourier transform and its inverse, using complex arithmetic
! for isign = (-1,1), input: all except t, output: f, t
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = nx/2
! f = input and output data
! t = complex scratch array
! indx = power of 2 which determines length of transform, nx = 2**indx
! if isign = -1, an inverse fourier transform is performed
! f(n) = (1/nx)*sum(f(j)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(j) = sum(f(n)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! fourier coefficients are stored as follows:
! f(1) = real part of mode 0, f(2) = real part of mode nx/2
! f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
! written by viktor k. decyk, ucla
! scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(nxd), t(nxhd), mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2
      real ani
      complex t1, t2
      if (isign.eq.0) return
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 70
! inverse fourier transform
! bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(j) = cmplx(f(2*j1-1),f(2*j1))
   10 continue
! transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      ani = 1./real(2*nx)
      do 50 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),-real(sct(j)))
      t(j) = ani*(t1 + t2)
      t(nxh2-j) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.*ani
      t(nxhh+1) = ani*conjg(t(nxhh+1))
      t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),                        &
     &                 real(t(1)) - aimag(t(1)))
! move to real destination
      do 60 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
   60 continue
      return
! forward fourier transform
! move to complex temporary
   70 do 80 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
   80 continue
! scramble coefficients
      do 90 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),real(sct(j)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
   90 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
! bit-reverse array elements to real destination
      do 100 j = 1, nxh
      j1 = mixup(j)
      f(2*j-1) = real(t(j1))
      f(2*j) = aimag(t(j1))
  100 continue
! move back to complex temporary
      do 110 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
  110 continue
! transform
      do 140 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 130 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 120 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  120 continue
  130 continue
  140 continue
! move to real destination
      do 150 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
  150 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
! this subroutine performs two one dimensional real to complex fast
! fourier transforms and their inverses, using complex arithmetic
! for isign = (-1,1), input: all except t, output: f, t
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = nx/2
! f = input and output data
! t = complex scratch array
! indx = power of 2 which determines length of transform, nx = 2**indx
! if isign = 0, the fft tables are prepared
! if isign = -1, an inverse fourier transform is performed
! f(1:2,n) = (1/nx)*sum(f(1:2,j)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(1:2,j) = sum(f(1:2,n)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! fourier coefficients are stored as follows:
! f(1:2,j) = real,imaginary part of mode j-1, 1 <= j <= nx/2
! except for
! real(f(1:2,1)) = real part of mode 0,
! aimag(f(1:2,1)) = real part of mode nx/2
! written by viktor k. decyk, ucla
! scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(2,nxd), t(2,nxhd), mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 90
! inverse fourier transform
! bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
   10 continue
! transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      ani = 1.0/real(2*nx)
      do 60 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 50 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
   50 continue
   60 continue
      ani = 2.*ani
      do 70 jj = 1, 2
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),               &
     &                    real(t(jj,1)) - aimag(t(jj,1)))
   70 continue
! move to complex destination
      do  80 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(1,2*j) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
   80 continue
      return
! forward fourier transform
! move complex source to complex temporary
   90 do 100 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(1,2*j),f(2,2*j))
  100 continue
! scramble coefficients
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 110 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 2
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),                   &
     &                real(t(jj,1)) - aimag(t(jj,1)))
  130 continue
! bit-reverse array elements to real destination
      do 140 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
  140 continue
! move back to complex temporary
      do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
  150 continue
! transform
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
  160 continue
  170 continue
  180 continue
! move to real destination
      do 190 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
  190 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
! this subroutine performs three one dimensional real to complex fast
! fourier transforms and their inverses, using complex arithmetic
! for isign = (-1,1), input: all except t, output: f, t
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = nx/2
! f = input and output data
! t = complex scratch array
! indx = power of 2 which determines length of transform, nx = 2**indx
! if isign = 0, the fft tables are prepared
! if isign = -1, an inverse fourier transform is performed
! f(1:3,n) = (1/nx)*sum(f(1:3,j)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(1:3,j) = sum(f(1:3,n)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! fourier coefficients are stored as follows:
! f(1) = real part of mode 0, f(2) = real part of mode nx/2
! f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
! written by viktor k. decyk, ucla
! scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(3,nxd), t(3,nxhd), mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real ani
      complex t1, t2, t3, t4
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 90
! inverse fourier transform
! bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
      t(3,j) = cmplx(f(3,2*j1-1),f(3,2*j1))
   10 continue
! transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize result
      ani = 1./real(2*nx)
      do 60 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 50 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
   50 continue
   60 continue
      ani = 2.0*ani
      do 70 jj = 1, 3
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),               &
     &                    real(t(jj,1)) - aimag(t(jj,1)))
   70 continue
! move to complex destination
      do 80 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(3,2*j-1) = real(t(2,j))
      f(1,2*j) = aimag(t(2,j))
      f(2,2*j) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
   80 continue
      return
! forward fourier transform
! move complex source to complex temporary
   90 do 100 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(3,2*j-1),f(1,2*j))
      t(3,j) = cmplx(f(2,2*j),f(3,2*j))
  100 continue
! scramble coefficients
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 110 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 3
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),                   &
     &                real(t(jj,1)) - aimag(t(jj,1)))
  130 continue
! bit-reverse array elements to real destination
      do 140 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
      f(3,2*j-1) = real(t(3,j1))
      f(3,2*j) = aimag(t(3,j1))
  140 continue
! move back to complex temporary
      do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
      t(3,j) = cmplx(f(3,2*j-1),f(3,2*j))
  150 continue
! transform
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
  160 continue
  170 continue
  180 continue
! move to real destination
      do 190 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
      f(3,2*j-1) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
  190 continue
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
