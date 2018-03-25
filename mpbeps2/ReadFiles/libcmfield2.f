!-----------------------------------------------------------------------
! Library for processing complex periodic 2d scalar data with OpenMP
! WRMODES2 places selected 2d fourier components into complex scalar
!          array.
! WRVMODES2 places selected 2d fourier components into complex vector
!           array.
! CSPECT2 performs frequency analysis of complex time series
! ICSPECT1 performs incremental frequency analysis of complex time
!          series for one time step
! IVCSPECT2 performs incremental frequency analysis of complex vector
!           time series for one time step
! WFFT2RINIT calculates tables needed by a two dimensional real to
!            complex fast fourier transform and its inverse.
! WFFT2RMX performs real to complex fft and its inverse for scalar array,
!          with packed data.
! WFFT2RM2 performs real to complex fft and its inverse for 2 component
!          vector array, with packed data.
! WFFT2RM3 performs real to complex fft and its inverse for 3 component
!          vector array, with packed data.
! FFT2RMXX performs x part of scalar 2d real/complex FFT
! FFT2RMXY performs y part of scalar 2d real/complex FFT
! FFT2RM2X performs x part of 2 component vector 2d real/complex FFT
! FFT2RM2Y performs y part of 2 component vector 2d real/complex FFT
! FFT2RM3X performs x part of 3 component vector 2d real/complex FFT
! FFT2RM3Y performs y part of 3 component vector 2d real/complex FFT
! written by Viktor K. Decyk, UCLA
!-----------------------------------------------------------------------
      subroutine WRMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,modesxd,&
     &modesyd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex array pott and stores them into a packed complex
! array pot
! modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! nxvh = first dimension of input array pot, nxvh >= nx/2
! nyv = second dimension of input array pot, nyv >= ny
! modesxd = first dimension of output array pott, modesxd >= modesx
! modesyd = second dimension of output array pott,
! where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, nxvh, nyv, modesxd, modesyd
      complex pot, pott
      dimension pot(nxvh,nyv), pott(modesxd,modesyd)
! local data
      integer nxh, nyh, kmax, jmax, ny2, j, k, j1, k1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1)
      do 30 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pot(j,k) = pott(j,2*k-2)
      pot(j,k1) = pott(j,2*k-1)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j,k) = zero
      pot(j,k1) = zero
   20 continue
! mode numbers kx = 0, nx/2
      pot(1,k) = pott(1,2*k-2)
      pot(1,k1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1) = conjg(pott(j1,2*k-2))
      endif
   30 continue
!$OMP END PARALLEL DO
      do 50 k = kmax+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(j,k) = zero
      pot(j,k1) = zero
   40 continue
   50 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, jmax
      pot(j,1) = pott(j,1)
      pot(j,k1) = zero
   60 continue
      do 70 j = jmax+1, nxh
      pot(j,1) = zero
      pot(j,k1) = zero
   70 continue
      pot(1,1) = cmplx(real(pott(1,1)),0.0)
      pot(1,k1) = zero
      if (modesx.gt.nxh) then
         pot(1,1) = cmplx(real(pot(1,1)),real(pott(j1,1)))
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 80 j = 2, jmax
         pot(j,k1) = pott(j,ny)
   80    continue
         pot(1,k1) = cmplx(real(pott(1,ny)),0.0)
         if (modesx.gt.nxh) then
            pot(1,k1) = cmplx(real(pot(1,k1)),real(pott(j1,ny)))
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh,nyv,&
     &modesxd,modesyd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex vector array vpott and stores them into a packed
! complex vector array vpot
! modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
! nx/ny = system length in x/y direction
! modesx/modesy = number of modes to store in x/y direction,
! where modesx <= nx/2+1, modesy <= ny/2+1
! ndim = number of field arrays, must be >= 1
! nxvh = second dimension of input array vpot, nxvh >= nx/2
! nyv = third dimension of input array vpot, nyv >= ny
! modesxd = second dimension of output array vpott, modesxd >= modesx
! modesyd = third dimension of output array vpott,
! where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, ndim, nxvh, nyv, modesxd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv), vpott(ndim,modesxd,modesyd)
! local data
      integer nxh, nyh, kmax, jmax, ny2, i, j, k, j1, k1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 60 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j,k) = vpott(i,j,2*k-2)
      vpot(i,j,k1) = vpott(i,j,2*k-1)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   30 continue
   40 continue
! mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1,k) = vpott(i,1,2*k-2)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,k1) = conjg(vpott(i,j1,2*k-2))
      endif
   50 continue
   60 continue
!$OMP END PARALLEL DO
      do 90 k = kmax+1, nyh
      k1 = ny2 - k
      do 80 j = 1, nxh
      do 70 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   70 continue
   80 continue
   90 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 110 j = 2, jmax
      do 100 i = 1, ndim
      vpot(i,j,1) = vpott(i,j,1)
      vpot(i,j,k1) = zero
  100 continue
  110 continue
      do 130 j = jmax+1, nxh
      do 120 i = 1, ndim
      vpot(i,j,1) = zero
      vpot(i,j,k1) = zero
  120 continue
  130 continue
      do 140 i = 1, ndim
      vpot(i,1,1) = cmplx(real(vpott(i,1,1)),0.0)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,1) = cmplx(real(vpot(i,1,1)),real(vpott(i,j1,1)))
      endif
  140 continue
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 160 j = 2, jmax
         do 150 i = 1, ndim
         vpot(i,j,k1) = vpott(i,j,ny)
  150    continue
  160    continue
         do 170 i = 1, ndim
         vpot(i,1,k1) = cmplx(real(vpott(i,1,ny)),0.0)
         if (modesx.gt.nxh) then
            vpot(i,1,k1) = cmplx(real(vpot(i,1,k1)),                    &
     &                           real(vpott(i,j1,ny)))
         endif
  170    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine CSPECT2(fc,wm,pkw,t0,dt,nt,iw,modesx,modesy2,ntd,iwd,  &
     &modesxd,modesyd)
! this subroutine performs frequency analysis of complex time series,
! pkw(j,k,w,1:2) = |(1/nt)*sum {fc(t,j,k)*exp(sqrt(-1)*w*(t-t0))}|**2
! where j, k represent the x, y modes,
! and wm(m) stores the positive frequency values w,
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed, real and imaginary
! parts stored adjacent, and wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,:,1) contains result for positive frequencies
! pkw(:,:,:,2) the negative frequencies.
! t0 = starting time value
! dt = time step
! nt = number of input data points
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! modesy2 = number of modes in y direction = 2*modesy - 1
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = dimension of input array fc, modesxd >= modesx
! modesyd = dimension of input array fc, modesyd >= modesy2
      implicit none
      integer nt, iw, modesx, modesy2, ntd, iwd, modesxd, modesyd
      real t0, dt
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,modesyd,iwd,2)
      complex fc
      dimension fc(ntd,modesxd,modesyd)
! local data
      integer i, j, k, m
      real anl, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4, cwdt, swdt
      anl = 1.0/real(nt)
! loop over frequencies
      do 40 m = 1, iw
      at3 = wm(m)*dt
      cwdt = dcos(at3)
      swdt = dsin(at3)
! loop over modes
      do 30 k = 1, modesy2
      do 20 j = 1, modesx
      at3 = wm(m)*t0
      at1 = dcos(at3)
      at2 = -dsin(at3)
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      sum4 = 0.0d0
! loop over time
      do 10 i = 1, nt
      fr = real(fc(i,j,k))
      fi = aimag(fc(i,j,k))
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
      pkw(j,k,m,1) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      pkw(j,k,m,2) = at1*at1 + at2*at2
   20 continue
   30 continue
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ICSPECT2(fc,wm,pkw,pks,time,t0,nt,iw,modesx,modesy,iwd,&
     &modesxd,modesyd)
! this subroutine performs incremental frequency analysis of complex
! time series for one time step
! pkw(j,k,w,1:2) = |(1/nt)*sum {fc(j,k)*exp(sqrt(-1)*w*(time-t0)}|**2
! where j, k represent the x, y modes,
! and wm(m) stores the positive frequency values w
! k modes are reordered in increasing mode number, negative to positive
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent
! on output, pkw(:,:,:,1) contains result for positive frequencies,
! pkw(:,:,:,2) the negative frequencies.
! pks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx/modesy = number of modes in x/y direction
! iwd = dimension of frequency array iwd >= iw
! modesxd = first dimension of input array fc, modesxd >= modesx
! modesyd = second dimension of input array fc, modesyd >= 2*modesy - 1
      implicit none
      integer nt, iw, modesx, modesy, iwd, modesxd, modesyd
      real time, t0
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,modesyd,iwd,2)
      complex fc
      dimension fc(modesxd,modesyd)
      double precision pks
      dimension pks(4,modesxd,modesyd,iwd)
! local data
      integer j, k, m, modesy2, kk, ks
      real anl, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4
      modesy2 = 2*modesy - 1
      anl = 1.0/real(nt)
! loop over frequencies
      do 30 m = 1, iw
! loop over modes
      do 20 k = 1, modesy2
      kk = k/2
      ks = k - 2*kk
      do 10 j = 1, modesx
      at3 = wm(m)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
! add contribution for current time
      fr = real(fc(j,k))
      fi = aimag(fc(j,k))
      sum1 = pks(1,j,k,m) + fr*at1
      sum2 = pks(2,j,k,m) + fi*at1
      sum3 = pks(3,j,k,m) + fi*at2
      sum4 = pks(4,j,k,m) + fr*at2
! save accumulation for next time
      pks(1,j,k,m) = sum1
      pks(2,j,k,m) = sum2
      pks(3,j,k,m) = sum3
      pks(4,j,k,m) = sum4
! calculate spectrum for accumulated data
! reorder y modes in increasing order of mode number
! positive frequencies
      at1 = anl*(sum1 - sum3)
      at2 = anl*(sum2 + sum4)
! move positive modes
      if (ks==0) then
         pkw(j,modesy+kk,m,1) = at1*at1 + at2*at2
! move negative modes
      else
         pkw(j,modesy-kk,m,1) = at1*at1 + at2*at2
      endif
! negative frequencies
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
! move positive modes
      if (ks==0) then
         pkw(j,modesy+kk,m,2) = at1*at1 + at2*at2
! move negative modes
      else
         pkw(j,modesy-kk,m,2) = at1*at1 + at2*at2
      endif
   10 continue
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine IVCSPECT2(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,modesy,&
     &ndim,iwd,modesxd,modesyd)
! this subroutine performs incremental frequency analysis of complex
! vector time series for one time step
! vpkw(1:ndim,j,k,w,1:2) = |(1/nt)*sum {fvc(1:ndim,j,k)*
!                                      exp(sqrt(-1)*w*(time-t0)}|**2
! where j, k represent the x, y modes,
! and wm(m) stores the positive frequency values w
! k modes are reordered in increasing mode number, negative to positive
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fvc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent
! on output, vpkw(:,:,:,:,1) contains result for positive frequencies,
! vpkw(:,:,:,:,2) the negative frequencies.
! vpks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! modesx/modesy = number of modes in x/y direction
! ndim = first dimension of input array fvc
! iwd = dimension of frequency array iwd >= iw
! modesxd = second dimension of input array fvc, modesxd >= modesx
! modesyd = third dimension of input array fvc, modesyd >= 2*modesy - 1
      implicit none
      integer nt, iw, modesx, modesy, ndim, iwd, modesxd, modesyd
      real time, t0
      real wm, vpkw
      dimension wm(iwd), vpkw(ndim,modesxd,modesyd,iwd,2)
      complex fvc
      dimension fvc(ndim,modesxd,modesyd)
      double precision vpks
      dimension vpks(ndim,4,modesxd,modesyd,iwd)
! local data
      integer i, j, k, m, modesy2, kk, ks
      real anl, fr, fi
      double precision at1, at2, at3, at4, sum1, sum2, sum3, sum4
      modesy2 = 2*modesy - 1
      anl = 1.0/real(nt)
! loop over frequencies
      do 40 m = 1, iw
! loop over modes
      do 30 k = 1, modesy2
      kk = k/2
      ks = k - 2*kk
      do 20 j = 1, modesx
      at3 = wm(m)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
      do 10 i = 1, ndim
! add contribution for current time
      fr = real(fvc(i,j,k))
      fi = aimag(fvc(i,j,k))
      sum1 = vpks(i,1,j,k,m) + fr*at1
      sum2 = vpks(i,2,j,k,m) + fi*at1
      sum3 = vpks(i,3,j,k,m) + fi*at2
      sum4 = vpks(i,4,j,k,m) + fr*at2
! save accumulation for next time
      vpks(i,1,j,k,m) = sum1
      vpks(i,2,j,k,m) = sum2
      vpks(i,3,j,k,m) = sum3
      vpks(i,4,j,k,m) = sum4
! calculate spectrum for accumulated data
! reorder y modes in increasing order of mode number
! positive frequencies
      at3 = anl*(sum1 - sum3)
      at4 = anl*(sum2 + sum4)
! move positive modes
      if (ks==0) then
         vpkw(i,j,modesy+kk,m,1) = at3*at3 + at4*at4
! move negative modes
      else
         vpkw(i,j,modesy-kk,m,1) = at3*at3 + at4*at4
      endif
! negative frequencies
      at3 = anl*(sum1 + sum3)
      at4 = anl*(sum2 - sum4)
! move positive modes
      if (ks==0) then
         vpkw(i,j,modesy+kk,m,2) = at3*at3 + at4*at4
! move negative modes
      else
         vpkw(i,j,modesy-kk,m,2) = at3*at3 + at4*at4
      endif
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
! this subroutine calculates tables needed by a two dimensional
! real to complex fast fourier transform and its inverse.
! input: indx, indy, nxhyd, nxyhd
! output: mixup, sct
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nxy, nxhy, nxyh
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
! bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
! sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/real(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   &
     &nxyhd)
! wrapper function for real to complex fft, with packed data
! parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
! perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RM2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   &
     &nxyhd)
! wrapper function for 2 2d real to complex ffts, with packed data
! parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
! perform y fft
         call FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform y fft
         call FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! perform x fft
         call FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   &
     &nxyhd)
! wrapper function for 3 2d real to complex ffts
! parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
! perform y fft
         call FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform y fft
         call FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    &
     &nxhyd,nxyhd)
! perform x fft
         call FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd&
     &,nxyhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the x part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of y,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform in x is performed
! f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform in x is performed
! f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = first dimension of f >= nx/2
! nyd = second dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1,1)) = real part of mode nx/2,0 and
! imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nrxb
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 70
! inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,ani,t1,t2,t3)
      do 60 i = nyi, nyt
! bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
      endif
   10 continue
! then transform in x
      do 40 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 50 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = ani*(t1 + t2)
      f(nxh2-j,i) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.0*ani
      f(nxhh+1,i) = ani*conjg(f(nxhh+1,i))
      f(1,i) = ani*cmplx(real(f(1,i)) + aimag(f(1,i)),                  &
     &                   real(f(1,i)) - aimag(f(1,i)))
   60 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
   70 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2,t3)
      do 130 i = nyi, nyt
! scramble coefficients
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = t1 + t2
      f(nxh2-j,i) = conjg(t1 - t2)
   80 continue
      f(nxhh+1,i) = 2.0*conjg(f(nxhh+1,i))
      f(1,i) = cmplx(real(f(1,i)) + aimag(f(1,i)),                      &
     &               real(f(1,i)) - aimag(f(1,i)))
! bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
      endif
   90 continue
! then transform in x
      do 120 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the y part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of x,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform in y is performed
! f(n,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform in y is performed
! f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxi = initial x index used
! nxp = number of x indices used
! nxhd = first dimension of f >= nx/2
! nyd = second dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1,1)) = real part of mode nx/2,0 and
! imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nryb
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 70
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 50 i = nxi, nxt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
      endif
   10 continue
! then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 60 k = 2, nyh
         t1 = f(1,ny2-k)
         f(1,ny2-k) = 0.5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = 0.5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
   60    continue
      endif
      return
! forward fourier transform
   70 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 80 k = 2, nyh
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
   80    continue
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 130 i = nxi, nxt
! bit-reverse array elements in y
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
      endif
   90 continue
! then transform in y
      do 120 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the x part of 2 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! y, using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, two inverse fourier transforms in x are performed
! f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, two forward fourier transforms in x are performed
! f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:2,j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1:2,1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1:2,1,1)) = real part of mode nx/2,0 and
! imag(f(1:2,1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nrxb
      real at1, ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
! inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,ani,t1,t2,t3)
      do 90 i = nyi, nyt
! swap complex components
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,i))
      f(1,j,i) = cmplx(real(f(1,j,i)),real(f(2,j,i)))
      f(2,j,i) = cmplx(at1,aimag(f(2,j,i)))
   10 continue
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
      endif
   20 continue
! then transform in x
      do 50 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = ani*(t1 + t2)
      f(jj,nxh2-j,i) = ani*conjg(t1 - t2)
   60 continue
   70 continue
      ani = 2.0*ani
      do 80 jj = 1, 2
      f(jj,nxhh+1,i) = ani*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = ani*cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),         &
     &                      real(f(jj,1,i)) - aimag(f(jj,1,i)))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,t1,t2,t3)
      do 190 i = nyi, nyt
! scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = t1 + t2
      f(jj,nxh2-j,i) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 2
      f(jj,nxhh+1,i) = 2.0*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),             &
     &                  real(f(jj,1,i)) - aimag(f(jj,1,i)))
  130 continue
! bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
      endif
  140 continue
! then transform in x
      do 170 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  150 continue
  160 continue
  170 continue
! swap complex components
      do 180 j = 1, nxh
      at1 = aimag(f(1,j,i))
      f(1,j,i) = cmplx(real(f(1,j,i)),real(f(2,j,i)))
      f(2,j,i) = cmplx(at1,aimag(f(2,j,i)))
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the y part of 2 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! x, using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, two inverse fourier transforms in y are performed
! f(1:2,n,m) = sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, two forward fourier transforms in y are performed
! f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxi = initial x index used
! nxp = number of x indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:2,j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1:2,1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1:2,1,1)) = real part of mode nx/2,0 and
! imag(f(1:2,1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nryb
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3)
      do 50 i = nxi, nxt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
      endif
   10 continue
! then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 70 k = 2, nyh
         do 60 jj = 1, 2
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               &
     &                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    &
     &                   aimag(f(jj,1,k) - t1))
   60    continue
   70    continue
      endif
      return
! forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 100 k = 2, nyh
         do 90 jj = 1, 2
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3)
      do 150 i = nxi, nxt
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
      endif
  110 continue
! then transform in y
      do 140 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the x part of 3 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! y, using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, three inverse fourier transforms are performed
! f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, two forward fourier transforms are performed
! f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:3,j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1:3,1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1:3,1,1)) = real part of mode nx/2,0 and
! imag(f(1:3,1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nrxb
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
! inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,ani,t1,t2,t3
!$OMP& ,t4)
      do 90 i = nyi, nyt
! swap complex components
      do 10 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   10 continue
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         t3 = f(3,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(3,j1,i) = f(3,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
         f(3,j,i) = t3
      endif
   20 continue
! then transform in x
      do 50 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = ani*(t1 + t2)
      f(jj,nxh2-j,i) = ani*conjg(t1 - t2)
   60 continue
   70 continue
      ani = 2.0*ani
      do 80 jj = 1, 3
      f(jj,nxhh+1,i) = ani*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = ani*cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),         &
     &                      real(f(jj,1,i)) - aimag(f(jj,1,i)))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,t1,t2,t3,t4)
      do 190 i = nyi, nyt
! scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = t1 + t2
      f(jj,nxh2-j,i) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 3
      f(jj,nxhh+1,i) = 2.0*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),             &
     &                  real(f(jj,1,i)) - aimag(f(jj,1,i)))
  130 continue
! bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         t3 = f(3,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(3,j1,i) = f(3,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
         f(3,j,i) = t3
      endif
  140 continue
! then transform in x
      do 170 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  150 continue
  160 continue
  170 continue
! swap complex components
      do 180 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, &
     &nxhyd,nxyhd)
! this subroutine performs the y part of 3 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! x, using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, three inverse fourier transforms are performed
! f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, two forward fourier transforms are performed
! f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxi = initial x index used
! nxp = number of x indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:3,j,k) = real, imaginary part of mode j-1,k-1, where
! 1 <= j <= nx/2 and 1 <= k <= ny, except for
! f(1:3,1,k) = real, imaginary part of mode nx/2,k-1, where
! ny/2+2 <= k <= ny, and
! imag(f(1:3,1,1)) = real part of mode nx/2,0 and
! imag(f(1:3,1,ny/2+1) ) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nryb
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3,t4)
      do 50 i = nxi, nxt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         t3 = f(3,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(3,i,k1) = f(3,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
         f(3,i,k) = t3
      endif
   10 continue
! then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 70 k = 2, nyh
         do 60 jj = 1, 3
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               &
     &                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    &
     &                         aimag(f(jj,1,k) - t1))
   60    continue
   70    continue
      endif
      return
! forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 100 k = 2, nyh
         do 90 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3,t4)
      do 150 i = nxi, nxt
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         t3 = f(3,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(3,i,k1) = f(3,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
         f(3,i,k) = t3
      endif
  110 continue
! first transform in y
      do 140 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
