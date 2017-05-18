!-----------------------------------------------------------------------
! Fortran Library for performing Fast Fourier Transforms
! 2D MPI/OpenMP PIC Code:
! WPFFT2RINIT calculates tables needed by 2d FFTs
! WPPFFT2RM wrapper function for scalar 2d real/complex FFT
! WPPFFT2RM2 wrapper function for 2 component vector 2d real/complex FFT
! WPPFFT2RM3 wrapper function for 3 component vector 2d real/complex FFT
! WPPFFT2RMN wrapper function for n component vector 2d real/complex FFT
! PPFFT2RMXX performs x part of scalar 2d real/complex FFT
! PPFFT2RMXY performs y part of scalar 2d real/complex FFT
! PPFFT2RM2XX performs x part of 2 component vector 2d real/complex FFT
! PPFFT2RM2XY performs y part of 2 component vector 2d real/complex FFT
! PPFFT2RM3XX performs x part of 3 component vector 2d real/complex FFT
! PPFFT2RM3XY performs y part of 3 component vector 2d real/complex FFT
! PPFFT2RMNXX performs x part of n component vector 2d real/complex FFT
! PPFFT2RMNXY performs y part of n component vector 2d real/complex FFT
! MPPSWAPC2N swaps components for multiple ffts
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: march 23, 2017
!-----------------------------------------------------------------------
      subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
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
      subroutine WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,   &
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
! wrapper function for parallel real to complex fft
! parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd), g(nyv,kxp)
      dimension bs(kxp,kyp), br(kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   &
     &nxvh,kypd,nxhyd,nxyhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,  &
     &kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv&
     &,kxp,nxhyd,nxyhd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,   &
     &kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp&
     &,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y fft
         call PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv&
     &,kxp,nxhyd,nxyhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd, &
     &kxp)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   &
     &nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT2RM2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,  &
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
! wrapper function for 2 component parallel real to complex fft
! parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(2,nxvh,kypd), g(2,nyv,kxp)
      dimension bs(2,kxp,kyp), br(2,kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT2RM2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp&
     &,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT2RM2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,nxhyd,nxyhd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,&
     &kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,&
     &kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y fft
         call PPFFT2RM2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,nxhyd,nxyhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,   &
     &kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT2RM2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT2RM3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,  &
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
! wrapper function for 3 component parallel real to complex fft
! parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd), g(3,nyv,kxp)
      dimension bs(3,kxp,kyp), br(3,kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp&
     &,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,nxhyd,nxyhd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,&
     &kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,&
     &kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y fft
         call PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,nxhyd,nxyhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,   &
     &kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT2RMN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,ndim,nxhyd,nxyhd)
! wrapper function for n component parallel real to complex fft
! parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, ndim, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, ss, sct
      dimension f(ndim,nxvh,kypd), g(ndim,nyv,kxp)
      dimension bs(ndim,kxp,kyp), br(ndim,kxp,kyp)
      dimension ss(ndim*nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT2RMNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,ndim,nxhyd,nxyhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,ndim,nxvh,nyv,&
     &kxp,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT2RMNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,ndim,nxhyd,nxyhd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,ndim,nyv,  &
     &nxvh,kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,ndim,nxvh, &
     &nyv,kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y fft
         call PPFFT2RMNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  &
     &nyv,kxp,ndim,nxhyd,nxyhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,ndim,nyv,nxvh,&
     &kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT2RMNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,ndim,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyhd)
! this subroutine performs the x part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of y,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
! if isign = 1, a forward fourier transform is performed
! f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kypd = second dimension of f
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! f(j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
! 1 <= j <= nx/2 and 1 <= kk <= ny, except for
! f(1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
! imaginary part of f(1,1) = real part of mode nx/2,0 on mode kstrt=1
! imaginary part of f(1,1) = real part of mode nx/2,ny/2
! on mode kstrt=(ny/2)/kyp
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 70
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 60 i = kypi, kypt
! bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t
      endif
   10 continue
! then transform in x
      do 40 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i))
      s = f(j,i) + t
      t = (f(j,i) - t)*t1
      f(j,i) = ani*(s + t)
      f(nxh2-j,i) = ani*conjg(s - t)
   50 continue
      f(1,i) = 2.0*ani*cmplx(real(f(1,i)) + aimag(f(1,i)),              &
     &                       real(f(1,i)) - aimag(f(1,i)))
      if (nxhh.gt.0) f(nxhh+1,i) = 2.0*ani*conjg(f(nxhh+1,i))
   60 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
   70 nrxb = nxhy/nxh
      nrx = nxhy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 130 i = kypi, kypt
! scramble coefficients
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i))
      s = f(j,i) + t
      t = (f(j,i) - t)*t1
      f(j,i) = s + t
      f(nxh2-j,i) = conjg(s - t)
   80 continue
      f(1,i) = cmplx(real(f(1,i)) + aimag(f(1,i)),                      &
     &               real(f(1,i)) - aimag(f(1,i)))
      if (nxhh.gt.0) f(nxhh+1,i) = 2.0*conjg(f(nxhh+1,i))
! bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t
      endif
   90 continue
! then transform in x
      do 120 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,&
     &nyv,kxp,nxhyd,nxyhd)
! this subroutine performs the y part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of x,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! g(m,n) = sum(g(k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(k,j) = sum(g(m,n)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! kxp = number of x indices per block
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g
! kxp = number of data values per block in x
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! g(k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
! 1 <= jj <= nx/2 and 1 <= k <= ny, except for
! g(k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! imaginary part of g(1,1) = real part of mode nx/2,0 and
! imaginary part of g(ny/2+1,1) = real part of mode nx/2,ny/2
! on node kstrt=1
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
      integer nxhyd, nxyhd, mixup
      complex g, sct
      dimension g(nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      complex s, t
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 70
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 50 i = kxpi, kxpt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i)
         g(k1,i) = g(k,i)
         g(k,i) = t
      endif
   10 continue
! then transform in y
      do 40 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 60 k = 2, nyh
         s = g(ny2-k,1)
         g(ny2-k,1) = 0.5*cmplx(aimag(g(k,1) + s),real(g(k,1) - s))
         g(k,1) = 0.5*cmplx(real(g(k,1) + s),aimag(g(k,1) - s))
   60    continue
      endif
      return
! forward fourier transform
   70 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 80 k = 2, nyh
         s = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - s)
         g(k,1) = g(k,1) + s
   80    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 130 i = kxpi, kxpt
! bit-reverse array elements in y
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i)
         g(k1,i) = g(k,i)
         g(k,i) = t
      endif
   90 continue
! then transform in y
      do 120 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RM2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyhd)
! this subroutine performs the x part of 2 two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of y,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
! if isign = 1, a forward fourier transform is performed
! f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kypd = second dimension of f
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! f(1:2,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
! 1 <= j <= nx/2 and 1 <= kk <= ny, except for
! f(1:2,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
! imaginary part of f(1:2,1,1) = real part of mode nx/2,0
! on mode kstrt=1
! imaginary part of f(1:2,1,1) = real part of mode nx/2,ny/2
! on mode kstrt=(ny/2)/kyp
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(2,nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      real ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 100
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,s,t,t1,t2)
      do 90 i = kypi, kypt
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
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 k = 1, 2
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = ani*(s + t)
      f(k,nxh2-j,i) = ani*conjg(s - t)
   60 continue
   70 continue
      do 80 k = 1, 2
      f(k,1,i) = 2.0*ani*cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),        &
     &                         real(f(k,1,i)) - aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*ani*conjg(f(k,nxhh+1,i))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,s,t,t1,t2)
      do 190 i = kypi, kypt
! scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = 1, 2
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = s + t
      f(k,nxh2-j,i) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = 1, 2
      f(k,1,i) = cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),                &
     &                 real(f(k,1,i)) - aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*conjg(f(k,nxhh+1,i))
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
      do 170 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
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
      subroutine PPFFT2RM2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxp,nxhyd,nxyhd)
! this subroutine performs the y part of 2 two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of x,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! g(1:2,m,n) = sum(g(1:2,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(1:2,k,j) = sum(g(1:2,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g
! kxp = number of data values per block in x
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! g(1:2,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
! 1 <= jj <= nx/2 and 1 <= k <= ny, except for
! g(1:2,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! imaginary part of g(1:2,1,1) = real part of mode nx/2,0 and
! imaginary part of g(1:2,ny/2+1,1) = real part of mode nx/2,ny/2
! on node kstrt=1
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
      integer nxhyd, nxyhd, mixup
      complex g, sct
      dimension g(2,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      complex s, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 80
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2)
      do 50 i = kxpi, kxpt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
      endif
   10 continue
! then transform in y
      do 40 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 70 k = 2, nyh
         do 60 j = 1, 2
         s = g(j,ny2-k,1)
         g(j,ny2-k,1) = 0.5*cmplx(aimag(g(j,k,1) + s),                  &
     &                            real(g(j,k,1) - s))
         g(j,k,1) = 0.5*cmplx(real(g(j,k,1) + s),aimag(g(j,k,1) - s))
   60    continue
   70    continue
      endif
      return
! forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 100 k = 2, nyh
         do 90 j = 1, 2
         s = cmplx(aimag(g(j,ny2-k,1)),real(g(j,ny2-k,1)))
         g(j,ny2-k,1) = conjg(g(j,k,1) - s)
         g(j,k,1) = g(j,k,1) + s
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2)
      do 150 i = kxpi, kxpt
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
      endif
  110 continue
! then transform in y
      do 140 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyhd)
! this subroutine performs the x part of 3 two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of y,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
! if isign = 1, a forward fourier transform is performed
! f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kypd = second dimension of f
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! f(1:3,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
! 1 <= j <= nx/2 and 1 <= kk <= ny, except for
! f(1:3,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
! imaginary part of f(1:3,1,1) = real part of mode nx/2,0
! on mode kstrt=1
! imaginary part of f(1:3,1,1) = real part of mode nx/2,ny/2
! on mode kstrt=(ny/2)/kyp
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, nxhyd, nxyhd
      complex f, sct
      dimension f(3,nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      real ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 100
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,s,t,t1,t2,t3)
      do 90 i = kypi, kypt
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
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 k = 1, 3
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = ani*(s + t)
      f(k,nxh2-j,i) = ani*conjg(s - t)
   60 continue
   70 continue
      do 80 k = 1, 3
      f(k,1,i) = 2.0*ani*cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),        &
     &                         real(f(k,1,i)) - aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*ani*conjg(f(k,nxhh+1,i))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,s,t,t1,t2,t3)
      do 190 i = kypi, kypt
! scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = 1, 3
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = s + t
      f(k,nxh2-j,i) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = 1, 3
      f(k,1,i) = cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),                &
     &                 real(f(k,1,i)) -aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*conjg(f(k,nxhh+1,i))
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
      do 170 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
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
      subroutine PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxp,nxhyd,nxyhd)
! this subroutine performs the y part of 3 two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of x,
! using complex arithmetic, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
! where N = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! g(1:3,m,n) = sum(g(1:3,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(1:3,k,j) = sum(g(1:3,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g
! kxp = number of data values per block in x
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! g(1:3,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
! 1 <= jj <= nx/2 and 1 <= k <= ny, except for
! g(1:3,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! imaginary part of g(1:3,1,1) = real part of mode nx/2,0 and
! imaginary part of g(1:3,ny/2+1,1) = real part of mode nx/2,ny/2
! on node kstrt=1
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, nxhyd, nxyhd
      complex g, sct
      dimension g(3,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      complex s, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 80
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 50 i = kxpi, kxpt
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         t3 = g(3,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(3,k1,i) = g(3,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
         g(3,k,i) = t3
      endif
   10 continue
! then transform in y
      do 40 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 70 k = 2, nyh
         do 60 j = 1, 3
         s = g(j,ny2-k,1)
         g(j,ny2-k,1) = 0.5*cmplx(aimag(g(j,k,1) + s),                  &
     &                            real(g(j,k,1) - s))
         g(j,k,1) = 0.5*cmplx(real(g(j,k,1) + s),aimag(g(j,k,1) - s))
   60    continue
   70    continue
      endif
      return
! forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 100 k = 2, nyh
         do 90 j = 1, 3
         s = cmplx(aimag(g(j,ny2-k,1)),real(g(j,ny2-k,1)))
         g(j,ny2-k,1) = conjg(g(j,k,1) - s)
         g(j,k,1) = g(j,k,1) + s
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 150 i = kxpi, kxpt
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         t3 = g(3,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(3,k1,i) = g(3,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
         g(3,k,i) = t3
      endif
  110 continue
! then transform in y
      do 140 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RMNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,ndim,nxhyd,nxyhd)
! this subroutine performs the x part of N two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of y,
! using complex arithmetic, where N = ndim, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: M*(5*log2(M) + 10)/nvp
! for isign = 1,  approximate flop count: M*(5*log2(M) + 8)/nvp
! where M = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(1:N,n,m) = (1/nx*ny)*sum(f(1:N,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
! if isign = 1, a forward fourier transform is performed
! f(1:N,j,k) = sum(f(1:N,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = second dimension of f
! kypd = third dimension of f
! ndim = leading dimension of arrays f and g
! ss = scratch array
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! f(1:N,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
! 1 <= j <= nx/2 and 1 <= kk <= ny, except for
! f(1:N,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
! imaginary part of f(1:N,1,1) = real part of mode nx/2,0
! on mode kstrt=1
! imaginary part of f(1:N,1,1) = real part of mode nx/2,ny/2
! on mode kstrt=(ny/2)/kyp
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, ndim, nxhyd, nxyhd
      complex f, ss, sct
      dimension f(ndim,nxvh,kypd), ss(ndim*nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrxb
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 110
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
! swap complex components
      call MPPSWAPC2N(f,ss,isign,nxh,kypi,kypt,nxvh,kypd,ndim)
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t,t1)
      do 100 i = kypi, kypt
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 10 jj = 1, ndim
         t1 = f(jj,j1,i)
         f(jj,j1,i) = f(jj,j,i)
         f(jj,j,i) = t1
   10    continue
      endif
   20 continue
! then transform in x
      do 60 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 jj = 1, ndim
      t1 = s*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t1
      f(jj,j1,i) = f(jj,j1,i) + t1
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = 1, ndim
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = ani*(s + t)
      f(k,nxh2-j,i) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 k = 1, ndim
      f(k,1,i) = 2.0*ani*cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),        &
     &                         real(f(k,1,i)) - aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*ani*conjg(f(k,nxhh+1,i))
   90 continue
  100 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  110 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t,t1)
      do 210 i = kypi, kypt
! scramble coefficients
      kmr = nxy/nx
      do 130 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 120 k = 1, ndim
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = s + t
      f(k,nxh2-j,i) = conjg(s - t)
  120 continue
  130 continue
      do 140 k = 1, ndim
      f(k,1,i) = cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),                &
     &                 real(f(k,1,i)) -aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*conjg(f(k,nxhh+1,i))
  140 continue
! bit-reverse array elements in x
      do 160 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 150 jj = 1, ndim
         t1 = f(jj,j1,i)
         f(jj,j1,i) = f(jj,j,i)
         f(jj,j,i) = t1
  150    continue
      endif
  160 continue
! then transform in x
      do 200 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 170 jj = 1, ndim
      t1 = s*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t1
      f(jj,j1,i) = f(jj,j1,i) + t1
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
!$OMP END PARALLEL DO
! swap complex components
      call MPPSWAPC2N(f,ss,isign,nxh,kypi,kypt,nxvh,kypd,ndim)
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT2RMNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxp,ndim,nxhyd,nxyhd)
! this subroutine performs the y part of N two dimensional real to
! complex fast fourier transforms and their inverses, for a subset of x,
! using complex arithmetic, where N = ndim, with OpenMP,
! for data which is distributed in blocks
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: M*(5*log2(M) + 10)/nvp
! for isign = 1,  approximate flop count: M*(5*log2(M) + 8)/nvp
! where M = (nx/2)*ny, and nvp = number of procs
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! g(1:N,m,n) = sum(g(1:N,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(1:N,k,j) = sum(g(1:N,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = second dimension of g
! kxp = number of data values per block in x
! ndim = leading dimension of arrays f and g
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! the real data is stored in a complex array of length nx/2, ny
! with the odd/even x points stored in the real/imaginary parts.
! in complex notation, fourier coefficients are stored as follows:
! g(1:N,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
! 1 <= jj <= nx/2 and 1 <= k <= ny, except for
! g(1:N,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! imaginary part of g(1:N,1,1) = real part of mode nx/2,0 and
! imaginary part of g(1:N,ny/2+1,1) = real part of mode nx/2,ny/2
! on node kstrt=1
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, ndim, nxhyd, nxyhd
      complex g, sct
      dimension g(ndim,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nryb
      complex s, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 100
! inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t1)
      do 70 i = kxpi, kxpt
! bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 10 jj = 1, ndim
         t1 = g(jj,k1,i)
         g(jj,k1,i) = g(jj,k,i)
         g(jj,k,i) = t1
   10    continue
      endif
   20 continue
! then transform in y
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 jj = 1, ndim
      t1 = s*g(jj,j2,i)
      g(jj,j2,i) = g(jj,j1,i) - t1
      g(jj,j1,i) = g(jj,j1,i) + t1
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 90 k = 2, nyh
         do 80 j = 1, ndim
         s = g(j,ny2-k,1)
         g(j,ny2-k,1) = 0.5*cmplx(aimag(g(j,k,1) + s),                  &
     &                            real(g(j,k,1) - s))
         g(j,k,1) = 0.5*cmplx(real(g(j,k,1) + s),aimag(g(j,k,1) - s))
   80    continue
   90    continue
      endif
      return
! forward fourier transform
  100 nryb = nxhy/ny
      nry = nxy/ny
! scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 120 k = 2, nyh
         do 110 j = 1, ndim
         s = cmplx(aimag(g(j,ny2-k,1)),real(g(j,ny2-k,1)))
         g(j,ny2-k,1) = conjg(g(j,k,1) - s)
         g(j,k,1) = g(j,k,1) + s
  110    continue
  120    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t1)
      do 190 i = kxpi, kxpt
! bit-reverse array elements in y
      do 140 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 130 jj = 1, ndim
         t1 = g(jj,k1,i)
         g(jj,k1,i) = g(jj,k,i)
         g(jj,k,i) = t1
  130    continue
      endif
  140 continue
! then transform in y
      do 180 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 150 jj = 1, ndim
      t1 = s*g(jj,j2,i)
      g(jj,j2,i) = g(jj,j1,i) - t1
      g(jj,j1,i) = g(jj,j1,i) + t1
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSWAPC2N(f,s,isign,nxh,kypi,kypt,nxvh,kypd,ndim)
! this subroutine swaps components for multiple ffts
! f = input  array
! s = scratch array
! isign = (-1,1) = swap (real-to-complex,complex-to-real)
! nxh = complex dimension in x direction
! kypi/kypt = initial/final y index used
! nxvh = half of the second dimension of f
! kypd = third dimension of f
! ndim = leading dimension of array f
      implicit none
      integer isign, nxh, kypi, kypt, nxvh, kypd, ndim
      real f, s
      dimension f(ndim,2*nxvh,kypd), s(2*ndim*nxvh,kypd)
! local data
      integer i, j, k, ioff
! swap complex components
! real to complex
      if (isign.lt.0) then
!$OMP PARALLEL DO PRIVATE(i,j,k,ioff)
         do 60 k = kypi, kypt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1,k) = f(i,2*j-1,k)
         s(2*i+ioff,k) = f(i,2*j,k)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k) = s(i+ioff,k)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k) = s(i+ioff,k)
   40    continue
   50    continue
   60    continue
!$OMP END PARALLEL DO
      else if (isign.gt.0) then
! complex to real
!$OMP PARALLEL DO PRIVATE(i,j,k,ioff)
         do 120 k = kypi, kypt
         do 90 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 70 i = 1, ndim
         s(i+ioff,k) = f(i,2*j-1,k)
   70    continue
         ioff = ioff + ndim
         do 80 i = 1, ndim
         s(i+ioff,k) = f(i,2*j,k)
   80    continue
   90    continue
         do 110 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 100 i = 1, ndim
         f(i,2*j-1,k) = s(2*i+ioff-1,k)
         f(i,2*j,k) = s(2*i+ioff,k)
  100    continue
  110    continue
  120    continue
!$OMP END PARALLEL DO
      endif
      return
      end
