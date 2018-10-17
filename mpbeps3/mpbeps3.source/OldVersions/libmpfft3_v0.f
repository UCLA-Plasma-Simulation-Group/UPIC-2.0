!-----------------------------------------------------------------------
! Fortran Library for performing Fast Fourier Transforms
! 3D MPI/OpenMP PIC Code:
! WPFFT32RINIT calculates tables needed by 3d FFTs
! WPPFFT32RM wrapper function for scalar 3d real/complex FFT
! WPPFFT32RM3 wrapper function for 3 component vector 3d real/complex
!             FFT
! WPPFFT32RMN wrapper function for n component vector 3d real/complex
!             FFT
! PPFFT32RMXX performs x part of scalar 3d real/complex FFT
! PPFFT32RMXY performs y part of scalar 3d real/complex FFT
! PPFFT32RMXZ performs z part of scalar 3d real/complex FFT
! PPFFT32RM3XX performs x part of 3 component vector 3d real/complex FFT
! PPFFT32RM3XY performs y part of 3 component vector 3d real/complex FFT
! PPFFT32RM3XZ performs z part of 3 component vector 3d real/complex FFT
! PPFFT32RMNXX performs x part of n component vector 3d real/complex FFT
! PPFFT32RMNXY performs y part of n component vector 3d real/complex FFT
! PPFFT32RMNXZ performs z part of n component vector 3d real/complex FFT
! MPPSWAPC32N swaps components for multiple ffts
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 15, 2018
!-----------------------------------------------------------------------
      subroutine WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
! this subroutine calculates tables needed by a three dimensional
! real to complex fast fourier transform and its inverse.
! input: indx, indy, indz, nxhyzd, nxyzhd
! output: mixup, sct
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, indz, nxhyzd, nxyzhd
      integer mixup
      complex sct
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh
      integer j, k, lb, ll, jb, it
      real dnxyz, arg
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
! bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 10 k = 1, ndx1yz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
! sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/real(nxyz)
      do 30 j = 1, nxyzh
      arg = dnxyz*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT32RM(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,   &
     &kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
! wrapper function for 3d real to complex fft, with packed data
! parallelized with MPI/OpenMP
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nvpy, nvpz
      integer nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
      integer kxypd, kypd, kyzpd, kzpd, kzyp, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(nxvh,kypd,kzpd), g(nyv,kxypd,kzpd), h(nzv,kxypd,kyzpd)
      dimension bs(kxyp*kzyp,kzp), br(kxyp*kzyp,kzp)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer nxh, ny, nz, kypi, kxypi, js, ks, kxypp, kypp, kzpp, nvp
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxypp = min(kxyp,max(0,nxh-kxyp*js))
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
      nvp = nvpy*nvpz
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,   &
     &kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxvh,&
     &nyv,kxypd,kypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PPTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz&
     &,nyv,nzv,kxypd,kyzpd,kzpd)
         call PWTIMERA(1,tp,dtime)
! perform z fft
         call PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy, &
     &nvpz,nzv,nyv,kxypd,kzpd,kyzpd)
            call PPTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,  &
     &nyv,nxvh,kypd,kxypd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,  &
     &nxvh,nyv,kxypd,kypd,kzpd)
            call PPTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy, &
     &nvpz,nyv,nzv,kxypd,kyzpd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform z fft
         call PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PPTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,nvpz&
     &,nzv,nyv,kxypd,kzpd,kyzpd)
         call PWTIMERA(1,tp,dtime)
! perform y fft
         call PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,nyv, &
     &nxvh,kypd,kxypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,   &
     &kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT32RM3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,  &
     &kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
! wrapper function for 3 3d real to complex ffts, with packed data
! c parallelized with MPI/OpenMP
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nvpy, nvpz
      integer nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
      integer kxypd, kypd, kyzpd, kzpd, kzyp, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(3,nxvh,kypd,kzpd), g(3,nyv,kxypd,kzpd)
      dimension h(3,nzv,kxypd,kyzpd)
      dimension bs(3,kxyp*kzyp,kzp), br(3,kxyp*kzyp,kzp)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer nxh, ny, nz, kypi, kxypi, js, ks, kxypp, kypp, kzpp, nvp
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxypp = min(kxyp,max(0,nxh-kxyp*js))
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
      nvp = nvpy*nvpz
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,  &
     &kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,3,  &
     &nxvh,nyv,kxypd,kypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   &
     &nvpz,3,nyv,nzv,kxypd,kyzpd,kzpd)
         call PWTIMERA(1,tp,dtime)
! perform z fft
         call PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,&
     &nvpz,3,nzv,nyv,kxypd,kzpd,kyzpd)
            call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,3&
     &,nyv,nxvh,kypd,kxypd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,3&
     &,nxvh,nyv,kxypd,kypd,kzpd)
            call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,&
     &nvpz,3,nyv,nzv,kxypd,kyzpd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform z fft
         call PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,   &
     &nvpz,3,nzv,nyv,kxypd,kzpd,kyzpd)
         call PWTIMERA(1,tp,dtime)
! perform y fft
         call PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,3,  &
     &nyv,nxvh,kypd,kxypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,  &
     &kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFFT32RMN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,ttp, &
     &indx,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,    &
     &kxypd,kypd,kyzpd,kzpd,kzyp,ndim,nxhyzd,nxyzhd)
! wrapper function for multiple 3d real to complex ffts,
! with packed data, parallelized with MPI/OpenMP
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nvpy, nvpz
      integer nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
      integer kxypd, kypd, kyzpd, kzpd, kzyp, ndim, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd), g(ndim,nyv,kxypd,kzpd)
      dimension h(ndim,nzv,kxypd,kyzpd)
      dimension bs(ndim,kxyp*kzyp,kzp), br(ndim,kxyp*kzyp,kzp)
      dimension ss(ndim,nxvh,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer nxh, ny, nz, kypi, kxypi, js, ks, kxypp, kypp, kzpp, nvp
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxypp = min(kxyp,max(0,nxh-kxyp*js))
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
      nvp = nvpy*nvpz
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call PPFFT32RMNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,nvp&
     &,kypi,kypp,nxvh,kzpp,kypd,kzpd,ndim,nxhyzd,nxyzhd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,ndim&
     &,nxvh,nyv,kxypd,kypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform y fft
         call PPFFT32RMNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,ndim,nxhyzd,nxyzhd)
! transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   &
     &nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
         call PWTIMERA(1,tp,dtime)
! perform z fft
         call PPFFT32RMNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,ndim,nxhyzd,nxyzhd)
! transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,&
     &nvpz,ndim,nzv,nyv,kxypd,kzpd,kyzpd)
            call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy, &
     &ndim,nyv,nxvh,kypd,kxypd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy, &
     &ndim,nxvh,nyv,kxypd,kypd,kzpd)
            call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,&
     &nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform z fft
         call PPFFT32RMNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,ndim,nxhyzd,nxyzhd)
! transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,   &
     &nvpz,ndim,nzv,nyv,kxypd,kzpd,kyzpd)
         call PWTIMERA(1,tp,dtime)
! perform y fft
         call PPFFT32RMNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, &
     &nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,ndim,nxhyzd,nxyzhd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,ndim&
     &,nyv,nxvh,kypd,kxypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
! perform x fft
         call PPFFT32RMNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,nvp&
     &,kypi,kypp,nxvh,kzpp,kypd,kzpd,ndim,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,&
     &kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
! this subroutine performs the x part of a three dimensional real to
! complex fast fourier transform and its inverse for a subset of y and z
! using complex arithmetic, for data which is distributed in blocks,
! with 2D spatial decomposition and OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! f(n,k,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(n,k,i) = sum(f(j,k,i)*exp(sqrt(-1)*2pi*n*j/nx))
! kstrt = starting data block number
! nvp = number of real or virtual processors
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kzpp = number of z indices used
! kypd = second dimension of f
! kzpd = third dimension of f
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! final fourier coefficients are stored as follows:
! h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
! and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1,1,1)) = real part of mode nx/2,0,0
! imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvp, kypi, kypp, nxvh
      integer kzpp, kypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex f, sct
      dimension f(nxvh,kypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      integer nrx, kypt, nn
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kypt = kypi + kypp - 1
      if (kstrt.gt.nvp) return
      if (isign.gt.0) go to 70
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 60 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t
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
      s = sct(1+kmr*(j-1))
      t = s*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t
      f(j1,i,n) = f(j1,i,n) + t
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      kmr = nxyz/nx
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i,n))
      s = f(j,i,n) + t
      t = (f(j,i,n) - t)*t1
      f(j,i,n) = ani*(s + t)
      f(nxh2-j,i,n) = ani*conjg(s - t)
   50 continue
      f(1,i,n) = 2.0*ani*cmplx(real(f(1,i,n)) + aimag(f(1,i,n)),        &
     &                         real(f(1,i,n)) - aimag(f(1,i,n)))
      if (nxhh.gt.0) f(nxhh+1,i,n) = 2.0*ani*conjg(f(nxhh+1,i,n))
   60 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
   70 nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 130 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! scramble coefficients
      kmr = nxyz/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i,n))
      s = f(j,i,n) + t
      t = (f(j,i,n) - t)*t1
      f(j,i,n) = s + t
      f(nxh2-j,i,n) = conjg(s - t)
   80 continue
      f(1,i,n) = cmplx(real(f(1,i,n)) + aimag(f(1,i,n)),                &
     &                 real(f(1,i,n)) - aimag(f(1,i,n)))
      if (nxhh.gt.0) f(nxhh+1,i,n) = 2.0*conjg(f(nxhh+1,i,n))
! bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t
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
      s = conjg(sct(1+kmr*(j-1)))
      t = s*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t
      f(j1,i,n) = f(j1,i,n) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy&
     &,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! this subroutine performs the y part of a three dimensional real to
! complex fast fourier transform and its inverse for a subset of x and z
! using complex arithmetic, for data which is distributed in blocks,
! with 2D spatial decomposition and OpenMP
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! g(m,j,i) = sum(g(k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(m,j,i) = sum(g(k,j,i)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nyv = first dimension of g
! kzpp = number of z indices used
! kxypd = second dimension of g
! kzpd = third dimension of g
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! final fourier coefficients are stored as follows:
! h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
! and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1,1,1)) = real part of mode nx/2,0,0
! imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nyv, kzpp, kxypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex g, sct
      dimension g(nyv,kxypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      integer js, ks, nry, kxypt, nn
      complex s, t
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 80
! inverse fourier transform
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 50 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i,n)
         g(k1,i,n) = g(k,i,n)
         g(k,i,n) = t
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
      s = sct(1+kmr*(j-1))
      t = s*g(j2,i,n)
      g(j2,i,n) = g(j1,i,n) - t
      g(j1,i,n) = g(j1,i,n) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,s)
         do 70 n = 1, kzpp
         do 60 k = 2, nyh
         s = g(ny2-k,1,n)
         g(ny2-k,1,n) = 0.5*cmplx(aimag(g(k,1,n) + s),                  &
     &                            real(g(k,1,n) - s))
         g(k,1,n) = 0.5*cmplx(real(g(k,1,n) + s),aimag(g(k,1,n) - s))
   60    continue
   70    continue
!$OMP END PARALLEL DO
      endif
      return
! forward fourier transform
   80 nryb = nxhyz/ny
      nry = nxyz/ny
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,s)
         do 100 n = 1, kzpp
         do 90 k = 2, nyh
         s = cmplx(aimag(g(ny2-k,1,n)),real(g(ny2-k,1,n)))
         g(ny2-k,1,n) = conjg(g(k,1,n) - s)
         g(k,1,n) = g(k,1,n) + s
   90    continue
  100    continue
!$OMP END PARALLEL DO
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 150 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i,n)
         g(k1,i,n) = g(k,i,n)
         g(k,i,n) = t
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
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i,n)
      g(j2,i,n) = g(j1,i,n) - t
      g(j1,i,n) = g(j1,i,n) + t
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy&
     &,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! this subroutine performs the z part of a three dimensional real to
! complex fast fourier transform and its inverse for a subset of x and y
! using complex arithmetic, for data which is distributed in blocks,
! with 2D spatial decomposition and OpenMP
! for isign = (-1,1), input: all, output: h
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! h(l,n,m) = sum(h(i,j,k)*exp(-sqrt(-1)*2pi*l*i/nz))
! if isign = 1, a forward fourier transform is performed
! h(l,n,m) = sum(h(i,j,k)*exp(sqrt(-1)*2pi*ll*ii/nz))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nzv = first dimension of h
! kyzpp = number of y indices used
! kxypd = second dimension of h
! kyzpd = third dimension of h
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! final fourier coefficients are stored as follows:
! h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
! and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1,1,1)) = real part of mode nx/2,0,0
! imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nzv, kyzp, kxypd, kyzpd, nxhyzd, nxyzhd
      integer mixup
      complex h, sct
      dimension h(nzv,kxypd,kyzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nrzb
      integer l1, js, ks, nrz, kxypt, kyzpp, kyzb, nn
      complex s, t
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      nz = 2**indz
      nzh = max(1,nz/2)
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyzpp = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 80
! inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t)
      do 50 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 10 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t = h(l1,i,n)
         h(l1,i,n) = h(l,i,n)
         h(l,i,n) = t
      endif
   10 continue
! finally transform in z
      do 40 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*h(j2,i,n)
      h(j2,i,n) = h(j1,i,n) - t
      h(j1,i,n) = h(j1,i,n) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 60 n = 2, nzh
            s = h(nz2-n,1,1)
            h(nz2-n,1,1) = 0.5*cmplx(aimag(h(n,1,1) + s),               &
     &                               real(h(n,1,1) - s))
            h(n,1,1) = 0.5*cmplx(real(h(n,1,1) + s),aimag(h(n,1,1) - s))
   60       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 70 n = 2, nzh
            s = h(nz2-n,1,k1)
            h(nz2-n,1,k1) = 0.5*cmplx(aimag(h(n,1,k1) + s),             &
     &                                real(h(n,1,k1) - s))
            h(n,1,k1) = 0.5*cmplx(real(h(n,1,k1) + s),                  &
     &                            aimag(h(n,1,k1) - s))
   70       continue
        endif
      endif
      return
! forward fourier transform
   80 nrzb = nxhyz/nz
      nrz = nxyz/nz
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 90 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,1)),real(h(nz2-n,1,1)))
            h(nz2-n,1,1) = conjg(h(n,1,1) - s)
            h(n,1,1) = h(n,1,1) + s
   90       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 100 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,k1)),real(h(nz2-n,1,k1)))
            h(nz2-n,1,k1) = conjg(h(n,1,k1) - s)
            h(n,1,k1) = h(n,1,k1) + s
  100       continue
         endif
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t)
      do 150 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 110 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t = h(l1,i,n)
         h(l1,i,n) = h(l,i,n)
         h(l,i,n) = t
      endif
  110 continue
! first transform in z
      do 140 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*h(j2,i,n)
      h(j2,i,n) = h(j1,i,n) - t
      h(j1,i,n) = h(j1,i,n) + t
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp&
     &,kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
! this subroutine performs the x part of 3 three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! y and z, using complex arithmetic,
! for data which is distributed in blocks, with 2D spatial decomposition
! with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! f(1:3,n,k,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*
!                                 exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(1:3,n,k,i) = sum(f(1:3,j,k,i)*exp(sqrt(-1)*2pi*n*j/nx))
! kstrt = starting data block number
! nvp = number of real or virtual processors
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kzpp = number of z indices used
! kypd = second dimension of f
! kzpd = third dimension of f
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! the real data is stored in a complex array of length nx/2, ny, nz
! with the odd/even x points stored in the real/imaginary parts.
! final fourier coefficients are stored as follows:
! h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvp, kypi, kypp, nxvh
      integer kzpp, kypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex f, sct
      dimension f(3,nxvh,kypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrxb
      integer nrx, kypt, nn
      real ani, at1, at2
      complex s, t, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kypt = kypi + kypp - 1
      if (kstrt.gt.nvp) return
      if (isign.gt.0) go to 100
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,s,t,t1,&
!$OMP& t2,t3)
      do 90 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! swap complex components
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = aimag(f(2,j,i,n))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
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
      s = sct(1+kmr*(j-1))
      t1 = s*f(1,j2,i,n)
      t2 = s*f(2,j2,i,n)
      t3 = s*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t1
      f(2,j2,i,n) = f(2,j1,i,n) - t2
      f(3,j2,i,n) = f(3,j1,i,n) - t3
      f(1,j1,i,n) = f(1,j1,i,n) + t1
      f(2,j1,i,n) = f(2,j1,i,n) + t2
      f(3,j1,i,n) = f(3,j1,i,n) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
      kmr = nxyz/nx
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 3
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = ani*(s + t)
      f(jj,nxh2-j,i,n) = ani*conjg(s - t)
   60 continue
   70 continue
      do 80 jj = 1, 3
      f(jj,1,i,n) =                                                     &
     &             2.0*ani*cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),&
     &                           real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*ani*conjg(f(jj,nxhh+1,i,n))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  100 nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,s,t,t1,&
!$OMP& t2,t3)
      do 190 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! scramble coefficients
      kmr = nxyz/nx
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 3
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = s + t
      f(jj,nxh2-j,i,n) = conjg(s - t)
  110 continue
  120 continue
      do 130 jj = 1, 3
      f(jj,1,i,n) = cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),       &
     &                    real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*conjg(f(jj,nxhh+1,i,n))
  130 continue
! bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
      endif
  140 continue
! finally transform in x
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
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*f(1,j2,i,n)
      t2 = s*f(2,j2,i,n)
      t3 = s*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t1
      f(2,j2,i,n) = f(2,j1,i,n) - t2
      f(3,j2,i,n) = f(3,j1,i,n) - t3
      f(1,j1,i,n) = f(1,j1,i,n) + t1
      f(2,j1,i,n) = f(2,j1,i,n) + t2
      f(3,j1,i,n) = f(3,j1,i,n) + t3
  150 continue
  160 continue
  170 continue
! swap complex components
      do 180 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,aimag(f(1,j,i,n)))
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
! this subroutine performs the y part of 3 three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! x and z, using complex arithmetic,
! for data which is distributed in blocks, with 2D spatial decomposition
! with OpenMP
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! g(1:3,m,j,i) = sum(g(1:3,k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(1:3,m,j,i) = sum(g(1:3,k,j,i)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nyv = first dimension of g
! kzpp = number of z indices used
! kxypd = second dimension of g
! kzpd = third dimension of g
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! final fourier coefficients are stored as follows:
! h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nyv, kzpp, kxypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex g, sct
      dimension g(3,nyv,kxypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nryb
      integer js, ks, nry, kxypt, nn
      complex s, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 90
! inverse fourier transform
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 50 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i,n)
         t2 = g(2,k1,i,n)
         t3 = g(3,k1,i,n)
         g(1,k1,i,n) = g(1,k,i,n)
         g(2,k1,i,n) = g(2,k,i,n)
         g(3,k1,i,n) = g(3,k,i,n)
         g(1,k,i,n) = t1
         g(2,k,i,n) = t2
         g(3,k,i,n) = t3
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
      s = sct(1+kmr*(j-1))
      t1 = s*g(1,j2,i,n)
      t2 = s*g(2,j2,i,n)
      t3 = s*g(3,j2,i,n)
      g(1,j2,i,n) = g(1,j1,i,n) - t1
      g(2,j2,i,n) = g(2,j1,i,n) - t2
      g(3,j2,i,n) = g(3,j1,i,n) - t3
      g(1,j1,i,n) = g(1,j1,i,n) + t1
      g(2,j1,i,n) = g(2,j1,i,n) + t2
      g(3,j1,i,n) = g(3,j1,i,n) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 80 n = 1, kzpp
         do 70 k = 2, nyh
         do 60 jj = 1, 3
         s = g(jj,ny2-k,1,n)
         g(jj,ny2-k,1,n) = 0.5*cmplx(aimag(g(jj,k,1,n) + s),            &
     &                               real(g(jj,k,1,n) - s))
         g(jj,k,1,n) = 0.5*cmplx(real(g(jj,k,1,n) + s),                 &
     &                           aimag(g(jj,k,1,n) - s))
   60    continue
   70    continue
   80    continue
!$OMP END PARALLEL DO
      endif
      return
! forward fourier transform
  90  nryb = nxhyz/ny
      nry = nxyz/ny
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 120 n = 1, kzpp
         do 110 k = 2, nyh
         do 100 jj = 1, 3
         s = cmplx(aimag(g(jj,ny2-k,1,n)),real(g(jj,ny2-k,1,n)))
         g(jj,ny2-k,1,n) = conjg(g(jj,k,1,n) - s)
         g(jj,k,1,n) = g(jj,k,1,n) + s
  100    continue
  110    continue
  120    continue
!$OMP END PARALLEL DO
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 170 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i,n)
         t2 = g(2,k1,i,n)
         t3 = g(3,k1,i,n)
         g(1,k1,i,n) = g(1,k,i,n)
         g(2,k1,i,n) = g(2,k,i,n)
         g(3,k1,i,n) = g(3,k,i,n)
         g(1,k,i,n) = t1
         g(2,k,i,n) = t2
         g(3,k,i,n) = t3
      endif
  130 continue
! then transform in y
      do 160 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*g(1,j2,i,n)
      t2 = s*g(2,j2,i,n)
      t3 = s*g(3,j2,i,n)
      g(1,j2,i,n) = g(1,j1,i,n) - t1
      g(2,j2,i,n) = g(2,j1,i,n) - t2
      g(3,j2,i,n) = g(3,j1,i,n) - t3
      g(1,j1,i,n) = g(1,j1,i,n) + t1
      g(2,j1,i,n) = g(2,j1,i,n) + t2
      g(3,j1,i,n) = g(3,j1,i,n) + t3
  140 continue
  150 continue
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
! this subroutine performs the z part of 3 three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! x and y, using complex arithmetic,
! for data which is distributed in blocks, with 2D spatial decomposition
! with OpenMP
! for isign = (-1,1), input: all, output: h
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
! where N = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! h(1:3,l,n,m) = sum(h(1:3,i,j,k)*exp(-sqrt(-1)*2pi*l*i/nz))
! if isign = 1, a forward fourier transform is performed
! h(1:3,l,n,m) = sum(h(1:3,i,j,k)*exp(sqrt(-1)*2pi*ll*ii/nz))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nzv = first dimension of h
! kyzpp = number of y indices used
! kxypd = second dimension of h
! kyzpd = third dimension of h
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nzv, kyzp, kxypd, kyzpd, nxhyzd, nxyzhd
      integer mixup
      complex h, sct
      dimension h(3,nzv,kxypd,kyzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrzb
      integer l1, js, ks, nrz, kxypt, kyzpp, kyzb, nn
      complex s, t1, t2, t3
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      nz = 2**indz
      nzh = max(1,nz/2)
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyzpp = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 100
! inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1,t2,t3)
      do 50 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 10 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t1 = h(1,l1,i,n)
         t2 = h(2,l1,i,n)
         t3 = h(3,l1,i,n)
         h(1,l1,i,n) = h(1,l,i,n)
         h(2,l1,i,n) = h(2,l,i,n)
         h(3,l1,i,n) = h(3,l,i,n)
         h(1,l,i,n) = t1
         h(2,l,i,n) = t2
         h(3,l,i,n) = t3
      endif
   10 continue
! finally transform in z
      do 40 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*h(1,j2,i,n)
      t2 = s*h(2,j2,i,n)
      t3 = s*h(3,j2,i,n)
      h(1,j2,i,n) = h(1,j1,i,n) - t1
      h(2,j2,i,n) = h(2,j1,i,n) - t2
      h(3,j2,i,n) = h(3,j1,i,n) - t3
      h(1,j1,i,n) = h(1,j1,i,n) + t1
      h(2,j1,i,n) = h(2,j1,i,n) + t2
      h(3,j1,i,n) = h(3,j1,i,n) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 70 n = 2, nzh
            do 60 jj = 1, 3
            s = h(jj,nz2-n,1,1)
            h(jj,nz2-n,1,1) = 0.5*cmplx(aimag(h(jj,n,1,1) + s),         &
     &                                  real(h(jj,n,1,1) - s))
            h(jj,n,1,1) = 0.5*cmplx(real(h(jj,n,1,1) + s),              &
     &                              aimag(h(jj,n,1,1) - s))
   60       continue
   70       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 90 n = 2, nzh
            do 80 jj = 1, 3
            s = h(jj,nz2-n,1,k1)
            h(jj,nz2-n,1,k1) = 0.5*cmplx(aimag(h(jj,n,1,k1) + s),       &
     &                                   real(h(jj,n,1,k1) - s))
            h(jj,n,1,k1) = 0.5*cmplx(real(h(jj,n,1,k1) + s),            &
     &                               aimag(h(jj,n,1,k1) - s))
   80       continue
   90       continue
        endif
      endif
      return
! forward fourier transform
  100 nrzb = nxhyz/nz
      nrz = nxyz/nz
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 120 n = 2, nzh
            do 110 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,1)),real(h(jj,nz2-n,1,1)))
            h(jj,nz2-n,1,1) = conjg(h(jj,n,1,1) - s)
            h(jj,n,1,1) = h(jj,n,1,1) + s
  110       continue
  120       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 140 n = 2, nzh
            do 130 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,k1)),real(h(jj,nz2-n,1,k1)))
            h(jj,nz2-n,1,k1) = conjg(h(jj,n,1,k1) - s)
            h(jj,n,1,k1) = h(jj,n,1,k1) + s
  130       continue
  140       continue
         endif
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1,t2,t3)
      do 190 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 150 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t1 = h(1,l1,i,n)
         t2 = h(2,l1,i,n)
         t3 = h(3,l1,i,n)
         h(1,l1,i,n) = h(1,l,i,n)
         h(2,l1,i,n) = h(2,l,i,n)
         h(3,l1,i,n) = h(3,l,i,n)
         h(1,l,i,n) = t1
         h(2,l,i,n) = t2
         h(3,l,i,n) = t3
      endif
  150 continue
! first transform in z
      do 180 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*h(1,j2,i,n)
      t2 = s*h(2,j2,i,n)
      t3 = s*h(3,j2,i,n)
      h(1,j2,i,n) = h(1,j1,i,n) - t1
      h(2,j2,i,n) = h(2,j1,i,n) - t2
      h(3,j2,i,n) = h(3,j1,i,n) - t3
      h(1,j1,i,n) = h(1,j1,i,n) + t1
      h(2,j1,i,n) = h(2,j1,i,n) + t2
      h(3,j1,i,n) = h(3,j1,i,n) + t3
  160 continue
  170 continue
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end

!-----------------------------------------------------------------------
      subroutine PPFFT32RMNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,&
     &nvp,kypi,kypp,nxvh,kzpp,kypd,kzpd,ndim,nxhyzd,nxyzhd)
! this subroutine performs the x part of N three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! y and z, using complex arithmetic, where N = ndim, with OpenMP
! for data which is distributed in blocks, with 2D spatial decomposition
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)/nvp
! for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)/nvp
! where M = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! f(1:N,n,k,i) = (1/nx*ny*nz)*sum(f(1:N,j,k,i)*
!                                 exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(1:N,n,k,i) = sum(f(1:N,j,k,i)*exp(sqrt(-1)*2pi*n*j/nx))
! kstrt = starting data block number
! nvp = number of real or virtual processors
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f
! kzpp = number of z indices used
! kypd = second dimension of f
! kzpd = third dimension of f
! ss = scratch array
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! ndim = leading dimension of array f
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! the real data is stored in a complex array of length nx/2, ny, nz
! with the odd/even x points stored in the real/imaginary parts.
! final fourier coefficients are stored as follows:
! h(1:N,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:N,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:N,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:N,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:N,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:N,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:N,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:N,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvp, kypi, kypp, nxvh
      integer kzpp, kypd, kzpd, ndim, nxhyzd, nxyzhd
      integer mixup
      complex f, ss, sct
      dimension ss(ndim,nxvh,kzpd)
      dimension f(ndim,nxvh,kypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrxb
      integer nrx, kypt, nn
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kypt = kypi + kypp - 1
      if (kstrt.gt.nvp) return
      if (isign.gt.0) go to 110
! inverse fourier transform
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
! swap complex components
      call MPPSWAPC32N(f,ss,isign,nxh,kypi,kypt,nxvh,kzpp,kypd,kzpd,ndim&
     &)
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t,t1)
      do 100 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 10 jj = 1, ndim
         t1 = f(jj,j1,i,n)
         f(jj,j1,i,n) = f(jj,j,i,n)
         f(jj,j,i,n) = t1
   10    continue
      endif
   20 continue
! then transform in x
      do 60 l = 1, indx1
      ns = 2**(l - 1)
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
      t1 = s*f(jj,j2,i,n)
      f(jj,j2,i,n) = f(jj,j1,i,n) - t1
      f(jj,j1,i,n) = f(jj,j1,i,n) + t1
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble coefficients and normalize
      kmr = nxyz/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = ani*(s + t)
      f(jj,nxh2-j,i,n) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 jj = 1, ndim
      f(jj,1,i,n) =                                                     &
     &             2.0*ani*cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),&
     &                           real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*ani*conjg(f(jj,nxhh+1,i,n))
   90 continue
  100 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  110 nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,s,t,t1)
      do 210 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
! scramble coefficients
      kmr = nxyz/nx
      do 130 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 120 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = s + t
      f(jj,nxh2-j,i,n) = conjg(s - t)
  120 continue
  130 continue
      do 140 jj = 1, ndim
      f(jj,1,i,n) = cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),       &
     &                    real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*conjg(f(jj,nxhh+1,i,n))
  140 continue
! bit-reverse array elements in x
      do 160 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 150 jj = 1, ndim
         t1 = f(jj,j1,i,n)
         f(jj,j1,i,n) = f(jj,j,i,n)
         f(jj,j,i,n) = t1
  150    continue
      endif
  160 continue
! finally transform in x
      do 200 l = 1, indx1
      ns = 2**(l - 1)
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
      t1 = s*f(jj,j2,i,n)
      f(jj,j2,i,n) = f(jj,j1,i,n) - t1
      f(jj,j1,i,n) = f(jj,j1,i,n) + t1
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
!$OMP END PARALLEL DO
! swap complex components
      call MPPSWAPC32N(f,ss,isign,nxh,kypi,kypt,nxvh,kzpp,kypd,kzpd,ndim&
     &)
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RMNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,ndim,nxhyzd,nxyzhd)
! this subroutine performs the y part of N three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! x and z, using complex arithmetic, where N = ndim, with OpenMP
! for data which is distributed in blocks, with 2D spatial decomposition
! for isign = (-1,1), input: all, output: g
! for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)/nvp
! for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)/nvp
! where M = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! g(1:N,m,j,i) = sum(g(1:N,k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! g(1:N,m,j,i) = sum(g(1:N,k,j,i)*exp(sqrt(-1)*2pi*m*k/ny))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nyv = first dimension of g
! kzpp = number of z indices used
! kxypd = second dimension of g
! kzpd = third dimension of g
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! ndim = leading dimension of array g
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! final fourier coefficients are stored as follows:
! h(1:N,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:N,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:N,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:N,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:N,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:N,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:N,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:N,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nyv, kzpp, kxypd, kzpd, ndim, nxhyzd, nxyzhd
      integer mixup
      complex g, sct
      dimension g(ndim,nyv,kxypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nryb
      integer js, ks, nry, kxypt, nn
      complex s, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 110
! inverse fourier transform
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1)
      do 70 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 10 jj = 1, ndim
         t1 = g(jj,k1,i,n)
         g(jj,k1,i,n) = g(jj,k,i,n)
         g(jj,k,i,n) = t1
   10    continue
      endif
   20 continue
! then transform in y
      do 60 l = 1, indy
      ns = 2**(l - 1)
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
      t1 = s*g(jj,j2,i,n)
      g(jj,j2,i,n) = g(jj,j1,i,n) - t1
      g(jj,j1,i,n) = g(jj,j1,i,n) + t1
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 100 n = 1, kzpp
         do 90 k = 2, nyh
         do 80 jj = 1, ndim
         s = g(jj,ny2-k,1,n)
         g(jj,ny2-k,1,n) = 0.5*cmplx(aimag(g(jj,k,1,n) + s),            &
     &                               real(g(jj,k,1,n) - s))
         g(jj,k,1,n) = 0.5*cmplx(real(g(jj,k,1,n) + s),                 &
     &                           aimag(g(jj,k,1,n) - s))
   80    continue
   90    continue
  100    continue
!$OMP END PARALLEL DO
      endif
      return
! forward fourier transform
  110 nryb = nxhyz/ny
      nry = nxyz/ny
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 140 n = 1, kzpp
         do 130 k = 2, nyh
         do 120 jj = 1, ndim
         s = cmplx(aimag(g(jj,ny2-k,1,n)),real(g(jj,ny2-k,1,n)))
         g(jj,ny2-k,1,n) = conjg(g(jj,k,1,n) - s)
         g(jj,k,1,n) = g(jj,k,1,n) + s
  120    continue
  130    continue
  140    continue
!$OMP END PARALLEL DO
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1)
      do 210 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
! bit-reverse array elements in y
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 150 jj = 1, ndim
         t1 = g(jj,k1,i,n)
         g(jj,k1,i,n) = g(jj,k,i,n)
         g(jj,k,i,n) = t1
  150    continue
      endif
  160 continue
! then transform in y
      do 200 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 170 jj = 1, ndim
      t1 = s*g(jj,j2,i,n)
      g(jj,j2,i,n) = g(jj,j1,i,n) - t1
      g(jj,j1,i,n) = g(jj,j1,i,n) + t1
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFFT32RMNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,   &
     &nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,ndim,nxhyzd,nxyzhd)
! this subroutine performs the z part of N three dimensional real to
! complex fast fourier transforms and their inverses for a subset of
! x and y, using complex arithmetic, where N = ndim, with OpenMP
! for data which is distributed in blocks, with 2D spatial decomposition
! for isign = (-1,1), input: all, output: h
! for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)/nvp
! for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)/nvp
! where M = (nx/2)*ny*nz, and nvp = number of procs
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform is performed
! h(1:N,l,n,m) = sum(h(1:N,i,j,k)*exp(-sqrt(-1)*2pi*l*i/nz))
! if isign = 1, a forward fourier transform is performed
! h(1:N,l,n,m) = sum(h(1:N,i,j,k)*exp(sqrt(-1)*2pi*ll*ii/nz))
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! kxypi = initial x index used
! kxypp = number of x indices used
! nzv = first dimension of h
! kyzpp = number of y indices used
! kxypd = second dimension of h
! kyzpd = third dimension of h
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! ndim = leading dimension of array h
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = one half of maximum of (nx,ny,nz)
! h(1:N,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
! kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
! 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
! h(1:N,l,1,k) = mode nx/2,kk-1,l-1,
! where ny/2+2 <= kk <= ny, 1 <= l <= nz,
! the following are located on node js = 0 and ks = 0:
! h(1:N,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
! imag(h(1:N,1,1,1)) = real part of mode nx/2,0,0
! imag(h(1:N,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
! the following are located on node js = 0 and ks = nyh/kyzp:
! h(1:N,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
! imag(h(1:N,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
! imag(h(1:N,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
! parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nzv, kyzp, kxypd, kyzpd, ndim, nxhyzd, nxyzhd
      integer mixup
      complex h, sct
      dimension h(ndim,nzv,kxypd,kyzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrzb
      integer l1, js, ks, nrz, kxypt, kyzpp, kyzb, nn
      complex s, t1
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      nz = 2**indz
      nzh = max(1,nz/2)
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyzpp = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 120
! inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1)
      do 70 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 10 jj = 1, ndim
         t1 = h(jj,l1,i,n)
         h(jj,l1,i,n) = h(jj,l,i,n)
         h(jj,l,i,n) = t1
   10    continue
      endif
   20 continue
! finally transform in z
      do 60 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 jj = 1, ndim
      t1 = s*h(jj,j2,i,n)
      h(jj,j2,i,n) = h(jj,j1,i,n) - t1
      h(jj,j1,i,n) = h(jj,j1,i,n) + t1
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 90 n = 2, nzh
            do 80 jj = 1, ndim
            s = h(jj,nz2-n,1,1)
            h(jj,nz2-n,1,1) = 0.5*cmplx(aimag(h(jj,n,1,1) + s),         &
     &                                  real(h(jj,n,1,1) - s))
            h(jj,n,1,1) = 0.5*cmplx(real(h(jj,n,1,1) + s),              &
     &                              aimag(h(jj,n,1,1) - s))
   80       continue
   90       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 110 n = 2, nzh
            do 100 jj = 1, ndim
            s = h(jj,nz2-n,1,k1)
            h(jj,nz2-n,1,k1) = 0.5*cmplx(aimag(h(jj,n,1,k1) + s),       &
     &                                   real(h(jj,n,1,k1) - s))
            h(jj,n,1,k1) = 0.5*cmplx(real(h(jj,n,1,k1) + s),            &
     &                               aimag(h(jj,n,1,k1) - s))
  100       continue
  110       continue
         endif
      endif
      return
! forward fourier transform
  120 nrzb = nxhyz/nz
      nrz = nxyz/nz
! scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 140 n = 2, nzh
            do 130 jj = 1, ndim
            s = cmplx(aimag(h(jj,nz2-n,1,1)),real(h(jj,nz2-n,1,1)))
            h(jj,nz2-n,1,1) = conjg(h(jj,n,1,1) - s)
            h(jj,n,1,1) = h(jj,n,1,1) + s
  130       continue
  140       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 160 n = 2, nzh
            do 150 jj = 1, ndim
            s = cmplx(aimag(h(jj,nz2-n,1,k1)),real(h(jj,nz2-n,1,k1)))
            h(jj,nz2-n,1,k1) = conjg(h(jj,n,1,k1) - s)
            h(jj,n,1,k1) = h(jj,n,1,k1) + s
  150       continue
  160       continue
         endif
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1)
      do 230 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
! bit-reverse array elements in z
      do 180 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 170 jj = 1, ndim
         t1 = h(jj,l1,i,n)
         h(jj,l1,i,n) = h(jj,l,i,n)
         h(jj,l,i,n) = t1
  170    continue
      endif
  180 continue
! first transform in z
      do 220 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 190 jj = 1, ndim
      t1 = s*h(jj,j2,i,n)
      h(jj,j2,i,n) = h(jj,j1,i,n) - t1
      h(jj,j1,i,n) = h(jj,j1,i,n) + t1
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSWAPC32N(f,s,isign,nxh,kypi,kypt,nxvh,kzpp,kypd,kzpd&
     &,ndim)
! this subroutine swaps components for multiple ffts
! f = input  array
! s = scratch array
! isign = (-1,1) = swap (real-to-complex,complex-to-real)
! nxh = complex dimension in x direction
! kypi/kypt = initial/final y index used
! nxvh = half of the second dimension of f
! kzpp = number of z indices used
! kypd, kzpd = third and fourth dimension of f
! ndim = leading dimension of array f
      implicit none
      integer isign, nxh, kypi, kypt, nxvh, kzpp, kypd, kzpd, ndim
      real f, s
      dimension f(ndim,2*nxvh,kypd,kzpd), s(2*ndim*nxvh,kzpd)
! local data
      integer i, j, k, l, ioff
! swap complex components
! real to complex
      if (isign.lt.0) then
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ioff)
         do 70 l = 1, kzpp
         do 60 k = kypi, kypt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1,l) = f(i,2*j-1,k,l)
         s(2*i+ioff,l) = f(i,2*j,k,l)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k,l) = s(i+ioff,l)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k,l) = s(i+ioff,l)
   40    continue
   50    continue
   60    continue
   70    continue
!$OMP END PARALLEL DO
      else if (isign.gt.0) then
! swap complex components
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ioff)
         do 140 l = 1, kzpp
         do 130 k = kypi, kypt
         do 100 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 80 i = 1, ndim
         s(i+ioff,l) = f(i,2*j-1,k,l)
   80    continue
         ioff = ioff + ndim
         do 90 i = 1, ndim
         s(i+ioff,l) = f(i,2*j,k,l)
   90    continue
  100    continue
         do 120 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 110 i = 1, ndim
         f(i,2*j-1,k,l) = s(2*i+ioff-1,l)
         f(i,2*j,k,l) = s(2*i+ioff,l)
  110    continue
  120    continue
  130    continue
  140    continue
!$OMP END PARALLEL DO
      endif
      return
      end
