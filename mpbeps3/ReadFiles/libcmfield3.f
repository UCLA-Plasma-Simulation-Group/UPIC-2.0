!-----------------------------------------------------------------------
! Library for processing complex periodic 3d scalar data with OpenMP
! WRMODES3 places selected 3d fourier components into complex scalar
!          array.
! WRVMODES3 places selected 3d fourier components into complex vector
!           array.
! CSPECT2 performs frequency analysis of complex time series
! ICSPECT1 performs incremental frequency analysis of complex time
!          series for one time step
! IVCSPECT2 performs incremental frequency analysis of complex vector
!           time series for one time step
! WFFT3RINIT calculates tables needed by a two dimensional real to
!            complex fast fourier transform and its inverse.
! WFFT3RMX performs real to complex fft and its inverse for scalar array,
!          with packed data.
! WFFT3RM3 performs real to complex fft and its inverse for 3 component
!          vector array, with packed data.
! FFT3RMXY performs x part of scalar 2d real/complex FFT
! FFT3RMXZ performs y part of scalar 2d real/complex FFT
! FFT3RM3XY performs x part of 3 component vector 2d real/complex FFT
! FFT3RM3Z performs y part of 3 component vector 2d real/complex FFT
! written by Viktor K. Decyk, UCLA
!-----------------------------------------------------------------------
      subroutine WRMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh,  &
     &nyv,nzv,modesxd,modesyd,modeszd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex array pott and stores them into a packed complex
! array pot
! modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
! kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! nxvh = first dimension of output array pot, nxvh >= nx/2
! nyv = second dimension of output array pot, nyv >= ny
! nzv = third dimension of output array pot, nzv >= nz
! modesxd = first dimension of input array pott, modesxd >= modesx
! modesyd = second dimension of output array pott,
! where modesyd  >= min(2*modesy-1,ny)
! modeszd = third dimension of output array pott,
! where modeszd  >= min(2*modesz-1,nz)
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, nxvh, nyv, nzv
      integer modesxd, modesyd, modeszd
      complex pot, pott
      dimension pot(nxvh,nyv,nzv), pott(modesxd,modesyd,modeszd)
! local data
      integer nxh, nyh, nzh, lmax, kmax, jmax, ny2, nz2, j, k, l
      integer j1, k1, l1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      ny2 = ny + 2
      nz2 = nz + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      lmax = min0(modesz,nzh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 90 l = 2, lmax
      l1 = nz2 - l
      do 30 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pot(j,k,l) = pott(j,2*k-2,2*l-2)
      pot(j,k1,l) = pott(j,2*k-1,2*l-2)
      pot(j,k,l1) = pott(j,2*k-2,2*l-1)
      pot(j,k1,l1) = pott(j,2*k-1,2*l-1)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
   20 continue
! mode numbers kx = 0, nx/2
      pot(1,k,l) = pott(1,2*k-2,2*l-2)
      pot(1,k1,l) = zero
      pot(1,k,l1) = pott(1,2*k-2,2*l-1)
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1,l) = conjg(pott(j1,2*k-2,2*l-1))
         pot(1,k1,l1) = conjg(pott(j1,2*k-2,2*l-2))
      endif
   30 continue
      do 50 k = kmax+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
   40 continue
   50 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, jmax
      pot(j,1,l) = pott(j,1,2*l-2)
      pot(j,k1,l) = zero
      pot(j,1,l1) = pott(j,1,2*l-1)
      pot(j,k1,l1) = zero
   60 continue
! mode numbers kx = 0, nx/2
      pot(1,1,l) = pott(1,1,2*l-2)
      pot(1,k1,l) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      do 70 j = jmax+1, nxh
      pot(j,1,l) = zero
      pot(j,k1,l) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
   70 continue
      if (modesx.gt.nxh) then
         pot(1,1,l1) = conjg(pott(j1,1,2*l-2))
      endif
! mode numbers ky = ny/2
      if (modesy.gt.nyh) then
         do 80 j = 2, jmax
         pot(j,k1,l) = pott(j,ny,2*l-2)
         pot(j,k1,l1) = pott(j,ny,2*l-1)
   80    continue
         pot(1,k1,l) = pott(1,ny,2*l-2)
         if (modesx.gt.nxh) then
            pot(1,k1,l1) = conjg(pott(j1,ny,2*l-2))
         endif
      endif
   90 continue
!$OMP END PARALLEL DO
      do 130 l = modesz+1, nzh
      l1 = nz2 - l
      do 110 k = 2, nyh
      k1 = ny2 - k
      do 100 j = 1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  100 continue
  110 continue
      k1 = nyh + 1
      do 120 j = 1, nxh
      pot(j,1,l) = zero
      pot(j,k1,l) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  120 continue
  130 continue
! mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 160 k = 2, kmax
      k1 = ny2 - k
      do 140 j = 2, jmax
      pot(j,k,1) = pott(j,2*k-2,1)
      pot(j,k1,1) = pott(j,2*k-1,1)
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  140 continue
      do 150 j = jmax+1, nxh
      pot(j,k,1) = zero
      pot(j,k1,1) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  150 continue
! mode numbers kx = 0, nx/2
      pot(1,k,1) = pott(1,2*k-2,1)
      pot(1,k1,1) = zero
      pot(1,k,l1) = zero
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1,1) = conjg(pott(j1,2*k-2,1))
      endif
  160 continue
      do 180 k = modesy+1, nyh
      k1 = ny2 - k
      do 170 j = 1, nxh
      pot(j,k,1) = zero
      pot(j,k1,1) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  170 continue
  180 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 190 j = 2, jmax
      pot(j,1,1) = pott(j,1,1)
      pot(j,k1,1) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  190 continue
      do 200 j = jmax+1, nxh
      pot(j,1,1) = zero
      pot(j,k1,1) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  200 continue
      pot(1,1,1) = cmplx(real(pott(1,1,1)),0.0)
      pot(1,k1,1) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,1,1) = cmplx(real(pot(1,1,1)),real(pott(j1,1,1)))
      endif
      if (modesy.gt.nyh) then
         do 210 j = 2, jmax
         pot(j,k1,1) = pott(j,ny,1)
  210    continue
         pot(1,k1,1) = cmplx(real(pott(1,ny,1)),0.0)
         if (modesx.gt.nxh) then
            pot(1,k1,1) = cmplx(real(pot(1,k1,1)),real(pott(j1,ny,1)))
         endif
      endif
! mode numbers kz = nz/2
      if (modesz.gt.nzh) then
         do 230 k = 2, kmax
         k1 = ny2 - k
         do 220 j = 2, jmax
         pot(j,k,l1) = pott(j,2*k-2,nz)
         pot(j,k1,l1) = pott(j,2*k-1,nz)
  220    continue
! mode numbers kx = 0, nx/2
         pot(1,k,l1) = pott(1,2*k-2,nz)
         if (modesx.gt.nxh) then
            pot(1,k1,l1) = conjg(pott(j1,2*k-2,nz))
         endif
  230 continue
! mode numbers ky = 0
         do 240 j = 2, jmax
         pot(j,1,l1) = pott(j,1,nz)
  240    continue
         pot(1,1,l1) = cmplx(real(pott(1,1,nz)),0.0)
         if (modesx.gt.nxh) then
            pot(1,1,l1) = cmplx(real(pot(1,1,l1)),real(pott(j1,1,nz)))
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            do 250 j = 2, jmax
            pot(j,k1,l1) = pott(j,ny,nz)
  250       continue
            pot(1,k1,l1) = cmplx(real(pott(1,ny,nz)),0.0)
            if (modesx.gt.nxh) then
               pot(1,k1,l1) = cmplx(real(pot(1,k1,l1)),                 &
     &                              real(pott(j1,ny,nz)))
            endif
         endif
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WRVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz,ndim&
     &,nxvh,nyv,nzv,modesxd,modesyd,modeszd)
! this subroutine extracts lowest order modes from a location in an
! unpacked complex vector array vpott and stores them into a packed
! complex vector array vpot
! modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
! kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
! nx/ny/nz = system length in x/y/z direction
! modesx/modesy/modesz = number of modes to store in x/y/z direction,
! where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
! ndim = number of field arrays, must be >= 1
! nxvh = second dimension of input array vpot, nxvh >= nx/2
! nyv = third dimension of input array vpot, nyv >= ny
! nzv = fourth dimension of input array vpot, nzv >= nz
! modesxd = second dimension of output array vpott, modesxd >= modesx
! modesyd = third dimension of output array vpott,
! where modesyd  >= min(2*modesy-1,ny)
! modeszd = fourth dimension of output array vpott,
! where modeszd  >= min(2*modesz-1,nz)
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, ndim, nxvh, nyv, nzv
      integer modesxd, modesyd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv,nzv)
      dimension vpott(ndim,modesxd,modesyd,modeszd)
! local data
      integer nxh, nyh, nzh, lmax, kmax, jmax, ny2, nz2, i, j, k, l
      integer j1, k1, l1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      ny2 = ny + 2
      nz2 = nz + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      lmax = min0(modesz,nzh)
      zero = cmplx(0.0,0.0)
      j1 = nxh + 1
! mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 190 l = 2, lmax
      l1 = nz2 - l
      do 60 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j,k,l) = vpott(i,j,2*k-2,2*l-2)
      vpot(i,j,k1,l) = vpott(i,j,2*k-1,2*l-2)
      vpot(i,j,k,l1) = vpott(i,j,2*k-2,2*l-1)
      vpot(i,j,k1,l1) = vpott(i,j,2*k-1,2*l-1)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j,k,l) = zero
      vpot(i,j,k1,l) = zero
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
   30 continue
   40 continue
! mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1,k,l) = vpott(i,1,2*k-2,2*l-2)
      vpot(i,1,k1,l) = zero
      vpot(i,1,k,l1) = vpott(i,1,2*k-2,2*l-1)
      vpot(i,1,k1,l1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,k1,l) = conjg(vpott(i,j1,2*k-2,2*l-1))
         vpot(i,1,k1,l1) = conjg(vpott(i,j1,2*k-2,2*l-2))
      endif
   50 continue
   60 continue
      do 90 k = kmax+1, nyh
      k1 = ny2 - k
      do 80 j = 1, nxh
      do 70 i = 1, ndim
      vpot(i,j,k,l) = zero
      vpot(i,j,k1,l) = zero
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
   70 continue
   80 continue
   90 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 110 j = 2, jmax
      do 100 i = 1, ndim
      vpot(i,j,1,l) = vpott(i,j,1,2*l-2)
      vpot(i,j,k1,l) = zero
      vpot(i,j,1,l1) = vpott(i,j,1,2*l-1)
      vpot(i,j,k1,l1) = zero
  100 continue
  110 continue
! mode numbers kx = 0, nx/2
      do 120 i = 1, ndim
      vpot(i,1,1,l) = vpott(i,1,1,2*l-2)
      vpot(i,1,k1,l) = zero
      vpot(i,1,1,l1) = zero
      vpot(i,1,k1,l1) = zero
  120 continue
      do 140 j = jmax+1, nxh
      do 130 i = 1, ndim
      vpot(i,j,1,l) = zero
      vpot(i,j,k1,l) = zero
      vpot(i,j,1,l1) = zero
      vpot(i,j,k1,l1) = zero
  130 continue
  140 continue
      do 150 i = 1, ndim
      if (modesx.gt.nxh) then
         vpot(i,1,1,l1) = conjg(vpott(i,j1,1,2*l-2))
      endif
  150 continue
! mode numbers ky = ny/2
      if (modesy.gt.nyh) then
         do 170 j = 2, jmax
         do 160 i = 1, ndim
         vpot(i,j,k1,l) = vpott(i,j,ny,2*l-2)
         vpot(i,j,k1,l1) = vpott(i,j,ny,2*l-1)
  160    continue
  170    continue
         do 180 i = 1, ndim
         vpot(i,1,k1,l) = vpott(i,1,ny,2*l-2)
         if (modesx.gt.nxh) then
            vpot(i,1,k1,l1) = conjg(vpott(i,j1,ny,2*l-2))
         endif
  180    continue
      endif
  190 continue
!$OMP END PARALLEL DO
      do 250 l = modesz+1, nzh
      l1 = nz2 - l
      do 220 k = 2, nyh
      k1 = ny2 - k
      do 210 j = 1, nxh
      do 200 i = 1, ndim
      vpot(i,j,k,l) = zero
      vpot(i,j,k1,l) = zero
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
  200 continue
  210 continue
  220 continue
      k1 = nyh + 1
      do 240 j = 1, nxh
      do 230 i = 1, ndim
      vpot(i,j,1,l) = zero
      vpot(i,j,k1,l) = zero
      vpot(i,j,1,l1) = zero
      vpot(i,j,k1,l1) = zero
  230 continue
  240 continue
  250 continue
! mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 310 k = 2, kmax
      k1 = ny2 - k
      do 270 j = 2, jmax
      do 260 i = 1, ndim
      vpot(i,j,k,1) = vpott(i,j,2*k-2,1)
      vpot(i,j,k1,1) = vpott(i,j,2*k-1,1)
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
  260 continue
  270 continue
      do 290 j = jmax+1, nxh
      do 280 i = 1, ndim
      vpot(i,j,k,1) = zero
      vpot(i,j,k1,1) = zero
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
  280 continue
  290 continue
! mode numbers kx = 0, nx/2
      do 300 i = 1, ndim
      vpot(i,1,k,1) = vpott(i,1,2*k-2,1)
      vpot(i,1,k1,1) = zero
      vpot(i,1,k,l1) = zero
      vpot(i,1,k1,l1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,k1,1) = conjg(vpott(i,j1,2*k-2,1))
      endif
  300 continue
  310 continue
      do 340 k = modesy+1, nyh
      k1 = ny2 - k
      do 330 j = 1, nxh
      do 320 i = 1, ndim
      vpot(i,j,k,1) = zero
      vpot(i,j,k1,1) = zero
      vpot(i,j,k,l1) = zero
      vpot(i,j,k1,l1) = zero
  320 continue
  330 continue
  340 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 360 j = 2, jmax
      do 350 i = 1, ndim
      vpot(i,j,1,1) = vpott(i,j,1,1)
      vpot(i,j,k1,1) = zero
      vpot(i,j,1,l1) = zero
      vpot(i,j,k1,l1) = zero
  350 continue
  360 continue
      do 380 j = jmax+1, nxh
      do 370 i = 1, ndim
      vpot(i,j,1,1) = zero
      vpot(i,j,k1,1) = zero
      vpot(i,j,1,l1) = zero
      vpot(i,j,k1,l1) = zero
  370 continue
  380 continue
      do 390 i = 1, ndim
      vpot(i,1,1,1) = cmplx(real(vpott(i,1,1,1)),0.0)
      vpot(i,1,k1,1) = zero
      vpot(i,1,1,l1) = zero
      vpot(i,1,k1,l1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,1,1) = cmplx(real(vpot(i,1,1,1)),                     &
     &                         real(vpott(i,j1,1,1)))
      endif
  390 continue
      if (modesy.gt.nyh) then
         do 410 j = 2, jmax
         do 400 i = 1, ndim
         vpot(i,j,k1,1) = vpott(i,j,ny,1)
  400    continue
  410    continue
         do 420 i = 1, ndim
         vpot(i,1,k1,1) = cmplx(real(vpott(i,1,ny,1)),0.0)
         if (modesx.gt.nxh) then
            vpot(i,1,k1,1) = cmplx(real(vpot(i,1,k1,1)),                &
     &                             real(vpott(i,j1,ny,1)))
         endif
  420    continue
      endif
! mode numbers kz = nz/2
      if (modesz.gt.nzh) then
         do 460 k = 2, kmax
         k1 = ny2 - k
         do 440 j = 2, jmax
         do 430 i = 1, ndim
         vpot(i,j,k,l1) = vpott(i,j,2*k-2,nz)
         vpot(i,j,k1,l1) = vpott(i,j,2*k-1,nz)
  430    continue
  440    continue
! mode numbers kx = 0, nx/2
         do 450 i = 1, ndim
         vpot(i,1,k,l1) = vpott(i,1,2*k-2,nz)
         if (modesx.gt.nxh) then
            vpot(i,1,k1,l1) = conjg(vpott(i,j1,2*k-2,nz))
         endif
  450    continue
  460    continue
! mode numbers ky = 0
         do 480 j = 2, jmax
         do 470 i = 1, ndim
         vpot(i,j,1,l1) = vpott(i,j,1,nz)
  470    continue
  480    continue
         do 490 i = 1, ndim
         vpot(i,1,1,l1) = cmplx(real(vpott(i,1,1,nz)),0.0)
         if (modesx.gt.nxh) then
            vpot(i,1,1,l1) = cmplx(real(vpot(i,1,1,l1)),                &
     &                             real(vpott(i,j1,1,nz)))
         endif
  490    continue
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            do 510 j = 2, jmax
            do 500 i = 1, ndim
            vpot(i,j,k1,l1) = vpott(i,j,ny,nz)
  500       continue
  510       continue
            do 520 i = 1, ndim
            vpot(i,1,k1,l1) = cmplx(real(vpott(i,1,ny,nz)),0.0)
            if (modesx.gt.nxh) then
               vpot(i,1,k1,l1) = cmplx(real(vpot(i,1,k1,l1)),           &
     &                                 real(vpott(i,j1,ny,nz)))
            endif
  520       continue
         endif
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
      subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
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
      subroutine WFFT3RMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,&
     &nxhyzd,nxyzhd)
! wrapper function for real to complex fft, with packed data
! parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
! calculate range of indices
      ny = 2**indy
      nz = 2**indz
! inverse fourier transform
      if (isign.lt.0) then
! perform xy fft
         call FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
! perform z fft
         call FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform z fft
         call FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
! perform xy fft
         call FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT3RM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,&
     &nxhyzd,nxyzhd)
! wrapper function for 3 2d real to complex ffts, with packed data
! parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
! calculate range of indices
      ny = 2**indy
      nz = 2**indz
! inverse fourier transform
      if (isign.lt.0) then
! perform xy fft
         call FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd&
     &,nzd,nxhyzd,nxyzhd)
! perform z fft
         call FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform z fft
         call FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
! perform xy fft
         call FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd&
     &,nzd,nxhyzd,nxyzhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd,&
     &nyd,nzd,nxhyzd,nxyzhd)
! this subroutine performs the x-y part of a three dimensional real to
! complex fast fourier transform and its inverse, for a subset of z,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny*nz
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform in x and y is performed
! f(n,m,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)*
!       exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform in x and y is performed
! f(j,k,l) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nzi = initial z index used
! nzp = number of z indices used
! nxhd = first dimension of f
! nyd,nzd = second and third dimensions of f
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = maximum of (nx,ny,nz)/2
! fourier coefficients are stored as follows:
! f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
! where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
! f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
! where ny/2+2 <= k <= ny and 1 <= l <= nz, and
! f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
! f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
! where nz/2+2 <= l <= nz, and
! imag(f(1,1,1)) = real part of mode nx/2,0,0
! imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
! imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
! imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb
      integer i, j, k, l, n, j1, j2, k1, k2, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 180
! inverse fourier transform
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ani,t1,t2,t3)
      do 170 n = nzi, nzt
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 10 i = 1, ny
         t1 = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t1
   10    continue
      endif
   20 continue
! first transform in x
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
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = ani*(t1 + t2)
      f(nxh2-j,k,n) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = 1, ny
      f(nxhh+1,k,n) = ani*conjg(f(nxhh+1,k,n))
      f(1,k,n) = ani*cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),            &
     &                     real(f(1,k,n)) - aimag(f(1,k,n)))
   90 continue
! bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 100 i = 1, nxh
         t1 = f(i,k1,n)
         f(i,k1,n) = f(i,k,n)
         f(i,k,n) = t1
  100    continue
      endif
  110 continue
! then transform in y
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 120 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  120 continue
  130 continue
  140 continue
  150 continue
! unscramble modes kx = 0, nx/2
      do 160 k = 2, nyh
      t1 = f(1,ny2-k,n)
      f(1,ny2-k,n) = 0.5*cmplx(aimag(f(1,k,n) + t1),real(f(1,k,n) - t1))
      f(1,k,n) = 0.5*cmplx(real(f(1,k,n) + t1),aimag(f(1,k,n) - t1))
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  180 nryb = nxhyz/ny
      nry = nxyz/ny
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2,t3)
      do 350 n = nzi, nzt
! scramble modes kx = 0, nx/2
      do 190 k = 2, nyh
      t1 = cmplx(aimag(f(1,ny2-k,n)),real(f(1,ny2-k,n)))
      f(1,ny2-k,n) = conjg(f(1,k,n) - t1)
      f(1,k,n) = f(1,k,n) + t1
  190 continue
! bit-reverse array elements in y
      do 210 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 200 i = 1, nxh
         t1 = f(i,k1,n)
         f(i,k1,n) = f(i,k,n)
         f(i,k,n) = t1
  200    continue
      endif
  210 continue
! then transform in y
      do 250 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  220 continue
  230 continue
  240 continue
  250 continue
! scramble coefficients
      kmr = nxyz/nx
      do 270 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 260 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = t1 + t2
      f(nxh2-j,k,n) = conjg(t1 - t2)
  260 continue
  270 continue
      do 280 k = 1, ny
      f(nxhh+1,k,n) = 2.0*conjg(f(nxhh+1,k,n))
      f(1,k,n) = cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),                &
     &                 real(f(1,k,n)) - aimag(f(1,k,n)))
  280 continue
! bit-reverse array elements in x
      do 300 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 290 i = 1, ny
         t1 = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t1
  290    continue
      endif
  300 continue
! finally transform in x
      do 340 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 330 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 320 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 310 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd,&
     &nyd,nzd,nxhyzd,nxyzhd)
! this subroutine performs the z part of a three dimensional real to
! complex fast fourier transform and its inverse, for a subset of y,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny*nz
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, an inverse fourier transform in z is performed
! f(j,k,l) = sum(f(j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
! if isign = 1, a forward fourier transform in z is performed
! f(n,m,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = first dimension of f
! nyd,nzd = second and third dimensions of f
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = maximum of (nx,ny,nz)/2
! fourier coefficients are stored as follows:
! f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
! where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
! f(1,k,l), = real, imaginary part of mode nx/2,k-1,l-1,
! where ny/2+2 <= k <= ny and 1 <= l <= nz, and
! f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
! f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
! where nz/2+2 <= l <= nz, and
! imag(f(1,1,1)) = real part of mode nx/2,0,0
! imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
! imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
! imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz, nrzb
      integer i, j, k, l, n, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 90
! inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,t1,t2)
      do 70 n = nyi, nyt
! bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 10 i = 1, nxh
         t1 = f(i,n,l1)
         f(i,n,l1) = f(i,n,l)
         f(i,n,l) = t1
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
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      do 80 n = 2, nzh
      if (nyi.eq.1) then
         t1 = f(1,1,nz2-n)
         f(1,1,nz2-n) = 0.5*cmplx(aimag(f(1,1,n) + t1),                 &
     &                            real(f(1,1,n) - t1))
         f(1,1,n) = 0.5*cmplx(real(f(1,1,n) + t1),aimag(f(1,1,n) - t1))
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = f(1,nyh+1,nz2-n)
         f(1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(1,nyh+1,n) + t1),         &
     &                                real(f(1,nyh+1,n) - t1))
         f(1,nyh+1,n) = 0.5*cmplx(real(f(1,nyh+1,n) + t1),              &
     &                            aimag(f(1,nyh+1,n) - t1))
      endif
   80 continue
      return
! forward fourier transform
   90 nrzb = nxhyz/nz
      nrz = nxyz/nz
! scramble modes kx = 0, nx/2
      do 100 n = 2, nzh
      if (nyi.eq.1) then
         t1 = cmplx(aimag(f(1,1,nz2-n)),real(f(1,1,nz2-n)))
         f(1,1,nz2-n) = conjg(f(1,1,n) - t1)
         f(1,1,n) = f(1,1,n) + t1
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = cmplx(aimag(f(1,nyh+1,nz2-n)),real(f(1,nyh+1,nz2-n)))
         f(1,nyh+1,nz2-n) = conjg(f(1,nyh+1,n) - t1)
         f(1,nyh+1,n) = f(1,nyh+1,n) + t1
      endif
  100 continue
! bit-reverse array elements in z
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,t1,t2)
      do 170 n = nyi, nyt
      do 120 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 110 i = 1, nxh
         t1 = f(i,n,l1)
         f(i,n,l1) = f(i,n,l)
         f(i,n,l) = t1
  110    continue
      endif
  120 continue
! first transform in z
      do 160 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 130 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
  130 continue
  140 continue
  150 continue
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd&
     &,nyd,nzd,nxhyzd,nxyzhd)
! this subroutine performs the x-y part of 3 three dimensional complex
! to real fast fourier transforms and their inverses, for a subset of z,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny*nz
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, three inverse fourier transforms in x and y are
! performed
! f(1:3,n,m,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
!       *exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, three forward fourier transforms in x and y are
! performed
! f(1:3,j,k,l) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nzi = initial z index used
! nzp = number of z indices used
! nxhd = second dimension of f
! nyd,nzd = third and fourth dimensions of f
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = maximum of (nx,ny,nz)/2
! fourier coefficients are stored as follows:
! f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
! where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
! f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
! where ny/2+2 <= k <= ny and 1 <= l <= nz, and
! f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
! f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
! where nz/2+2 <= l <= nz, and
! imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
! imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
! imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
! imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd,nxyzhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb
      integer i, j, k, l, n, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 230
! inverse fourier transform
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,ani,t1,t2,&
!$OMP& t3,t4)
      do 220 n = nzi, nzt
! swap complex components
      do 20 i = 1, ny
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = aimag(f(2,j,i,n))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
   20 continue
! bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
      do 30 i = 1, ny
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
   30    continue
      endif
   40 continue
! first transform in x
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
   50 continue
   60 continue
   70 continue
   80 continue
! unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = 1, ny
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = ani*(t1 + t2)
      f(jj,nxh2-j,k,n) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.0*ani
      do 130 k = 1, ny
      do 120 jj = 1, 3
      f(jj,nxhh+1,k,n) = ani*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = ani*cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),   &
     &                        real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  120 continue
  130 continue
! bit-reverse array elements in y
      do 150 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 140 i = 1, nxh
         t1 = f(1,i,k1,n)
         t2 = f(2,i,k1,n)
         t3 = f(3,i,k1,n)
         f(1,i,k1,n) = f(1,i,k,n)
         f(2,i,k1,n) = f(2,i,k,n)
         f(3,i,k1,n) = f(3,i,k,n)
         f(1,i,k,n) = t1
         f(2,i,k,n) = t2
         f(3,i,k,n) = t3
  140    continue
      endif
  150 continue
! then transform in y
      do 190 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 160 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  160 continue
  170 continue
  180 continue
  190 continue
! unscramble modes kx = 0, nx/2
      do 210 k = 2, nyh
      do 200 jj = 1, 3
      t1 = f(jj,1,ny2-k,n)
      f(jj,1,ny2-k,n) = 0.5*cmplx(aimag(f(jj,1,k,n) + t1),              &
     &                            real(f(jj,1,k,n) - t1))
      f(jj,1,k,n) = 0.5*cmplx(real(f(jj,1,k,n) + t1),                   &
     &                        aimag(f(jj,1,k,n) - t1))
  200 continue
  210 continue
  220 continue
!$OMP END PARALLEL DO
      return
! forward fourier transform
  230 nryb = nxhyz/ny
      nry = nxyz/ny
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,t1,t2,t3, &
!$OMP& t4)
      do 450 n = nzi, nzt
! scramble modes kx = 0, nx/2
      do 250 k = 2, nyh
      do 240 jj = 1, 3
      t1 = cmplx(aimag(f(jj,1,ny2-k,n)),real(f(jj,1,ny2-k,n)))
      f(jj,1,ny2-k,n) = conjg(f(jj,1,k,n) - t1)
      f(jj,1,k,n) = f(jj,1,k,n) + t1
  240 continue
  250 continue
! bit-reverse array elements in y
      do 270 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 260 i = 1, nxh
         t1 = f(1,i,k1,n)
         t2 = f(2,i,k1,n)
         t3 = f(3,i,k1,n)
         f(1,i,k1,n) = f(1,i,k,n)
         f(2,i,k1,n) = f(2,i,k,n)
         f(3,i,k1,n) = f(3,i,k,n)
         f(1,i,k,n) = t1
         f(2,i,k,n) = t2
         f(3,i,k,n) = t3
  260 continue
      endif
  270 continue
! then transform in y
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 280 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  280 continue
  290 continue
  300 continue
  310 continue
! scramble coefficients
      kmr = nxyz/nx
      do 340 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 330 k = 1, ny
      do 320 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = t1 + t2
      f(jj,nxh2-j,k,n) = conjg(t1 - t2)
  320 continue
  330 continue
  340 continue
      do 360 k = 1, ny
      do 350 jj = 1, 3
      f(jj,nxhh+1,k,n) = 2.0*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),       &
     &                    real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  350 continue
  360 continue
! bit-reverse array elements in x
      do 380 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
      do 370 i = 1, ny
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
  370 continue
      endif
  380 continue
! finally transform in x
      do 420 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 410 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 400 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 390 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
  390 continue
  400 continue
  410 continue
  420 continue
! swap complex components
      do 440 i = 1, ny
      do 430 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,aimag(f(1,j,i,n)))
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  430 continue
  440 continue
  450 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd,&
     &nyd,nzd,nxhyzd,nxyzhd)
! this subroutine performs the z part of 3 three dimensional complex to
! real fast fourier transforms and their inverses, for a subset of y,
! using complex arithmetic, with OpenMP
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny*nz
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
! if isign = -1, three inverse fourier transforms in z are performed
! f(1:3,j,k,l) = sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
! if isign = 1, three forward fourier transforms in z are performed
! f(1:3,n,m,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = second dimension of f
! nyd,nzd = third and fourth dimensions of f
! nxhyzd = maximum of (nx/2,ny,nz)
! nxyzhd = maximum of (nx,ny,nz)/2
! fourier coefficients are stored as follows:
! f(1:3,2*j-1,k,l),f(2*j,k,l) = real, imaginary part of mode j-1,k-1,l-1
! where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
! f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
! where ny/2+2 <= k <= ny and 1 <= l <= nz, and
! f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
! f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
! where nz/2+2 <= l <= nz, and
! imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
! imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
! imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
! imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
! using jpl storage convention, as described in:
! E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
! Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
! Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
! December 1993.
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
! local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz, nrzb
      integer i, j, k, l, n, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 110
! inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,t1,t2,t3,t4)
      do 70 n = nyi, nyt
! bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
      do 10 i = 1, nxh
         t1 = f(1,i,n,l1)
         t2 = f(2,i,n,l1)
         t3 = f(3,i,n,l1)
         f(1,i,n,l1) = f(1,i,n,l)
         f(2,i,n,l1) = f(2,i,n,l)
         f(3,i,n,l1) = f(3,i,n,l)
         f(1,i,n,l) = t1
         f(2,i,n,l) = t2
         f(3,i,n,l) = t3
   10 continue
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
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
! unscramble modes kx = 0, nx/2
      do 100 n = 2, nzh
      if (nyi.eq.1) then
         do 80 jj = 1, 3
         t1 = f(jj,1,1,nz2-n)
         f(jj,1,1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,1,n) + t1),           &
     &                               real(f(jj,1,1,n) - t1))
         f(jj,1,1,n) = 0.5*cmplx(real(f(jj,1,1,n) + t1),                &
     &                           aimag(f(jj,1,1,n) - t1))
   80    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 90 jj = 1, 3
         t1 = f(jj,1,nyh+1,nz2-n)
         f(jj,1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,nyh+1,n) + t1),   &
     &                                  real(f(jj,1,nyh+1,n) - t1))
         f(jj,1,nyh+1,n) = 0.5*cmplx(real(f(jj,1,nyh+1,n) + t1),        &
     &                              aimag(f(jj,1,nyh+1,n) - t1))
   90    continue
      endif
  100 continue
      return
! forward fourier transform
  110 nrzb = nxhyz/nz
      nrz = nxyz/nz
! scramble modes kx = 0, nx/2
      do 140 n = 2, nzh
      if (nyi.eq.1) then
         do 120 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,1,nz2-n)),real(f(jj,1,1,nz2-n)))
         f(jj,1,1,nz2-n) = conjg(f(jj,1,1,n) - t1)
         f(jj,1,1,n) = f(jj,1,1,n) + t1
  120    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 130 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,nyh+1,nz2-n)),                         &
     &              real(f(jj,1,nyh+1,nz2-n)))
         f(jj,1,nyh+1,nz2-n) = conjg(f(jj,1,nyh+1,n) - t1)
         f(jj,1,nyh+1,n) = f(jj,1,nyh+1,n) + t1
  130    continue
      endif
  140 continue
! bit-reverse array elements in z
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,t1,t2,t3,t4)
      do 210 n = nyi, nyt
      do 160 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 150 i = 1, nxh
         t1 = f(1,i,n,l1)
         t2 = f(2,i,n,l1)
         t3 = f(3,i,n,l1)
         f(1,i,n,l1) = f(1,i,n,l)
         f(2,i,n,l1) = f(2,i,n,l)
         f(3,i,n,l1) = f(3,i,n,l)
         f(1,i,n,l) = t1
         f(2,i,n,l) = t2
         f(3,i,n,l) = t3
  150    continue
      endif
  160 continue
! first transform in z
      do 200 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 170 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
!$OMP END PARALLEL DO
      return
      end
