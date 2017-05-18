!-----------------------------------------------------------------------
! Fortran Library for performing Fast Sine/Cosine Transforms
! 2D MPI/OpenMP PIC Code:
! WPFST2RINIT calculates tables needed by 2d Sine/Cosine Transforms
! WPPFSST2RM wrapper function for scalar 2d real sine/sine transform
! WPPFSCT2RM wrapper function for scalar 2d real sine/cosine transform
! WPPFCST2RM wrapper function for scalar 2d real cosine/sine transform
! WPPFCCT2RM wrapper function for scalar 2d real cosine/cosine transform
! PPFST2RMXX performs x part of scalar 2d sine Transform
! PPFCT2RMXX performs x part of scalar 2d cosine Transform
! PPFST2RMXY performs y part of scalar 2d sine Transform
! PPFCT2RMXY performs y part of scalar 2d cosine Transform
! WPPFCST2RM2 wrapper function for 2 component vector 2d real
!             cosine/sine transforms for electric field with dirichlet
!             or magnetic field with neumann boundary conditions
! WPPFSCT2RM2 wrapper function for 2 component vector 2d real
!             sine/cosine transforms for magnetic field with dirichlet
!             or electric field with neumann boundary conditions
! PPFCST2RM2X performs x part of 2 component vector 2d real cosine/sine
!             transform
! PPFSCT2RM2X performs x part of 2 component vector 2d real sine/cosine
!             transform
! PPFSCT2RM2Y performs y part of 2 component vector 2d real sine/cosine
!             transform
! PPFCST2RM2Y performs y part of 2 component vector 2d real cosine/sine
!             transform
! WPPFCST2RM3 wrapper function for 3 component vector 2d real
!             cosine/sine transforms for electric field with dirichlet
!             or magnetic field with neumann boundary conditions
! WPPFSCT2RM3 wrapper function for 3 component vector 2d real
!             sine/cosine transforms for magnetic field with dirichlet
!             or electric field with neumann boundary conditions
! PPFCSST2RM3X performs x part of 3 component vector 2d real 
!              cosine/sine/sine transform
! PPFSCCT2RM3X performs x part of 3 component vector 2d real
!              sine/cosine/cosine transform
! PPFSCST2RM3Y performs y part of 3 component vector 2d real
!              sine/cosine/sine transform
! PPFCSCT2RM3Y performs y part of 3 component vector 2d real
!              cosine/sine/cosine transform
! WPPFSCT2RM4 wrapper function for 4 component tensor 2d real
!             sine/cosine transforms
! PPFSCCST2RM4X performs x part of 4 component vector 2d real 
!               sine/cosine/cosine/sine transform
! PPFSCSCT2RM4Y performs y part of 4 component vector 2d real
!               sine/cosine/sine/cosine transform
! WPPFSCT2RM22 wrapper function for 2 component tensor 2d real
!              sine/cosine transforms
! PPFSCCST2RM22X performs x part of 2 component vector 2d real 
!                sine/cosine/cosine/sine transform
! PPFSCSCT2RM22Y performs y part of 2 component vector 2d real
!                sine/cosine/sine/cosine transform
! WPPFSST2RM23 wrapper function for 3 component vector 2d real
!             sine/cosine transforms
! PPFSSCT2RM23X performs x part of 3 component vector 2d real
!               sine/sine/cosine transform
! PPFSSCT2RM23Y performs y part of 3 component vector 2d real
!               sine/sine/cosine transform
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: april 2, 2017
!-----------------------------------------------------------------------
      subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
! this subroutine calculates tables needed by a two dimensional
! fast real sine and cosine transforms and their inverses.
! input: indx, indy, nxhyd, nxyd
! output: mixup, sctd
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyd
      integer mixup
      complex sctd
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, ny, nxy, nxhy
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
! sine/cosine table for the angles n*pi/nxy
      dnxy = 0.5*6.28318530717959/real(nxy)
      do 30 j = 1, nxy
      arg = dnxy*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real sine/sine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real sine/cosine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real cosine/sine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real cosine/cosine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of a two dimensional fast real
! sine transform and its inverse, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse sine transform is performed
! f(n,k) = (1/nx*ny)*sum(f(j,k)*sin(pi*n*j/nx))
! if isign = 1, a forward sine transform is performed
! f(j,k) = sum(f(n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp, nxvh, kypd
      integer nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb
      real at1, at2, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,i)
      at1 = f(j,i) + at2
      at2 = f(j,i) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = 0.5*at2
      f(j,i) = at1 + at2
      f(nx+2-j,i) = at1 - at2
   10 continue
      f(1,i) = 0.0
      f(nxh+1,i) = 2.0*f(nxh+1,i)
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t2 = f(2*j1-1,i)
         t3 = f(2*j1,i)
         f(2*j1-1,i) = f(2*j-1,i)
         f(2*j1,i) = f(2*j,i)
         f(2*j-1,i) = t2
         f(2*j,i) = t3
      endif
   20 continue
! then transform in x
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 60 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = ani*(t2 + t4)
         f(2*j,i) = ani*(t3 + t5)
         f(nx3-2*j,i) = ani*(t2 - t4)
         f(nx3-2*j+1,i) = ani*(t5 - t3)
   60    continue
         f(nxh+1,i) = 2.0*ani*f(nxh+1,i)
         f(nxh+2,i) = -2.0*ani*f(nxh+2,i)
         t2 = 2.0*ani*(f(1,i) + f(2,i))
         f(2,i) = 2.0*ani*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*ani*f(nx+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 70 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = t2 + t4
         f(2*j,i) = t3 + t5
         f(nx3-2*j,i) = t2 - t4
         f(nx3-2*j+1,i) = t5 - t3
   70    continue
         f(nxh+1,i) = 2.0*f(nxh+1,i)
         f(nxh+2,i) = -2.0*f(nxh+2,i)
         t2 = 2.0*(f(1,i) + f(2,i))
         f(2,i) = 2.0*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*f(nx+1,i)
      endif
! perform recursion for sine transform
      sum1 = 0.5*f(1,i)
      f(1,i) = 0.0
      f(2,i) = sum1
      do 80 j = 2, nxh
      sum1 = sum1 + f(2*j-1,i)
      f(2*j-1,i) = -f(2*j,i)
      f(2*j,i) = sum1
   80 continue
      f(nx+1,i) = 0.0
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of a two dimensional fast real
! cosine transform and its inverse, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse cosine transform is performed
! f(n,k) = (1/nx*ny)*(.5*f(1,k) + ((-1)**n)*f(nx+1,k)
!            + sum(f(j,k)*cos(pi*n*j/nx)))
! if isign = 1, a forward cosine transform is performed
! f(j,k) = 2*(.5*f(1,k) + ((-1)**j)*f(n+1,k) + sum(f(n,k)*
!            cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp+1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp, nxvh, kypd
      integer nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb
      real at1, at2, t2, t3, t4, t5, t6, ani
      double precision sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,i) - f(nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,i)
      at1 = f(j,i) + at2
      at2 = f(j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = 0.5*at1
      f(j,i) = at1 - at2
      f(nx+2-j,i) = at1 + at2
   10 continue
      f(1,i) = 0.5*(f(1,i) + f(nx+1,i))
      f(nx+1,i) = sum1
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t2 = f(2*j1-1,i)
         t3 = f(2*j1,i)
         f(2*j1-1,i) = f(2*j-1,i)
         f(2*j1,i) = f(2*j,i)
         f(2*j-1,i) = t2
         f(2*j,i) = t3
      endif
   20 continue
! then transform in x
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 60 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = ani*(t2 + t4)
         f(2*j,i) = ani*(t3 + t5)
         f(nx3-2*j,i) = ani*(t2 - t4)
         f(nx3-2*j+1,i) = ani*(t5 - t3)
   60    continue
         f(nxh+1,i) = 2.0*ani*f(nxh+1,i)
         f(nxh+2,i) = -2.0*ani*f(nxh+2,i)
         t2 = 2.0*ani*(f(1,i) + f(2,i))
         f(2,i) = 2.0*ani*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*ani*f(nx+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 70 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = t2 + t4
         f(2*j,i) = t3 + t5
         f(nx3-2*j,i) = t2 - t4
         f(nx3-2*j+1,i) = t5 - t3
   70    continue
         f(nxh+1,i) = 2.0*f(nxh+1,i)
         f(nxh+2,i) = -2.0*f(nxh+2,i)
         t2 = 2.0*(f(1,i) + f(2,i))
         f(2,i) = 2.0*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*f(nx+1,i)
      endif
! perform recursion for cosine transform
      sum1 = f(nx+1,i)
      f(nx+1,i) = f(2,i)
      f(2,i) = sum1
      do 80 j = 2, nxh
      sum1 = sum1 - f(2*j,i)
      f(2*j,i) = sum1
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of a two dimensional fast real
! sine transform and its inverse, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse sine transform is performed
! g(m,n) = sum(g(k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward sine transform is performed
! g(k,n) = sum(g(m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxpd
      integer nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb
      real at1, at2, t2, t3, t4, t5, t6
      complex t1
      double precision sum1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,i)
      at1 = g(k,i) + at2
      at2 = g(k,i) - at2
      at1 = -aimag(sctd(k1))*at1
      at2 = 0.5*at2
      g(k,i) = at1 + at2
      g(ny+2-k,i) = at1 - at2
   10 continue
      g(1,i) = 0.0
      g(nyh+1,i) = 2.0*g(nyh+1,i)
! bit-reverse array elements in y
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t2 = g(2*k1-1,i)
         t3 = g(2*k1,i)
         g(2*k1-1,i) = g(2*k-1,i)
         g(2*k1,i) = g(2*k,i)
         g(2*k-1,i) = t2
         g(2*k,i) = t3
      endif
   20 continue
! then transform in y
      do 50 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*g(2*j2-1,i) - aimag(t1)*g(2*j2,i)
      t3 = aimag(t1)*g(2*j2-1,i) + real(t1)*g(2*j2,i)
      g(2*j2-1,i) = g(2*j1-1,i) - t2
      g(2*j2,i) = g(2*j1,i) - t3
      g(2*j1-1,i) = g(2*j1-1,i) + t2
      g(2*j1,i) = g(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 60 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = 0.5*(t2 + t4)
         g(2*k,i) = 0.5*(t3 + t5)
         g(ny3-2*k,i) = 0.5*(t2 - t4)
         g(ny3-2*k+1,i) = 0.5*(t5 - t3)
   60    continue
         g(nyh+1,i) = g(nyh+1,i)
         g(nyh+2,i) = -g(nyh+2,i)
         t2 = g(1,i) + g(2,i)
         g(2,i) = g(1,i) - g(2,i)
         g(1,i) = t2
         g(ny+1,i) = g(ny+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 70 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = t2 + t4
         g(2*k,i) = t3 + t5
         g(ny3-2*k,i) = t2 - t4
         g(ny3-2*k+1,i) = t5 - t3
   70    continue
         g(nyh+1,i) = 2.0*g(nyh+1,i)
         g(nyh+2,i) = -2.0*g(nyh+2,i)
         t2 = 2.0*(g(1,i) + g(2,i))
         g(2,i) = 2.0*(g(1,i) - g(2,i))
         g(1,i) = t2
         g(ny+1,i) = 2.0*g(ny+1,i)
      endif
! perform recursion for sine transform
      sum1 = 0.5*g(1,i)
      g(1,i) = 0.0
      g(2,i) = sum1
      do 80 k = 2, nyh
      sum1 = sum1 + g(2*k-1,i)
      g(2*k-1,i) = -g(2*k,i)
      g(2*k,i) = sum1
   80 continue
      g(ny+1,i) = 0.0
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of a two dimensional fast real
! cosine transform and its inverse, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse cosine transform is performed
! g(m,n) = (.5*g(1,n) + ((-1)**m)*g(ny+1,n)
!            + sum(g(k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward cosine transform is performed
! g(k,n) = 2*(.5*g(1,n) + ((-1)**m)*g(ny+1,n) + sum(g(m,n)*
!            cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxpd
      integer nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb
      real at1, at2, t2, t3, t4, t5, t6
      complex t1
      double precision sum1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,i) - g(ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,i)
      at1 = g(k,i) + at2
      at2 = g(k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = -aimag(sctd(k1))*at2
      at1 = 0.5*at1
      g(k,i) = at1 - at2
      g(ny+2-k,i) = at1 + at2
   10 continue
      g(1,i) = 0.5*(g(1,i) + g(ny+1,i))
      g(ny+1,i) = sum1
! bit-reverse array elements in y
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t2 = g(2*k1-1,i)
         t3 = g(2*k1,i)
         g(2*k1-1,i) = g(2*k-1,i)
         g(2*k1,i) = g(2*k,i)
         g(2*k-1,i) = t2
         g(2*k,i) = t3
      endif
   20 continue
! then transform in y
      do 50 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*g(2*j2-1,i) - aimag(t1)*g(2*j2,i)
      t3 = aimag(t1)*g(2*j2-1,i) + real(t1)*g(2*j2,i)
      g(2*j2-1,i) = g(2*j1-1,i) - t2
      g(2*j2,i) = g(2*j1,i) - t3
      g(2*j1-1,i) = g(2*j1-1,i) + t2
      g(2*j1,i) = g(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 60 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = 0.5*(t2 + t4)
         g(2*k,i) = 0.5*(t3 + t5)
         g(ny3-2*k,i) = 0.5*(t2 - t4)
         g(ny3-2*k+1,i) = 0.5*(t5 - t3)
   60    continue
         g(nyh+1,i) = g(nyh+1,i)
         g(nyh+2,i) = -g(nyh+2,i)
         t2 = g(1,i) + g(2,i)
         g(2,i) = g(1,i) - g(2,i)
         g(1,i) = t2
         g(ny+1,i) = g(ny+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 70 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = t2 + t4
         g(2*k,i) = t3 + t5
         g(ny3-2*k,i) = t2 - t4
         g(ny3-2*k+1,i) = t5 - t3
   70    continue
         g(nyh+1,i) = 2.0*g(nyh+1,i)
         g(nyh+2,i) = -2.0*g(nyh+2,i)
         t2 = 2.0*(g(1,i) + g(2,i))
         g(2,i) = 2.0*(g(1,i) - g(2,i))
         g(1,i) = t2
         g(ny+1,i) = 2.0*g(ny+1,i)
      endif
! perform recursion for cosine transform
      sum1 = g(ny+1,i)
      g(ny+1,i) = g(2,i)
      g(2,i) = sum1
      do 80 k = 2, nyh
      sum1 = sum1 - g(2*k,i)
      g(2*k,i) = sum1
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real cosine/sine transforms
! for the electric field with dirichlet or magnetic field with neumann
! boundary conditions
! x component has a cosine/sine transform in x and y, respectively
! y component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine-sine transform
         call PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine transform
         call PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine-sine transform
         call PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real sine/cosine transforms
! for the magnetic field with dirichlet or electric field with neumann
! boundary conditions
! x component has a sine/cosine transform in x and y, respectively
! y component has a cosine/sine transform in x and y, respectively
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine transform
         call PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine-sine transform
         call PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine transform
         call PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine-cosine transform
         call PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine transforms are performed
! f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! f(2,n,k) = (1/nx*ny)*sum(f(2,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward cosine-sine transforms are performed
! f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! f(2,j,k) = sum(f(2,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = second dimension of f >= nx/2 + 1
! kypd = third dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,1,i) - f(1,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(1,j,i) = at1 - at2
      f(1,nx+2-j,i) = at1 + at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.5*(f(1,1,i) + f(1,nx+1,i))
      f(1,nx+1,i) = sum1
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = f(1,nx+1,i)
      f(1,nx+1,i) = f(1,2,i)
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 - f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(2,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transforms are performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
!              + sum(f(2,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
!              + sum(f(2,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = second dimension of f >= nx/2 + 1
! kypd = third dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum1
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward sine-cosine transforms are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = second dimension of g >= ny + 1
! kxpd = third dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum1
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(1,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine transform are performed
! g(1,m,n) = (.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,k,n)*cos(pi*m*k/ny))
! g(2,m,n) = sum(g(2,k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward cosine-sine transforms are performed
! g(1,k,n) = 2*(.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,m,n)*cos(pi*m*k/ny))
! g(2,k,n) = sum(g(2,m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = second dimension of g >= ny + 1
! kxpd = third dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,1,i) - g(1,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(1,k,i) = at1 - at2
      g(1,ny+2-k,i) = at1 + at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
   10 continue
      g(1,1,i) = 0.5*(g(1,1,i) + g(1,ny+1,i))
      g(1,ny+1,i) = sum1
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = g(1,ny+1,i)
      g(1,ny+1,i) = g(1,2,i)
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 - g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real cosine/sine transforms
! for the electric field with dirichlet or magnetic field with neumann
! boundary conditions
! x component has a cosine/sine transform in x and y, respectively
! y/z component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine-sine-sine transform
         call PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine transform
         call PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine-sine transform
         call PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine-sine-sine transform
         call PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real sine/cosine transforms
! for the magnetic field with dirichlet or electric field with neumann
! boundary conditions
! x component has a sine/cosine transform in x and y, respectively
! y component has a cosine/sine transform in x and y, respectively
! z component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine-cosine transform
         call PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine-sine-cosine transform
         call PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine-cosine transform
         call PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine-cosine-cosine transform
         call PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y/z component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine-sine transforms are performed
! f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! f(2:3,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward cosine-sine-sine transforms are performed
! f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! f(2:3,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,1,i) - f(1,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(1,j,i) = at1 - at2
      f(1,nx+2-j,i) = at1 + at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(3,j,i) = at1 + at2
      f(3,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.5*(f(1,1,i) + f(1,nx+1,i))
      f(1,nx+1,i) = sum1
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
      f(3,1,i) = 0.0
      f(3,nxh+1,i) = 2.0*f(3,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = f(1,nx+1,i)
      f(1,nx+1,i) = f(1,2,i)
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      sum3 = 0.5*f(3,1,i)
      f(3,1,i) = 0.0
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 - f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 + f(3,2*j-1,i)
      f(3,2*j-1,i) = -f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(2,nx+1,i) = 0.0
      f(3,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y/z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-cosine transforms are performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
!              + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      sum2 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum1
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum2
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-sine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! g(3,m,n) = sum(g(3,k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward sine-cosine-sine transforms are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! g(3,k,n) = sum(g(3,m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(3,k,i) = at1 + at2
      g(3,ny+2-k,i) = at1 - at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum1
      g(3,1,i) = 0.0
      g(3,nyh+1,i) = 2.0*g(3,nyh+1,i)
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      sum3 = 0.5*g(3,1,i)
      g(3,1,i) = 0.0
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 + g(3,2*k-1,i)
      g(3,2*k-1,i) = -g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(1,ny+1,i) = 0.0
      g(3,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine-cosine transform are performed
! g(1,m,n) = (.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,k,n)*cos(pi*m*k/ny))
! g(2,m,n) = sum(g(2,k,n)*sin(pi*m*k/ny))
! g(3,m,n) = (.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward cosine-sine-cosine transforms are performed
! g(1,k,n) = 2*(.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,m,n)*cos(pi*m*k/ny))
! g(2,k,n) = sum(g(2,m,n)*sin(pi*m*k/ny))
! g(3,k,n) = 2*(.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,1,i) - g(1,ny+1,i))
      sum2 = 0.5*(g(3,1,i) - g(3,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(1,k,i) = at1 - at2
      g(1,ny+2-k,i) = at1 + at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(3,k,i) = at1 - at2
      g(3,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.5*(g(1,1,i) + g(1,ny+1,i))
      g(1,ny+1,i) = sum1
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
      g(3,1,i) = 0.5*(g(3,1,i) + g(3,ny+1,i))
      g(3,ny+1,i) = sum2
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = g(1,ny+1,i)
      g(1,ny+1,i) = g(1,2,i)
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      sum3 = g(3,ny+1,i)
      g(3,ny+1,i) = g(3,2,i)
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 - g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 - g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM4(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 4 parallel real sine/cosine transforms
! for the momentum flux with dirichlet boundary conditions
! x component has a sine/sine transform in x and y, respectively
! y component has a cosine/cosine transform in x and y, respectively
! z component has a cosine/sine transform in x and y, respectively
! w component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(4,2*nxvh,kypd), g(4,nyv,kxp2d)
      dimension bs(4,kxp2+1,kyp+1), br(4,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine-cosine-sine transform
         call PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,4,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine-cosine transform
         call PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,4,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,4,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine-sine-cosine transform
         call PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,4,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine-cosine transform
         call PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 4 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/w component has a sine transform, y/z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-cosine-sine transforms are
! performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
!              + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
! f(4,n,k) = (1/nx*ny)*sum(f(4,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward sine-cosine-cosine-sine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! f(4,j,k) = sum(f(4,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(4,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3, sum4
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3,sum4)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum2 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      sum3 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum3 = sum3 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
      at2 = f(4,nx+2-j,i)
      at1 = f(4,j,i) + at2
      at2 = f(4,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(4,j,i) = at1 + at2
      f(4,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum2
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum3
      f(4,1,i) = 0.0
      f(4,nxh+1,i) = 2.0*f(4,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 4
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 4
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 4
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 4
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 4
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 4
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      sum4 = 0.5*f(4,1,i)
      f(4,1,i) = 0.0
      f(4,2,i) = sum4
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
      sum4 = sum4 + f(4,2*j-1,i)
      f(4,2*j-1,i) = -f(4,2*j,i)
      f(4,2*j,i) = sum4
  140 continue
      f(1,nx+1,i) = 0.0
      f(4,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 4 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a sine transform, y/w component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! g(3,m,n) = sum(g(3,k,n)*sin(pi*m*k/ny))
! g(4,m,n) = (.5*g(4,1,n) + ((-1)**m)*g(4,ny+1,n)
!              + sum(g(4,k,n)*cos(pi*m*k/ny))
! if isign = 1, forward sine-cosine-sine-cosine transform are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! g(3,k,n) = sum(g(3,m,n)*sin(pi*m*k/ny))
! g(4,k,n) = 2*(.5*g(4,1,n) + ((-1)**m)*g(4,ny+1,n)
!              + sum(g(4,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(4,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3, sum4
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3,sum4)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum2 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      sum4 = 0.5*(g(4,1,i) - g(4,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(3,k,i) = at1 + at2
      g(3,ny+2-k,i) = at1 - at2
      at2 = g(4,ny+2-k,i)
      at1 = g(4,k,i) + at2
      at2 = g(4,k,i) - at2
      sum4 = sum4 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(4,k,i) = at1 - at2
      g(4,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum2
      g(3,1,i) = 0.0
      g(3,nyh+1,i) = 2.0*g(3,nyh+1,i)
      g(4,1,i) = 0.5*(g(4,1,i) + g(4,ny+1,i))
      g(4,ny+1,i) = sum4
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 4
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 4
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 4
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 4
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 4
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 4
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      sum3 = 0.5*g(3,1,i)
      g(3,1,i) = 0.0
      g(3,2,i) = sum3
      sum4 = g(4,ny+1,i)
      g(4,ny+1,i) = g(4,2,i)
      g(4,2,i) = sum4
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 + g(3,2*k-1,i)
      g(3,2*k-1,i) = -g(3,2*k,i)
      g(3,2*k,i) = sum3
      sum4 = sum4 - g(4,2*k,i)
      g(4,2*k,i) = sum4
  140 continue
      g(1,ny+1,i) = 0.0
      g(3,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM22(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real sine/cosine transforms
! for the momentum flux with dirichlet boundary conditions
! x component has a sine/sine transform in x and y, respectively
! y component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine transform
         call PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine transform
         call PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transforms are
! performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
!              + sum(f(2,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum2 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum2
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! if isign = 1, forward sine-cosine transform are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum2 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum2
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(1,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real sine transforms
! x/y component has a sine/sine transform in x and y, respectively
! z component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transforms
         call PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transforms
         call PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transforms
         call PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transforms
         call PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/y component has a sine transform, z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-sine-cosine transforms are performed
! f(1:2,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
! f(3,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-sine-cosine transforms are performed
! f(1:2,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
! f(3,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum3 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum3 = sum3 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum3
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for sine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(1,nx+1,i) = 0.0
      f(2,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/y component has a sine transform, z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-sine-cosine transform are performed
! g(1:2,m,n) = sum(g(1:2,k,n)*sin(pi*m*k/ny))
! g(3,m,n) = (.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward sine-sine-cosine transforms are performed
! g(1:2,k,n) = sum(g(1:2,m,n)*sin(pi*m*k/ny))
! g(3,k,n) = 2*(.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum3 = 0.5*(g(3,1,i) - g(3,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      sum3 = sum3 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(3,k,i) = at1 - at2
      g(3,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
      g(3,1,i) = 0.5*(g(3,1,i) + g(3,ny+1,i))
      g(3,ny+1,i) = sum3
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      sum3 = g(3,ny+1,i)
      g(3,ny+1,i) = g(3,2,i)
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 - g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(1,ny+1,i) = 0.0
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
