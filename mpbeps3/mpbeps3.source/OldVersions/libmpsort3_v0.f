!-----------------------------------------------------------------------
! Fortran Library for reordering particles
! 3D MPI/OpenMP PIC Codes:
! PPPORDER32LA performs first part of particle reordering into tiles,
!              creates list of particles which are leaving tile, and
!              buffers outgoing particles
! PPPORDERF32LA performs first part of particle reordering into tiles,
!              buffers outgoing particles
! PPPORDER2LB performs second part of particle reordering into tiles,
!             copies buffers into particle array
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: february 11, 2016
!-----------------------------------------------------------------------
      subroutine PPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,  &
     &ncll,nclr,noff,nyzp,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,  &
     &mxzyp1,npbmx,ntmax,nbmax,idds,irc)
! this subroutine performs first part of a particle sort by x,y,z grid
! in tiles of mx, my, mz
! linear interpolation, with periodic boundary conditions
! for distributed data, with 2d domain decomposition in y/z.
! tiles are assumed to be arranged in 3D linear memory
! this part of the algorithm has 3 steps.  first, one finds particles
! leaving tile and stores their number in each directon, location, and
! destination in ncl and ihole.  then, a prefix scan of ncl is performed
! and departing particles are buffered in ppbuff in direction order.
! finally, we buffer particles leaving the processor in y/z direction in
! sbufl and sbufr, and store particle number offsets in ncll and nclr.
! input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = position z of particle n in tile m
! ppbuff(i,n,l) = i co-ordinate of particle n in tile l
! sbufl = buffer for particles being sent to lower/back processor
! sbufr = buffer for particles being sent to upper/forward processor
! kpic(l) = number of particles in tile l
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower/back processor
! nclr = number offset being sent to upper/forward processor
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! nx/ny/nz = system length in x/y/z direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! mxzyp1 = mx1*max(myp1,mzp1)
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! idds = dimensionality of domain decomposition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, myp1, mzp1
      integer mxzyp1, npbmx, ntmax, nbmax, idds, irc
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, noff, nyzp, ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension sbufl(idimp,nbmax,2), sbufr(idimp,nbmax,2)
      dimension kpic(mx1*myp1*mzp1), noff(idds), nyzp(idds)
      dimension ncl(26,mx1*myp1*mzp1), ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension ncll(3,mxzyp1,3,2), nclr(3,mxzyp1,3,2)
! local data
      integer mxyp1, mxzp1, mxyzp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, n, ii, ih, nh, ist, nn, mm, ll, isum
      integer ip, j1, j2, k1, k2, kk, nr, nl
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dx, dy, dz
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,ih,nh,ist,dx,dy,dz,
!$OMP& edgelx,edgely,edgelz,edgerx,edgery,edgerz)
      do 30 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      ih = 0
      nh = 0
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
! clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
   10 continue
! loop over particles in tile
      do 20 j = 1, nppp
      dx = ppart(1,j,l)
      dy = ppart(2,j,l)
      dz = ppart(3,j,l)
! find particles going out of bounds
      ist = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,l) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,l) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,l) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,l) = dy
         else
            ist = ist + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) ppart(3,j,l) = dz - anz
         ist = ist + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               ist = ist + 9
            else
               dz = 0.0
            endif
            ppart(3,j,l) = dz
         else
            ist = ist + 9
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,l) = ncl(ist,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = ist
         else
            nh = 1
         endif
      endif
   20 continue
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) return
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 70 l = 1, mxyzp1
! find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,l)
      ip = 0
! loop over particles leaving tile
      do 60 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,l) = ppart(i,j1,l)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   60 continue
! set error
      if (ip.gt.0) irc = ncl(26,l)
   70 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node up or down:
! update sbufl(:,1), sbufr(:,1), ncll(:,1), nclr(:,1)
      mxzp1 = mx1*mzp1
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,l,ll)
      do 80 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
! going straight up and down
      ncll(1,ll,1,1) = ncl(5,k) - ncl(2,k)
      nclr(1,ll,1,1) = ncl(8,k+kk) - ncl(5,k+kk)
! particles going up and back
      ncll(1,ll,2,1) = ncl(14,k) - ncl(11,k)
! particles going down and back
      nclr(1,ll,2,1) = ncl(17,k+kk) - ncl(14,k+kk)
! particles going up and forward
      ncll(1,ll,3,1) = ncl(23,k) - ncl(20,k)
! particles going down and forward
      nclr(1,ll,3,1) = ncl(26,k+kk) - ncl(23,k+kk)
   80 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   90 if (kk.ge.mxzp1) go to 110
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 100 ll = 1, mxzp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxzp1) then
         ncll(1,nn,1,1) = ncll(1,nn,1,1) + ncll(1,mm+1,1,1)
         nclr(1,nn,1,1) = nclr(1,nn,1,1) + nclr(1,mm+1,1,1)
         ncll(1,nn,2,1) = ncll(1,nn,2,1) + ncll(1,mm+1,2,1)
         nclr(1,nn,2,1) = nclr(1,nn,2,1) + nclr(1,mm+1,2,1)
         ncll(1,nn,3,1) = ncll(1,nn,3,1) + ncll(1,mm+1,3,1)
         nclr(1,nn,3,1) = nclr(1,nn,3,1) + nclr(1,mm+1,3,1)
      endif
  100 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 90
  110 kk = mx1*(myp1 - 1)
      j1 = ncll(1,mxzp1,1,1)
      k1 = j1 + ncll(1,mxzp1,2,1)
      j2 = nclr(1,mxzp1,1,1)
      k2 = j2 + nclr(1,mxzp1,2,1)
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,ii,nn,mm)
      do 300 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
! particles going straight up
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,ll,1,1) - ii
      do 130 j = 1, min(ii,nbmax-nn)
      do 120 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(2,k),k)
  120 continue
  130 continue
      do 140 i = 1, 3
      ncll(i,ll,1,1) = ncl(i+2,k) - ncl(2,k) + nn
  140 continue
! particles going straight down
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,ll,1,1) - ii
      do 160 j = 1, min(ii,nbmax-mm)
      do 150 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  150 continue
  160 continue
      do 170 i = 1, 3
      nclr(i,ll,1,1) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  170 continue
! particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,1) - ii
      do 190 j = 1, min(ii,nbmax-nn)
      do 180 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(11,k),k)
  180 continue
  190 continue
      do 200 i = 1, 3
      ncll(i,ll,2,1) = ncl(i+11,k) - ncl(11,k) + nn
  200 continue
! particles going down and back
      ii = ncl(17,k+kk) - ncl(14,k+kk)
      mm = j2 + nclr(1,ll,2,1) - ii
      do 220 j = 1, min(ii,nbmax-mm)
      do 210 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(14,k+kk),k+kk)
  210 continue
  220 continue
      do 230 i = 1, 3
      nclr(i,ll,2,1) = ncl(i+14,k+kk) - ncl(14,k+kk) + mm
  230 continue
! particles going up and forward
      ii = ncl(23,k) - ncl(20,k)
      nn = k1 + ncll(1,ll,3,1) - ii
      do 250 j = 1, min(ii,nbmax-nn)
      do 240 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(20,k),k)
  240 continue
  250 continue
      do 260 i = 1, 3
      ncll(i,ll,3,1) = ncl(i+20,k) - ncl(20,k) + nn
  260 continue
! particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,1) - ii
      do 280 j = 1, min(ii,nbmax-mm)
      do 270 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  270 continue
  280 continue
      do 290 i = 1, 3
      nclr(i,ll,3,1) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  290 continue
  300 continue
!$OMP END PARALLEL DO
!
! buffer particles and their number leaving the node back or forward:
! update sbufl(:,2), sbufr(:,2), ncll(:,2), nclr(:,2)
      kk = mxyp1*(mzp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,ll)
      do 310 ll = 1, mxyp1
      k = ll
! going straight back or forward
      ncll(1,ll,1,2) = ncl(11,k) - ncl(8,k)
      nclr(1,ll,1,2) = ncl(20,k+kk) - ncl(17,k+kk)
! particles going back and up
      ncll(1,ll,2,2) = ncl(14,k) - ncl(11,k)
! particles going forward and up
      nclr(1,ll,2,2) = ncl(23,k+kk) - ncl(20,k+kk)
! particles going back and down
      ncll(1,ll,3,2) = ncl(17,k) - ncl(14,k)
! particles going forward and down
      nclr(1,ll,3,2) = ncl(26,k+kk) - ncl(23,k+kk)
  310 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
  320 if (kk.ge.mxyp1) go to 340
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 330 ll = 1, mxyp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxyp1) then
         ncll(1,nn,1,2) = ncll(1,nn,1,2) + ncll(1,mm+1,1,2)
         nclr(1,nn,1,2) = nclr(1,nn,1,2) + nclr(1,mm+1,1,2)
         ncll(1,nn,2,2) = ncll(1,nn,2,2) + ncll(1,mm+1,2,2)
         nclr(1,nn,2,2) = nclr(1,nn,2,2) + nclr(1,mm+1,2,2)
         ncll(1,nn,3,2) = ncll(1,nn,3,2) + ncll(1,mm+1,3,2)
         nclr(1,nn,3,2) = nclr(1,nn,3,2) + nclr(1,mm+1,3,2)
      endif
  330 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 320
  340 kk = mxyp1*(mzp1 - 1)
      j1 = ncll(1,mxyp1,1,2)
      k1 = j1 + ncll(1,mxyp1,2,2)
      j2 = nclr(1,mxyp1,1,2)
      k2 = j2 + nclr(1,mxyp1,2,2)
!$OMP PARALLEL DO PRIVATE(i,j,k,ll,ii,nn,mm)
      do 530 ll = 1, mxyp1
      k = ll
! particles going straight up
      ii = ncl(11,k) - ncl(8,k)
      nn = ncll(1,ll,1,2) - ii
      do 360 j = 1, min(ii,nbmax-nn)
      do 350 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(8,k),k)
  350 continue
  360 continue
      do 370 i = 1, 3
      ncll(i,ll,1,2) = ncl(i+8,k) - ncl(8,k) + nn
  370 continue
! particles going straight down
      ii = ncl(20,k+kk) - ncl(17,k+kk)
      mm = nclr(1,ll,1,2) - ii
      do 390 j = 1, min(ii,nbmax-mm)
      do 380 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(17,k+kk),k+kk)
  380 continue
  390 continue
      do 400 i = 1, 3
      nclr(i,ll,1,2) = ncl(i+17,k+kk) - ncl(17,k+kk) + mm
  400 continue
! particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,2) - ii
      do 420 j = 1, min(ii,nbmax-nn)
      do 410 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(11,k),k)
  410 continue
  420 continue
      do 430 i = 1, 3
      ncll(i,ll,2,2) = ncl(i+11,k) - ncl(11,k) + nn
  430 continue
! particles going down and back
      ii = ncl(23,k+kk) - ncl(20,k+kk)
      mm = j2 + nclr(1,ll,2,2) - ii
      do 450 j = 1, min(ii,nbmax-mm)
      do 440 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(20,k+kk),k+kk)
  440 continue
  450 continue
      do 460 i = 1, 3
      nclr(i,ll,2,2) = ncl(i+20,k+kk) - ncl(20,k+kk) + mm
  460 continue
! particles going up and forward
      ii = ncl(17,k) - ncl(14,k)
      nn = k1 + ncll(1,ll,3,2) - ii
      do 480 j = 1, min(ii,nbmax-nn)
      do 470 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(14,k),k)
  470 continue
  480 continue
      do 490 i = 1, 3
      ncll(i,ll,3,2) = ncl(i+14,k) - ncl(14,k) + nn
  490 continue
! particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,2) - ii
      do 510 j = 1, min(ii,nbmax-mm)
      do 500 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  500 continue
  510 continue
      do 520 i = 1, 3
      nclr(i,ll,3,2) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  520 continue
  530 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      kk = max(ncll(3,mxzp1,3,1),nclr(3,mxzp1,3,1))
      ll = max(ncll(3,mxyp1,3,2),nclr(3,mxyp1,3,2))
      ii = max(kk,ll)
! corners overflow
      nn = nclr(3,mx1,2,1) - nclr(3,mxzp1,1,1)
      mm = ncll(3,mx1,2,1) - ncll(3,mxzp1,1,1)
      n = mx1*(mzp1 - 1)
      if (n.gt.0) then
         nr = nclr(3,n,3,1)
         nl = ncll(3,n,3,1)
      else
         nr = nclr(3,mxzp1,2,1)
         nl = ncll(3,mxzp1,2,1)
      endif
      kk = nclr(3,n+mx1,3,1) - nr
      ll = ncll(3,n+mx1,3,1) - nl
! total overflow: result valid only for one processor case
      ii = ii + max(nn+kk,mm+ll)
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll, &
     &nclr,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,irc)
! this subroutine performs first part of a particle sort by x,y,z grid
! in tiles of mx, my, mz
! linear interpolation, with periodic boundary conditions
! for distributed data, with 2d domain decomposition in y/z.
! tiles are assumed to be arranged in 3D linear memory
! this part of the algorithm has 2 steps.  first, a prefix scan of ncl
! is performed and departing particles are buffered in ppbuff in
! direction order. then, we buffer particles leaving the processor in
! sbufl and sbufr, and store particle number offsets in ncll and nclr.
! it assumes that the number, location, and destination of particles 
! leaving a tile have been previously stored in ncl and ihole by the
! PPGPPUSHF32L subroutine.
! input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = position z of particle n in tile m
! ppbuff(i,n,l) = i co-ordinate of particle n in tile l
! sbufl = buffer for particles being sent to lower/back processor
! sbufr = buffer for particles being sent to upper/forward processor
! ncl(i,l) = number of particles going to destination i, tile l
! ihole(1,:,l) = location of hole in array left by departing particle
! ihole(2,:,l) = direction destination of particle leaving hole
! all for tile l
! ihole(1,1,l) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower/back processor
! nclr = number offset being sent to upper/forward processor
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! mxzyp1 = mx1*max(myp1,mzp1)
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, mzp1, mxzyp1, npbmx, ntmax, nbmax
      integer irc
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension sbufl(idimp,nbmax,2), sbufr(idimp,nbmax,2)
      dimension ncl(26,mx1*myp1*mzp1), ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension ncll(3,mxzyp1,3,2), nclr(3,mxzyp1,3,2)
! local data
      integer mxyp1, mxzp1, mxyzp1
      integer i, j, k, l, n, ii, nh, ist, nn, mm, ll, isum
      integer ip, j1, j2, k1, k2, kk, nr, nl
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 40 l = 1, mxyzp1
! find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,l)
      ip = 0
! loop over particles leaving tile
      do 30 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,l) = ppart(i,j1,l)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   30 continue
! set error
      if (ip.gt.0) irc = ncl(26,l)
   40 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node up or down:
! update sbufl(:,1), sbufr(:,1), ncll(:,1), nclr(:,1)
      mxzp1 = mx1*mzp1
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,l,ll)
      do 50 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
! going straight up and down
      ncll(1,ll,1,1) = ncl(5,k) - ncl(2,k)
      nclr(1,ll,1,1) = ncl(8,k+kk) - ncl(5,k+kk)
! particles going up and back
      ncll(1,ll,2,1) = ncl(14,k) - ncl(11,k)
! particles going down and back
      nclr(1,ll,2,1) = ncl(17,k+kk) - ncl(14,k+kk)
! particles going up and forward
      ncll(1,ll,3,1) = ncl(23,k) - ncl(20,k)
! particles going down and forward
      nclr(1,ll,3,1) = ncl(26,k+kk) - ncl(23,k+kk)
   50 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   60 if (kk.ge.mxzp1) go to 80
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 70 ll = 1, mxzp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxzp1) then
         ncll(1,nn,1,1) = ncll(1,nn,1,1) + ncll(1,mm+1,1,1)
         nclr(1,nn,1,1) = nclr(1,nn,1,1) + nclr(1,mm+1,1,1)
         ncll(1,nn,2,1) = ncll(1,nn,2,1) + ncll(1,mm+1,2,1)
         nclr(1,nn,2,1) = nclr(1,nn,2,1) + nclr(1,mm+1,2,1)
         ncll(1,nn,3,1) = ncll(1,nn,3,1) + ncll(1,mm+1,3,1)
         nclr(1,nn,3,1) = nclr(1,nn,3,1) + nclr(1,mm+1,3,1)
      endif
   70 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 60
   80 kk = mx1*(myp1 - 1)
      j1 = ncll(1,mxzp1,1,1)
      k1 = j1 + ncll(1,mxzp1,2,1)
      j2 = nclr(1,mxzp1,1,1)
      k2 = j2 + nclr(1,mxzp1,2,1)
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,ii,nn,mm)
      do 270 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
! particles going straight up
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,ll,1,1) - ii
      do 100 j = 1, min(ii,nbmax-nn)
      do 90 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(2,k),k)
   90 continue
  100 continue
      do 110 i = 1, 3
      ncll(i,ll,1,1) = ncl(i+2,k) - ncl(2,k) + nn
  110 continue
! particles going straight down
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,ll,1,1) - ii
      do 130 j = 1, min(ii,nbmax-mm)
      do 120 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  120 continue
  130 continue
      do 140 i = 1, 3
      nclr(i,ll,1,1) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  140 continue
! particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,1) - ii
      do 160 j = 1, min(ii,nbmax-nn)
      do 150 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(11,k),k)
  150 continue
  160 continue
      do 170 i = 1, 3
      ncll(i,ll,2,1) = ncl(i+11,k) - ncl(11,k) + nn
  170 continue
! particles going down and back
      ii = ncl(17,k+kk) - ncl(14,k+kk)
      mm = j2 + nclr(1,ll,2,1) - ii
      do 190 j = 1, min(ii,nbmax-mm)
      do 180 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(14,k+kk),k+kk)
  180 continue
  190 continue
      do 200 i = 1, 3
      nclr(i,ll,2,1) = ncl(i+14,k+kk) - ncl(14,k+kk) + mm
  200 continue
! particles going up and forward
      ii = ncl(23,k) - ncl(20,k)
      nn = k1 + ncll(1,ll,3,1) - ii
      do 220 j = 1, min(ii,nbmax-nn)
      do 210 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(20,k),k)
  210 continue
  220 continue
      do 230 i = 1, 3
      ncll(i,ll,3,1) = ncl(i+20,k) - ncl(20,k) + nn
  230 continue
! particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,1) - ii
      do 250 j = 1, min(ii,nbmax-mm)
      do 240 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  240 continue
  250 continue
      do 260 i = 1, 3
      nclr(i,ll,3,1) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  260 continue
  270 continue
!$OMP END PARALLEL DO
!
! buffer particles and their number leaving the node back or forward:
! update sbufl(:,2), sbufr(:,2), ncll(:,2), nclr(:,2)
      kk = mxyp1*(mzp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,ll)
      do 280 ll = 1, mxyp1
      k = ll
! going straight back or forward
      ncll(1,ll,1,2) = ncl(11,k) - ncl(8,k)
      nclr(1,ll,1,2) = ncl(20,k+kk) - ncl(17,k+kk)
! particles going back and up
      ncll(1,ll,2,2) = ncl(14,k) - ncl(11,k)
! particles going forward and up
      nclr(1,ll,2,2) = ncl(23,k+kk) - ncl(20,k+kk)
! particles going back and down
      ncll(1,ll,3,2) = ncl(17,k) - ncl(14,k)
! particles going forward and down
      nclr(1,ll,3,2) = ncl(26,k+kk) - ncl(23,k+kk)
  280 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
  290 if (kk.ge.mxyp1) go to 310
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 300 ll = 1, mxyp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxyp1) then
         ncll(1,nn,1,2) = ncll(1,nn,1,2) + ncll(1,mm+1,1,2)
         nclr(1,nn,1,2) = nclr(1,nn,1,2) + nclr(1,mm+1,1,2)
         ncll(1,nn,2,2) = ncll(1,nn,2,2) + ncll(1,mm+1,2,2)
         nclr(1,nn,2,2) = nclr(1,nn,2,2) + nclr(1,mm+1,2,2)
         ncll(1,nn,3,2) = ncll(1,nn,3,2) + ncll(1,mm+1,3,2)
         nclr(1,nn,3,2) = nclr(1,nn,3,2) + nclr(1,mm+1,3,2)
      endif
  300 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 290
  310 kk = mxyp1*(mzp1 - 1)
      j1 = ncll(1,mxyp1,1,2)
      k1 = j1 + ncll(1,mxyp1,2,2)
      j2 = nclr(1,mxyp1,1,2)
      k2 = j2 + nclr(1,mxyp1,2,2)
!$OMP PARALLEL DO PRIVATE(i,j,k,ll,ii,nn,mm)
      do 500 ll = 1, mxyp1
      k = ll
! particles going straight up
      ii = ncl(11,k) - ncl(8,k)
      nn = ncll(1,ll,1,2) - ii
      do 330 j = 1, min(ii,nbmax-nn)
      do 320 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(8,k),k)
  320 continue
  330 continue
      do 340 i = 1, 3
      ncll(i,ll,1,2) = ncl(i+8,k) - ncl(8,k) + nn
  340 continue
! particles going straight down
      ii = ncl(20,k+kk) - ncl(17,k+kk)
      mm = nclr(1,ll,1,2) - ii
      do 360 j = 1, min(ii,nbmax-mm)
      do 350 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(17,k+kk),k+kk)
  350 continue
  360 continue
      do 370 i = 1, 3
      nclr(i,ll,1,2) = ncl(i+17,k+kk) - ncl(17,k+kk) + mm
  370 continue
! particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,2) - ii
      do 390 j = 1, min(ii,nbmax-nn)
      do 380 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(11,k),k)
  380 continue
  390 continue
      do 400 i = 1, 3
      ncll(i,ll,2,2) = ncl(i+11,k) - ncl(11,k) + nn
  400 continue
! particles going down and back
      ii = ncl(23,k+kk) - ncl(20,k+kk)
      mm = j2 + nclr(1,ll,2,2) - ii
      do 420 j = 1, min(ii,nbmax-mm)
      do 410 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(20,k+kk),k+kk)
  410 continue
  420 continue
      do 430 i = 1, 3
      nclr(i,ll,2,2) = ncl(i+20,k+kk) - ncl(20,k+kk) + mm
  430 continue
! particles going up and forward
      ii = ncl(17,k) - ncl(14,k)
      nn = k1 + ncll(1,ll,3,2) - ii
      do 450 j = 1, min(ii,nbmax-nn)
      do 440 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(14,k),k)
  440 continue
  450 continue
      do 460 i = 1, 3
      ncll(i,ll,3,2) = ncl(i+14,k) - ncl(14,k) + nn
  460 continue
! particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,2) - ii
      do 480 j = 1, min(ii,nbmax-mm)
      do 470 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  470 continue
  480 continue
      do 490 i = 1, 3
      nclr(i,ll,3,2) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  490 continue
  500 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      kk = max(ncll(3,mxzp1,3,1),nclr(3,mxzp1,3,1))
      ll = max(ncll(3,mxyp1,3,2),nclr(3,mxyp1,3,2))
      ii = max(kk,ll)
! corners overflow
      nn = nclr(3,mx1,2,1) - nclr(3,mxzp1,1,1)
      mm = ncll(3,mx1,2,1) - ncll(3,mxzp1,1,1)
      n = mx1*(mzp1 - 1)
      if (n.gt.0) then
         nr = nclr(3,n,3,1)
         nl = ncll(3,n,3,1)
      else
         nr = nclr(3,mxzp1,2,1)
         nl = ncll(3,mxzp1,2,1)
      endif
      kk = nclr(3,n+mx1,3,1) - nr
      ll = ncll(3,n+mx1,3,1) - nl
! total overflow: result valid only for one processor case
      ii = ii + max(nn+kk,mm+ll)
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,  &
     &mcll,mclr,mcls,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,&
     &irc)
! this subroutine performs second part of a particle sort by x,y,z grid
! in tiles of mx, my, mz
! linear interpolation, with periodic boundary conditions
! for distributed data, with 2d domain decomposition in y/z.
! tiles are assumed to be arranged in 3D linear memory
! incoming particles from other tiles are copied from ppbuff, rbufl, and
! rbufr into ppart
! input: all except ppart, kpic, irc
! output: ppart, kpic, irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppart(3,n,k) = position z of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! rbufl = buffer for particles being received from lower/back processor
! rbufr = buffer for particles being received from upper/forward
! processor
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! mcll = number offset being received from lower/back processor
! mclr = number offset being received from upper/forward processor
! mcls = number ofsets received from corner processors
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! mxzyp1 = mx1*max(myp1,mzp1)
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, mzp1, mxzyp1, npbmx, ntmax, nbmax
      integer irc
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr, mcls
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension rbufl(idimp,nbmax,2), rbufr(idimp,nbmax,2)
      dimension kpic(mx1*myp1*mzp1), ncl(26,mx1*myp1*mzp1)
      dimension ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension mcll(3,mxzyp1,3,2), mclr(3,mxzyp1,3,2), mcls(3,mx1+1,4)
! local data
      integer mxyp1, mxyzp1, nppp, ncoff, joff, koff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, lorr
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, ll, lk, lr, mm, kzs
      logical inside
      integer ks
      dimension ks(26)
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ii,kk,nppp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,mm,kzs
!$OMP& ,ih,nh,ncoff,joff,koff,ist,j1,j2,ip,lorr,inside,ks)
      do 180 l = 1, mxyzp1
      nppp = kpic(l)
      kz = (l - 1)/mxyp1
      k = l - mxyp1*kz
      kzs = kz*mx1
      kz = kz + 1
! loop over tiles in z
      lk = (kz - 1)*mxyp1
! find tile behind
      ll = (kz - 2)*mxyp1
! find tile in front
      lr = kz*mxyp1
      ky = (k - 1)/mx1 + 1
! loop over tiles in y
      kk = (ky - 1)*mx1
! find tile above
      kl = (ky - 2)*mx1
! find tile below
      kr = ky*mx1
! loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
! find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      if (ky.eq.myp1) then
         ks(3) = -kx
         ks(4) = -kxr
         ks(5) = -kxl
      else
         ks(3) = kx + kr + lk
         ks(4) = kxr + kr + lk
         ks(5) = kxl + kr + lk
      endif
      if (ky.eq.1) then
         ks(6) = -kx
         ks(7) = -kxr 
         ks(8) = -kxl
      else
         ks(6) = kx + kl + lk
         ks(7) = kxr + kl + lk
         ks(8) = kxl + kl + lk
      endif
      if (kz.eq.mzp1) then
         ks(9) = -kx
         ks(10) = -kxr
         ks(11) = -kxl
      else
         ks(9) = kx + kk + lr
         ks(10) = kxr + kk + lr
         ks(11) = kxl + kk + lr
      endif
      if ((ky.eq.myp1).or.(kz.eq.mzp1)) then
         ks(12) = -kx
         ks(13) = -kxr
         ks(14) = -kxl
      else
         ks(12) = kx + kr + lr
         ks(13) = kxr + kr + lr
         ks(14) = kxl + kr + lr
      endif
      if ((ky.eq.1).or.(kz.eq.mzp1)) then
         ks(15) = -kx
         ks(16) = -kxr
         ks(17) = -kxl
      else
         ks(15) = kx + kl + lr
         ks(16) = kxr + kl + lr
         ks(17) = kxl + kl + lr
      endif
      if (kz.eq.1) then
         ks(18) = -kx
         ks(19) = -kxr 
         ks(20) = -kxl
      else
         ks(18) = kx + kk + ll
         ks(19) = kxr + kk + ll
         ks(20) = kxl + kk + ll
      endif
      if ((ky.eq.myp1).or.(kz.eq.1)) then
         ks(21) = -kx
         ks(22) = -kxr
         ks(23) = -kxl
      else
         ks(21) = kx + kr + ll
         ks(22) = kxr + kr + ll
         ks(23) = kxl + kr + ll
      endif
      if ((ky.eq.1).or.(kz.eq.1)) then
         ks(24) = -kx
         ks(25) = -kxr
         ks(26) = -kxl
      else
         ks(24) = kx + kl + ll
         ks(25) = kxr + kl + ll
         ks(26) = kxl + kl + ll
      endif
! identify interior
      if ((ky.gt.1).and.(ky.lt.myp1).and.(kz.gt.1).and.(kz.lt.mzp1))    &
     & then
         inside = .true.
      else
         inside = .false.
      endif
! loop over directions
      nh = ihole(1,1,l)
      joff = 0
      koff = 0
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 150 ii = 1, 26
      lorr = 0
      ip = -1
! ip = number of particles coming from direction ii
! interior layers
      if (inside) then
         if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
         ip = ncl(ii,ks(ii)) - ncoff
! edges
      else
! top layer
         if (ky.eq.1) then
            if ((ii.ge.6).and.(ii.le.8)) then
               lorr = -1
               joff = 0
               mm = kzs - ks(ii)
               if (ii.eq.6) then
                  if (mm.gt.1) joff = mcll(3,mm-1,1,1)
               else
                  joff = mcll(ii-6,mm,1,1)
               endif
               ip = mcll(ii-5,mm,1,1) - joff
            else if ((ii.ge.15).and.(ii.le.17)) then
               lorr = -2
               joff = mcll(3,mx1*mzp1,1,1)
               if (kz.lt.mzp1) then
                  mm = kzs + mx1 - ks(ii)
                  if (ii.eq.15) then
                     if (mm.gt.1) joff = mcll(3,mm-1,2,1)
                  else
                     joff = mcll(ii-15,mm,2,1)
                  endif
                  ip = mcll(ii-14,mm,2,1) - joff
! corner data, (ky=1,kz=mzp1)
               else
                  mm = -ks(ii)
                  if (ii.eq.15) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,1)
                     if (mm.gt.1) joff = mcls(3,mm-1,1)
                  else
                     joff = mcls(ii-15,mm,1)
                  endif
                  ip = mcls(ii-14,mm,1) - joff
               endif
            else if ((ii.ge.24).and.(ii.le.26)) then
               lorr = -3
               joff = mcll(3,mx1*mzp1,2,1)
               if (kz.gt.1) then
                  mm = kzs - mx1 - ks(ii)
                  if (ii.eq.24) then
                     if (mm.gt.1) joff = mcll(3,mm-1,3,1)
                  else
                     joff = mcll(ii-24,mm,3,1)
                  endif
                  ip = mcll(ii-23,mm,3,1) - joff
! corner data, (ky=1,kz=1)
               else
                  mm = -ks(ii)
                  if (ii.eq.24) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,2)
                     if (mm.gt.1) joff = mcls(3,mm-1,2)
                  else
                     joff = mcls(ii-24,mm,2)
                  endif
                  ip = mcls(ii-23,mm,2) - joff
               endif
! internal data
            else
               if (ks(ii).gt.0) then
                  if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                  ip = ncl(ii,ks(ii)) - ncoff
               endif
            endif
         endif
! bottom layer
         if (ky.eq.myp1) then
            if ((ii.ge.3).and.(ii.le.5)) then
               lorr = 1
               joff = 0
               mm = kzs - ks(ii)
               if (ii.eq.3) then
                  if (mm.gt.1) joff = mclr(3,mm-1,1,1)
               else
                  joff = mclr(ii-3,mm,1,1)
               endif
               ip = mclr(ii-2,mm,1,1) - joff
            else if ((ii.ge.12).and.(ii.le.14)) then
               lorr = 2
               joff = mclr(3,mx1*mzp1,1,1)
               if (kz.lt.mzp1) then
                  mm = kzs + mx1 - ks(ii)
                  if (ii.eq.12) then
                     if (mm.gt.1) joff = mclr(3,mm-1,2,1)
                  else
                     joff = mclr(ii-12,mm,2,1)
                  endif
                  ip = mclr(ii-11,mm,2,1) - joff
! corner data, (ky=myp1,kz=mzp1)
               else
                  mm = -ks(ii)
                  if (ii.eq.12) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,3)
                     if (mm.gt.1) joff = mcls(3,mm-1,3)
                  else
                     joff = mcls(ii-12,mm,3)
                  endif
                  ip = mcls(ii-11,mm,3) - joff
               endif
            else if ((ii.ge.21).and.(ii.le.23)) then
               lorr = 3
               joff = mclr(3,mx1*mzp1,2,1)
               if (kz.gt.1) then
                  mm = kzs - mx1 - ks(ii)
                  if (ii.eq.21) then
                    if (mm.gt.1) joff = mclr(3,mm-1,3,1)
                  else
                     joff = mclr(ii-21,mm,3,1)
                  endif
                  ip = mclr(ii-20,mm,3,1) - joff
! corner data, (ky=myp1,kz=1)
               else
                  mm = -ks(ii)
                  if (ii.eq.21) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,4)
                     if (mm.gt.1) joff = mcls(3,mm-1,4)
                  else
                     joff = mcls(ii-21,mm,4)
                  endif
                  ip = mcls(ii-20,mm,4) - joff
               endif
! internal data
            else
               if (ks(ii).gt.0) then
                  if (ky.gt.1) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
! front layer
         if (kz.eq.1) then
            if ((ii.ge.18).and.(ii.le.20)) then
               koff = 0
               mm = kk - ks(ii)
               if (ii.eq.18) then
                  if (mm.gt.1) koff = mcll(3,mm-1,1,2)
               else
                  koff = mcll(ii-18,mm,1,2)
               endif
               ip = mcll(ii-17,mm,1,2) - koff
            else if ((ii.ge.21).and.(ii.le.23)) then
               koff = mcll(3,mx1*myp1,1,2)
               if (ky.lt.myp1) then
                  mm = kr - ks(ii)
                  if (ii.eq.21) then
                     if (mm.gt.1) koff = mcll(3,mm-1,2,2)
                  else
                     koff = mcll(ii-21,mm,2,2)
                  endif
                  ip = mcll(ii-20,mm,2,2) - koff
! corner data, already done
!              else
               endif
            else if ((ii.ge.24).and.(ii.le.26)) then
               koff = mcll(3,mx1*myp1,2,2)
               if (ky.gt.1) then
                  mm = kl - ks(ii)
                  if (ii.eq.24) then
                     if (mm.gt.1) koff = mcll(3,mm-1,3,2)
                  else
                     koff = mcll(ii-24,mm,3,2)
                  endif
                  ip = mcll(ii-23,mm,3,2) - koff
! corner data, already done
!              else
               endif
! internal data
            else
               if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
! back layer
         if (kz.eq.mzp1) then
            if ((ii.ge.9).and.(ii.le.11)) then
               koff = 0
               mm = kk - ks(ii)
               if (ii.eq.9) then
                  if (mm.gt.1) koff = mclr(3,mm-1,1,2)
               else
                  koff = mclr(ii-9,mm,1,2)
               endif
               ip = mclr(ii-8,mm,1,2) - koff
            else if ((ii.ge.12).and.(ii.le.14)) then
               koff = mclr(3,mx1*myp1,1,2)
               if (ky.lt.myp1) then
                  mm = kr - ks(ii)
                  if (ii.eq.12) then
                     if (mm.gt.1) koff = mclr(3,mm-1,2,2)
                  else
                     koff = mclr(ii-12,mm,2,2)
                  endif
                  ip = mclr(ii-11,mm,2,2) - koff
! corner data, already done
!              else
               endif
            else if ((ii.ge.15).and.(ii.le.17)) then
               koff = mclr(3,mx1*myp1,2,2)
               if (ky.gt.1) then
                  mm = kl - ks(ii)
                  if (ii.eq.15) then
                     if (mm.gt.1) koff = mclr(3,mm-1,3,2)
                  else
                     koff = mclr(ii-15,mm,3,2)
                  endif
                  ip = mclr(ii-14,mm,3,2) - koff
! corner data, already done
!              else
               endif
! internal data
            else
               if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
      endif
!
      if (ip.lt.0) write (*,*) 'help, ip undefined:l,ii=',l,ii
! copy incoming particles
      do 140 j = 1, ip
      ih = ih + 1
! insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
! place overflow at end of array
      else
         j1 = nppp + 1
         nppp = j1
      endif
      if (j1.le.nppmx) then
! interior layers
         if (inside) then
            do 10 i = 1, idimp
            ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   10       continue
! edges
         else
! top layer
            if (ky.eq.1) then
! external data
               if (lorr.lt.0) then
                  do 20 i = 1, idimp
                  ppart(i,j1,l) = rbufl(i,j+joff,1)
   20             continue
! internal data
               else if (ks(ii).gt.0) then
                  do 30 i = 1, idimp
                  ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   30             continue
               endif
            endif
! bottom layer
            if (ky.eq.myp1) then
! external data
               if (lorr.gt.0) then
                  do 40 i = 1, idimp
                  ppart(i,j1,l) = rbufr(i,j+joff,1)
   40             continue
! internal data
               else if (ks(ii).gt.0) then
                  if (ky.gt.1) then
                     do 50 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   50                continue
                  endif
               endif
            endif
! front layer
            if (kz.eq.1) then
              if ((ii.ge.18).and.(ii.le.20)) then
                  do 60 i = 1, idimp
                  ppart(i,j1,l) = rbufl(i,j+koff,2)
   60             continue
               else if ((ii.ge.21).and.(ii.le.23)) then
                  if (ky.lt.myp1) then
                     do 70 i = 1, idimp
                     ppart(i,j1,l) = rbufl(i,j+koff,2)
   70                continue
                  endif
               else if ((ii.ge.24).and.(ii.le.26)) then
                  if (ky.gt.1) then
                     do 80 i = 1, idimp
                     ppart(i,j1,l) = rbufl(i,j+koff,2)
   80                continue
                  endif
! internal data
               else if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     do 90 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   90                continue
                  endif
               endif
            endif
! back layer
            if (kz.eq.mzp1) then
               if ((ii.ge.9).and.(ii.le.11)) then
                  do 100 i = 1, idimp
                  ppart(i,j1,l) = rbufr(i,j+koff,2)
  100             continue
               else if ((ii.ge.12).and.(ii.le.14)) then
                  if (ky.lt.myp1) then
                     do 110 i = 1, idimp
                     ppart(i,j1,l) = rbufr(i,j+koff,2)
  110                continue
                  endif
               else if ((ii.ge.15).and.(ii.le.17)) then
                  if (ky.gt.1) then
                     do 120 i = 1, idimp
                     ppart(i,j1,l) = rbufr(i,j+koff,2)
  120                continue
                  endif
! internal data
               else if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     do 130 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
  130                continue
                  endif
               endif
            endif
         endif
      else
         ist = 1
      endif
  140 continue
  150 continue
! set error
      if (ist.gt.0) irc = j1
! fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 170 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,l)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 160 i = 1, idimp
            ppart(i,j2,l) = ppart(i,j1,l)
  160       continue
         endif
  170    continue
         nppp = nppp - ip
      endif
      kpic(l) = nppp
  180 continue
!$OMP END PARALLEL DO
      return
      end
