!-----------------------------------------------------------------------
! Fortran Library for reordering particles
! 2D Vector/MPI/OpenMP PIC Codes:
! VPPPORDER2LA performs first part of particle reordering into tiles,
!             creates list of particles which are leaving tile, and
!             buffers outgoing particles
! VPPPORDERF2LA performs first part of particle reordering into tiles,
!              buffers outgoing particles
! VPPPORDER2LB performs second part of particle reordering into tiles,
!             copies buffers into particle array
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: october 17, 2018
!-----------------------------------------------------------------------
      subroutine VPPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,  &
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc2)
! this subroutine performs first part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! this part of the algorithm has 3 steps.  first, one finds particles
! leaving tile and stores their number in each directon, location, and
! destination in ncl and ihole.  then, a prefix scan of ncl is performed
! and departing particles are buffered in ppbuff in direction order.
! finally, we buffer particles leaving the processor in sbufl and sbufr,
! and store particle number offsets in ncll and nclr.
! input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc2
! output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc2
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower processor
! nclr = number offset being sent to upper processor
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc2 = error codes, returned only if error occurs, when irc2(1) > 0
! when (irc2(1).eq.1), ihole overflow, irc2(2) = new ntmax required
! when (irc2(1).eq.2), ppbuff overflow, irc2(2) = new npbmx required
! when (irc2(1).eq.3), sbufr/sbufl overflow, irc2(2)=new nbmax required
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, mx1, myp1, npbmx
      integer ntmax, nbmax
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, ncl, ihole, ncll, nclr, irc2
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
      dimension irc2(2)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=2)
      integer mxyp1, noffp, moffp, nppp, ipp, joff, nps
      integer i, j, k, m, ii, ih, nh, ist, nn, mm, ip, j1, kk, lb, kxs
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      integer sncl, ks
!dir$ attributes align : 64 :: sncl, ks
      dimension sncl(8), ks(8)
! scratch array
      integer n
!dir$ attributes align : 64 :: n
      dimension n(npblk,lvect)
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,m,noffp,moffp,nppp,ipp,joff,nps,nn,mm,ih,nh,ist,dx,dy&
!$OMP& ,edgelx,edgely,edgerx,edgery,n)
      do 60 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      ih = 0
      nh = 0
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
! clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
! loop over particles in tile
      ipp = nppp/npblk
! outer loop over number of full blocks
      do 40 m = 1, ipp
      joff = npblk*(m - 1)
! inner loop over particles in block
! !dir$ vector aligned
      do 20 j = 1, npblk
      dx = ppart(1,j+joff,k)
      dy = ppart(2,j+joff,k)
! find particles going out of bounds
      ist = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         ist = 2
      else if (dx.lt.edgelx) then
         ist = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j+joff,k) = 0.0
               ist = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         ist = ist + 6
      else if (dy.lt.edgely) then
         ist = ist + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               ppart(2,j+joff,k) = 0.0
               ist = ist - 3
            endif
         endif
      endif
      n(j,1) = ist
   20 continue
! store outgoing particle address and destination
      do 30 j = 1, npblk
      ist = n(j,1)
      if (ist.gt.0) then
         ncl(ist,k) = ncl(ist,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = ist
         else
            nh = 1
         endif
      endif
   30 continue
   40 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 50 j = nps, nppp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
! find particles going out of bounds
      ist = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         ist = 2
      else if (dx.lt.edgelx) then
         ist = 1
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.ge.anx) then
               ppart(1,j,k) = 0.0
               ist = 0
            endif
         endif
      endif
      if (dy.ge.edgery) then
         ist = ist + 6
      else if (dy.lt.edgely) then
         ist = ist + 3
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.ge.any) then
               ppart(2,j,k) = 0.0
               ist = ist - 3
            endif
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,k) = ncl(ist,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = ist
         else
            nh = 1
         endif
      endif
   50 continue
! set error and end of file flag
      if (nh.gt.0) irc2(1) = 1
      ihole(1,1,k) = ih
   60 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc2(1).gt.0) then
         ih = 0
         do 70 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   70    continue
         irc2(2) = ih
         return
      endif
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 200 k = 1, mxyp1
! find address offset for ordered ppbuff array
! !dir$ vector aligned
      do 80 j = 1, 8
      sncl(j) = ncl(j,k)
      ks(j) = j - 1
   80 continue
      kxs = 1
   90 if (kxs.lt.8) then
!dir$ ivdep
         do 100 j = 1, 4
         lb = kxs*ks(j)
!        if ((j+lb+kxs).le.8) then
            sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
!        endif
         ks(j) = ks(j)/2
  100    continue
         kxs = kxs + kxs
         go to 90
      endif
      do 110 j = 1, 8
      sncl(j) = sncl(j) - ncl(j,k)
  110 continue
      nh = ihole(1,1,k)
      ip = 0
! buffer particles that are leaving tile, in direction order
! loop over particles leaving tile
      ipp = nh/npblk
! outer loop over number of full blocks
      do 160 m = 1, ipp
      joff = npblk*(m - 1) + 1
! inner loop over particles in block
      do 120 j = 1, npblk
      n(j,1) = ihole(1,j+joff,k)
      n(j,2) = ihole(2,j+joff,k)
  120 continue
! calculate offsets
      do 130 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
  130 continue
! buffer particles that are leaving tile, in direction order
      do 150 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         do 140 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
  140    continue
      else
         ip = 1
      endif
  150 continue
  160 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 180 j = nps, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 170 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
  170    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  180 continue
      do 190 j = 1, 8
      ncl(j,k) = sncl(j)
  190 continue
! set error
      if (ip.gt.0) irc2(1) = 2
  200 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc2(1).gt.0) then
         ii = 0
         do 210 k = 1, mxyp1
         ii = max(ii,ncl(8,k))
  210    continue
         irc2(2) = ii
         return
      endif
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 220 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
  220 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
  230 if (kk.ge.mx1) go to 250
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 240 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  240 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 230
! load particle buffers
  250 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 320 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 270 j = 1, min(ii,nbmax-nn)
      do 260 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  260 continue
  270 continue
      do 280 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  280 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 300 j = 1, min(ii,nbmax-mm)
      do 290 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  290 continue
  300 continue
      do 310 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  310 continue
  320 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc2(1) = 3
         irc2(2) = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll, &
     &nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc2)
! this subroutine performs first part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! this part of the algorithm has 2 steps.  first, a prefix scan of ncl
! is performed and departing particles are buffered in ppbuff in
! direction order. then, we buffer particles leaving the processor in
! sbufl and sbufr, and store particle number offsets in ncll and nclr.
! it assumes that the number, location, and destination of particles 
! leaving a tile have been previously stored in ncl and ihole by the
! PPGPPUSHF2L subroutine.
! input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc2
! output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc2
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower processor
! nclr = number offset being sent to upper processor
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! when (irc2(1).eq.2), ppbuff overflow, irc2(2) = new npbmx required
! when (irc2(1).eq.3), sbufr/sbufl overflow, irc2(2)=new nbmax required
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr, irc2
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
      dimension irc2(2)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=2)
      integer mxyp1, ipp, joff, nps
      integer i, j, k, m, ii, nh, ist, nn, mm, ip, j1, kk, lb, kxs
      integer sncl, ks
!dir$ attributes align : 64 :: sncl, ks
      dimension sncl(8), ks(8)
! scratch array
      integer n
!dir$ attributes align : 64 :: n
      dimension n(npblk,lvect)
      mxyp1 = mx1*myp1
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 130 k = 1, mxyp1
! find address offset for ordered ppbuff array
! !dir$ vector aligned
      do 10 j = 1, 8
      sncl(j) = ncl(j,k)
      ks(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.8) then
!dir$ ivdep
         do 30 j = 1, 4
         lb = kxs*ks(j)
!        if ((j+lb+kxs).le.8) then
            sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
!        endif
         ks(j) = ks(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      do 40 j = 1, 8
      sncl(j) = sncl(j) - ncl(j,k)
   40 continue
      nh = ihole(1,1,k)
      ip = 0
! buffer particles that are leaving tile, in direction order
! loop over particles leaving tile
      ipp = nh/npblk
! outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1) + 1
! inner loop over particles in block
      do 50 j = 1, npblk
      n(j,1) = ihole(1,j+joff,k)
      n(j,2) = ihole(2,j+joff,k)
   50 continue
! calculate offsets
      do 60 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
   60 continue
! buffer particles that are leaving tile, in direction order
      do 80 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         do 70 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   70    continue
      else
         ip = 1
      endif
   80 continue
   90 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 110 j = nps, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 100 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
  100    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  110 continue
      do 120 j = 1, 8
      ncl(j,k) = sncl(j)
  120 continue
! set error
      if (ip.gt.0) irc2(1) = 2
  130 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc2(1).gt.0) then
         ii = 0
         do 140 k = 1, mxyp1
         ii = max(ii,ncl(8,k))
  140    continue
         irc2(2) = ii
         return
      endif
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 150 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
  150 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
  160 if (kk.ge.mx1) go to 180
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 170 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  170 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 160
  180 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 250 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 200 j = 1, min(ii,nbmax-nn)
      do 190 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  190 continue
  200 continue
      do 210 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  210 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 230 j = 1, min(ii,nbmax-mm)
      do 220 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  220 continue
  230 continue
      do 240 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  240 continue
  250 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc2(1) = 3
         irc2(2) = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine VPPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,  &
     &mcll,mclr,idimp,nppmx,nx,ny,mx1,myp1,npbmx,ntmax,nbmax,irc2)
! this subroutine performs second part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! incoming particles from other tiles are copied from ppbuff, rbufl, and
! rbufr into ppart
! input: all except ppart, kpic, irc2
! output: ppart, kpic, irc2
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! mcll = number offset being received from lower processor
! mclr = number offset being received from upper processor
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! when (irc2(1).eq.4), ppart overflow, irc2(2) = new nppmx required
      implicit none
      integer idimp, nppmx, nx, ny, mx1, myp1, npbmx, ntmax, nbmax
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr, irc2
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension mcll(3,mx1), mclr(3,mx1)
      dimension irc2(2)
! local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer mxyp1, nppp, ncoff, noff, moff, joff, ipp, nps
      integer i, j, k, m, nn, mm, ii, in, kx, ky, ih, nh, ist
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      real anx, any, dx, dy
      integer ks
!dir$ attributes align : 64 :: ks
      dimension ks(8)
! scratch integer array
      integer n
!dir$ attributes align : 64 :: n
      dimension n(npblk,lvect)
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,   &
!$OMP& moff,ist,j1,ip,dx,dy,ks)
      do 60 k = 1, mxyp1
      nppp = kpic(k)
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
      ks(1) = kxr + kk
      ks(2) = kxl + kk
      ks(3) = kx + kr
      ks(4) = kxr + kr
      ks(5) = kxl + kr
      ks(6) = kx + kl
      ks(7) = kxr + kl
      ks(8) = kxl + kl
! loop over directions
      nh = ihole(1,1,k)
      noff = 0
      moff = 0
      if (ky.eq.1) then
         if (kx.gt.1) noff = mcll(3,kx-1)
      endif
      if (ky.eq.myp1) then
         if (kx.gt.1) moff = mclr(3,kx-1)
      endif
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 50 ii = 1, 8
! ip = number of particles coming from direction ii
      if (ks(ii).le.0) then
         if (ii.gt.6) noff = mcll(ii-6,ks(ii)+mx1)
         ip = mcll(ii-5,ks(ii)+mx1) - noff
      else if (ks(ii).gt.mxyp1) then
         if (ii.gt.3) moff = mclr(ii-3,ks(ii)-mxyp1)
         ip = mclr(ii-2,ks(ii)-mxyp1) - moff
      else
         if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
         ip = ncl(ii,ks(ii)) - ncoff
      endif
      do 40 j = 1, ip
      ih = ih + 1
! insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
! place overflow at end of array
      else
         j1 = nppp + 1
         nppp = j1
      endif
      if (j1.le.nppmx) then
! particle coming from left node
         if (ks(ii).le.0) then
! check for periodic boundary conditions
            dx = rbufl(1,j+noff)
            if (dx.lt.0.0) dx = dx + anx
            if (dx.ge.anx) dx = dx - anx
            ppart(1,j1,k) = dx
            dy = rbufl(2,j+noff)
            if (dy.lt.0.0) dy = dy + any
            if (dy.ge.any) dy = dy - any
            ppart(2,j1,k) = dy
! copy remaining particle data
            do 10 i = 3, idimp
            ppart(i,j1,k) = rbufl(i,j+noff)
   10       continue
! particle coming from right node
         else if (ks(ii).gt.mxyp1) then
! check for periodic boundary conditions
            dx = rbufr(1,j+moff)
            if (dx.lt.0.0) dx = dx + anx
            if (dx.ge.anx) dx = dx - anx
            ppart(1,j1,k) = dx
            dy = rbufr(2,j+moff)
            if (dy.lt.0.0) dy = dy + any
            if (dy.ge.any) dy = dy - any
            ppart(2,j1,k) = dy
! copy remaining particle data
            do 20 i = 3, idimp
            ppart(i,j1,k) = rbufr(i,j+moff)
   20       continue
! particle coming from interior
         else
! check for periodic boundary conditions
            dx = ppbuff(1,j+ncoff,ks(ii))
            if (dx.lt.0.0) dx = dx + anx
            if (dx.ge.anx) dx = dx - anx
            ppart(1,j1,k) = dx
            dy = ppbuff(2,j+ncoff,ks(ii))
            if (dy.lt.0.0) dy = dy + any
            if (dy.ge.any) dy = dy - any
            ppart(2,j1,k) = dy
! copy remaining particle data
            do 30 i = 3, idimp
            ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   30       continue
         endif
      else
         ist = 1
      endif
   40 continue
   50 continue
! save parameters for next loop
      if (ih.lt.nh) then
         ihole(2,1,k) = -(ih+1)
      else
         ihole(2,1,k) = nppp
      endif
! set error
      if (ist.gt.0) irc2(1) = 4
   60 continue
!$OMP END PARALLEL DO
! ppart overflow
      if (irc2(1).gt.0) then
         j1 = 0
         do 70 k = 1, mxyp1
         j1 = max(j1,ihole(2,1,k))
   70    continue
         irc2(2) = j1
         return
      endif
! fill up remaining holes in particle array with particles from bottom
! holes with locations great than npp-ip do not need to be filled
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,m,ih,nppp,nh,ip,ii,ipp,joff,nps,in,nn,mm,j1,j2,n)
      do 150 k = 1, mxyp1
      ih = ihole(2,1,k)
      if (ih.lt.0) then
         nppp = kpic(k)
         nh = ihole(1,1,k)
         ip = nh + ih + 1
! move particles from end into remaining holes
! holes are processed in increasing order
         ii = nh + 1
         ih = 1 - ih
         ipp = ip/npblk
! outer loop over number of full blocks
         do 120 m = 1, ipp
         joff = npblk*(m - 1)
! inner loop over particles in block
         do 80 j = 1, npblk
         n(j,2) = ihole(1,ih+j-1,k)
         n(j,3) = ihole(1,ii-j+1,k)
   80    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 90 j = 1, npblk
         j1 = nppp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            n(j,1) = -1
            in = in + 1
            if (j.lt.npblk) nn = n(in,3)
         else
            mm = mm + 1
         endif
   90    continue
         do 110 j = 1, npblk
         j1 = nppp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) then
            do 100 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  100       continue
         endif
  110    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  120    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,k)
         j2 = ihole(1,ih,k)
! loop over remaining particles
         do 140 j = nps, ip
         j1 = nppp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,k)
         else
            do 130 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  130       continue
            ih = ih + 1
            j2 = ihole(1,ih,k)
         endif
  140    continue
         kpic(k) = nppp - ip
      else
         kpic(k) = ihole(2,1,k)
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
