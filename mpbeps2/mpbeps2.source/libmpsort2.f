!-----------------------------------------------------------------------
! Fortran Library for reordering particles
! 2D MPI/OpenMP PIC Codes:
! PPPORDER2LA performs first part of particle reordering into tiles,
!             creates list of particles which are leaving tile, and
!             buffers outgoing particles
! PPPORDERF2LA performs first part of particle reordering into tiles,
!              buffers outgoing particles
! PPPORDER2LB performs second part of particle reordering into tiles,
!             copies buffers into particle array
! PPPRSNCL2L restores initial values of ncl array
! PPPORDERF2LAF performs the final section of first part of a particle
!               reordering into tiles, buffers outgoing particles
! PPPRSTOR2L restores particle coordinates from ppbuff
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 26, 2017
!-----------------------------------------------------------------------
      subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,   &
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
! idimp = size of phase space = 4
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
      integer mxyp1, noffp, moffp, nppp
      integer i, j, k, ii, ih, nh, ist, nn, mm, isum, ip, j1, kk
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,
!$OMP& edgerx,edgery)
      do 30 k = 1, mxyp1
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
      do 20 j = 1, nppp
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
   20 continue
! set error and end of file flag
      if (nh.gt.0) irc2(1) = 1
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc2(1).gt.0) then
         ih = 0
         do 40 k = 1, mxyp1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc2(2) = ih
         return
      endif
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 80 k = 1, mxyp1
! find address offset for ordered ppbuff array
      isum = 0
      do 50 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   50 continue
      nh = ihole(1,1,k)
      ip = 0
! loop over particles leaving tile
      do 70 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 60 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   60    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   70 continue
! set error
      if (ip.gt.0) irc2(1) = 2
   80 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc2(1).gt.0) then
         ii = 0
         do 90 k = 1, mxyp1
         ii = max(ii,ncl(8,k))
   90    continue
         irc2(2) = ii
         return
      endif
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 100 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
  100 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
  110 if (kk.ge.mx1) go to 130
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 120 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  120 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 110
! load particle buffers
  130 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 200 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 150 j = 1, min(ii,nbmax-nn)
      do 140 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  140 continue
  150 continue
      do 160 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  160 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 180 j = 1, min(ii,nbmax-mm)
      do 170 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  170 continue
  180 continue
      do 190 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  190 continue
  200 continue
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
      subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,  &
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
! idimp = size of phase space = 4
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
      integer mxyp1
      integer i, j, k, ii, nh, ist, nn, mm, isum, ip, j1, kk
      mxyp1 = mx1*myp1
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mxyp1
! find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,k)
      ip = 0
! loop over particles leaving tile
      do 30 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   30 continue
! set error
      if (ip.gt.0) irc2(1) = 2
   40 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc2(1).gt.0) then
         ii = 0
         do 50 k = 1, mxyp1
         ii = max(ii,ncl(8,k))
   50    continue
         irc2(2) = ii
         return
      endif
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 60 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   60 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   70 if (kk.ge.mx1) go to 90
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 80 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
   80 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 70
   90 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 160 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 110 j = 1, min(ii,nbmax-nn)
      do 100 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  100 continue
  110 continue
      do 120 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  120 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 140 j = 1, min(ii,nbmax-mm)
      do 130 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  130 continue
  140 continue
      do 150 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  150 continue
  160 continue
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
      subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,   &
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
! idimp = size of phase space = 4
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
      integer mxyp1, nppp, ncoff, noff, moff
      integer i, j, k, ii, kx, ky, ih, nh, ist
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      real anx, any, dx, dy
      integer ks
      dimension ks(8)
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,
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
            if (dx.ge.anx) dx = dx - anx
            if (dx.lt.0.0) dx = dx + anx
            ppart(1,j1,k) = dx
            dy = rbufl(2,j+noff)
            if (dy.ge.any) dy = dy - any
            if (dy.lt.0.0) dy = dy + any
            ppart(2,j1,k) = dy
! copy remaining particle data
            do 10 i = 3, idimp
            ppart(i,j1,k) = rbufl(i,j+noff)
   10       continue
! particle coming from right node
         else if (ks(ii).gt.mxyp1) then
! check for periodic boundary conditions
            dx = rbufr(1,j+moff)
            if (dx.ge.anx) dx = dx - anx
            if (dx.lt.0.0) dx = dx + anx
            ppart(1,j1,k) = dx
            dy = rbufr(2,j+moff)
            if (dy.ge.any) dy = dy - any
            if (dy.lt.0.0) dy = dy + any
            ppart(2,j1,k) = dy
! copy remaining particle data
            do 20 i = 3, idimp
            ppart(i,j1,k) = rbufr(i,j+moff)
   20       continue
! particle coming from interior
         else
! check for periodic boundary conditions
            dx = ppbuff(1,j+ncoff,ks(ii))
            if (dx.ge.anx) dx = dx - anx
            if (dx.lt.0.0) dx = dx + anx
            ppart(1,j1,k) = dx
            dy = ppbuff(2,j+ncoff,ks(ii))
            if (dy.ge.any) dy = dy - any
            if (dy.lt.0.0) dy = dy + any
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
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,nppp,ih,nh,j1,j2,ip)
      do 100 k = 1, mxyp1
      ih = ihole(2,1,k)
      if (ih.lt.0) then
         nppp = kpic(k)
         nh = ihole(1,1,k)
         ip = nh + ih + 1
         do 90 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 80 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
   80       continue
         endif
   90    continue
         kpic(k) = nppp - ip
      else
         kpic(k) = ihole(2,1,k)
      endif
  100 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPRSNCL2L(ncl,mxyp1)
! this subroutine restores initial values of ncl array
! for distributed data, with 1d domain decomposition in y.
! input: all, output: ncl
! ncl(i,k) = number of particles going to destination i, tile k
! mxyp1 = total number of tiles
      implicit none
      integer mxyp1
      integer ncl
      dimension ncl(8,mxyp1)
! local data
      integer j, k, noff, ist
! restores address offset array: update ncl
! !$OMP PARALLEL DO PRIVATE(j,k,noff,ist)
      do 20 k = 1, mxyp1
! find restore ncl for ordered ppbuff array
      noff = 0
      do 10 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = ist - noff
      noff = ist
   10 continue
   20 continue
! !$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDERF2LAF(ppbuff,sbufl,sbufr,ncl,ncll,nclr,idimp,  &
     &mx1,myp1,npbmx,nbmax,irc2)
! this subroutine performs the final section of first part of a particle
! sort by x,y grid in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! particles leaving the processor are buffered in sbufl and sbufr, and
! particle number offsets are stored in ncll and nclr.
! input: all except sbufl, sbufr, ncll, nclr, irc2
! output: sbufl, sbufr, ncll, nclr, irc2
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! ncl(i,k) = number of particles going to destination i, tile k
! ncll = number offset being sent to lower processor
! nclr = number offset being sent to upper processor
! idimp = size of phase space = 4
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! nbmax =  size of buffers for passing particles between processors
! when (irc2(1).eq.3), sbufr/sbufl overflow, irc2(2)=new nbmax required
      implicit none
      integer idimp, mx1, myp1, npbmx, nbmax
      real ppbuff, sbufl, sbufr
      integer ncl, ncll, nclr, irc2
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension ncl(8,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
      dimension irc2(2)
! local data
      integer i, j, k, ii, nn, mm, kk
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 60 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   60 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   70 if (kk.ge.mx1) go to 90
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 80 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
   80 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 70
   90 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 160 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 110 j = 1, min(ii,nbmax-nn)
      do 100 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  100 continue
  110 continue
      do 120 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  120 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 140 j = 1, min(ii,nbmax-mm)
      do 130 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  130 continue
  140 continue
      do 150 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  150 continue
  160 continue
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
      subroutine PPPRSTOR2L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mxyp1,   &
     &npbmx,ntmax)
! this subroutine restores particle coordinates from ppbuff
! used in resizing segmented particle array ppart if overflow occurs
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! input: all, output: ppart, ncl
! ppart(i,n,k) = i co-ordinate of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = total number of tiles
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
      implicit none
      integer idimp, nppmx, mxyp1, npbmx, ntmax
      real ppart, ppbuff
      integer ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension ppbuff(idimp,npbmx,mxyp1)
      dimension ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
      integer i, j, k, nh, j1, ist, ii
! restores particles that are leaving tile: update ppart, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,nh,j1,ist,ii)
      do 40 k = 1, mxyp1
! find restore address offset for ordered ppbuff array
      do 10 j = 1, 7
      ncl(9-j,k) = ncl(8-j,k)
   10 continue
      ncl(1,k) = 0
      nh = ihole(1,1,k)
! loop over particles leaving tile
      do 30 j = 1, nh
! restore particles from buffer, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      do 20 i = 1, idimp
      ppart(i,j1,k) = ppbuff(i,ii,k)
   20 continue
      ncl(ist,k) = ii
   30 continue
   40 continue
!$OMP END PARALLEL DO
      return
      end
