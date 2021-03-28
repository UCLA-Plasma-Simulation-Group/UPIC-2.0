!-----------------------------------------------------------------------
! Fortran Library for reordering particles
! 1D OpenMP PIC Codes:
! PPORDER1L performs particle reordering into tiles, creates list of
!           particles which are leaving tiles, buffers outgoing
!           particles, then copies buffers into particle array
! PPORDERF1L performs particle reordering into tiles, buffers outgoing
!            particles, and copies buffers into particle array
! PPRSNCL1L restores initial values of ncl array
! PPRSTOR1L restores particle coordinates from ppbuff
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 30, 2017
!-----------------------------------------------------------------------
      subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx&
     &,mx1,npbmx,ntmax,irc2)
! this subroutine sorts particles by x grid in tiles of mx, 
! linear interpolation, with periodic boundary conditions
! tiles are assumed to be arranged in 1D linear memory
! algorithm has 3 steps.  first, one finds particles leaving tile and
! stores their number in each directon, location, and destination in ncl
! and ihole.  second, a prefix scan of ncl is performed and departing
! particles are buffered in ppbuff in direction order.  finally, we copy
! the incoming particles from other tiles into ppart.
! input: all except ppbuff, ncl, ihole, irc2
! output: ppart, ppbuff, kpic, ncl, ihole, irc2
! ppart(1,n,k) = position x of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! irc2 = error codes, returned only if error occurs, when irc2(1) > 0
! when (irc2(1).eq.1), ihole overflow, irc2(2) = new ntmax required
! when (irc2(1).eq.2), ppbuff overflow, irc2(2) = new npbmx required
! when (irc2(1).eq.3), ppart overflow, irc2(2) = new nppmx required
! when (irc2(1).eq.(-1)), incoming particles not in proper tile
      implicit none
      integer idimp, nppmx, nx, mx, mx1, npbmx, ntmax
      real ppart, ppbuff
      integer kpic, ncl, ihole, irc2
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
      dimension irc2(2)
! local data
      integer noff, npp, ncoff
      integer i, j, k, ii, ih, nh, ist, nn, isum
      integer ip, j1, j2, kxl, kxr
      real anx, edgelx, edgerx, dx, dxt
      integer ks
      dimension ks(2)
      anx = real(nx)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,ist,dx,edgelx,edgerx)              &
!$OMP& SCHEDULE(dynamic)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
! clear counters
      do 10 j = 1, 2
      ncl(j,k) = 0
   10 continue
! loop over particles in tile
      do 20 j = 1, npp
      dx = ppart(1,j,k)
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
         do 40 k = 1, mx1
         ih = max(ih,ihole(1,1,k))
   40    continue
         irc2(2) = ih
         return
      endif
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 80 k = 1, mx1
! find address offset for ordered ppbuff array
      isum = 0
      do 50 j = 1, 2
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
         do 90 k = 1, mx1
         ii = max(ii,ncl(2,k))
   90    continue
         irc2(2) = ii
         return
      endif
!
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noff,npp,nn,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,dx,&
!$OMP& dxt,edgelx,edgerx,ks)
      do 130 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
! loop over tiles in x, assume periodic boundary conditions
      kxl = k - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = k + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
! find tile number for different directions
      ks(1) = kxr
      ks(2) = kxl
! loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      j2 = 0
      do 120 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
! ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 110 j = 1, ip
      ih = ih + 1
! insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
! place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         dx = ppbuff(1,j+ncoff,ks(ii))
! check for out of bounds
         if ((dx.ge.edgerx).or.(dx.lt.edgelx)) then
            dxt = dx
! check for periodic boundary conditions
            if (dx.ge.anx) dx = dx - anx
            if (dx.lt.0.0) dx = dx + anx
! check again for out of bounds
            if ((dx.ge.edgerx).or.(dx.lt.edgelx)) then
               dx = dxt
               j2 = 1
            endif
         endif
         ppart(1,j1,k) = dx
! copy remaining particle data
         do 100 i = 2, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
  100    continue
      else
         ist = 1
      endif
  110 continue
  120 continue
! save parameters for next loop
      if (ih.lt.nh) then
         ihole(2,1,k) = -(ih+1)
      else
         ihole(2,1,k) = npp
      endif
! set overflow error
      if (ist.gt.0) irc2(1) = 3
! set particle out of bounds error
      if (j2.gt.0) irc2(2) = -j2
  130 continue
!$OMP END PARALLEL DO
! ppart overflow
      if (irc2(1).gt.0) then
         j1 = 0
         do 140 k = 1, mx1
         j1 = max(j1,ihole(2,1,k))
  140    continue
         irc2(2) = j1
         return
! particle out of bounds
      else if (irc2(2).lt.0) then
         irc2(1) = -1
         irc2(2) = -irc2(2)
      endif
! fill up remaining holes in particle array with particles from bottom
!$OMP PARALLEL DO PRIVATE(i,j,k,npp,ih,nh,j1,j2,ip)
      do 170 k = 1, mx1
      ih = ihole(2,1,k)
      if (ih.lt.0) then
         npp = kpic(k)
         nh = ihole(1,1,k)
         ip = nh + ih + 1
         do 160 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 150 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  150       continue
         endif
  160    continue
         kpic(k) = npp - ip
      else
         kpic(k) = ihole(2,1,k)
      endif
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx, &
     &mx,mx1,npbmx,ntmax,irc2)
! this subroutine sorts particles by x grid in tiles of mx, 
! linear interpolation, with periodic boundary conditions
! tiles are assumed to be arranged in 1D linear memory.
! the algorithm has 2 steps.  first, a prefix scan of ncl is performed
! and departing particles are buffered in ppbuff in direction order.
! then we copy the incoming particles from other tiles into ppart.
! it assumes that the number, location, and destination of particles 
! leaving a tile have been previously stored in ncl and ihole by the
! GPPUSHF1L subroutine.
! input: all except ppbuff, irc2
! output: ppart, ppbuff, kpic, ncl, irc2
! ppart(1,n,k) = position x of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! mx1 = (system length in x direction - 1)/mx + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! irc2 = error codes, returned only if error occurs, when irc2(1) > 0
! when (irc2(1).eq.2), ppbuff overflow, irc2(2) = new npbmx required
! when (irc2(1).eq.3), ppart overflow, irc2(2) = new nppmx required
! when (irc2(1).eq.(-1)), incoming particles not in proper tile
      implicit none
      integer idimp, nppmx, nx, mx, mx1, npbmx, ntmax
      real ppart, ppbuff
      integer kpic, ncl, ihole, irc2
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
      dimension irc2(2)
! local data
      integer noff, npp, ncoff
      integer i, j, k, ii, ih, nh, ist, nn, isum
      integer ip, j1, j2, kxl, kxr
      real anx, edgelx, edgerx, dx, dxt
      integer ks
      dimension ks(2)
      anx = real(nx)
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mx1
! find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 2
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
         do 50 k = 1, mx1
         ii = max(ii,ncl(2,k))
   50    continue
         irc2(2) = ii
         return
      endif
!
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noff,npp,nn,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,dx,&
!$OMP& dxt,edgelx,edgerx,ks)
      do 90 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
! loop over tiles in x, assume periodic boundary conditions
      kxl = k - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = k + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
! find tile number for different directions
      ks(1) = kxr
      ks(2) = kxl
! loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      j2 = 0
      do 80 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
! ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 70 j = 1, ip
      ih = ih + 1
! insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
! place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         dx = ppbuff(1,j+ncoff,ks(ii))
! check for out of bounds
         if ((dx.ge.edgerx).or.(dx.lt.edgelx)) then
            dxt = dx
! check for periodic boundary conditions
            if (dx.ge.anx) dx = dx - anx
            if (dx.lt.0.0) dx = dx + anx
! check again for out of bounds
            if ((dx.ge.edgerx).or.(dx.lt.edgelx)) then
               dx = dxt
               j2 = 1
            endif
         endif
         ppart(1,j1,k) = dx
! copy remaining particle data
         do 60 i = 2, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   60    continue
      else
         ist = 1
      endif
   70 continue
   80 continue
! save parameters for next loop
      if (ih.lt.nh) then
         ihole(2,1,k) = -(ih+1)
      else
         ihole(2,1,k) = npp
      endif
! set overflow error
      if (ist.gt.0) irc2(1) = 3
! set particle out of bounds error
      if (j2.gt.0) irc2(2) = -j2
   90 continue
!$OMP END PARALLEL DO
! ppart overflow
      if (irc2(1).gt.0) then
         j1 = 0
         do 100 k = 1, mx1
         j1 = max(j1,ihole(2,1,k))
  100    continue
         irc2(2) = j1
         return
! particle out of bounds
      else if (irc2(2).lt.0) then
         irc2(1) = -1
         irc2(2) = -irc2(2)
      endif
! fill up remaining holes in particle array with particles from bottom
!$OMP PARALLEL DO PRIVATE(i,j,k,npp,ih,nh,j1,j2,ip)
      do 130 k = 1, mx1
      ih = ihole(2,1,k)
      if (ih.lt.0) then
         npp = kpic(k)
         nh = ihole(1,1,k)
         ip = nh + ih + 1
         do 120 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 110 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  110       continue
         endif
  120    continue
         kpic(k) = npp - ip
      else
         kpic(k) = ihole(2,1,k)
      endif
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRSNCL1L(ncl,mx1)
! this subroutine restores initial values of ncl array
! input: all, output: ncl
! ncl(i,k) = number of particles going to destination i, tile k
! mx1 = total number of tiles
      implicit none
      integer mx1
      integer ncl
      dimension ncl(2,mx1)
! local data
      integer j, k, noff, ist
! restores address offset array: update ncl
! !$OMP PARALLEL DO PRIVATE(j,k,noff,ist)
      do 20 k = 1, mx1
! find restore ncl for ordered ppbuff array
      noff = 0
      do 10 j = 1, 2
      ist = ncl(j,k)
      ncl(j,k) = ist - noff
      noff = ist
   10 continue
   20 continue
! !$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPRSTOR1L(ppart,ppbuff,ncl,ihole,idimp,nppmx,mx1,npbmx,&
     &ntmax)
! this subroutine restores particle coordinates from ppbuff
! used in resizing segmented particle array ppart if overflow occurs
! tiles are assumed to be arranged in 1D linear memory.
! input: all, output: ppart, ncl
! ppart(i,n,k) = i co-ordinate of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
      implicit none
      integer idimp, nppmx, mx1, npbmx, ntmax
      real ppart, ppbuff
      integer ncl, ihole
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
! restores particles that are leaving tile: update ppart, ncl
! loop over tiles
      integer i, j, k, nh, j1, ist, ii
!$OMP PARALLEL DO PRIVATE(i,j,k,nh,j1,ist,ii)
      do 30 k = 1, mx1
! find restore address offset for ordered ppbuff array
      ncl(2,k) = ncl(1,k)
      ncl(1,k) = 0
      nh = ihole(1,1,k)
! loop over particles leaving tile
      do 20 j = 1, nh
! restore particles from buffer, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      do 10 i = 1, idimp
      ppart(i,j1,k) = ppbuff(i,ii,k)
   10 continue
      ncl(ist,k) = ii
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
