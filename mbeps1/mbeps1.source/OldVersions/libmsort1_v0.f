!-----------------------------------------------------------------------
! Fortran Library for reordering particles
! 1D OpenMP PIC Codes:
! PPORDER1L performs particle reordering into tiles, creates list of
!           particles which are leaving tiles, buffers outgoing
!           particles, then copies buffers into particle array
! PPORDERF1L performs particle reordering into tiles, buffers outgoing
!            particles, and copies buffers into particle array
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: july 23, 2016
!-----------------------------------------------------------------------
      subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx&
     &,mx1,npbmx,ntmax,irc)
! this subroutine sorts particles by x grid in tiles of mx, 
! linear interpolation, with periodic boundary conditions
! tiles are assumed to be arranged in 1D linear memory
! algorithm has 3 steps.  first, one finds particles leaving tile and
! stores their number in each directon, location, and destination in ncl
! and ihole.  second, a prefix scan of ncl is performed and departing
! particles are buffered in ppbuff in direction order.  finally, we copy
! the incoming particles from other tiles into ppart.
! input: all except ppbuff, ncl, ihole, irc
! output: ppart, ppbuff, kpic, ncl, ihole, irc
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
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, mx1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
! local data
      integer noff, npp, ncoff
      integer i, j, k, ii, ih, nh, ist, nn, isum
      integer ip, j1, j2, kxl, kxr
      real anx, edgelx, edgerx, dx
      integer ks
      dimension ks(2)
      anx = real(nx)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,ist,dx,edgelx,edgerx)
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
! use periodic boundary conditions and check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,k) = dx
         else
            ist = 1
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
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) return
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 70 k = 1, mx1
! find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 2
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,k)
      ip = 0
! loop over particles leaving tile
      do 60 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   60 continue
! set error
      if (ip.gt.0) irc = ncl(2,k)
   70 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,npp,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
      do 130 k = 1, mx1
      npp = kpic(k)
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
      do 100 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
! ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 90 j = 1, ip
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
         do 80 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   80    continue
      else
         ist = 1
      endif
   90 continue
  100 continue
! set error
      if (ist.gt.0) irc = j1
! fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
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
         npp = npp - ip
      endif
      kpic(k) = npp
  130 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1,&
     &npbmx,ntmax,irc)
! this subroutine sorts particles by x grid in tiles of mx, 
! linear interpolation, with periodic boundary conditions
! tiles are assumed to be arranged in 1D linear memory.
! the algorithm has 2 steps.  first, a prefix scan of ncl is performed
! and departing particles are buffered in ppbuff in direction order.
! then we copy the incoming particles from other tiles into ppart.
! it assumes that the number, location, and destination of particles 
! leaving a tile have been previously stored in ncl and ihole by the
! GPPUSHF1L subroutine.
! input: all except ppbuff, irc
! output: ppart, ppbuff, kpic, ncl, irc
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
! mx1 = (system length in x direction - 1)/mx + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
! local data
      integer npp, ncoff
      integer i, j, k, ii, ih, nh, ist, isum
      integer ip, j1, j2, kxl, kxr
      integer ks
      dimension ks(2)
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
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
      if (ip.gt.0) irc = ncl(2,k)
   40 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,npp,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
      do 100 k = 1, mx1
      npp = kpic(k)
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
      do 70 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
! ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 60 j = 1, ip
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
         do 50 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   50    continue
      else
         ist = 1
      endif
   60 continue
   70 continue
! set error
      if (ist.gt.0) irc = j1
! fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 90 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 80 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
   80       continue
         endif
   90    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  100 continue
!$OMP END PARALLEL DO
      return
      end
