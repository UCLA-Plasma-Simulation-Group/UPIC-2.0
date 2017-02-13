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
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 14, 2016
!-----------------------------------------------------------------------
      subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,   &
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
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
! input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
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
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
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
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,k) = dy
         else
            ist = ist + 3
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
      do 70 k = 1, mxyp1
! find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 8
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
      if (ip.gt.0) irc = ncl(8,k)
   70 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 80 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   80 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   90 if (kk.ge.mx1) go to 110
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 100 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  100 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 90
  110 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 180 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 130 j = 1, min(ii,nbmax-nn)
      do 120 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  120 continue
  130 continue
      do 140 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  140 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 160 j = 1, min(ii,nbmax-mm)
      do 150 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  150 continue
  160 continue
      do 170 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  170 continue
  180 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,  &
     &nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
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
! input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
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
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
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
      if (ip.gt.0) irc = ncl(8,k)
   40 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 50 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   50 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   60 if (kk.ge.mx1) go to 80
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 70 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
   70 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 60
   80 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 150 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 100 j = 1, min(ii,nbmax-nn)
      do 90 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
   90 continue
  100 continue
      do 110 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  110 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 130 j = 1, min(ii,nbmax-mm)
      do 120 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  120 continue
  130 continue
      do 140 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  140 continue
  150 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,   &
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
! this subroutine performs second part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! incoming particles from other tiles are copied from ppbuff, rbufl, and
! rbufr into ppart
! input: all except ppart, kpic, irc
! output: ppart, kpic, irc
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
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension mcll(3,mx1), mclr(3,mx1)
! local data
      integer mxyp1, nppp, ncoff, noff, moff
      integer i, j, k, ii, kx, ky, ih, nh, ist
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer ks
      dimension ks(8)
      mxyp1 = mx1*myp1
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,
!$OMP& moff,ist,j1,j2,ip,ks)
      do 200 k = 1, mxyp1
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
      do 170 ii = 1, 8
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
      do 160 j = 1, ip
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
         if (ks(ii).le.0) then
            do 130 i = 1, idimp
            ppart(i,j1,k) = rbufl(i,j+noff)
  130       continue
         else if (ks(ii).gt.mxyp1) then
            do 140 i = 1, idimp
            ppart(i,j1,k) = rbufr(i,j+moff)
  140       continue
         else
            do 150 i = 1, idimp
            ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
  150       continue
         endif
      else
         ist = 1
      endif
  160 continue
  170 continue
! set error
      if (ist.gt.0) irc = j1
! fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 190 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 180 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  180       continue
         endif
  190    continue
         nppp = nppp - ip
      endif
      kpic(k) = nppp
  200 continue
!$OMP END PARALLEL DO
      return
      end
