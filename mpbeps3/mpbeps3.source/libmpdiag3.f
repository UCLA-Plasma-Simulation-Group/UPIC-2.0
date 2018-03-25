!-----------------------------------------------------------------------
! Fortran Library for diagnostics
! 3D MPI/OpenMP PIC Codes:
! PPVDIST32 calculates 3 component velocity distribution, velocity
!           moments and entropy for segmented particle array,
!           for distributed data.
! PERPDIST32 calculates 3d energy distribution for relativistic
!            particles
! PPVBDIST32 for 3d code, calculates 3d velocity distribution and
!            velocity moments for magnetized plasma with segmented
!            particle array
! PPVSDIST32 for 3d code, calculates 3d velocity distribution, in
!            different regions of space, for segmented particle array
!            with distributed data.
! PVDIST32 calculates 3 component velocity distribution, velocity
!          moments, and entropy for standard particle array,
!          for distributed data.
! PVBDIST32 for 3d code, calculates 3d velocity distribution and
!           velocity moments for magnetized plasma
! PPROFX32L calculates fluid moments from particle quantities: density,
!           momentum, momentum flux, energy, energy flux, assumes
!           particle positions and velocities at same time level
! PRPROFX32L calculates fluid moments from relativistic particle
!            quantities: density, velocity, velocity flux, energy,
!            energy flux, assumes particle positions and velocities at 
!            same time level
! PGPROFX32L calculates fluid moments from particle quantities:
!            density, momentum, momentum flux, energy, energy flux,
!            assumes particle positions and velocities not at same time
!            levels and electrostatic fields
! PGRPROFX32L calculates fluid moments from relativistic particle
!             quantities: density, velocity, velocity flux, energy,
!             energy flux, assumes particle positions and velocities
!             not at same time levels and electrostatic fields
! PGBPROFX32L calculates fluid moments from particle quantities:
!             density, momentum, momentum flux, energy, energy flux,
!             assumes particle positions and velocities not at same time
!             levels and electromagnetic fields
! PGRBPROFX32L calculates fluid moments from relativistic particle
!              quantities: density, velocity, velocity flux, energy,
!              energy flux
!              assumes particle positions and velocities not at same
!              time levels and electromagnetic fields
! FLUIDQS3 calculates fluid quantities from fluid moments:
!          density, velocity field, pressure tensor, energy, heat flux
! PSTPTRAJ3 sets test charge distribution by setting a particle id
!           in particle location 7 for 3d code
! PPTRAJ3 copies tagged particles in ppart to array partt for 3 code
! PORDTRAJ3 reorders tagged particles in partt to array spartt for 3d
! PCPYTRAJ2 copies tagged particles in partt to array part
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: march 21, 2018
!-----------------------------------------------------------------------
      subroutine PPVDIST32(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyzp1,&
     &nmv,nmvf)
! for 3d code, this subroutine calculates 3d velocity distribution,
! velocity moments, and entropy, with 2D domain decomposition
! particles stored in segmented array
! input: all except fvm, output: fv, sfv, fvm
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = velocity vy of particle n in tile m
! ppart(6,n,m) = velocity vz of particle n in tile m
! fv = distribution function, number of particles in each
! velocity range, summed over tiles. maximum velocity (used for scaling)
! is contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! velocity range
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! entropy for i-th dimension is contained in fvm(i,3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v*delta_x)).
! Assumes that distributions in each dimension are independent.
! nvp = number of real or virtual processors
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx+1 and 
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nvp, idimp, nppmx, mxyzp1, nmv, nmvf
      real ppart, fv, sfv, fvm
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fv(nmvf,3), sfv(nmvf,3,mxyzp1), fvm(3,3)
      integer kpic
      dimension kpic(mxyzp1)
! local data
      integer j, k, nppp, nmv21, nvx, nvy, nvz
      real anmv, svx, svy, svz, svxx, svyx, svzx, vx, vy, vz
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      double precision ssumvx, ssumvy, ssumvz, ssumvx2, ssumvy2, ssumvz2
      double precision sanp
      double precision sum1, sum2, sum3
      double precision sum7, work7
      dimension sum7(7), work7(7)
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svy = anmv/fv(nmv21+1,2)
      svz = anmv/fv(nmv21+1,3)
! normalization constant for entropy
      svxx = svx*real(mxyzp1*nvp)
      svyx = svy*real(mxyzp1*nvp)
      svzx = svz*real(mxyzp1*nvp)
! zero out distribution
      do 20 k = 1, mxyzp1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
      sfv(j,2,k) = 0.0
      sfv(j,3,k) = 0.0
   10 continue
      sfv(nmv21+1,1,k) = fv(nmv21+1,1)
      sfv(nmv21+1,2,k) = fv(nmv21+1,2)
      sfv(nmv21+1,3,k) = fv(nmv21+1,3)
   20 continue
! count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,nppp,nvx,nvy,nvz,vx,vy,vz,sanp,ssumvx,ssumvy,ssumvz, &
!$OMP& ssumvx2,ssumvy2,ssumvz2)                                         &
!$OMP& REDUCTION(+:anp) REDUCTION(+:sumvx) REDUCTION(+:sumvy)           &
!$OMP& REDUCTION(+:sumvz) REDUCTION(+:sumvx2) REDUCTION(+:sumvy2)       &
!$OMP& REDUCTION(+:sumvz2)
      do 40 k = 1, mxyzp1
      nppp = kpic(k)
      sanp = 0.0d0
      ssumvx = 0.0d0
      ssumvy = 0.0d0
      ssumvz = 0.0d0
      ssumvx2 = 0.0d0
      ssumvy2 = 0.0d0
      ssumvz2 = 0.0d0
      do 30 j = 1, nppp
      sanp = sanp + 1.0d0
      vx = ppart(4,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      vy = ppart(5,j,k)
      nvy = vy*svy + anmv
      ssumvy = ssumvy + vy
      ssumvy2 = ssumvy2 + vy*vy
      if ((nvy.ge.1).and.(nvy.le.nmv21)) sfv(nvy,2,k) = sfv(nvy,2,k)+1.0
      vz = ppart(6,j,k)
      nvz = vz*svz + anmv
      ssumvz = ssumvz + vz
      ssumvz2 = ssumvz2 + vz*vz
      if ((nvz.ge.1).and.(nvz.le.nmv21)) sfv(nvz,3,k) = sfv(nvz,3,k)+1.0
   30 continue
! calculate global sums
      anp = anp + sanp
      sumvx = sumvx + ssumvx
      sumvy = sumvy + ssumvy
      sumvz = sumvz + ssumvz
      sumvx2 = sumvx2 + ssumvx2
      sumvy2 = sumvy2 + ssumvy2
      sumvz2 = sumvz2 + ssumvz2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 50 k = 1, mxyzp1
      sum1 = sum1 + sfv(j,1,k)
      sum2 = sum2 + sfv(j,2,k)
      sum3 = sum3 + sfv(j,3,k)
   50 continue
      fv(j,1) = sum1
      fv(j,2) = sum2
      fv(j,3) = sum3
   60 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      sum7(7) = anp
      call PPDSUM(sum7,work7,7)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      anp = sum7(7)
! calculate global velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(2,1) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(3,1) = sumvz
      fvm(3,2) = dsqrt(sumvz2*anp - sumvz**2)
! count number of particles in global distribution
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      do 70 j = 1, nmv21
      sumvx = sumvx + fv(j,1)
      sumvy = sumvy + fv(j,2)
      sumvz = sumvz + fv(j,3)
   70 continue
! calculate entropy
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 90 k = 1, mxyzp1
      do 80 j = 1, nmv21
      if (sfv(j,1,k).gt.0.0) then
         sumvx2 = sumvx2 + sfv(j,1,k)*dlog(dble(sfv(j,1,k)*svxx))
      endif
      if (sfv(j,2,k).gt.0.0) then
         sumvy2 = sumvy2 + sfv(j,2,k)*dlog(dble(sfv(j,2,k)*svyx))
      endif
      if (sfv(j,3,k).gt.0.0) then
         sumvz2 = sumvz2 + sfv(j,3,k)*dlog(dble(sfv(j,3,k)*svzx))
      endif
   80 continue
   90 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      call PPDSUM(sum7,work7,6)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(1,3) = sumvx
      fvm(2,3) = sumvy
      fvm(3,3) = sumvz
      return
      end
!-----------------------------------------------------------------------
      subroutine PERPDIST32(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyzp1, &
     &nmv,nmvf)
! for 3d code, this subroutine calculates 3d energy distribution,
! for relativistic particles, with 2D domain decomposition
! the function calculated is of the form g*exp(-e/vth**2), where
! e = (gamma-1)*(c*c) is the kinetic energy per mass, and where
! vth = sqrt(KT/m).  Note vth is a momentum/mass and can be > c
! gamma = sqrt(1 + (p*p)/(c*c)), where p = is the momentum per mass
! g = p*p/(de/dp) = p*gamma
! one can express this quantity g as a function of e as follows:
! e = (p*p)/(gamma+1) => p = sqrt((gamma+1)*e), and gamma = 1 + e/(c*c)
! particles stored in segmented array
! input: all except wk, output: fv, sfv, wk
! ppart(4,n,m) = momentum px of particle n in tile m
! ppart(5,n,m) = momentum py of particle n in tile m
! ppart(6,n,m) = momentum pz of particle n in tile m
! fv = global distribution function, number of particles in each energy
! range, summed over tiles.  maximum energy (used for scaling) is
! contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! energy range
! ci = reciprocal of velocity of light
! wk = total energy contained in distribution
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx+1 and 
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! nmvf = first dimension of fv
! the number of energy bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mxyzp1, nmv, nmvf
      real ci, wk
      real ppart, fv, sfv
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fv(nmvf,3), sfv(nmvf,3,mxyzp1)
      integer kpic
      dimension kpic(mxyzp1)
! local data
      integer j, k, nppp, nmv21, nvx
      real ci2, anmv, svx, px, py, pz, p2
      double precision sumpx, ssumpx, sum1
      double precision dsum1, dwork1
      dimension dsum1(1), dwork1(1)
      ci2 = ci*ci
! energy scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
! zero out distribution
      do 20 k = 1, mxyzp1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
   10 continue
      sfv(nmv21+1,1,k) = fv(nmv21+1,1)
   20 continue
! count particles in each energy region
      anmv = 1.0
      sumpx = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,nppp,nvx,px,py,pz,p2,ssumpx) REDUCTION(+:sumpx) 
      do 40 k = 1, mxyzp1
      nppp = kpic(k)
      ssumpx = 0.0d0
      do 30 j = 1, nppp
      px = ppart(4,j,k)
      py = ppart(5,j,k)
      pz = ppart(6,j,k)
      p2 = px*px + py*py + pz*pz
      px = p2/(1.0 + sqrt(1.0 + p2*ci2))
      nvx = px*svx + anmv
      ssumpx = ssumpx + px
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
   30 continue
! calculate global sums
      sumpx = sumpx + ssumpx
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      do 50 k = 1, mxyzp1
      sum1 = sum1 + sfv(j,1,k)
   50 continue
      fv(j,1) = sum1
   60 continue
      dsum1(1) = sumpx
      call PPDSUM(dsum1,dwork1,1)
      sumpx = dsum1(1)
! return energy
      wk = sumpx
      return
      end
!-----------------------------------------------------------------------
      subroutine PPVBDIST32(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp,    &
     &nppmx,mxyzp1,nmv,nmvf)
! for 3d code, this subroutine calculates 3d velocity distribution
! and velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! with 2D domain decomposition
! particles stored in segmented array
! input: all except fvm, output: fv, sfv, fvm
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = velocity vy of particle n in tile m
! ppart(6,n,m) = velocity vz of particle n in tile m
! fv = global distribution function, number of particles in each
! velocity range, summed over tiles. maximum velocity (used for scaling)
! is contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! velocity range
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx+1 and 
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mxyzp1, nmv, nmvf
      real omx, omy, omz
      real ppart, fv, sfv, fvm
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fv(nmvf,3), sfv(nmvf,3,mxyzp1), fvm(3,3)
      integer kpic
      dimension kpic(mxyzp1)
! local data
      integer j, k, ndir, nppp, nmv21, nvx, nvz
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real anmv, svx, svz
      double precision sumvx, sumvz, sumvx2, sumvz2, anp
      double precision ssumvx, ssumvz, ssumvx2, ssumvz2, sanp
      double precision sum1, sum2
      double precision sum5, work5
      dimension sum5(5), work5(5)
! find rotation to convert to cylindrical co-ordinates
      at1 = sqrt(omx*omx + omy*omy + omz*omz)
! no rotation if zero B field
      if (at1.eq.0.0) then
         ox = 0.0
         oy = 0.0
         oz = 1.0
! create rotation vectors
      else
! first create unit vector in B direction
         at1 = 1.0/at1
         ox = omx*at1
         oy = omy*at1
         oz = omz*at1
      endif
! then create unit vector in first perpendicular direction
! find direction with smallest component of B
      ndir = 1
      at1 = abs(omx)
      at2 = abs(omy)
      if (at2.le.at1) then
         ndir = 2
         at1 = at2
      endif
      if (abs(omz).lt.at1) ndir = 3
! take the cross product of that direction with B
! vpr1 = x cross B
      if (ndir.eq.1) then
         at1 = 1.0/sqrt(oy*oy + oz*oz)
         px = 0.0
         py = -oz*at1
         pz = oy*at1
! vpr1 = y cross B
      else if (ndir.eq.2) then
         at1 = 1.0/sqrt(ox*ox + oz*oz)
         px = oz*at1
         py = 0.0
         pz = -ox*at1
! vpr1 = z cross B
      else if (ndir.eq.3) then
         at1 = 1.0/sqrt(ox*ox + oy*oy)
         px = -oy*at1
         py = ox*at1
         pz = 0.0
      endif
! finally create unit vector in second perpendicular direction
! vpr2 = B cross vpr1
      qx = oy*pz - oz*py
      qy = oz*px - ox*pz
      qz = ox*py - oy*px
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svz = anmv/fv(nmv21+1,2)
! zero out distribution
      do 20 k = 1, mxyzp1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
      sfv(j,2,k) = 0.0
   10 continue
      sfv(nmv21+1,1,k) = fv(nmv21+1,1)
      sfv(nmv21+1,2,k) = fv(nmv21+1,2)
   20 continue
! count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvz2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,nppp,nvx,nvz,vx,vy,vz,at1,at2,sanp,ssumvx,ssumvz,    &
!$OMP& ssumvx2,ssumvz2)                                                 &
!$OMP& REDUCTION(+:anp) REDUCTION(+:sumvx) REDUCTION(+:sumvz)           &
!$OMP& REDUCTION(+:sumvx2) REDUCTION(+:sumvz2)
      do 40 k = 1, mxyzp1
      nppp = kpic(k)
      sanp = 0.0d0
      ssumvx = 0.0d0
      ssumvz = 0.0d0
      ssumvx2 = 0.0d0
      ssumvz2 = 0.0d0
      do 30 j = 1, nppp
      sanp = sanp + 1.0d0
      vx = ppart(4,j,k)
      vy = ppart(5,j,k)
      vz = ppart(6,j,k)
! vperp1 co-ordinate
      at1 = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
      at2 = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
      vz = vx*ox + vy*oy + vz*oz
      vx = sqrt(at1*at1 + at2*at2)
      nvx = vx*svx + anmv
      sumvx = sumvx + vx
      sumvx2 = sumvx2 + vx**2
      nvz = vz*svz + anmv
      sumvz = sumvz + vz
      sumvz2 = sumvz2 + vz**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      if ((nvz.ge.1).and.(nvz.le.nmv21)) sfv(nvz,2,k) = sfv(nvz,2,k)+1.0
   30 continue
! calculate global sums
      anp = anp + sanp
      sumvx = sumvx + ssumvx
      sumvz = sumvz + ssumvz
      sumvx2 = sumvx2 + ssumvx2
      sumvz2 = sumvz2 + ssumvz2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 50 k = 1, mxyzp1
      sum1 = sum1 + sfv(j,1,k)
      sum2 = sum2 + sfv(j,2,k)
   50 continue
      fv(j,1) = sum1
      fv(j,2) = sum2
      fv(j,3) = 0.0
   60 continue
      sum5(1) = sumvx
      sum5(2) = sumvz
      sum5(3) = sumvx2
      sum5(4) = sumvz2
      sum5(5) = anp
      call PPDSUM(sum5,work5,5)
      sumvx = sum5(1)
      sumvz = sum5(2)
      sumvx2 = sum5(3)
      sumvz2 = sum5(4)
      anp = sum5(5)
! calculate global velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      fvm(1,3) = 0.0
      sumvz = sumvz*anp
      fvm(2,1) = sumvz
      fvm(2,2) = dsqrt(sumvz2*anp - sumvz**2)
      fvm(2,3) = 0.0
      fvm(3,1) = 0.0
      fvm(3,2) = 0.0
      fvm(3,3) = 0.0
      return
      end
!-----------------------------------------------------------------------
      subroutine PPVSDIST32(ppart,kpic,fvs,noff,nmv,mvx,mvy,mvz,nxb,nyb,&
     &nzb,idimp,nppmx,mxyzp1,nmvf,nybmx,idds)
! for 3d code, this subroutine calculates 3d velocity distribution,
! in different regions of space, with 2D domain decomposition
! particles stored in segmented array
! input: all except fvs, output: fvs
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = position z of particle n in tile m
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = velocity vy of particle n in tile m
! ppart(6,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! fvs = spatially resolved distribution function, number of particles in
! each velocity and spatial range.  maximum velocity (used for scaling)
! is contained in last element of first dimension of fvs
! noff(1:2) = lowermost global gridpoint in y/z in particle partition
! nmv = number of segments in v for velocity distribution
! mvx/mvy/mvz = number of grids in x/y/z for phase space aggregation
! nxb/nyb/nzb = number of segments in x/y/z for velocity distribution
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx+1 and 
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! nmvf = first dimension of fvs
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
! nybmx = maximum number of segments in y for velocity distribution
! idds = dimensionality of domain decomposition = 2
      implicit none
      integer nmv, mvx, mvy, mvz, nxb, nyb, nzb, idimp, nppmx, mxyzp1
      integer nmvf, nybmx, idds
      real ppart, fvs
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fvs(nmvf,3,nxb,nybmx+1,nzb+1)
      integer kpic, noff
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer j, k, ns, ms, nppp, nmv21, nn, mm, ll, nvx, nvy, nvz
      real anmv, anoff, amoff, svx, svy, svz, at1, at2, at3
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      ns = noff(1) - mvy*(noff(1)/mvy)
      ms = noff(2) - mvz*(noff(2)/mvz)
      anmv = real(nmv)
      anoff = real(noff(1) - ns)
      amoff = real(noff(2) - ms)
      svx = anmv/fvs(nmv21+1,1,1,1,1)
      svy = anmv/fvs(nmv21+1,2,1,1,1)
      svz = anmv/fvs(nmv21+1,3,1,1,1)
! spatial scaling
      at1 = 1.0/real(mvx)
      at2 = 1.0/real(mvy)
      at3 = 1.0/real(mvz)
! zero out distribution
      do 40 ll = 1, nzb+1
      do 30 mm = 1, nyb+1
      do 20 nn = 1, nxb
      do 10 j = 1, nmv21
      fvs(j,1,nn,mm,ll) = 0.0
      fvs(j,2,nn,mm,ll) = 0.0
      fvs(j,3,nn,mm,ll) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
! count particles in each velocity region
      anmv = anmv + 1.5
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,nn,mm,ll,nvx,nvy,nvz)
      do 60 k = 1, mxyzp1
      nppp = kpic(k)
      do 50 j = 1, nppp
      nn = ppart(1,j,k)*at1 + 1.0
      mm = (ppart(2,j,k) - anoff)*at2 + 1.0
      ll = (ppart(3,j,k) - amoff)*at3 + 1.0
      nvx = ppart(4,j,k)*svx + anmv
      nvy = ppart(5,j,k)*svy + anmv
      nvz = ppart(6,j,k)*svz + anmv
      if ((nvx.ge.1).and.(nvx.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvx,1,nn,mm,ll) = fvs(nvx,1,nn,mm,ll) + 1.0
      endif
      if ((nvy.ge.1).and.(nvy.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvy,2,nn,mm,ll) = fvs(nvy,2,nn,mm,ll) + 1.0
      endif
      if ((nvz.ge.1).and.(nvz.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvz,3,nn,mm,ll) = fvs(nvz,3,nn,mm,ll) + 1.0
      endif
   50 continue
   60 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDIST32(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
! for 3d code, this subroutine calculates 3d velocity distribution
! and velocity moments, with 2D domain decomposition
! input: all, output: fv, fvm
! part(4,n) = velocity vx of particle n
! part(5,n) = velocity vy of particle n
! part(6,n) = velocity vz of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! entropy for i-th dimension is contained in fvm(i,3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
! is uniform in space and distributions in each dimension are
! independent.
! npp = number of particles in partition
! nvp = number of real or virtual processors
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer npp, nvp, idimp, npmax, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,npmax)
      dimension fv(nmvf,3), fvm(3,3)
! local data
      integer j, nmv21, nvx, nvy, nvz
      real anmv, svx, svy, svz, svxx, svyx, svzx
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      double precision sum7, work7
      dimension sum7(7), work7(7)
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svy = anmv/fv(nmv21+1,2)
      svz = anmv/fv(nmv21+1,3)
! normalization constant for entropy
      svxx = svx*real(nvp)
      svyx = svy*real(nvp)
      svzx = svz*real(nvp)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j,1) = 0.0
      fv(j,2) = 0.0
      fv(j,3) = 0.0
   10 continue
! count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, npp
      anp = anp + 1.0d0
      nvx = part(4,j)*svx + anmv
      sumvx = sumvx + part(4,j)
      sumvx2 = sumvx2 + part(4,j)**2
      nvy = part(5,j)*svy + anmv
      sumvy = sumvy + part(5,j)
      sumvy2 = sumvy2 + part(5,j)**2
      nvz = part(6,j)*svz + anmv
      sumvz = sumvz + part(6,j)
      sumvz2 = sumvz2 + part(6,j)**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx,1) = fv(nvx,1) + 1.0
      if ((nvy.ge.1).and.(nvy.le.nmv21)) fv(nvy,2) = fv(nvy,2) + 1.0
      if ((nvz.ge.1).and.(nvz.le.nmv21)) fv(nvz,3) = fv(nvz,3) + 1.0
   20 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      sum7(7) = anp
      call PPDSUM(sum7,work7,7)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      anp = sum7(7)
! calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(2,1) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(3,1) = sumvz
      fvm(3,2) = dsqrt(sumvz2*anp - sumvz**2)
! calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 30 j = 1, nmv21
      if (fv(j,1).gt.0.0) then
         sumvx = sumvx + fv(j,1)
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svxx))
      endif
      if (fv(j,2).gt.0.0) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svyx))
      endif
      if (fv(j,3).gt.0.0) then
         sumvz = sumvz + fv(j,3)
         sumvz2 = sumvz2 + fv(j,3)*dlog(dble(fv(j,3)*svzx))
      endif
   30 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      call PPDSUM(sum7,work7,6)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(1,3) = sumvx
      fvm(2,3) = sumvy
      fvm(3,3) = sumvz
      return
      end
!-----------------------------------------------------------------------
      subroutine PVBDIST32(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,nmv, &
     &nmvf)
! for 3d code, this subroutine calculates 3d velocity distribution and
! velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! with 2D domain decomposition
! input: all except fvm, output: fv, fvm
! part(4,n) = velocity vx of particle n
! part(5,n) = velocity vy of particle n
! part(6,n) = velocity vz of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! npp = number of particles in partition
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer npp, idimp, npmax, nmv, nmvf
      real omx, omy, omz
      real part, fv, fvm
      dimension part(idimp,npmax)
      dimension fv(nmvf,2), fvm(3,2)
! local data
      integer j, ndir, nmv21, nvx, nvz
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real anmv, svx, svz
      double precision sumvx,sumvz, sumvx2, sumvz2, anp
      double precision sum5, work5
      dimension sum5(5), work5(5)
! find rotation to convert to cylindrical co-ordinates
      at1 = sqrt(omx*omx + omy*omy + omz*omz)
! no rotation if zero B field
      if (at1.eq.0.0) then
         ox = 0.0
         oy = 0.0
         oz = 1.0
! create rotation vectors
      else
! first create unit vector in B direction
         at1 = 1.0/at1
         ox = omx*at1
         oy = omy*at1
         oz = omz*at1
      endif
! then create unit vector in first perpendicular direction
! find direction with smallest component of B
      ndir = 1
      at1 = abs(omx)
      at2 = abs(omy)
      if (at2.le.at1) then
         ndir = 2
         at1 = at2
      endif
      if (abs(omz).lt.at1) ndir = 3
! take the cross product of that direction with B
! vpr1 = x cross B
      if (ndir.eq.1) then
         at1 = 1.0/sqrt(oy*oy + oz*oz)
         px = 0.0
         py = -oz*at1
         pz = oy*at1
! vpr1 = y cross B
      else if (ndir.eq.2) then
         at1 = 1.0/sqrt(ox*ox + oz*oz)
         px = oz*at1
         py = 0.0
         pz = -ox*at1
! vpr1 = z cross B
      else if (ndir.eq.3) then
         at1 = 1.0/sqrt(ox*ox + oy*oy)
         px = -oy*at1
         py = ox*at1
         pz = 0.0
      endif
! finally create unit vector in second perpendicular direction
! vpr2 = B cross vpr1
      qx = oy*pz - oz*py
      qy = oz*px - ox*pz
      qz = ox*py - oy*px
! velocity scaling
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svz = anmv/fv(nmv21+1,2)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j,1) = 0.0
      fv(j,2) = 0.0
   10 continue
! count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, npp
      anp = anp + 1.0d0
      vx = part(4,j)
      vy = part(5,j)
      vz = part(6,j)
! vperp1 co-ordinate
      at1 = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
      at2 = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
      vz = vx*ox + vy*oy + vz*oz
      vx = sqrt(at1*at1 + at2*at2)
      nvx = vx*svx + anmv
      sumvx = sumvx + vx
      sumvx2 = sumvx2 + vx**2
      nvz = vz*svz + anmv
      sumvz = sumvz + vz
      sumvz2 = sumvz2 + vz**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx,1) = fv(nvx,1) + 1.0
      if ((nvz.ge.1).and.(nvz.le.nmv21)) fv(nvz,2) = fv(nvz,2) + 1.0
   20 continue
      sum5(1) = sumvx
      sum5(2) = sumvz
      sum5(3) = sumvx2
      sum5(4) = sumvz2
      sum5(5) = anp
      call PPDSUM(sum5,work5,5)
      sumvx = sum5(1)
      sumvz = sum5(2)
      sumvx2 = sum5(3)
      sumvz2 = sum5(4)
      anp = sum5(5)
! calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvz = sumvz*anp
      fvm(2,1) = sumvz
      fvm(2,2) = dsqrt(sumvz2*anp - sumvz**2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PPROFX32L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,mz&
     &,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 243 flops/particle, 114 loads, 108 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, npro, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real ppart, fms
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z, w
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfms(nprd,mx+1,my+1,mz+1), sg(14)
      mxyp1 = mx1*myp1
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sfms,&
!$OMP& sg)
      do 280 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 40 k = 1, mz+1
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
! loop over particles in tile
      do 60 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit fluid moments
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         x = vx*vx
         sg(5) = x
         sg(6) = vx*vy
         sg(7) = vx*vz
         y = vy*vy
         w = x + y
         sg(8) = y
         sg(9) = vy*vz
         z = vz*vz
         sg(10) = z
         w = 0.5*(w + z)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 50 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
   50 continue
   60 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 100 k = 2, ll
      do 90 j = 2, mm
      do 80 i = 2, nn
      do 70 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
   70 continue
   80 continue
   90 continue
  100 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 140 j = 2, mm
      do 130 i = 2, nn
      do 110 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  110 continue
      if (lm > mz) then
         do 120 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  120    continue
      endif
  130 continue
  140 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 210 k = 1, ll
      do 170 i = 2, nn
      do 150 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  150 continue
      if (mm > my) then
         do 160 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  160    continue
      endif
  170 continue
      do 200 j = 1, mm
      do 180 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  180 continue
      if (nm > mx) then
         do 190 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  190    continue
      endif
  200 continue
  210 continue
      if (lm > mz) then
         do 240 i = 2, nn
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  220    continue
         if (mm > my) then
            do 230 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  230       continue
         endif
  240    continue
         do 270 j = 1, mm
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  250    continue
         if (nm > mx) then
            do 260 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  260       continue
         endif
  270    continue
      endif
  280 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PRPROFX32L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx, &
     &my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux
! for relativistic particles
! assumes particle positions and momenta are at the same time level
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 253 flops/particle, 2 divides, 1 sqrt, 114 loads, 108 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci =  gami*p2/(1.0 + gami),
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = position z of particle n in partition in tile m
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 6
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, npro, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real ci
      real ppart, fms
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real ci2, x, y, z, w, p2, gami
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfms(nprd,mx+1,my+1,mz+1), sg(14)
      mxyp1 = mx1*myp1
      ci2 = ci*ci
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,p2,  &
!$OMP& gami,sfms,sg)
      do 280 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! zero out local accumulator
      do 40 k = 1, mz+1
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
! loop over particles in tile
      do 60 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
! find inverse gamma
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! deposit fluid moments
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         sg(5) = vx*vx
         sg(6) = vx*vy
         sg(7) = vx*vz
         sg(8) = vy*vy
         sg(9) = vy*vz
         sg(10) = vz*vz
         w = gami*p2/(1.0 + gami)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 50 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
   50 continue
   60 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 100 k = 2, ll
      do 90 j = 2, mm
      do 80 i = 2, nn
      do 70 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
   70 continue
   80 continue
   90 continue
  100 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 140 j = 2, mm
      do 130 i = 2, nn
      do 110 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  110 continue
      if (lm > mz) then
         do 120 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  120    continue
      endif
  130 continue
  140 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 210 k = 1, ll
      do 170 i = 2, nn
      do 150 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  150 continue
      if (mm > my) then
         do 160 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  160    continue
      endif
  170 continue
      do 200 j = 1, mm
      do 180 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  180 continue
      if (nm > mx) then
         do 190 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  190    continue
      endif
  200 continue
  210 continue
      if (lm > mz) then
         do 240 i = 2, nn
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  220    continue
         if (mm > my) then
            do 230 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  230       continue
         endif
  240    continue
         do 270 j = 1, mm
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  250    continue
         if (nm > mx) then
            do 260 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  260       continue
         endif
  270    continue
      endif
  280 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,idimp, &
     &nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 320 flops/particle, 144 loads, 108 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
! z(t+dt) = z(t) + vz(t+dt/2)*dt
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, npro, nx, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real qbm, dt
      real ppart, fxyz, fms
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real x, y, z, w, vx, vy, vz
      real sfxyz, sfms, sg
!     dimension sfxyz(3,MXV,MYV,MZV), sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfxyz(3,mx+1,my+1,mz+1), sfms(nprd,mx+1,my+1,mz+1)
      dimension sg(14)
      mxyp1 = mx1*myp1
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,     &
!$OMP& sfxyz,sfms,sg)
      do 310 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
! zero out local accumulators
      do 70 k = 1, mz+1
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      do 40 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   40 continue
   50 continue
   60 continue
   70 continue
! loop over particles in tile
      do 90 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(4,j,l) + x
      vy = ppart(5,j,l) + y
      vz = ppart(6,j,l) + z
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         x = vx*vx
         sg(5) = x
         sg(6) = vx*vy
         sg(7) = vx*vz
         y = vy*vy
         w = x + y
         sg(8) = y
         sg(9) = vy*vz
         z = vz*vz
         sg(10) = z
         w = 0.5*(w + z)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 80 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
   80 continue
   90 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
      do 100 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
  100 continue
  110 continue
  120 continue
  130 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 170 j = 2, mm
      do 160 i = 2, nn
      do 140 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  140 continue
      if (lm > mz) then
         do 150 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  150    continue
      endif
  160 continue
  170 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 200 i = 2, nn
      do 180 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  180 continue
      if (mm > my) then
         do 190 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  190    continue
      endif
  200 continue
      do 230 j = 1, mm
      do 210 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  210 continue
      if (nm > mx) then
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  220    continue
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 270 i = 2, nn
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  250    continue
         if (mm > my) then
            do 260 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  260       continue
         endif
  270    continue
         do 300 j = 1, mm
         do 280 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  280    continue
         if (nm > mx) then
            do 290 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  290       continue
         endif
  300    continue
      endif
  310 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGRPROFX32L(ppart,fxyz,fms,kpic,noff,nyzp,qbm,dt,ci,   &
     &idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux
! for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 328 flops/particle, 2 divides, 1 sqrt, 144 loads, 108 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci =  gami*p2/(1.0 + gami),
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
! pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
!                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
!                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
! fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
!                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
!                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, npro, nx, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real qbm, dt, ci
      real ppart, fxyz, fms
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, ci2, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real x, y, z, w, vx, vy, vz, p2, gami
      real sfxyz, sfms, sg
!     dimension sfxyz(3,MXV,MYV,MZV), sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfxyz(3,mx+1,my+1,mz+1), sfms(nprd,mx+1,my+1,mz+1)
      dimension sg(14)
      mxyp1 = mx1*myp1
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,p2,  &
!$OMP& gami,sfxyz,sfms,sg)
      do 310 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
! zero out local accumulators
      do 70 k = 1, mz+1
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      do 40 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   40 continue
   50 continue
   60 continue
   70 continue
! loop over particles in tile
      do 90 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     &
     &                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(4,j,l) + x
      vy = ppart(5,j,l) + y
      vz = ppart(6,j,l) + z
! update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate average velocity
      vx = gami*vx
      vy = gami*vy
      vz = gami*vz
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         sg(5) = vx*vx
         sg(6) = vx*vy
         sg(7) = vx*vz
         sg(8) = vy*vy
         sg(9) = vy*vz
         sg(10) = vz*vz
         w = gami*p2/(1.0 + gami)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 80 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
   80 continue
   90 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 130 k = 2, ll
      do 120 j = 2, mm
      do 110 i = 2, nn
      do 100 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
  100 continue
  110 continue
  120 continue
  130 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 170 j = 2, mm
      do 160 i = 2, nn
      do 140 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  140 continue
      if (lm > mz) then
         do 150 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  150    continue
      endif
  160 continue
  170 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 240 k = 1, ll
      do 200 i = 2, nn
      do 180 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  180 continue
      if (mm > my) then
         do 190 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  190    continue
      endif
  200 continue
      do 230 j = 1, mm
      do 210 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  210 continue
      if (nm > mx) then
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  220    continue
      endif
  230 continue
  240 continue
      if (lm > mz) then
         do 270 i = 2, nn
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  250    continue
         if (mm > my) then
            do 260 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  260       continue
         endif
  270    continue
         do 300 j = 1, mm
         do 280 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  280    continue
         if (nm > mx) then
            do 290 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  290       continue
         endif
  300    continue
      endif
  310 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt, &
     &idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,&
     &idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 421 flops/particle, 1 divide, 162 loads, 108 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
! omz = (q/m)*bz(x(t),y(t),z(t)).
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, npro, nx, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real qbm, dt
      real ppart, fxyz, bxyz, fms
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, w, vx, vy, vz
      real sfxyz, sbxyz, sfms, sg
!     dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
!     dimension sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
      dimension sfms(nprd,mx+1,my+1,mz+1), sg(14)
      mxyp1 = mx1*myp1
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,  &
!$OMP& oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,  &
!$OMP& rot5,rot6,rot7,rot8,rot9,sfxyz,sbxyz,sfms,sg)
      do 340 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 100 k = 1, mz+1
      do 90 j = 1, my+1
      do 80 i = 1, mx+1
      do 70 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
  100 continue
! loop over particles in tile
      do 120 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      acx = vx + x
      acy = vy + y
      acz = vz + z
! calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      x = (rot1*acx + rot2*acy + rot3*acz)*anorm + x
      y = (rot4*acx + rot5*acy + rot6*acz)*anorm + y
      z = (rot7*acx + rot8*acy + rot9*acz)*anorm + z
! calculate average velocity
      vx = 0.5*(x + vx)
      vy = 0.5*(y + vy)
      vz = 0.5*(z + vz)
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         x = vx*vx
         sg(5) = x
         sg(6) = vx*vy
         sg(7) = vx*vz
         y = vy*vy
         w = x + y
         sg(8) = y
         sg(9) = vy*vz
         z = vz*vz
         sg(10) = z
         w = 0.5*(w + z)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 110 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
  110 continue
  120 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 160 k = 2, ll
      do 150 j = 2, mm
      do 140 i = 2, nn
      do 130 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
  130 continue
  140 continue
  150 continue
  160 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 200 j = 2, mm
      do 190 i = 2, nn
      do 170 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  170 continue
      if (lm > mz) then
         do 180 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  180    continue
      endif
  190 continue
  200 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 270 k = 1, ll
      do 230 i = 2, nn
      do 210 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  210 continue
      if (mm > my) then
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  220    continue
      endif
  230 continue
      do 260 j = 1, mm
      do 240 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  240 continue
      if (nm > mx) then
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  250    continue
      endif
  260 continue
  270 continue
      if (lm > mz) then
         do 300 i = 2, nn
         do 280 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  280    continue
         if (mm > my) then
            do 290 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  290       continue
         endif
  300    continue
         do 330 j = 1, mm
         do 310 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  310    continue
         if (nm > mx) then
            do 320 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  320       continue
         endif
  330    continue
      endif
  340 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGRBPROFX32L(ppart,fxyz,bxyz,fms,kpic,noff,nyzp,qbm,dt,&
     &ci,idimp,nppmx,npro,nx,mx,my,mz,nprd,nxv,nypmx,nzpmx,mx1,myp1,    &
     &mxyzp1,idds)
! for 3d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux
! for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 2D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 449 flops/particle, 4 divides, 2 sqrt, 162 loads, 108 stores,
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m,l)=mci*(1.-dx)*(1.-dy)*(1.-dz)
! fms(i,n+1,m,l)=mci*dx*(1.-dy)*(1.-dz)
! fms(i,n,m+1,l)=mci*(1.-dx)*dy*(1.-dz)
! fms(i,n+1,m+1,l)=mci*dx*dy*(1.-dz)
! fms(i,n,m,l+1)=mci*(1.-dx)*(1.-dy)*dz
! fms(i,n+1,m,l+1)=mci*dx*(1.-dy)*dz
! fms(i,n,m+1,l+1)=mci*(1.-dx)*dy*dz
! fms(i,n+1,m+1,l+1)=mci*dx*dy*dz
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci =  gami*p2/(1.0 + gami),
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami.
! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
! bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
!                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
!           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
!                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
! where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
! similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = position z of particle n in partition in tile m at t
! ppart(4,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(6,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
! fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,lll)
! fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
! that is, convolution of electric field over particle shape
! bxyz(1,j,k,l) = x component of magnetic field at grid (j,kk,ll)
! bxyz(2,j,k,l) = y component of magnetic field at grid (j,kk,ll)
! bxyz(3,j,k,l) = z component of magnetic field at grid (j,kk,ll)
! that is, the convolution of magnetic field over particle shape
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! fms(i,j,k,l) = ith component of fluid moments at grid point j,kk,ll
! where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
! kpic(l) = number of particles in tile l
! noff(1) = lowermost global gridpoint in y in particle partition
! noff(2) = backmost global gridpoint in z in particle partition
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! qbm = particle charge/mass ratio
! dt = time interval between successive force calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my/mz = number of grids in sorting cell in x/y/z
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mxyzp1 = mx1*myp1*mzp1,
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idds = dimensionality of domain decomposition
      implicit none
      integer idimp, nppmx, npro, nx, mx, my, mz
      integer nprd, nxv, nypmx, nzpmx, mx1, myp1, mxyzp1, idds
      real qbm, dt, ci
      real ppart, fxyz, bxyz, fms
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension fxyz(3,nxv,nypmx,nzpmx), bxyz(3,nxv,nypmx,nzpmx)
      dimension fms(nprd,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
! local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp, npr
      integer i, j, k, l, ii, mnoff, lnoff, nn, mm, ll, nm, lm
      real qtmh, ci2, gami, qtmg, dxp, dyp, dzp
      real amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, w, vx, vy, vz, p2
      real sfxyz, sbxyz, sfms, sg
!     dimension sfxyz(3,MXV,MYV,MZV), sbxyz(3,MXV,MYV,MZV)
!     dimension sfms(nprd,MXV,MYV,MZV), sg(14)
      dimension sfxyz(3,mx+1,my+1,mz+1), sbxyz(3,mx+1,my+1,mz+1)
      dimension sfms(nprd,mx+1,my+1,mz+1), sg(14)
      mxyp1 = mx1*myp1
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,l,ii,noffp,moffp,loffp,nppp,nn,mm,ll,mnoff,lnoff,  &
!$OMP& nm,lm,x,y,z,w,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,  &
!$OMP& oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,  &
!$OMP& rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,sfxyz,sbxyz,sfms,sg)
      do 340 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
! load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      do 60 k = 1, min(mz,nyzp(2)-loffp)+1
      do 50 j = 1, min(my,nyzp(1)-moffp)+1
      do 40 i = 1, min(mx,nx-noffp)+1
      sbxyz(1,i,j,k) = bxyz(1,i+noffp,j+moffp,k+loffp)
      sbxyz(2,i,j,k) = bxyz(2,i+noffp,j+moffp,k+loffp)
      sbxyz(3,i,j,k) = bxyz(3,i+noffp,j+moffp,k+loffp)
   40 continue
   50 continue
   60 continue
! zero out local accumulators
      do 100 k = 1, mz+1
      do 90 j = 1, my+1
      do 80 i = 1, mx+1
      do 70 ii = 1, npr
      sfms(ii,i,j,k) = 0.0
   70 continue
   80 continue
   90 continue
  100 continue
! loop over particles in tile
      do 120 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
! find electric field
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sfxyz(3,nn+1,mm+1,ll))
      acx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      acy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      acz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(acx + dyp*sfxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(acy + dyp*sfxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(acz + dyp*sfxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sfxyz(3,nn+1,mm+1,ll+1))
! find magnetic field
      ox = amx*sbxyz(1,nn,mm,ll) + amy*sbxyz(1,nn+1,mm,ll)
      oy = amx*sbxyz(2,nn,mm,ll) + amy*sbxyz(2,nn+1,mm,ll)
      oz = amx*sbxyz(3,nn,mm,ll) + amy*sbxyz(3,nn+1,mm,ll)
      ox = amz*(ox + dyp*sbxyz(1,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(1,nn+1,mm+1,ll))
      oy = amz*(oy + dyp*sbxyz(2,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(2,nn+1,mm+1,ll))
      oz = amz*(oz + dyp*sbxyz(3,nn,mm+1,ll)                            &
     &             + dx1*sbxyz(3,nn+1,mm+1,ll))
      acx = amx*sbxyz(1,nn,mm,ll+1) + amy*sbxyz(1,nn+1,mm,ll+1)
      acy = amx*sbxyz(2,nn,mm,ll+1) + amy*sbxyz(2,nn+1,mm,ll+1)
      acz = amx*sbxyz(3,nn,mm,ll+1) + amy*sbxyz(3,nn+1,mm,ll+1)
      ox = ox + dzp*(acx + dyp*sbxyz(1,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(1,nn+1,mm+1,ll+1))
      oy = oy + dzp*(acy + dyp*sbxyz(2,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(2,nn+1,mm+1,ll+1))
      oz = oz + dzp*(acz + dyp*sbxyz(3,nn,mm+1,ll+1)                    &
     &                   + dx1*sbxyz(3,nn+1,mm+1,ll+1))
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      acx = vx + x
      acy = vy + y
      acz = vz + z
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new momentum
      x = (rot1*acx + rot2*acy + rot3*acz)*anorm + x
      y = (rot4*acx + rot5*acy + rot6*acz)*anorm + y
      z = (rot7*acx + rot8*acy + rot9*acz)*anorm + z
! calculate average momentum
      vx = 0.5*(x + vx)
      vy = 0.5*(y + vy)
      vz = 0.5*(z + vz)
! find inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate average velocity
      vx = gami*vx
      vy = gami*vy
      vz = gami*vz
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         sg(5) = vx*vx
         sg(6) = vx*vy
         sg(7) = vx*vz
         sg(8) = vy*vy
         sg(9) = vy*vz
         sg(10) = vz*vz
         w = gami*p2/(1.0 + gami)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amz
      y = amy*amz
      z = dyp*amz
      w = dx1*amz
      dx = amx*dzp
      dy = amy*dzp
      dz = dyp*dzp
      dx1 = dx1*dzp
      do 110 ii = 1, npr
      sfms(ii,nn,mm,ll) = sfms(ii,nn,mm,ll) + sg(ii)*x
      sfms(ii,nn+1,mm,ll) = sfms(ii,nn+1,mm,ll) + sg(ii)*y
      sfms(ii,nn,mm+1,ll) = sfms(ii,nn,mm+1,ll) + sg(ii)*z
      sfms(ii,nn+1,mm+1,ll) = sfms(ii,nn+1,mm+1,ll) + sg(ii)*w
      sfms(ii,nn,mm,ll+1) = sfms(ii,nn,mm,ll+1) + sg(ii)*dx
      sfms(ii,nn+1,mm,ll+1) = sfms(ii,nn+1,mm,ll+1) + sg(ii)*dy
      sfms(ii,nn,mm+1,ll+1) = sfms(ii,nn,mm+1,ll+1) + sg(ii)*dz
      sfms(ii,nn+1,mm+1,ll+1) = sfms(ii,nn+1,mm+1,ll+1) + sg(ii)*dx1
  110 continue
  120 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 160 k = 2, ll
      do 150 j = 2, mm
      do 140 i = 2, nn
      do 130 ii = 1, npr
      fms(ii,i+noffp,j+moffp,k+loffp) = fms(ii,i+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,i,j,k)
  130 continue
  140 continue
  150 continue
  160 continue
! deposit fluid moments to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 200 j = 2, mm
      do 190 i = 2, nn
      do 170 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,j+moffp,1+loffp) = fms(ii,i+noffp,j+moffp,1+loffp) &
     &+ sfms(ii,i,j,1)
  170 continue
      if (lm > mz) then
         do 180 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,i+noffp,j+moffp,lm+loffp) + sfms(ii,i,j,lm)
  180    continue
      endif
  190 continue
  200 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 270 k = 1, ll
      do 230 i = 2, nn
      do 210 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp,k+loffp) = fms(ii,i+noffp,1+moffp,k+loffp) &
     &+ sfms(ii,i,1,k)
  210 continue
      if (mm > my) then
         do 220 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp,k+loffp) =                             &
     &   fms(ii,i+noffp,mm+moffp,k+loffp) + sfms(ii,i,mm,k)
  220    continue
      endif
  230 continue
      do 260 j = 1, mm
      do 240 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp,k+loffp) = fms(ii,1+noffp,j+moffp,k+loffp) &
     &+ sfms(ii,1,j,k)
  240 continue
      if (nm > mx) then
         do 250 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nm+noffp,j+moffp,k+loffp) =                             &
     &   fms(ii,nm+noffp,j+moffp,k+loffp) + sfms(ii,nm,j,k)
  250    continue
      endif
  260 continue
  270 continue
      if (lm > mz) then
         do 300 i = 2, nn
         do 280 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,1+moffp,lm+loffp) =                             &
     &fms(ii,i+noffp,1+moffp,lm+loffp) + sfms(ii,i,1,lm)
  280    continue
         if (mm > my) then
            do 290 ii = 1, npr
!$OMP ATOMIC
            fms(ii,i+noffp,mm+moffp,lm+loffp) =                         &
     &      fms(ii,i+noffp,mm+moffp,lm+loffp) + sfms(ii,i,mm,lm)
  290       continue
         endif
  300    continue
         do 330 j = 1, mm
         do 310 ii = 1, npr
!$OMP ATOMIC
         fms(ii,1+noffp,j+moffp,lm+loffp) =                             &
     &   fms(ii,1+noffp,j+moffp,lm+loffp)+ sfms(ii,1,j,lm)
  310    continue
         if (nm > mx) then
            do 320 ii = 1, npr
!$OMP ATOMIC
            fms(ii,nm+noffp,j+moffp,lm+loffp) =                         &
     &      fms(ii,nm+noffp,j+moffp,lm+loffp) + sfms(ii,nm,j,lm)
  320       continue
         endif
  330    continue
      endif
  340 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FLUIDQS3(fms,npro,nx,ny,nz,kstrt,nvpy,kyp,kzp,nprd,nxv,&
     &nypmx,nzpmx)
! for 3d code, this subroutine calculates fluid quantities from fluid
! moments: density, velocity field, pressure tensor, energy, heat flux
! assumes guard cells have been added
! for distributed data, with 2D spatial decomposition
! OpenMP version
! 42 flops/grid, 1 divide, 26 loads, 12 stores,
! if all profiles calculated
! input: all, output: fms
! fluid quantities are calculated as follows:
! fms(1,:) = density n is unchanged
! fms(2:4,:) = velocity field u, calculated from momentum and density n:
! u = fms(2:4,:)/n where n /= 0.0
! fms(5:10) = pressure tensor P, calculated from momentum flux, velocity
! field u and density n: P = fms(5:10,:) - n*[u,u]
! fms(11,:) = energy density U is unchanged
! fms(12:14,:) = heat flux Q, calculated from energy U, energy flux,
! velocity field u and pressure P: Q = fms(12:14,:) - u.P - u*U
! on entry:
! fms(i,j,k,l) = ith component of fluid moments at grid point j,k,l
! on exit
! fms(i,j,k,l) = ith component of fluid quantities at grid point j,k,l
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx/ny/nz = system length in x/y/z direction
! kstrt = starting data block number
! nvpy = number of real or virtual processors in y
! kyp/kzp = number of real grids in each field partition in y/z direction
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
      implicit none
      integer npro, nx, ny, nz, kstrt, nvpy, kyp, kzp
      integer nprd, nxv, nypmx, nzpmx
      real fms
      dimension fms(nprd,nxv,nypmx,nzpmx)
! local data
      integer j, k, l, n, js, ks, kypp, kzpp, npr
      real at1, at2
      double precision dt1, dt2, dtx, dty, dtz
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 4
      else if (npro==3) then
         npr = 10
      else if (npro==4) then
         npr = 14
      endif
      if (npr > nprd) npr = nprd
! exit if error
      if (npr < 4) return
!$OMP PARALLEL DO PRIVATE(j,k,l,n,at1,at2,dt1,dt2,dtx,dty,dtz)
      do 20 n = 1, kypp*kzpp
      l = (n - 1)/kypp
      k = n - kypp*l
      l = l + 1
      do 10 j = 1, nx
      at1 = fms(1,j,k,l)
! calculate velocity field
      if (at1.gt.0.0) then
         at2 = 1.0/at1
         fms(2,j,k,l) = fms(2,j,k,l)*at2
         fms(3,j,k,l) = fms(3,j,k,l)*at2
         fms(4,j,k,l) = fms(4,j,k,l)*at2
      else
         fms(2,j,k,l) = 0.0
         fms(3,j,k,l) = 0.0
         fms(4,j,k,l) = 0.0
      endif
      if (npr < 10) go to 10
! calculate pressure tensor
      dt1 = dble(at1)
      dtx = dble(fms(2,j,k,l))
      dty = dble(fms(3,j,k,l))
      dtz = dble(fms(4,j,k,l))
      dt2 = dble(fms(5,j,k,l))
      fms(5,j,k,l) = dt2 - dt1*dtx*dtx
      dt2 = dble(fms(6,j,k,l))
      fms(6,j,k,l) = dt2 - dt1*dtx*dty
      dt2 = dble(fms(7,j,k,l))
      fms(7,j,k,l) = dt2 - dt1*dtx*dtz
      dt2 = dble(fms(8,j,k,l))
      fms(8,j,k,l) = dt2 - dt1*dty*dty
      dt2 = dble(fms(9,j,k,l))
      fms(9,j,k,l) = dt2 - dt1*dty*dtz
      dt2 = dble(fms(10,j,k,l))
      fms(10,j,k,l) = dt2 - dt1*dtz*dtz
! calculate heat flux
      if (npr < 14) go to 10
      dt1 = fms(11,j,k,l)
      dt2 = dtx*fms(5,j,k,l) + dty*fms(6,j,k,l) + dtz*fms(7,j,k,l)
      dt2 = fms(12,j,k,l) - dt2 - dtx*dt1
      fms(12,j,k,l) = dt2
      dt2 = dtx*fms(6,j,k,l) + dty*fms(8,j,k,l) + dtz*fms(9,j,k,l)
      dt2 = fms(13,j,k,l) - dt2 - dty*dt1
      fms(13,j,k,l) = dt2
      dt2 = dtx*fms(7,j,k,l) + dty*fms(9,j,k,l) + dtz*fms(10,j,k,l)
      dt2 = fms(14,j,k,l) - dt2 - dtz*dt1
      fms(14,j,k,l) = dt2
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PSTPTRAJ3(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx, &
     &dvtx,nvpy,nvpz,idimp,nppmx,mxyzp1,idps,np,nprobt)
! for 3d code, this procedure sets test charge distribution by
! setting a particle id in particle location 7, whose values are between
! 1 and nprobt
! particles stored in segmented array
! input: all, output: iprobt, nprobt
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(7,n,m) = particle id of tagged particle n in tile m
! tedges(1:2) = back:front z boundaries of particle tags
! tedges(3:4) = lower:upper y boundaries of particle tags
! kpic = number of particles per tile
! iprobt = scratch array of size nprobt, used for nst = 2
! kstrt = starting data block number (processor id + 1)
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 7
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx + 1 and
! myp1 = (partition length in y direction - 1)/my + 1 and
! where mzp1 = (partition length in z direction - 1)/mz + 1
! idps = number of partition boundaries
! np = total number of particles in part
! nprobt = total number of test charges whose trajectories are stored.
! particle id should be <= 16777215
      implicit none
      integer kstrt, nst, nvpy, nvpz, idimp, nppmx, mxyzp1, idps, nprobt
      real vtx, vtsx, dvtx
      double precision np
      real ppart, tedges
      dimension ppart(idimp,nppmx,mxyzp1)
      dimension tedges(idps)
      integer kpic, iprobt
      dimension kpic(mxyzp1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, nppp, noff, joff, it, nt, itt, itcom
      real st, at
      double precision dit, ditt
      integer nty
      dimension nty(1)
      double precision dnpl, dnp
      dimension dnpl(1), dnp(1)
! itcom = (0,1) = (no,yes) transposed communicator created
      data itcom /0/
      save itcom
      if (idimp < 7) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         dit = aint((np-1.0d0)/dble(nprobt)) + 1.0d0
         ditt = 1.0
! find how many particles on global left processor
         npp = 0
         do 10 k = 1, mxyzp1
         npp = npp + kpic(k)
   10    continue
         dnp(1) = dble(npp)
         call PPDSCAN(dnp,dnpl,1)
! dnpl(1) = integrated number of particles on left processor
         dnpl(1) = dnp(1) - dble(npp)
! noff = integrated number of test particles on left processor
         noff = int((dnpl(1) - 1.0d0)/dit + 1.0d0)
! ditt = global address of next test particle
         if (dnpl(1).gt.0.0d0) then
            ditt = dble(noff)*dit + 1.0d0
         endif
! itt = local address of next test particle
         itt = int(ditt - dnpl(1))
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         st = 4.0*vtx
         it = nprobt/2
         at = real(it)
         st = at/st
         at = at + 1.0 + 0.5*(nprobt - 2*it)
         do 20 j = 1, nprobt
         iprobt(j) = 0
   20    continue
! velocity slice in x direction
      else if (nst.eq.3) then
         st = 1.0/dvtx
         itt = vtsx*st + 0.5
      endif
      joff = 0
      nt = 0
! loop over tiles
      do 40 k = 1, mxyzp1
      nppp = kpic(k)
! loop over particles in tile
      do 30 j = 1, nppp
! clear tag
      ppart(7,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(7,j,k) = real(nt)
            ditt = ditt + dit
            itt = int(ditt - dnpl(1))
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         if (kstrt==1) then
            it = ppart(4,j,k)*st + at
            if ((it.gt.0).and.(it.le.nprobt)) then
               if (iprobt(it).eq.0) then
                  nt = nt + 1
                  iprobt(it) = j + joff
                  ppart(7,j,k) = real(nt)
               endif
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(4,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(7,j,k) = real(nt)
         endif
      endif
   30 continue
      joff = joff + nppp
   40 continue
! create transposed communicator to sum in z direction
      if (itcom.eq.0) then
         call PPCOMM_T(nvpy,nvpz)
         itcom = 1
      endif
! local count of test particles on this processor
      nprobt = nt
! find how many test particles on lower processor
      dnp(1) = dble(nt)
      call PPDSCAN_T(dnp,dnpl,1)
      st = real(dnp(1)) - real(nt)
! tedges(3) = integrated number of test particles on lower processor
      tedges(3) = st + 1.0
! tedges(4) = integrated number of test particles on current processor
      tedges(4) = real(dnp(1)) + 1.0
      nty(1) = st
      call PPBICASTZ(nty,1,nvpy,nvpz)
! tedges(1) = number of test particles in a z row on processor below
      tedges(1) = real(nty(1)) + 1.0
      nty(1) = real(dnp(1))
      call PPBICASTZR(nty,1,nvpy,nvpz)
! tedges(2) = number of test particles in a z row on current processor
      tedges(2) = real(nty(1)) + 1.0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,at)
      do 60 k = 1, mxyzp1
      nppp = kpic(k)
! add offset to test particle tags
      do 50 j = 1, nppp
      at = ppart(7,j,k)
      if (at.gt.0.0) then
         ppart(7,j,k) = at + st
      endif
   50 continue
   60 continue
!$OMP END PARALLEL DO
! find total number of test particles on each node
      dnp(1) = dble(nprobt)
      call PPDSUM(dnp,dnpl,1)
      nprobt = min(16777216.0d0,dnp(1))
      return
      end
!-----------------------------------------------------------------------
      subroutine PPTRAJ3(ppart,kpic,partt,numtp,idimp,nppmx,mxyzp1,     &
     &nprobt)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, numtp, output: partt, numtp
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = position y of particle n in tile m
! ppart(4,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = velocity vy of particle n in tile m
! ppart(6,n,m) = velocity vz of particle n in tile m
! ppart(7,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! partt = tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 7
! nppmx = maximum number of particles in tile
! mxyzp1 = mx1*myp1*mzp1, where
! mx1 = (system length in x direction - 1)/mx + 1 and
! myp1 = (partition length in y direction - 1)/my + 1 and
! where mzp1 = (partition length in z direction - 1)/mz + 1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer numtp, idimp, nppmx, mxyzp1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mxyzp1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mxyzp1)
! local data
      integer i, j, k, nppp, nt
      if (idimp < 7) return
      nt = 0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,nppp)
      do 30 k = 1, mxyzp1
      nppp = kpic(k)
! loop over particles in tile
      do 20 j = 1, nppp
      if (ppart(7,j,k).gt.0.0) then
!$OMP ATOMIC
         nt = nt + 1
         if (nt.le.nprobt) then
            do 10 i = 1, idimp
            partt(i,nt) = ppart(i,j,k)
   10       continue
         endif
      endif
   20 continue
   30 continue
!$OMP END PARALLEL DO
      numtp = nt
      return
      end
!-----------------------------------------------------------------------
      subroutine PORDTRAJ3(partt,spartt,tedges,numtp,idimp,idps,nprobt)
! this procedure reorders tagged particles in partt to array spartt
! tagged particle coordinates stored in tag order
! input: all except spartt, output: spartt
! partt = tagged particle coordinates
! spartt = reordered tagged particle coordinates
! tedges(1:2) = back:front z boundaries of particle tags
! tedges(3:4) = lower:upper y boundaries of particle tags
! numtp = number of test particles found on this node
! idimp = size of phase space = 7
! idps = number of partition boundaries
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer numtp, idimp, idps, nprobt
      real partt, spartt, tedges
      dimension partt(idimp,nprobt), spartt(idimp,nprobt)
      dimension tedges(idps)
! local data
      integer i, j, nt
      real tn, st
      if (idimp < 7) return
      st = tedges(3) - 1.0
! loop over tagged particles
!$OMP PARALLEL DO PRIVATE(i,j,nt,tn)
      do 20 j = 1, numtp
      tn = partt(7,j)
      if (tn.gt.0.0) then
         nt = tn - st
         do 10 i = 1, idimp
         spartt(i,nt) = partt(i,j)
   10    continue
      endif
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PCPYTRAJ3(partt,part,numtp,idimp,nprobt)
! this procedure copies tagged particles in partt to array part
! tagged particle coordinates stored in tag order
! input: all except part, output: part
! spartt = tagged particle coordinates
! spartt = reordered tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 6
! nprobt = number of test charges whose trajectories will be stored.
      implicit none
      integer numtp, idimp, nprobt
      real partt, part
      dimension partt(idimp,nprobt), part(idimp,nprobt)
! local data
      integer i, j, np
      np = min(numtp,nprobt)
! loop over tagged particles
!$OMP PARALLEL DO PRIVATE(i,j)
      do 20 j = 1, np
      do 10 i = 1, idimp
      part(i,j) = partt(i,j)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end