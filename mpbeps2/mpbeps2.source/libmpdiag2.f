!-----------------------------------------------------------------------
! Fortran Library for diagnostics
! 2-1/2D MPI/OpenMP PIC Codes:
! PPVDIST2 calculates 2 component velocity distribution, velocity
!          moments and entropy for segmented particle array,
!          with distributed data.
! PPVDIST23 calculates 3 component velocity distribution, velocity
!           moments and entropy for segmented particle array,
!           with distributed data.
! PERPDIST2 calculates 2d energy distribution for relativistic particles
! PERPDIST23 calculates 2-1/2d energy distribution for relativistic
!            particles
! PPVBDIST23 for 2-1/2d code, calculates 3d velocity distribution and
!            velocity moments for magnetized plasma with segmented
!            particle array
! PPVSDIST2 for 2d code, calculates 2d velocity distribution, in
!           different regions of space, for segmented particle array
!           with distributed data.
! PPVSDIST23 for 2-1/2d code, calculates 3d velocity distribution, in
!            different regions of space, for segmented particle array
!            with distributed data.
! PVDIST2 calculates 2 component velocity distribution, velocity moments
!         and entropy for standard particle array, for distributed data.
! PVDIST23 calculates 3 component velocity distribution, velocity
!          moments, and entropy for standard particle array,
!          for distributed data.
! PVBDIST23 for 2-1/2d code, calculates 3d velocity distribution and
!           velocity moments for magnetized plasma with standard
!           particle array
! PPROFX23L calculates fluid moments from particle quantities: density,
!           momentum, momentum flux, energy, energy flux, assumes
!           particle positions and velocities at same time level
!           for 2-1/2d code
! PRPROFX23L calculates fluid moments from relativistic particle
!            quantities: density, velocity, velocity flux, energy,
!            energy flux, assumes particle positions and velocities at 
!            same time level, for 2-1/2d code
! PPROFX22L calculates fluid moments from particle quantities: density,
!           momentum, momentum flux, energy, energy flux, assumes
!           particle positions and velocities at same time level
!           for 2d code
! PRPROFX22L calculates fluid moments from relativistic particle
!            quantities: density, velocity, velocity flux, energy,
!            energy flux, assumes particle positions and velocities at 
!            same time level, for 2d code
! PGPROFX2L calculates fluid moments from particle quantities:
!           density, momentum, momentum flux, energy, energy flux,
!           assumes particle positions and velocities not at same time
!           levels and electrostatic fields
! PGRPROFX2L calculates fluid moments from relativistic particle
!            quantities: density, velocity, velocity flux, energy,
!            energy flux, assumes particle positions and velocities
!            not at same time levels and electrostatic fields
! PGBPROFX23L calculates fluid moments from particle quantities:
!             density, momentum, momentum flux, energy, energy flux,
!             assumes particle positions and velocities not at same time
!             levels and electromagnetic fields
! PGRBPROFX23L calculates fluid moments from relativistic particle
!              quantities: density, velocity, velocity flux, energy,
!              energy flux
!              assumes particle positions and velocities not at same
!              time levels and electromagnetic fields
! FLUIDQS23 calculates fluid quantities from fluid moments:
!           density, velocity field, pressure tensor, energy, heat flux
!           for 2-1/2d code
! FLUIDQS22 calculates fluid quantities from fluid moments:
!           density, velocity field, pressure tensor, energy, heat flux
!           for 2d code
! PSTPTRAJ2 sets test charge distribution by setting a particle id
!           in particle location 5 for 2d code
! PSTPTRAJ23 sets test charge distribution by setting a particle id
!            in particle location 6 for 2-1/2d code
! PPTRAJ2 copies tagged particles in ppart to array partt for 2d code
! PPTRAJ23 copies tagged particles in ppart to array partt for 2-1/2d
!          code
! PORDTRAJ2 reorders tagged particles in partt to array spartt for 2d
!           code
! PORDTRAJ23 reorders tagged particles in partt to array spartt for
!            2-1/2d
! PCPYTRAJ2 copies tagged particles in partt to array part
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: march 13, 2018
!-----------------------------------------------------------------------
      subroutine PPVDIST2(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyp1,  &
     &nmv,nmvf)
! for 2d code, this subroutine calculates 2d velocity distribution,
! velocity moments, and entropy, with 1D domain decomposition
! particles stored in segmented array
! input: all except fvm, output: fv, sfv, fvm
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! kpic = number of particles per tile
! fv = global distribution function, number of particles in each
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
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nvp, idimp, nppmx, mxyp1, nmv, nmvf
      real ppart, fv, sfv, fvm
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fv(nmvf,2), sfv(nmvf,2,mxyp1), fvm(2,3)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer j, k, nmv21, nppp, nvx, nvy
      real anmv, svx, svy, svxx, svyx, vx, vy
      double precision sumvx, sumvy, sumvx2, sumvy2, anp
      double precision ssumvx, ssumvy, ssumvx2, ssumvy2, sanp
      double precision sum1, sum2
      double precision sum5, work5
      dimension sum5(5), work5(5)
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svy = anmv/fv(nmv21+1,2)
! normalization constant for entropy
      svxx = svx*real(mxyp1*nvp)
      svyx = svy*real(mxyp1*nvp)
! zero out distribution
      do 20 k = 1, mxyp1
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
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,nppp,nvx,nvy,vx,vy,sanp,ssumvx,ssumvy,ssumvx2,       &
!$OMP& ssumvy2) REDUCTION(+:anp) REDUCTION(+:sumvx) REDUCTION(+:sumvy)  &
!$OMP& REDUCTION(+:sumvx2) REDUCTION(+:sumvy2)
      do 40 k = 1, mxyp1
      nppp = kpic(k)
      sanp = 0.0d0
      ssumvx = 0.0d0
      ssumvy = 0.0d0
      ssumvx2 = 0.0d0
      ssumvy2 = 0.0d0
      do 30 j = 1, nppp
      sanp = sanp + 1.0d0
      vx = ppart(3,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      vy = ppart(4,j,k)
      nvy = vy*svy + anmv
      ssumvy = ssumvy + vy
      ssumvy2 = ssumvy2 + vy*vy
      if ((nvy.ge.1).and.(nvy.le.nmv21)) sfv(nvy,2,k) = sfv(nvy,2,k)+1.0
   30 continue
! calculate global sums
      anp = anp + sanp
      sumvx = sumvx + ssumvx
      sumvy = sumvy + ssumvy
      sumvx2 = sumvx2 + ssumvx2
      sumvy2 = sumvy2 + ssumvy2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 50 k = 1, mxyp1
      sum1 = sum1 + sfv(j,1,k)
      sum2 = sum2 + sfv(j,2,k)
   50 continue
      fv(j,1) = sum1
      fv(j,2) = sum2
   60 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      sum5(5) = anp
      call PPDSUM(sum5,work5,5)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      anp = sum5(5)
! calculate global velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(2,1) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
! count number of particles in global distribution
      sumvx = 0.0d0
      sumvy = 0.0d0
      do 70 j = 1, nmv21
      sumvx = sumvx + fv(j,1)
      sumvy = sumvy + fv(j,2)
   70 continue
! calculate entropy
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 90 k = 1, mxyp1
      do 80 j = 1, nmv21
      if (sfv(j,1,k).gt.0.0) then
         sumvx2 = sumvx2 + sfv(j,1,k)*dlog(dble(sfv(j,1,k)*svxx))
      endif
      if (sfv(j,2,k).gt.0.0) then
         sumvy2 = sumvy2 + sfv(j,2,k)*dlog(dble(sfv(j,2,k)*svyx))
      endif
   80 continue
   90 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      call PPDSUM(sum5,work5,4)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(1,3) = sumvx
      fvm(2,3) = sumvy
      return
      end
!-----------------------------------------------------------------------
      subroutine PPVDIST23(ppart,kpic,fv,sfv,fvm,nvp,idimp,nppmx,mxyp1, &
     &nmv,nmvf)
! for 2-1/2d code, this subroutine calculates 3d velocity distribution,
! velocity moments, and entropy, with 1D domain decomposition
! particles stored in segmented array
! input: all except fvm, output: fv, sfv, fvm
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! ppart(5,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! fv = global distribution function, number of particles in each
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
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nvp, idimp, nppmx, mxyp1, nmv, nmvf
      real ppart, fv, sfv, fvm
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fv(nmvf,3), sfv(nmvf,3,mxyp1), fvm(3,3)
      integer kpic
      dimension kpic(mxyp1)
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
      svxx = svx*real(mxyp1*nvp)
      svyx = svy*real(mxyp1*nvp)
      svzx = svz*real(mxyp1*nvp)
! zero out distribution
      do 20 k = 1, mxyp1
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
!$OMP& ssumvx2,ssumvy2,ssumvz2) REDUCTION(+:anp) REDUCTION(+:sumvx)     &
!$OMP& REDUCTION(+:sumvy) REDUCTION(+:sumvz) REDUCTION(+:sumvx2)        &
!$OMP& REDUCTION(+:sumvy2) REDUCTION(+:sumvz2)
      do 40 k = 1, mxyp1
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
      vx = ppart(3,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      vy = ppart(4,j,k)
      nvy = vy*svy + anmv
      ssumvy = ssumvy + vy
      ssumvy2 = ssumvy2 + vy*vy
      if ((nvy.ge.1).and.(nvy.le.nmv21)) sfv(nvy,2,k) = sfv(nvy,2,k)+1.0
      vz = ppart(5,j,k)
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
      do 50 k = 1, mxyp1
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
      do 90 k = 1, mxyp1
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
      subroutine PERPDIST2(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,nmv&
     &,nmvf)
! for 2d code, this subroutine calculates 2d energy distribution,
! for relativistic particles, with 1D domain decomposition
! the function calculated is of the form g*exp(-e/vth**2), where
! e = (gamma-1)*(c*c) is the kinetic energy per mass, and where
! vth = sqrt(KT/m).  Note vth is a momentum/mass and can be > c
! gamma = sqrt(1 + (p*p)/(c*c)), where p = is the momentum per mass
! g = p/(de/dp) = gamma
! one can express this quantity g as a function of e as follows:
! e = (p*p)/(gamma+1) => p = sqrt((gamma+1)*e), and gamma = 1 + e/(c*c)
! particles stored in segmented array
! input: all except wk, output: fv, sfv, wk
! ppart(3,n,m) = momentum px of particle n in tile m
! ppart(4,n,m) = momentum py of particle n in tile m
! kpic = number of particles per tile
! fv = global distribution function, number of particles in each energy
! range, summed over tiles.  maximum energy (used for scaling) is
! contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! energy range
! ci = reciprocal of velocity of light
! wk = total energy contained in distribution
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fv
! the number of energy bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mxyp1, nmv, nmvf
      real ci, wk
      real ppart, fv, sfv
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fv(nmvf,1), sfv(nmvf,2,mxyp1)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer j, k, nppp, nmv21, nvx
      real ci2, anmv, svx, px, py, p2
      double precision sumpx, ssumpx, sum1
      double precision dsum1, dwork1
      dimension dsum1(1), dwork1(1)
      ci2 = ci*ci
! energy scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
! zero out distribution
      do 20 k = 1, mxyp1
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
!$OMP& PRIVATE(j,k,nppp,nvx,px,py,p2,ssumpx) REDUCTION(+:sumpx) 
      do 40 k = 1, mxyp1
      nppp = kpic(k)
      ssumpx = 0.0d0
      do 30 j = 1, nppp
      px = ppart(3,j,k)
      py = ppart(4,j,k)
      p2 = px*px + py*py
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
      do 50 k = 1, mxyp1
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
      subroutine PERPDIST23(ppart,kpic,fv,sfv,ci,wk,idimp,nppmx,mxyp1,  &
     &nmv,nmvf)
! for 2-1/2d code, this subroutine calculates 3d energy distribution,
! for relativistic particles, with 1D domain decomposition
! the function calculated is of the form g*exp(-e/vth**2), where
! e = (gamma-1)*(c*c) is the kinetic energy per mass, and where
! vth = sqrt(KT/m).  Note vth is a momentum/mass and can be > c
! gamma = sqrt(1 + (p*p)/(c*c)), where p = is the momentum per mass
! g = p*p/(de/dp) = p*gamma
! one can express this quantity g as a function of e as follows:
! e = (p*p)/(gamma+1) => p = sqrt((gamma+1)*e), and gamma = 1 + e/(c*c)
! particles stored in segmented array
! input: all except wk, output: fv, sfv, wk
! ppart(3,n,m) = momentum px of particle n in tile m
! ppart(4,n,m) = momentum py of particle n in tile m
! ppart(5,n,m) = momentum pz of particle n in tile m
! kpic = number of particles per tile
! fv = global distribution function, number of particles in each energy
! range, summed over tiles.  maximum energy (used for scaling) is
! contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! energy range
! ci = reciprocal of velocity of light
! wk = total energy contained in distribution
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fv
! the number of energy bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mxyp1, nmv, nmvf
      real ci, wk
      real ppart, fv, sfv
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fv(nmvf,1), sfv(nmvf,3,mxyp1)
      integer kpic
      dimension kpic(mxyp1)
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
      do 20 k = 1, mxyp1
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
      do 40 k = 1, mxyp1
      nppp = kpic(k)
      ssumpx = 0.0d0
      do 30 j = 1, nppp
      px = ppart(3,j,k)
      py = ppart(4,j,k)
      pz = ppart(5,j,k)
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
      do 50 k = 1, mxyp1
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
      subroutine PPVBDIST23(ppart,kpic,fv,sfv,fvm,omx,omy,omz,idimp,    &
     &nppmx,mxyp1,nmv,nmvf)
! for 2-1/2d code, this subroutine calculates 3d velocity distribution,
! and velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! with 1D domain decomposition
! particles stored in segmented array
! input: all except fvm, output: fv, sfv, fvm
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! ppart(5,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! fv = global distribution function, number of particles in each
! velocity range, summed over tiles. maximum velocity (used for scaling)
! is contained in last element of fv
! sfv = distribution function in tile, number of particles in each
! velocity range
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mxyp1, nmv, nmvf
      real omx, omy, omz
      real ppart, fv, sfv, fvm
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fv(nmvf,2), sfv(nmvf,3,mxyp1), fvm(3,3)
      integer kpic
      dimension kpic(mxyp1)
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
      do 20 k = 1, mxyp1
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
!$OMP& ssumvx2,ssumvz2) REDUCTION(+:anp) REDUCTION(+:sumvx)             &
!$OMP& REDUCTION(+:sumvz) REDUCTION(+:sumvx2) REDUCTION(+:sumvz2)
      do 40 k = 1, mxyp1
      nppp = kpic(k)
      sanp = 0.0d0
      ssumvx = 0.0d0
      ssumvz = 0.0d0
      ssumvx2 = 0.0d0
      ssumvz2 = 0.0d0
      do 30 j = 1, nppp
      sanp = sanp + 1.0d0
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
      do 50 k = 1, mxyp1
      sum1 = sum1 + sfv(j,1,k)
      sum2 = sum2 + sfv(j,2,k)
   50 continue
      fv(j,1) = sum1
      fv(j,2) = sum2
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
      subroutine PPVSDIST2(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb,idimp&
     &,nppmx,mxyp1,nmvf)
! for 2d code, this subroutine calculates 2d velocity distribution, in
! different regions of space, with 1D domain decomposition
! particles stored in segmented array
! input: all except fvs, output: fvs
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! kpic = number of particles per tile
! fvs = spatially resolved distribution function, number of particles in
! each velocity and spatial range.  maximum velocity (used for scaling)
! is contained in last element of first dimension of fvs
! noff = lowermost global gridpoint in particle partition.
! nmv = number of segments in v for velocity distribution
! mvx/mvy = number of grids in x/y for phase space aggregation
! nxb/nyb = number of segments in x/y for velocity distribution
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fvs
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer noff, nmv, mvx, mvy, nxb, nyb, idimp, nppmx, mxyp1, nmvf
      real ppart, fvs
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fvs(nmvf,2,nxb,nyb+1)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer j, k, ns, nppp, nmv21, nn, mm, nvx, nvy
      real anmv, anoff, svx, svy, at1, at2
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      ns = noff - mvy*(noff/mvy)
      anmv = real(nmv)
      anoff = real(noff - ns)
      svx = anmv/fvs(nmv21+1,1,1,1)
      svy = anmv/fvs(nmv21+1,2,1,1)
! spatial scaling
      at1 = 1.0/real(mvx)
      at2 = 1.0/real(mvy)
! zero out distribution
      do 30 mm = 1, nyb+1
      do 20 nn = 1, nxb
      do 10 j = 1, nmv21
      fvs(j,1,nn,mm) = 0.0
      fvs(j,2,nn,mm) = 0.0
   10 continue
   20 continue
   30 continue
! count particles in each velocity region
      anmv = anmv + 1.5
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,nn,mm,nvx,nvy)
      do 50 k = 1, mxyp1
      nppp = kpic(k)
      do 40 j = 1, nppp
      nn = ppart(1,j,k)*at1 + 1.0
      mm = (ppart(2,j,k) - anoff)*at2 + 1.0
      nvx = ppart(3,j,k)*svx + anmv
      nvy = ppart(4,j,k)*svy + anmv
      if ((nvx.ge.1).and.(nvx.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvx,1,nn,mm) = fvs(nvx,1,nn,mm) + 1.0
      endif
      if ((nvy.ge.1).and.(nvy.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvy,2,nn,mm) = fvs(nvy,2,nn,mm) + 1.0
      endif
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPVSDIST23(ppart,kpic,fvs,noff,nmv,mvx,mvy,nxb,nyb,    &
     &idimp,nppmx,mxyp1,nmvf)
! for 2-1/2d code, this subroutine calculates 3d velocity distribution,
! in different regions of space, with 1D domain decomposition
! particles stored in segmented array
! input: all except fvs, output: fvs
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! ppart(5,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! fvs = spatially resolved distribution function, number of particles in
! each velocity and spatial range.  maximum velocity (used for scaling)
! is contained in last element of first dimension of fvs
! noff = lowermost global gridpoint in particle partition.
! nmv = number of segments in v for velocity distribution
! mvx/mvy = number of grids in x/y for phase space aggregation
! nxb/nyb = number of segments in x/y for velocity distribution
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nmvf = first dimension of fvs
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer noff, nmv, mvx, mvy, nxb, nyb, idimp, nppmx, mxyp1, nmvf
      real ppart, fvs
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fvs(nmvf,3,nxb,nyb+1)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer j, k, ns, nppp, nmv21, nn, mm, nvx, nvy, nvz
      real anmv, anoff, svx, svy, svz, at1, at2
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      ns = noff - mvy*(noff/mvy)
      anmv = real(nmv)
      anoff = real(noff - ns)
      svx = anmv/fvs(nmv21+1,1,1,1)
      svy = anmv/fvs(nmv21+1,2,1,1)
      svz = anmv/fvs(nmv21+1,3,1,1)
! spatial scaling
      at1 = 1.0/real(mvx)
      at2 = 1.0/real(mvy)
! zero out distribution
      do 30 mm = 1, nyb+1
      do 20 nn = 1, nxb
      do 10 j = 1, nmv21
      fvs(j,1,nn,mm) = 0.0
      fvs(j,2,nn,mm) = 0.0
      fvs(j,3,nn,mm) = 0.0
   10 continue
   20 continue
   30 continue
! count particles in each velocity region
      anmv = anmv + 1.5
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,nn,mm,nvx,nvy,nvz)
      do 50 k = 1, mxyp1
      nppp = kpic(k)
      do 40 j = 1, nppp
      nn = ppart(1,j,k)*at1 + 1.0
      mm = (ppart(2,j,k) - anoff)*at2 + 1.0
      nvx = ppart(3,j,k)*svx + anmv
      nvy = ppart(4,j,k)*svy + anmv
      nvz = ppart(5,j,k)*svz + anmv
      if ((nvx.ge.1).and.(nvx.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvx,1,nn,mm) = fvs(nvx,1,nn,mm) + 1.0
      endif
      if ((nvy.ge.1).and.(nvy.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvy,2,nn,mm) = fvs(nvy,2,nn,mm) + 1.0
      endif
      if ((nvz.ge.1).and.(nvz.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvz,3,nn,mm) = fvs(nvz,3,nn,mm) + 1.0
      endif
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDIST2(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
! for 2d code, this subroutine calculates 2d velocity distribution,
! velocity moments, and entropy, with 1D domain decomposition
! input: all except fvm, output: fv, fvm
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
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
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer npp, nvp, idimp, npmax, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,npmax)
      dimension fv(nmvf,2), fvm(2,3)
! local data
      integer j, nmv21, nvx, nvy
      real anmv, svx, svy, svxx, svyx
      double precision sumvx, sumvy, sumvx2, sumvy2, anp
      double precision sum5, work5
      dimension sum5(5), work5(5)
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svy = anmv/fv(nmv21+1,2)
! normalization constant for entropy
      svxx = svx*real(nvp)
      svyx = svy*real(nvp)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j,1) = 0.0
      fv(j,2) = 0.0
   10 continue
! count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 20 j = 1, npp
      anp = anp + 1.0d0
      nvx = part(3,j)*svx + anmv
      sumvx = sumvx + part(3,j)
      sumvx2 = sumvx2 + part(3,j)**2
      nvy = part(4,j)*svy + anmv
      sumvy = sumvy + part(4,j)
      sumvy2 = sumvy2 + part(4,j)**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx,1) = fv(nvx,1) + 1.0
      if ((nvy.ge.1).and.(nvy.le.nmv21)) fv(nvy,2) = fv(nvy,2) + 1.0
   20 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      sum5(5) = anp
      call PPDSUM(sum5,work5,5)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      anp = sum5(5)
! calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(2,1) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
! calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 30 j = 1, nmv21
      if (fv(j,1).gt.0.0) then
         sumvx = sumvx + fv(j,1)
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svxx))
      endif
      if (fv(j,2).gt.0.0) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svyx))
      endif
   30 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      call PPDSUM(sum5,work5,4)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(1,3) = sumvx
      fvm(2,3) = sumvy
      return
      end
!-----------------------------------------------------------------------
      subroutine PVDIST23(part,fv,fvm,npp,nvp,idimp,npmax,nmv,nmvf)
! for 2-1/2d code, this subroutine calculates 3d velocity distribution,
! velocity moments, and entropy, with 1D domain decomposition
! input: all except fvm, output: fv, fvm
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
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
! idimp = size of phase space = 5
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
      nvx = part(3,j)*svx + anmv
      sumvx = sumvx + part(3,j)
      sumvx2 = sumvx2 + part(3,j)**2
      nvy = part(4,j)*svy + anmv
      sumvy = sumvy + part(4,j)
      sumvy2 = sumvy2 + part(4,j)**2
      nvz = part(5,j)*svz + anmv
      sumvz = sumvz + part(5,j)
      sumvz2 = sumvz2 + part(5,j)**2
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
      subroutine PVBDIST23(part,fv,fvm,omx,omy,omz,npp,idimp,npmax,nmv, &
     &nmvf)
! for 2-1/2d code, this subroutine calculates 3d velocity distribution,
! and velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! with 1D domain decomposition
! input: all except fvm, output: fv, fvm
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! npp = number of particles in partition
! idimp = size of phase space = 5
! npmax = maximum number of particles in each partition
! nmvf = first dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer npp, idimp, npmax, nmv, nmvf
      real omx, omy, omz
      real part, fv, fvm
      dimension part(idimp,npmax)
      dimension fv(nmvf,2), fvm(3,3)
! local data
      integer j, ndir, nmv21, nvx, nvz
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real anmv, svx, svz
      double precision sumvx, sumvz, sumvx2, sumvz2, anp
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
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
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
      subroutine PPROFX23L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,  &
     &nprd,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates fluid moments from
! particle quantities: density, momentum, momentum flux, energy,
! energy flux
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 138 flops/particle, 61 loads, 56 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, npro, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real dxp, dyp, amx, amy, dx, dy, x, y, w, vx, vy, vz
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV), sg(14)
      dimension sfms(nprd,mx+1,my+1), sg(14)
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
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,vz,sfms,sg)
      do 150 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! deposit fluid moments
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
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
         dx = vz*vz
         sg(10) = dx
         w = 0.5*(w + dx)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 40 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
      do 70 i = 2, nn
      do 60 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   60 continue
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 110 i = 2, nn
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
   90 continue
      if (mm > my) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  100    continue
      endif
  110 continue
      nn = min(mx+1,nxv-noffp)
      do 140 j = 1, mm
      do 120 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  120 continue
      if (nn > mx) then
         do 130 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  130    continue
      endif
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PRPROFX23L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx, &
     &my,nprd,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates fluid moments from
! particle quantities: density, velocity, velocity flux, energy, energy
! flux, for relativistic particles
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 144 flops/particle, 2 divides, 1 sqrt, 61 loads, 56 stores
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci =  gami*p2/(1.0 + gami),
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! ppart(5,n,m) = momentum pz of particle n in partition in tile m
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, npro, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real ci
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real dxp, dyp, amx, amy, dx, dy, x, y, w, vx, vy, vz
      real ci2, p2, gami
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV), sg(14)
      dimension sfms(nprd,mx+1,my+1), sg(14)
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
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,vz,p2,gami,sfms,sg)
      do 150 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 40 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
      do 70 i = 2, nn
      do 60 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   60 continue
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 110 i = 2, nn
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
   90 continue
      if (mm > my) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  100    continue
      endif
  110 continue
      nn = min(mx+1,nxv-noffp)
      do 140 j = 1, mm
      do 120 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  120 continue
      if (nn > mx) then
         do 130 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  130    continue
      endif
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPROFX22L(ppart,fms,kpic,noff,nppmx,idimp,npro,mx,my,  &
     &nprd,nxv,nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 93 flops/particle, 41 loads, 36 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 3, mci = vi, where i = x,y
! where for i = 4, 6, mci = vj*vk, where jk = xx,xy,yy
! where for i = 7, mci = (vx**2+vy**2)
! where for i = 8, 9, mci = (vx**2+vy**2)*vi, where i = x,y
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 9
! nxv = second dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, npro, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real dxp, dyp, amx, amy, dx, dy, x, y, w, vx, vy
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV), sg(9)
      dimension sfms(nprd,mx+1,my+1), sg(9)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 3
      else if (npro==3) then
         npr = 6
      else if (npro==4) then
         npr = 9
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,sfms,sg)
      do 150 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! deposit fluid moments
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      if (npr > 3) then
         x = vx*vx
         sg(4) = x
         sg(5) = vx*vy
         y = vy*vy
         w = 0.5*(x + y)
         sg(6) = y
         sg(7) = w
         sg(8) = w*vx
         sg(9) = w*vy
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 40 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
      do 70 i = 2, nn
      do 60 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   60 continue
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 110 i = 2, nn
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
   90 continue
      if (mm > my) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  100    continue
      endif
  110 continue
      nn = min(mx+1,nxv-noffp)
      do 140 j = 1, mm
      do 120 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  120 continue
      if (nn > mx) then
         do 130 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  130    continue
      endif
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PRPROFX22L(ppart,fms,kpic,noff,ci,nppmx,idimp,npro,mx, &
     &my,nprd,nxv,nypmx,mx1,mxyp1)
! for d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux,
! for relativistic particles
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 97 flops/particle, 2 divides, 1 sqrt, 41 loads, 36 stores
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 3, mci = vi, where i = x,y
! where for i = 4, 6, mci = vj*vk, where jk = xx,xy,yy
! where for i = 7, mci = gami*p2/(1.0 + gami), where p2 = px*px + py*py
! where for i = 8, 9, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nppmx, idimp, npro, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real ci
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real dxp, dyp, amx, amy, dx, dy, x, y, w, vx, vy
      real ci2, p2, gami
      real sfms, sg
!     dimension sfms(nprd,MXV,MYV), sg(9)
      dimension sfms(nprd,mx+1,my+1), sg(9)
      ci2 = ci*ci
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 3
      else if (npro==3) then
         npr = 6
      else if (npro==4) then
         npr = 9
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,p2,gami,sfms,sg)
      do 150 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 30 j = 1, my+1
      do 20 i = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,i,j) = 0.0
   10 continue
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
! deposit fluid moments
      vx = vx*gami
      vy = vy*gami
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      if (npr > 3) then
         sg(4) = vx*vx
         sg(5) = vx*vy
         sg(6) = vy*vy
         w = gami*p2/(1.0 + gami)
         sg(7) = w
         sg(8) = w*vx
         sg(9) = w*vy
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 40 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
      do 70 i = 2, nn
      do 60 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   60 continue
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 110 i = 2, nn
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
   90 continue
      if (mm > my) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  100    continue
      endif
  110 continue
      nn = min(mx+1,nxv-noffp)
      do 140 j = 1, mm
      do 120 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  120 continue
      if (nn > mx) then
         do 130 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  130    continue
      endif
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,idimp,    &
     &nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 113 flops/particle, 48 loads, 36 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 3, mci = vi, where i = x,y
! where for i = 4, 6, mci = vj*vk, where jk = xx,xy,yy
! where for i = 7, mci = (vx**2+vy**2)
! where for i = 8, 9, mci = (vx**2+vy**2)*vi, where i = x,y
! velocity equations used are:
! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 9
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, npro, nx, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real qbm, dt
      real ppart, fxy, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(2,nxv,nypmx), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff,  nn, mm
      real qtmh, dxp, dyp, amx, amy, x, y, w, dx, dy, vx, vy
      real sfxy, sfms, sg
!     dimension sfxy(2,MXV,MYV), sfms(nprd,MXV,MYV), sg(9)
      dimension sfxy(2,mx+1,my+1), sfms(nprd,mx+1,my+1), sg(9)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 3
      else if (npro==3) then
         npr = 6
      else if (npro==4) then
         npr = 9
      endif
      qtmh = 0.5*qbm*dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,mnoff,nn,mm,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,sfxy,sfms,sg)
      do 170 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
! zero out local accumulator
      do 50 j = 1, my+1
      do 40 i = 1, mx+1
      do 30 ii = 1, npr
      sfms(ii,i,j) = 0.0
   30 continue
   40 continue
   50 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
! half acceleration
      vx = ppart(3,j,k) + x
      vy = ppart(4,j,k) + y
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      if (npr > 3) then
         x = vx*vx
         sg(4) = x
         sg(5) = vx*vy
         y = vy*vy
         w = 0.5*(x + y)
         sg(6) = y
         sg(7) = w
         sg(8) = w*vx
         sg(9) = w*vy
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 60 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   60 continue
   70 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 100 j = 2, mm
      do 90 i = 2, nn
      do 80 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   80 continue
   90 continue
  100 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 130 i = 2, nn
      do 110 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
  110 continue
      if (mm > my) then
         do 120 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  120    continue
      endif
  130 continue
      nn = min(mx+1,nxv-noffp)
      do 160 j = 1, mm
      do 140 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  140 continue
      if (nn > mx) then
         do 150 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  150    continue
      endif
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGRPROFX2L(ppart,fxy,fms,kpic,noff,nyp,qbm,dt,ci,idimp,&
     &nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! for 2d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux,
! for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 120 flops/particle, 2 divides, 1 sqrt, 48 loads, 36 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 3, mci = vi, where i = x,y
! where for i = 4, 6, mci = vj*vk, where jk = xx,xy,yy
! where for i = 7, mci = gami*p2/(1.0 + gami), where p2 = px*px + py*py
! where for i = 8, 9, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y
! momentum equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
! py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
! fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
! the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
!    + dx*fy(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! in other words, fxy are the convolutions of the electric field
! over the particle shape, where kk = k + noff - 1
! fms(i,j,k) = ith component of fluid moments at grid (j,kk)
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 9
! nxv = first dimension of field array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, npro, nx, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real qbm, dt, ci
      real ppart, fxy, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(2,nxv,nypmx), fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff,  nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy, x, y, w, dx, dy, vx, vy
      real p2, gami
      real sfxy, sfms, sg
!     dimension sfxy(2,MXV,MYV), sfms(nprd,MXV,MYV), sg(9)
      dimension sfxy(2,mx+1,my+1), sfms(nprd,mx+1,my+1), sg(9)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 3
      else if (npro==3) then
         npr = 6
      else if (npro==4) then
         npr = 9
      endif
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,mnoff,nn,mm,x,y,w,dxp,dyp,amx, &
!$OMP& amy,dx,dy,vx,vy,p2,gami,sfxy,sfms,sg)
      do 170 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global array
      do 20 j = 1, min(my,nyp-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
! zero out local accumulator
      do 50 j = 1, my+1
      do 40 i = 1, mx+1
      do 30 ii = 1, npr
      sfms(ii,i,j) = 0.0
   30 continue
   40 continue
   50 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
! half acceleration
      vx = ppart(3,j,k) + x
      vy = ppart(4,j,k) + y
! update inverse gamma
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate average velocity
      vx = gami*vx
      vy = gami*vy
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      if (npr > 3) then
         sg(4) = vx*vx
         sg(5) = vx*vy
         sg(6) = vy*vy
         w = gami*p2/(1.0 + gami)
         sg(7) = w
         sg(8) = w*vx
         sg(9) = w*vy
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 60 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   60 continue
   70 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 100 j = 2, mm
      do 90 i = 2, nn
      do 80 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
   80 continue
   90 continue
  100 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 130 i = 2, nn
      do 110 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
  110 continue
      if (mm > my) then
         do 120 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  120    continue
      endif
  130 continue
      nn = min(mx+1,nxv-noffp)
      do 160 j = 1, mm
      do 140 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  140 continue
      if (nn > mx) then
         do 150 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  150    continue
      endif
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,    &
     &idimp,nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates fluid moments from
! particle quantities: density, momentum, momentum flux, energy,
! energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 252 flops/particle, 1 divide, 141 loads, 112 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! fms(i,j,k) = ith component of fluid moments at grid point j,kk
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, npro, nx, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real qbm, dt
      real ppart, fxy, bxy, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, w, vx, vy, vz
      real sfxy, sbxy, sfms, sg
!     dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
!     dimension sfms(nprd,MXV,MYV), sg(14)
      dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      dimension sfms(nprd,mx+1,my+1), sg(14)
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
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,vx,vy,vz,dxp,&
!$OMP& dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,    &
!$OMP& sfms,sg)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 70 j = 1, my+1
      do 60 i = 1, mx+1
      do 50 ii = 1, npr
      sfms(ii,i,j) = 0.0
   50 continue
   60 continue
   70 continue
! loop over particles in tile
      do 90 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      w = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      acx = vx + x
      acy = vy + y
      acz = vz + w
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
      w = (rot7*acx + rot8*acy + rot9*acz)*anorm + w
! calculate average velocity
      vx = 0.5*(x + vx)
      vy = 0.5*(y + vy)
      vz = 0.5*(w + vz)
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
         dx = vz*vz
         sg(10) = dx
         w = 0.5*(w + dx)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 80 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   80 continue
   90 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 120 j = 2, mm
      do 110 i = 2, nn
      do 100 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
  100 continue
  110 continue
  120 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 150 i = 2, nn
      do 130 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
  130 continue
      if (mm > my) then
         do 140 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  140    continue
      endif
  150 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
      do 160 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  160 continue
      if (nn > mx) then
         do 170 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  170    continue
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PGRBPROFX23L(ppart,fxy,bxy,fms,kpic,noff,nyp,qbm,dt,ci,&
     &idimp,nppmx,npro,nx,mx,my,nprd,nxv,nypmx,mx1,mxyp1)
! for 2-1/2d code, this subroutine calculates fluid moments from
! particle quantities: density, velocity, velocity flux, energy,
! energy flux, for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! for distributed data, with 1D spatial decomposition
! OpenMP version using guard cells
! data deposited in tiles
! particles stored segmented array
! 269 flops/particle, 4 divides, 2 sqrt, 141 loads, 112 stores,
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n,m)=mci*(1.-dx)*(1.-dy)
! fms(i,n+1,m)=mci*dx*(1.-dy)
! fms(i,n,m+1)=mci*(1.-dx)*dy
! fms(i,n+1,m+1)=mci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci =  gami*p2/(1.0 + gami),
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
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
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! fms(i,j,k) = ith component of fluid moments at grid point j,kk
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, npro, nx, mx, my
      integer nprd, nxv, nypmx, mx1, mxyp1
      real qbm, dt, ci
      real ppart, fxy, bxy, fms
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension fms(nprd,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp, npr
      integer i, j, k, ii, mnoff, nn, mm
      real qtmh, ci2, gami, qtmg, dxp, dyp, amx, amy, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, w, vx, vy, vz, p2
      real sfxy, sbxy, sfms, sg
!     dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
!     dimension sfms(nprd,MXV,MYV), sg(14)
      dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      dimension sfms(nprd,mx+1,my+1), sg(14)
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
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(i,j,k,ii,noffp,moffp,nppp,nn,mm,mnoff,x,y,w,vx,vy,vz,dxp,&
!$OMP& dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    &
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg, &
!$OMP& sfxy,sbxy,sfms,sg)
      do 190 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! zero out local accumulators
      do 70 j = 1, my+1
      do 60 i = 1, mx+1
      do 50 ii = 1, npr
      sfms(ii,i,j) = 0.0
   50 continue
   60 continue
   70 continue
! loop over particles in tile
      do 90 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      w = qtmh*dz
! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      acx = vx + x
      acy = vy + y
      acz = vz + w
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
! new velocity
      x = (rot1*acx + rot2*acy + rot3*acz)*anorm + x
      y = (rot4*acx + rot5*acy + rot6*acz)*anorm + y
      w = (rot7*acx + rot8*acy + rot9*acz)*anorm + w
! calculate average momentum
      vx = 0.5*(x + vx)
      vy = 0.5*(y + vy)
      vz = 0.5*(w + vz)
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
      x = amx*amy
      y = dxp*amy
      dx = amx*dyp
      dy = dxp*dyp
      do 80 ii = 1, npr
      sfms(ii,nn,mm) = sfms(ii,nn,mm) + sg(ii)*x
      sfms(ii,nn+1,mm) = sfms(ii,nn+1,mm) + sg(ii)*y
      sfms(ii,nn,mm+1) = sfms(ii,nn,mm+1) + sg(ii)*dx
      sfms(ii,nn+1,mm+1) = sfms(ii,nn+1,mm+1) + sg(ii)*dy
   80 continue
   90 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 120 j = 2, mm
      do 110 i = 2, nn
      do 100 ii = 1, npr
      fms(ii,i+noffp,j+moffp) = fms(ii,i+noffp,j+moffp) + sfms(ii,i,j)
  100 continue
  110 continue
  120 continue
! deposit fluid moments to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 150 i = 2, nn
      do 130 ii = 1, npr
!$OMP ATOMIC
      fms(ii,i+noffp,1+moffp) = fms(ii,i+noffp,1+moffp) + sfms(ii,i,1)
  130 continue
      if (mm > my) then
         do 140 ii = 1, npr
!$OMP ATOMIC
         fms(ii,i+noffp,mm+moffp) = fms(ii,i+noffp,mm+moffp)            &
     & + sfms(ii,i,mm)
  140    continue
      endif
  150 continue
      nn = min(mx+1,nxv-noffp)
      do 180 j = 1, mm
      do 160 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noffp,j+moffp) = fms(ii,1+noffp,j+moffp) + sfms(ii,1,j)
  160 continue
      if (nn > mx) then
         do 170 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noffp,j+moffp) = fms(ii,nn+noffp,j+moffp)            &
     & + sfms(ii,nn,j)
  170    continue
      endif
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FLUIDQS23(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
! for 2-1/2d code, this subroutine calculates fluid quantities from
! fluid moments: density, velocity field, pressure tensor, energy,
! heat flux
! assumes guard cells have been added
! for distributed data, with 1D spatial decomposition
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
! fms(i,j,k) = ith component of fluid moments at grid point j,k
! on exit
! fms(i,j,k) = ith component of fluid quantities at grid point j,k
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! kyp = number of real grids in each field partition in y direction
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
      implicit none
      integer npro, nx, ny, kstrt, kyp, nprd, nxv, nypmx
      real fms
      dimension fms(nprd,nxv,nypmx)
! local data
      integer j, k, kypp, npr
      real at1, at2
      double precision dt1, dt2, dtx, dty, dtz
      kypp = min(kyp,max(0,ny-kyp*(kstrt-1)))
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
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,dt1,dt2,dtx,dty,dtz)
      do 20 k = 1, kypp
      do 10 j = 1, nx
      at1 = fms(1,j,k)
! calculate velocity field
      if (at1.gt.0.0) then
         at2 = 1.0/at1
         fms(2,j,k) = fms(2,j,k)*at2
         fms(3,j,k) = fms(3,j,k)*at2
         fms(4,j,k) = fms(4,j,k)*at2
      else
         fms(2,j,k) = 0.0
         fms(3,j,k) = 0.0
         fms(4,j,k) = 0.0
      endif
      if (npr < 10) go to 10
! calculate pressure tensor
      dt1 = dble(at1)
      dtx = dble(fms(2,j,k))
      dty = dble(fms(3,j,k))
      dtz = dble(fms(4,j,k))
      dt2 = dble(fms(5,j,k))
      fms(5,j,k) = dt2 - dt1*dtx*dtx
      dt2 = dble(fms(6,j,k))
      fms(6,j,k) = dt2 - dt1*dtx*dty
      dt2 = dble(fms(7,j,k))
      fms(7,j,k) = dt2 - dt1*dtx*dtz
      dt2 = dble(fms(8,j,k))
      fms(8,j,k) = dt2 - dt1*dty*dty
      dt2 = dble(fms(9,j,k))
      fms(9,j,k) = dt2 - dt1*dty*dtz
      dt2 = dble(fms(10,j,k))
      fms(10,j,k) = dt2 - dt1*dtz*dtz
! calculate heat flux
      if (npr < 14) go to 10
      dt1 = fms(11,j,k)
      dt2 = dtx*fms(5,j,k) + dty*fms(6,j,k) + dtz*fms(7,j,k)
      dt2 = fms(12,j,k) - dt2 - dtx*dt1
      fms(12,j,k) = dt2
      dt2 = dtx*fms(6,j,k) + dty*fms(8,j,k) + dtz*fms(9,j,k)
      dt2 = fms(13,j,k) - dt2 - dty*dt1
      fms(13,j,k) = dt2
      dt2 = dtx*fms(7,j,k) + dty*fms(9,j,k) + dtz*fms(10,j,k)
      dt2 = fms(14,j,k) - dt2 - dtz*dt1
      fms(14,j,k) = dt2
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FLUIDQS22(fms,npro,nx,ny,kstrt,kyp,nprd,nxv,nypmx)
! for 2d code, this subroutine calculates fluid quantities from fluid
! moments: density, velocity field, pressure tensor, energy, heat flux
! assumes guard cells have been added
! for distributed data, with 1D spatial decomposition
! OpenMP version
! 21 flops/grid, 1 divide, 15 loads, 7 stores
! if all profiles calculated
! input: all, output: fms
! fluid quantities are calculated as follows:
! fms(1,:) = density n is unchanged
! fms(2:3,:) = velocity field u, calculated from momentum and density n:
! u = fms(2:3,:)/n where n /= 0.0
! fms(4:6) = pressure tensor P, calculated from momentum flux, velocity
! field u and density n: P = fms(4:6,:) - n*[u,u]
! fms(7,:) = energy density U is unchanged
! fms(8:9,:) = heat flux Q, calculated from energy U, energy flux,
! velocity field u and pressure P: Q = fms(8:9,:) - u.P - u*U
! on entry:
! fms(i,j,k) = ith component of fluid moments at grid point j,k
! on exit
! fms(i,j,k) = ith component of fluid quantities at grid point j,k
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! kyp = number of real grids in each field partition in y direction
! nprd = maximum number of fluid components, nprd >= 9
! nxv = second dimension of current array, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
      implicit none
      integer npro, nx, ny, kstrt, kyp, nprd, nxv, nypmx
      real fms
      dimension fms(nprd,nxv,nypmx)
! local data
      integer j, k, kypp, npr
      real at1, at2
      double precision dt1, dt2, dtx, dty
      kypp = min(kyp,max(0,ny-kyp*(kstrt-1)))
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 3
      else if (npro==3) then
         npr = 6
      else if (npro==4) then
         npr = 9
      endif
      if (npr > nprd) npr = nprd
! exit if error
      if (npr < 4) return
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,dt1,dt2,dtx,dty)
      do 20 k = 1, kypp
      do 10 j = 1, nx
      at1 = fms(1,j,k)
! calculate velocity field
      if (at1.gt.0.0) then
         at2 = 1.0/at1
         fms(2,j,k) = fms(2,j,k)*at2
         fms(3,j,k) = fms(3,j,k)*at2
      else
         fms(2,j,k) = 0.0
         fms(3,j,k) = 0.0
      endif
      if (npr < 6) go to 10
! calculate pressure tensor
      dt1 = dble(at1)
      dtx = dble(fms(2,j,k))
      dty = dble(fms(3,j,k))
      dt2 = dble(fms(4,j,k))
      fms(4,j,k) = dt2 - dt1*dtx*dtx
      dt2 = dble(fms(5,j,k))
      fms(5,j,k) = dt2 - dt1*dtx*dty
      dt2 = dble(fms(6,j,k))
      fms(6,j,k) = dt2 - dt1*dty*dty
! calculate heat flux
      if (npr < 9) go to 10
      dt1 = fms(7,j,k)
      dt2 = dtx*fms(4,j,k) + dty*fms(5,j,k)
      dt2 = fms(8,j,k) - dt2 - dtx*dt1
      fms(8,j,k) = dt2
      dt2 = dtx*fms(5,j,k) + dty*fms(6,j,k)
      dt2 = fms(9,j,k) - dt2 - dty*dt1
      fms(9,j,k) = dt2
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PSTPTRAJ2(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx, &
     &dvtx,idimp,nppmx,mxyp1,idps,np,nprobt)
! for 2d code, this procedure sets test charge distribution by
! setting a particle id in particle location 5, whose values are between
! 1 and nprobt
! particles stored in segmented array
! input: all, output: iprobt, nprobt
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = particle id of tagged particle n in tile m
! tedges(1) = lower boundary of particle tags
! tedges(2) = upper boundary of particle tags
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
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! idps = number of partition boundaries
! np = total number of particles in part
! nprobt = total number of test charges whose trajectories are stored.
! particle id should be <= 16777215
      implicit none
      integer kstrt, nst, idimp, nppmx, mxyp1, idps, nprobt
      real vtx, vtsx, dvtx
      double precision np
      real ppart, tedges
      dimension ppart(idimp,nppmx,mxyp1)
      dimension tedges(idps)
      integer kpic, iprobt
      dimension kpic(mxyp1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, nppp, noff, joff, it, nt, itt
      real st, at
      double precision dit, ditt
      double precision dnpl, dnp
      dimension dnpl(1), dnp(1)
      if (idimp < 5) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         dit = aint((np-1.0d0)/dble(nprobt)) + 1.0d0
         ditt = 1.0
! find how many particles on left processor
         npp = 0
         do 10 k = 1, mxyp1
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
      do 40 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 30 j = 1, nppp
! clear tag
      ppart(5,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
            ditt = ditt + dit
            itt = int(ditt - dnpl(1))
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         if (kstrt==1) then
            it = ppart(3,j,k)*st + at
            if ((it.gt.0).and.(it.le.nprobt)) then
               if (iprobt(it).eq.0) then
                  nt = nt + 1
                  iprobt(it) = j + joff
                  ppart(5,j,k) = real(nt)
               endif
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(3,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
         endif
      endif
   30 continue
      joff = joff + nppp
   40 continue
! local count of test particles on this processor
      nprobt = nt
! find how many test particles on left processor
      dnp(1) = dble(nt)
      call PPDSCAN(dnp,dnpl,1)
      st = real(dnp(1)) - real(nt)
! tedges(1) = integrated number of test particles on left processor
      tedges(1) = st + 1.0
! tedges(1) = integrated number of test particles on current processor
      tedges(2) = real(dnp(1)) + 1.0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,at)
      do 60 k = 1, mxyp1
      nppp = kpic(k)
! add offset to test particle tags
      do 50 j = 1, nppp
      at = ppart(5,j,k)
      if (at.gt.0.0) then
         ppart(5,j,k) = at + st
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
      subroutine PSTPTRAJ23(ppart,tedges,kpic,iprobt,kstrt,nst,vtx,vtsx,&
     &dvtx,idimp,nppmx,mxyp1,idps,np,nprobt)
! for 2-1/2d code, this procedure sets test charge distribution by
! setting a particle id in particle location 6, whose values are between
! 1 and nprobt
! particles stored in segmented array
! input: all, output: iprobt, nprobt
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(6,n,m) = particle id of tagged particle n in tile m
! tedges(1) = lower boundary of particle tags
! tedges(2) = upper boundary of particle tags
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
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! idps = number of partition boundaries
! np = total number of particles in part
! nprobt = total number of test charges whose trajectories are stored.
! particle id should be <= 16777215
      implicit none
      integer kstrt, nst, idimp, nppmx, mxyp1, idps, nprobt
      real vtx, vtsx, dvtx
      double precision np
      real ppart, tedges
      dimension ppart(idimp,nppmx,mxyp1)
      dimension tedges(idps)
      integer kpic, iprobt
      dimension kpic(mxyp1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, nppp, noff, joff, it, nt, itt
      real st, at
      double precision dit, ditt
      double precision dnpl, dnp
      dimension dnpl(1), dnp(1)
      if (idimp < 6) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         dit = aint((np-1.0d0)/dble(nprobt)) + 1.0d0
         ditt = 1.0
! find how many particles on left processor
         npp = 0
         do 10 k = 1, mxyp1
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
      do 40 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 30 j = 1, nppp
! clear tag
      ppart(6,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(6,j,k) = real(nt)
            ditt = ditt + dit
            itt = int(ditt - dnpl(1))
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         if (kstrt==1) then
            it = ppart(3,j,k)*st + at
            if ((it.gt.0).and.(it.le.nprobt)) then
               if (iprobt(it).eq.0) then
                  nt = nt + 1
                  iprobt(it) = j + joff
                  ppart(6,j,k) = real(nt)
               endif
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(3,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(6,j,k) = real(nt)
         endif
      endif
   30 continue
      joff = joff + nppp
   40 continue
! local count of test particles on this processor
      nprobt = nt
! find how many test particles on left processor
      dnp(1) = dble(nt)
      call PPDSCAN(dnp,dnpl,1)
      st = real(dnp(1)) - real(nt)
! tedges(1) = integrated number of test particles on left processor
      tedges(1) = st + 1.0
! tedges(2) = integrated number of test particles on current processor
      tedges(2) = real(dnp(1)) + 1.0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,nppp,at)
      do 60 k = 1, mxyp1
      nppp = kpic(k)
! add offset to test particle tags
      do 50 j = 1, nppp
      at = ppart(6,j,k)
      if (at.gt.0.0) then
         ppart(6,j,k) = at + st
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
      subroutine PPTRAJ2(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,nprobt&
     &)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, numtp, output: partt, numtp
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! ppart(5,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! partt = tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer numtp, idimp, nppmx, mxyp1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mxyp1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer i, j, k, nppp, nt
      if (idimp < 5) return
      nt = 0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,nppp)
      do 30 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 20 j = 1, nppp
      if (ppart(5,j,k).gt.0.0) then
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
      subroutine PPTRAJ23(ppart,kpic,partt,numtp,idimp,nppmx,mxyp1,     &
     &nprobt)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, numtp, output: partt, numtp
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = position y of particle n in tile m
! ppart(3,n,m) = velocity vx of particle n in tile m
! ppart(4,n,m) = velocity vy of particle n in tile m
! ppart(5,n,m) = velocity vz of particle n in tile m
! ppart(6,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! partt = tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 6
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where mx1 = (system length in x direction - 1)/mx+1
! and where myp1=(partition length in y direction-1)/my+1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer numtp, idimp, nppmx, mxyp1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mxyp1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mxyp1)
! local data
      integer i, j, k, nppp, nt
      if (idimp < 6) return
      nt = 0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,nppp)
      do 30 k = 1, mxyp1
      nppp = kpic(k)
! loop over particles in tile
      do 20 j = 1, nppp
      if (ppart(6,j,k).gt.0.0) then
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
      subroutine PORDTRAJ2(partt,spartt,tedges,numtp,idimp,idps,nprobt)
! this procedure reorders tagged particles in partt to array spartt
! tagged particle coordinates stored in tag order
! input: all except spartt, output: spartt
! partt = tagged particle coordinates
! spartt = reordered tagged particle coordinates
! tedges(1) = lower boundary of particle tags
! tedges(2) = upper boundary of particle tags
! numtp = number of test particles found on this node
! idimp = size of phase space = 5
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
      if (idimp < 5) return
      st = tedges(1) - 1.0
! loop over tagged particles
!$OMP PARALLEL DO PRIVATE(i,j,nt,tn)
      do 20 j = 1, numtp
      tn = partt(5,j)
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
      subroutine PORDTRAJ23(partt,spartt,tedges,numtp,idimp,idps,nprobt)
! this procedure reorders tagged particles in partt to array spartt
! tagged particle coordinates stored in tag order
! input: all except spartt, output: spartt
! partt = tagged particle coordinates
! spartt = reordered tagged particle coordinates
! tedges(1) = lower boundary of particle tags
! tedges(2) = upper boundary of particle tags
! numtp = number of test particles found on this node
! idimp = size of phase space = 6
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
      if (idimp < 6) return
      st = tedges(1) - 1.0
! loop over tagged particles
!$OMP PARALLEL DO PRIVATE(i,j,nt,tn)
      do 20 j = 1, numtp
      tn = partt(6,j)
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
      subroutine PCPYTRAJ2(partt,part,numtp,idimp,nprobt)
! this procedure copies tagged particles in partt to array part
! tagged particle coordinates stored in tag order
! input: all except part, output: part
! spartt = tagged particle coordinates
! spartt = reordered tagged particle coordinates
! numtp = number of test particles found on this node
! idimp = size of phase space = 5
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