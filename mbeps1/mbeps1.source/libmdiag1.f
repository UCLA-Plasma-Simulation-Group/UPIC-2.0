!-----------------------------------------------------------------------
! Fortran Library for diagnostics
! 1D OpenMP PIC Codes:
! CSPECT1 performs frequency analysis of complex time series
! ICSPECT1 performs incremental frequency analysis of complex time
!          series for one time step
! IVCSPECT1 performs incremental frequency analysis of complex vector
!           time series for one time step
! VPDIST1 calculates 1 component velocity distribution, velocity moments,
!         and entropy, for segmented particle array
! VPDIST13 calculates 3 component velocity distribution, velocity moments,
!          and entropy, for segmented particle array
! VBPDIST13 calculates 3d velocity distribution and velocity moments for
!           magnetized plasma, for segmented particle array
! ERPDIST1 calculates 1d energy distribution for relativistic particles
! ERPDIST13 calculates 1-2/2d energy distribution for relativistic
!           particles
! PVSDIST1 for 1d code, calculates 1d velocity distribution, in
!          different regions of space, for segmented particle array
! PVSDIST13 for 1-2/2d code, calculates 1d velocity distribution, in
!           different regions of space, for segmented particle array
! VDIST1 calculates 1 component velocity distribution, velocity moments,
!        and entropy for standard particle array
! VDIST13 calculates 3 component velocity distribution, velocity
!         moments, and entropy.
! VBDIST13 for 1-2/2d code, calculates 3d velocity distribution,
!          and velocity moments for magnetized plasma
! PROFX13L calculates fluid moments from particle quantities: density,
!           momentum, momentum flux, energy, energy flux, assumes
!           particle positions and velocities at same time level
!           for 1-2/2d code
! RPROFX13L calculates fluid moments from relativistic particle
!           quantities: density, velocity, velocity flux, energy,
!           energy flux, assumes particle positions and velocities at 
!           same time level, for 1-2/2d code
! PROFX1L calculates fluid moments from particle quantities: density,
!         momentum, momentum flux, energy, energy flux, assumes
!         particle positions and velocities at same time level
!         for 1d code
! RPROFX1L calculates fluid moments from relativistic particle
!          quantities: density, velocity, velocity flux, energy,
!          energy flux, assumes particle positions and velocities at 
!          same time level, for 1d code
! GPROFX1L calculates fluid moments from particle quantities:
!          density, momentum, momentum flux, energy, energy flux,
!          assumes particle positions and velocities not at same time
!          levels and electrostatic fields
! GRPROFX1L calculates fluid moments from relativistic particle
!           quantities: density, velocity, velocity flux, energy,
!           energy flux, assumes particle positions and velocities
!           not at same time levels and electrostatic fields
! GBPROFX13L calculates fluid moments from particle quantities:
!            density, momentum, momentum flux, energy, energy flux,
!            assumes particle positions and velocities not at same time
!            levels and electromagnetic fields
! FLUIDQS13 calculates fluid quantities from fluid moments:
!           density, velocity field, pressure tensor, energy, heat flux
!           for 1-2/2d code
! STPTRAJ1 sets test charge distribution by setting a particle id
!          in particle location 3 for 1d code
! STPTRAJ13 sets test charge distribution by setting a particle id
!           in particle location 5 for 1-2/2d code
! PTRAJ1 copies tagged particles in ppart to array partt for 1d code
! PTRAJ13 copies tagged particles in ppart to array partt for 1-2/2d
!         code
! FNPTRAJ1 finds how many tagged particles are in ppart for 1d code
! FNPTRAJ13 finds how many tagged particles are in ppart for 1-2/2d code
! STPBEAM1 for 1d code, marks beam particles by setting a particle id in
!          particle location 3
! STPBEAM13 for 1-2/2d code, marks beam particles by setting a particle
!           id in particle location 5
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 4, 2021
!-----------------------------------------------------------------------
      subroutine CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
! this subroutine performs frequency analysis of complex time series,
! pkw(w,k,1:2) = |(1/nt)*sum {fc(t,k)*exp(sqrt(-1)*w*(t-t0))}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed, real and imaginary
! parts stored adjacent, and wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! t0 = starting time value
! dt = time step
! nt = number of input data points, 
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = dimension of input array fc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, ntd, iwd, modesxd
      real t0, dt
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(ntd,modesxd)
! local data
      integer i, j, k
      real anl, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4, cwdt, swdt
      anl = 1.0/real(nt)
! loop over frequencies
      do 30 j = 1, iw
      at3 = wm(j)*dt
      cwdt = dcos(at3)
      swdt = dsin(at3)
! loop over modes
      do 20 k = 1, modesx
      at3 = wm(j)*t0
      at1 = dcos(at3)
      at2 = -dsin(at3)
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      sum4 = 0.0d0
! loop over time
      do 10 i = 1, nt
      fr = real(fc(i,k))
      fi = aimag(fc(i,k))
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
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd&
     &,modesxd)
! this subroutine performs incremental frequency analysis of complex
! time series for one time step
! pkw(w,k,1:2) = |(anorm/nt)*sum {fc(k)*exp(sqrt(-1)*w*(time-t0)}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent, and
! wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! pks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! nx = system length in x direction
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
! norm = 1 for potential or norm = -1 for density gives spectrum as
!        electric field energy
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = dimension of input array fc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, nx, norm, iwd, modesxd
      real time, t0
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(modesxd)
      double precision pks
      dimension pks(4,modesxd,iwd)
! local data
      integer j, k
      real anl, dnx, anorm, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4
      anl = 1.0/real(nt)
      dnx = 6.28318530717959/real(nx)
      anorm = anl
! loop over frequencies
      do 20 j = 1, iw
! loop over modes
      do 10 k = 1, modesx
      at3 = wm(j)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
! add contribution for current time
      fr = real(fc(k))
      fi = aimag(fc(k))
      sum1 = pks(1,k,j) + fr*at1
      sum2 = pks(2,k,j) + fi*at1
      sum3 = pks(3,k,j) + fi*at2
      sum4 = pks(4,k,j) + fr*at2
! save accumulation for next time
      pks(1,k,j) = sum1
      pks(2,k,j) = sum2
      pks(3,k,j) = sum3
      pks(4,k,j) = sum4
! normalize
      if (k.gt.1) then
         if (norm==1) then
            anorm = anl*(dnx*real(k - 1))
         else if (norm==(-1)) then
            anorm = anl/(dnx*real(k - 1))
         endif
      endif
! calculate spectrum for accumulated data
      at1 = anorm*(sum1 - sum3)
      at2 = anorm*(sum2 + sum4)
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anorm*(sum1 + sum3)
      at2 = anorm*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   10 continue
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine IVCSPECT1(fvc,wm,vpkw,vpks,time,t0,nt,iw,modesx,nx,norm&
     &,iwd,modesxd)
! this subroutine performs incremental frequency analysis of complex
! vector time series for one time step
! vpkw(1:2,w,k,1:2) = |(anorm/nt)*sum {fvc(1:2,k)*
!                                      exp(sqrt(-1)*w*(time-t0)}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent, and
! wm(w) contains the (positive) frequencies.
! on output, vpkw(1:2,:,:,1) contains result for positive frequencies,
! vpkw(1:2,:,:,2) the negative frequencies.
! vpks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! nx = system length in x direction
! norm = (-1,0,1) = normalize with (inverse curl,null,curl) op
! norm = 1 for vector potential or norm = -1 for current density gives
!        spectrum as magnetic field energy
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iwd >= iw
! modesxd = second dimension of input array fvc, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, nx, norm, iwd, modesxd
      real time, t0
      real wm, vpkw
      dimension wm(iwd), vpkw(2,modesxd,iwd,2)
      complex fvc
      dimension fvc(2,modesxd)
      double precision vpks
      dimension vpks(2,4,modesxd,iwd)
! local data
      integer i, j, k
      real anl, dnx, anorm, fr, fi
      double precision at1, at2, at3, at4, sum1, sum2, sum3, sum4
      anl = 1.0/real(nt)
      dnx = 6.28318530717959/real(nx)
      anorm = anl
! loop over frequencies
      do 30 j = 1, iw
! loop over modes
      do 20 k = 1, modesx
      at3 = wm(j)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
      do 10 i = 1, 2
! add contribution for current time
      fr = real(fvc(i,k))
      fi = aimag(fvc(i,k))
      sum1 = vpks(i,1,k,j) + fr*at1
      sum2 = vpks(i,2,k,j) + fi*at1
      sum3 = vpks(i,3,k,j) + fi*at2
      sum4 = vpks(i,4,k,j) + fr*at2
! save accumulation for next time
      vpks(i,1,k,j) = sum1
      vpks(i,2,k,j) = sum2
      vpks(i,3,k,j) = sum3
      vpks(i,4,k,j) = sum4
! normalize
      if (k.gt.1) then
         if (norm==1) then
            anorm = anl*(dnx*real(k - 1))
         else if (norm==(-1)) then
            anorm = anl/(dnx*real(k - 1))
         endif
      endif
! calculate spectrum for accumulated data
      at3 = anorm*(sum1 - sum3)
      at4 = anorm*(sum2 + sum4)
      vpkw(i,k,j,1) = at3*at3 + at4*at4
      at3 = anorm*(sum1 + sum3)
      at4 = anorm*(sum2 - sum4)
      vpkw(i,k,j,2) = at3*at3 + at4*at4
   10 continue
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine VPDIST1(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution,
! velocity moments, and entropy
! particles stored in segmented array
! input: all except fvm, output: sfv, fvm
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum velocity (used for scaling) is contained in last element sfv.
! vdrift is contained in fvm(1)
! vth is contained in fvm(2)
! entropy is contained in fvm(3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v*delta_x)).
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = number of particles
! nmvf = dimension of sfv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, np, nmv, nmvf
      real ppart, sfv, fvm
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,mx1+1), fvm(3)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, nmv21, npp, nvx
      real anmv, svx, svxx, vx
      double precision sumvx, sumvx2, ssumvx, ssumvx2, anp, sum1
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/sfv(nmv21+1,mx1+1)
! normalization constant for entropy
      svxx = svx*real(mx1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 1, nmv21
      sfv(j,k) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,nvx,vx,ssumvx,ssumvx2)                           &
!$OMP& REDUCTION(+:sumvx) REDUCTION(+:sumvx2) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumvx = 0.0d0
      ssumvx2 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      vx = ppart(2,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,k) = sfv(nvx,k) + 1.0
   30 continue
! calculate global sums
      sumvx = sumvx + ssumvx
      sumvx2 = sumvx2 + ssumvx2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,k)
   50 continue
      sfv(j,mx1+1) = sum1
   60 continue
! calculate global velocity moments
      anp = 0.0d0
      if (np.gt.0) anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
! count number of particles in global distribution
      sumvx = 0.0d0
      do 70 j = 1, nmv21
      sumvx = sumvx + sfv(j,mx1+1)
   70 continue
! calculate entropy
      sumvx2 = 0.0d0
      do 90 k = 1, mx1
      do 80 j = 1, nmv21
      if (sfv(j,k).gt.0.0) then
         sumvx2 = sumvx2 + sfv(j,k)*dlog(dble(sfv(j,k)*svxx))
      endif
   80 continue
   90 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(3) = sumvx
      return
      end
!-----------------------------------------------------------------------
      subroutine VPDIST13(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf&
     &)
! for 1-2/2d code, this subroutine calculates 3d velocity distribution,
! velocity moments, and entropy
! particles stored in segmented array
! input: all except fvm, output: sfv, fvm
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum velocity (used for scaling) is contained in last element sfv.
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! entropy for i-th dimension is contained in fvm(i,3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v*delta_x)).
! Assumes that distributions in each dimension are independent.
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = number of particles
! nmvf = dimension of sfv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, np, nmv, nmvf
      real ppart, sfv, fvm
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,3,mx1+1), fvm(3,3)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, nmv21, npp, nvx, nvy, nvz
      real anmv, svx, svy, svz, svxx, svyx, svzx, vx, vy, vz
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      double precision ssumvx, ssumvy, ssumvz, ssumvx2, ssumvy2, ssumvz2
      double precision sum1, sum2, sum3
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/sfv(nmv21+1,1,mx1+1)
      svy = anmv/sfv(nmv21+1,2,mx1+1)
      svz = anmv/sfv(nmv21+1,3,mx1+1)
! normalization constant for entropy
      svxx = svx*real(mx1)
      svyx = svy*real(mx1)
      svzx = svz*real(mx1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
      sfv(j,2,k) = 0.0
      sfv(j,3,k) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,nvx,nvy,nvz,vx,vy,vz,ssumvx,ssumvy,ssumvz,       &
!$OMP& ssumvx2,ssumvy2,ssumvz2)                                         &
!$OMP& REDUCTION(+:sumvx) REDUCTION(+:sumvy) REDUCTION(+:sumvz)         &
!$OMP& REDUCTION(+:sumvx2) REDUCTION(+:sumvy2) REDUCTION(+:sumvz2)      &
!$OMP& SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumvx = 0.0d0
      ssumvy = 0.0d0
      ssumvz = 0.0d0
      ssumvx2 = 0.0d0
      ssumvy2 = 0.0d0
      ssumvz2 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      vx = ppart(2,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      vy = ppart(3,j,k)
      nvy = vy*svy + anmv
      ssumvy = ssumvy + vy
      ssumvy2 = ssumvy2 + vy*vy
      if ((nvy.ge.1).and.(nvy.le.nmv21)) sfv(nvy,2,k) = sfv(nvy,2,k)+1.0
      vz = ppart(4,j,k)
      nvz = vz*svz + anmv
      ssumvz = ssumvz + vz
      ssumvz2 = ssumvz2 + vz*vz
      if ((nvz.ge.1).and.(nvz.le.nmv21)) sfv(nvz,3,k) = sfv(nvz,3,k)+1.0
   30 continue
! calculate global sums
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
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,1,k)
      sum2 = sum2 + sfv(j,2,k)
      sum3 = sum3 + sfv(j,3,k)
   50 continue
      sfv(j,1,mx1+1) = sum1
      sfv(j,2,mx1+1) = sum2
      sfv(j,3,mx1+1) = sum3
   60 continue
! calculate global velocity moments
      anp = 0.0d0
      if (np.gt.0) anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      sumvy = sumvy*anp
      sumvz = sumvz*anp
      fvm(1,1) = sumvx
      fvm(2,1) = sumvy
      fvm(3,1) = sumvz
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
      fvm(3,2) = dsqrt(sumvz2*anp - sumvz**2)
! count number of particles in global distribution
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      do 70 j = 1, nmv21
      sumvx = sumvx + sfv(j,1,mx1+1)
      sumvy = sumvy + sfv(j,2,mx1+1)
      sumvz = sumvz + sfv(j,3,mx1+1)
   70 continue
! calculate entropy
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 90 k = 1, mx1
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
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(1,3) = sumvx
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(2,3) = sumvy
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,3) = sumvz
      return
      end
!-----------------------------------------------------------------------
      subroutine VBPDIST13(ppart,kpic,sfv,fvm,omx,omy,omz,idimp,nppmx,  &
     &mx1,np,nmv,nmvf)
! for 1-2/2d code, this subroutine calculates 3d velocity distribution,
! and velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! particles stored in segmented array
! input: all except fvm, output: sfv, fvm
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum velocity (used for scaling) is contained in last element sfv.
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = number of particles
! nmvf = dimension of sfv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, np, nmv, nmvf
      real omx, omy, omz
      real ppart, sfv, fvm
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,3,mx1+1), fvm(3,2)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, nmv21, ndir, npp, nvx, nvz
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real anmv, svx, svz
      double precision sumvx, sumvz, sumvx2, sumvz2, anp
      double precision ssumvx, ssumvz, ssumvx2, ssumvz2
      double precision sum1, sum3
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
      svx = anmv/sfv(nmv21+1,1,mx1+1)
      svz = anmv/sfv(nmv21+1,2,mx1+1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
      sfv(j,2,k) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvz2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,nvx,nvz,vx,vy,vz,at1,at2,ssumvx,ssumvz,ssumvx2,  &
!$OMP& ssumvz2) REDUCTION(+:sumvx) REDUCTION(+:sumvz)                   &
!$OMP& REDUCTION(+:sumvx2) REDUCTION(+:sumvz2) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumvx = 0.0d0
      ssumvz = 0.0d0
      ssumvx2 = 0.0d0
      ssumvz2 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
! vperp1 co-ordinate
      at1 = vx*px + vy*py + vz*pz
! vperp2 co-ordinate
      at2 = vx*qx + vy*qy + vz*qz
! vparallel co-ordinate
      vz = vx*ox + vy*oy + vz*oz
      vx = sqrt(at1*at1 + at2*at2)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,1,k) = sfv(nvx,1,k)+1.0
      nvz = vz*svz + anmv
      ssumvz = ssumvz + vz
      ssumvz2 = ssumvz2 + vz*vz
      if ((nvz.ge.1).and.(nvz.le.nmv21)) sfv(nvz,2,k) = sfv(nvz,2,k)+1.0
   30 continue
! calculate global sums
      sumvx = sumvx + ssumvx
      sumvz = sumvz + ssumvz
      sumvx2 = sumvx2 + ssumvx2
      sumvz2 = sumvz2 + ssumvz2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      sum3 = 0.0d0
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,1,k)
      sum3 = sum3 + sfv(j,2,k)
   50 continue
      sfv(j,1,mx1+1) = sum1
      sfv(j,2,mx1+1) = sum3
   60 continue
! calculate global velocity moments
      anp = 0.0d0
      if (np.gt.0) anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      sumvz = sumvz*anp
      fvm(1,1) = sumvx
      fvm(2,1) = sumvz
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      fvm(2,2) = dsqrt(sumvz2*anp - sumvz**2)
      return
      end
!-----------------------------------------------------------------------
      subroutine ERPDIST1(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv,nmvf)
! for 1d code, this subroutine calculates 1d energy distribution
! for relativistic particles
! the function calculated is of the form g*exp(-e/vth**2), where
! e = (gamma-1)*(c*c) is the kinetic energy per mass, and where
! vth = sqrt(KT/m).  Note vth is not a physical velocity and can be > c
! gamma = sqrt(1 + (p*p)/(c*c)), where p = is the momentum per mass
! g = 1/(de/dp) = gamma/p
! one can express this quantity g as a function of e as follows:
! e = (p*p)/(gamma+1) => p = sqrt((gamma+1)*e), and gamma = 1 + e/(c*c)
! particles stored in segmented array
! input: all except wk, output: sfv, wk
! ppart(2,n,m) = momentum px of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum energy (used for scaling) is contained in last element sfv.
! ci = reciprocal of velocity of light
! wk = total energy is contained in distribution
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nmvf = dimension of sfv
! the number of energy bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, nmv, nmvf
      real ci, wk
      real ppart, sfv
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,mx1+1)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, nmv21, npp, nvx
      real ci2, anmv, svx, px, p2
      double precision sumpx, ssumpx, sum1
      ci2 = ci*ci
! energy scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/sfv(nmv21+1,mx1+1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 1, nmv21
      sfv(j,k) = 0.0
   10 continue
   20 continue
! count particles in each energy region
      anmv = 1.0
      sumpx = 0.0d0
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,npp,nvx,px,p2,ssumpx) REDUCTION(+:sumpx)  &
!$OMP& SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumpx = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      px = ppart(2,j,k)
      p2 = px*px
      px = p2/(1.0 + sqrt(1.0 + p2*ci2))
      nvx = px*svx + anmv
      ssumpx = ssumpx + px
      if ((nvx.ge.1).and.(nvx.le.nmv21)) sfv(nvx,k) = sfv(nvx,k) + 1.0
   30 continue
! calculate global sums
      sumpx = sumpx + ssumpx
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 1, nmv21
      sum1 = 0.0d0
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,k)
   50 continue
      sfv(j,mx1+1) = sum1
   60 continue
! return energy
      wk = sumpx
      return
      end
!-----------------------------------------------------------------------
      subroutine ERPDIST13(ppart,kpic,sfv,ci,wk,idimp,nppmx,mx1,nmv,nmvf&
     &)
! for 1-2/2d code, this subroutine calculates 3d energy distribution,
! for relativistic particles
! the function calculated is of the form g*exp(-e/vth**2), where
! e = (gamma-1)*(c*c) is the kinetic energy per mass, and where
! vth = sqrt(KT/m).  Note vth is not a physical velocity and can be > c
! gamma = sqrt(1 + (p*p)/(c*c)), where p = is the momentum per mass
! g = p*p/(de/dp) = p*gamma
! one can express this quantity g as a function of e as follows:
! e = (p*p)/(gamma+1) => p = sqrt((gamma+1)*e), and gamma = 1 + e/(c*c)
! particles stored in segmented array
! input: all except wk, output: sfv, wk
! ppart(2,n,m) = momentum px of particle n in tile m
! ppart(3,n,m) = momentum py of particle n in tile m
! ppart(4,n,m) = momentum pz of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum energy (used for scaling) is contained in last element sfv.
! ci = reciprocal of velocity of light
! wk = total energy is contained in distribution
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nmvf = dimension of sfv
! the number of energy bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, nmv, nmvf
      real ci, wk
      real ppart, sfv
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,3,mx1+1)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, nmv21, npp, nvx
      real ci2, anmv, svx, px, py, pz, p2
      double precision sumpx, ssumpx, sum1
      ci2 = ci*ci
! energy scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/sfv(nmv21+1,1,mx1+1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 1, nmv21
      sfv(j,1,k) = 0.0
   10 continue
   20 continue
! count particles in each energy region
      anmv = 1.0
      sumpx = 0.0d0
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,npp,nvx,px,py,pz,p2,ssumpx) REDUCTION(+:sumpx)       &
!$OMP& SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumpx = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      px = ppart(2,j,k)
      py = ppart(3,j,k)
      pz = ppart(4,j,k)
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
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,1,k)
   50 continue
      sfv(j,1,mx1+1) = sum1
   60 continue
! return energy
      wk = sumpx
      return
      end
!-----------------------------------------------------------------------
      subroutine PVSDIST1(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1,   &
     &nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution, in
! different regions of space, particles stored in segmented array
! input: all except fvs, output: fvs
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic = number of particles per tile
! fvs = spatially resolved distribution function, number of particles in
! each velocity and spatial range.  maximum velocity (used for scaling)
! is contained in last element of first dimension of fvs
! nmv = number of segments in v for velocity distribution
! mvx = number of grids in x for phase space aggregation
! nxb = number of segments in x for velocity distribution
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx+1
! nmvf = first dimension of fvs
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nmv, mvx, nxb, idimp, nppmx, mx1, nmvf
      real ppart, fvs
      dimension ppart(idimp,nppmx,mx1)
      dimension fvs(nmvf,nxb)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, npp, nmv21, nn, nvx
      real anmv, svx, at1
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fvs(nmv21+1,1)
! spatial scaling
      at1 = 1.0/real(mvx)
! zero out distribution
      do 20 nn = 1, nxb
      do 10 j = 1, nmv21
      fvs(j,nn) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 1.5
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,npp,nn,nvx) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      do 30 j = 1, npp
      nn = ppart(1,j,k)*at1 + 1.0
      nvx = ppart(2,j,k)*svx + anmv
      if ((nvx.ge.1).and.(nvx.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvx,nn) = fvs(nvx,nn) + 1.0
      endif
   30 continue
   40 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PVSDIST13(ppart,kpic,fvs,nmv,mvx,nxb,idimp,nppmx,mx1,  &
     &nmvf)
! for 1-2/2d code, this subroutine calculates 3d velocity distribution,
! in different regions of space, particles stored in segmented array
! input: all except fvs, output: fvs
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! kpic = number of particles per tile
! fvs = spatially resolved distribution function, number of particles in
! each velocity and spatial range.  maximum velocity (used for scaling)
! is contained in last element of last dimension of fvs
! nmv = number of segments in v for velocity distribution
! mvx = number of grids in x for phase space aggregation
! nxb = number of segments in x for velocity distribution
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx+1
! nmvf = first dimension of fvs
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nmv, mvx, nxb, idimp, nppmx, mx1, nmvf
      real ppart, fvs
      dimension ppart(idimp,nppmx,mx1)
      dimension fvs(nmvf,3,nxb)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, npp, nmv21, nn, nvx, nvy, nvz
      real anmv, svx, svy, svz, at1
! velocity scaling, same scaling used for all tiles
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fvs(nmv21+1,1,1)
      svy = anmv/fvs(nmv21+1,2,1)
      svz = anmv/fvs(nmv21+1,3,1)
! spatial scaling
      at1 = 1.0/real(mvx)
! zero out distribution
      do 20 nn = 1, nxb
      do 10 j = 1, nmv21
      fvs(j,1,nn) = 0.0
      fvs(j,2,nn) = 0.0
      fvs(j,3,nn) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 1.5
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,npp,nn,nvx,nvy,nvz) SCHEDULE(dynamic)
      do 40 k = 1, mx1
      npp = kpic(k)
      do 30 j = 1, npp
      nn = ppart(1,j,k)*at1 + 1.0
      nvx = ppart(2,j,k)*svx + anmv
      nvy = ppart(3,j,k)*svy + anmv
      nvz = ppart(4,j,k)*svz + anmv
      if ((nvx.ge.1).and.(nvx.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvx,1,nn) = fvs(nvx,1,nn) + 1.0
      endif
      if ((nvy.ge.1).and.(nvy.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvy,2,nn) = fvs(nvy,2,nn) + 1.0
      endif
      if ((nvz.ge.1).and.(nvz.le.nmv21)) then
!$OMP ATOMIC
         fvs(nvz,3,nn) = fvs(nvz,3,nn) + 1.0
      endif
   30 continue
   40 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution,
! velocity moments, and entropy
! input: all except fvm, output: fv, fvm
! part(2,n) = velocity vx of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv.
! vdrift is contained in fvm(1)
! vth is contained in fvm(2)
! entropy is contained in fvm(3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
! is uniform in space
! idimp = size of phase space = 2
! np = number of particles
! nmvf = dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf), fvm(3)
! local data
      double precision sumvx, sumvx2, anp
      real anmv, svx
      integer j, nmv21, nvx
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j) = 0.0
   10 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 20 j = 1, np
      nvx = part(2,j)*svx + anmv
      sumvx = sumvx + part(2,j)
      sumvx2 = sumvx2 + part(2,j)**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx) = fv(nvx) + 1.0
   20 continue
! calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
! calculate entropy
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 30 j = 1, nmv21
      if (fv(j).gt.0.0) then
         sumvx = sumvx + fv(j)
         sumvx2 = sumvx2 + fv(j)*dlog(dble(fv(j)*svx))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(3) = sumvx
      return
      end
!-----------------------------------------------------------------------
      subroutine VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
! for 1-2/2d code, this subroutine calculates 3d velocity distribution,
! velocity moments, and entropy
! input: all except fvm, output: fv, fvm
! part(2,n) = velocity vx of particle n
! part(3,n) = velocity vy of particle n
! part(4,n) = velocity vz of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv.
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! entropy for i-th dimension is contained in fvm(i,3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
! is uniform in space and distributions in each dimension are
! independent.
! idimp = size of phase space = 4
! np = number of particles
! nmvf = dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf,3), fvm(3,3)
! local data
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      real anmv, svx, svy, svz
      integer j, nmv21, nvx, nvy, nvz
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1,1)
      svy = anmv/fv(nmv21+1,2)
      svz = anmv/fv(nmv21+1,3)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j,1) = 0.0
      fv(j,2) = 0.0
      fv(j,3) = 0.0
   10 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, np
      nvx = part(2,j)*svx + anmv
      sumvx = sumvx + part(2,j)
      sumvx2 = sumvx2 + part(2,j)**2
      nvy = part(3,j)*svy + anmv
      sumvy = sumvy + part(3,j)
      sumvy2 = sumvy2 + part(3,j)**2
      nvz = part(4,j)*svz + anmv
      sumvz = sumvz + part(4,j)
      sumvz2 = sumvz2 + part(4,j)**2
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx,1) = fv(nvx,1) + 1.0
      if ((nvy.ge.1).and.(nvy.le.nmv21)) fv(nvy,2) = fv(nvy,2) + 1.0
      if ((nvz.ge.1).and.(nvz.le.nmv21)) fv(nvz,3) = fv(nvz,3) + 1.0
   20 continue
! calculate velocity moments
      anp = 1.0d0/dble(np)
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
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svx))
      endif
      if (fv(j,2).gt.0.0) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svy))
      endif
      if (fv(j,3).gt.0.0) then
         sumvz = sumvz + fv(j,3)
         sumvz2 = sumvz2 + fv(j,3)*dlog(dble(fv(j,3)*svz))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(1,3) = sumvx
      fvm(2,3) = sumvy
      fvm(3,3) = sumvz
      return
      end
!-----------------------------------------------------------------------
      subroutine VBDIST13(part,fv,fvm,omx,omy,omz,idimp,np,nmv,nmvf)
! for 1-2/2d code, this subroutine calculates 3d velocity distribution,
! and velocity moments for magnetized plasma
! rotating cartesian co-ordinates so that B points in the z direction.
! input: all except fvm, output: fv, fvm
! part(2,n) = velocity vx of particle n
! part(3,n) = velocity vy of particle n
! part(4,n) = velocity vz of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) contained in last element of fv.
! vdrift for i-th dimension is contained in fvm(i,1)
! vth for i-th dimension is contained in fvm(i,2)
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! idimp = size of phase space = 4
! np = number of particles
! nmvf = dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real omx, omy, omz
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf,2), fvm(3,2)
! local data
      integer j, nmv21, ndir, nvx, nvz
      real at1, at2, ox, oy, oz, px, py, pz, qx, qy, qz, vx, vy, vz
      real anmv, svx, svz
      double precision sumvx, sumvz, sumvx2, sumvz2, anp
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
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, np
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
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
! calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(1,2) = dsqrt(sumvx2*anp - sumvx**2)
      sumvz = sumvz*anp
      fvm(2,1) = sumvz
      fvm(2,2) = dsqrt(sumvz2*anp - sumvz**2)
      return
      end
!-----------------------------------------------------------------------
      subroutine PROFX13L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv,  &
     &mx1)
! for 1-2/2d code, this subroutine calculates fluid moments from
! particle quantities: density, momentum, momentum flux, energy,
! energy flux
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 72 flops/particle, 32 loads, 28 stores
! input: all, output: ppart, fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic = number of particles per tile
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, npro, mx, nprd, nxv, mx1
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real dxp, amx, vx, vy, vz, x, w, dx
      real sfms, sg
!     dimension sfms(nprd,MXV), sg(14)
      dimension sfms(nprd,mx+1), sg(14)
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
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,dx,vx,vy,vz,sfms,sg)      &
!$OMP& SCHEDULE(dynamic)
      do 90 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 20 j = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! deposit fluid moments
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      sg(1) = 1.0
      sg(2) = vx
      sg(3) = vy
      sg(4) = vz
      if (npr > 4) then
         x = vx*vx
         sg(5) = x
         sg(6) = vx*vy
         sg(7) = vx*vz
         dx = vy*vy
         w = x + dx
         sg(8) = dx
         sg(9) = vy*vz
         dx = vz*vz
         sg(10) = dx
         w = 0.5*(w + dx)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      do 30 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   30 continue
   40 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 60 j = 2, nn
      do 50 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   50 continue
   60 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 70 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   70 continue
      if (nn > mx) then
         do 80 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   80    continue
      endif
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine RPROFX13L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd,  &
     &nxv,mx1)
! for 1-2/2d code, this subroutine calculates fluid moments from
! particle quantities: density, velocity, velocity flux, energy, energy
! flux, for relativistic particles
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 81 flops/particle, 2 divides, 1 sqrt, 32 loads, 28 stores
! input: all, output: ppart, fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = gami*p2/(1.0 + gami)
! where p2 = px*px + py*py + pz*pz
! where for i = 12, 14, mci = (gami*p2/(1.0 + gami))*vi, where i = x,y,z
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! ppart(3,n,m) = momentum py of particle n in tile m
! ppart(4,n,m) = momentum pz of particle n in tile m
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic = number of particles per tile
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 4
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, npro, mx, nprd, nxv, mx1
      real ci
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real dxp, amx, vx, vy, vz, x, w
      real ci2, p2, gami
      real sfms, sg
!     dimension sfms(nprd,MXV), sg(14)
      dimension sfms(nprd,mx+1), sg(14)
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
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,vx,vy,vz,p2,gami,sfms,sg) &
!$OMP& SCHEDULE(dynamic)
      do 90 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 20 j = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
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
      do 30 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   30 continue
   40 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 60 j = 2, nn
      do 50 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   50 continue
   60 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 70 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   70 continue
      if (nn > mx) then
         do 80 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   80    continue
      endif
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PROFX1L(ppart,fms,kpic,nppmx,idimp,npro,mx,nprd,nxv,mx1&
     &)
! for 1d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 27 flops/particle, 12 loads, 10 stores
! input: all, output: ppart, fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, mci = vi, where i = x
! where for i = 3, mci = vj*vk, where jk = xx
! where for i = 4, mci = (vx**2)
! where for i = 5, mci = (vx**2+)*vi, where i = x
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic = number of particles per tile
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 2
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 5
! nxv = second dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, npro, mx, nprd, nxv, mx1
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real dxp, amx, vx, x, w
      real sfms, sg
!     dimension sfms(nprd,MXV), sg(5)
      dimension sfms(nprd,mx+1), sg(5)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 2
      else if (npro==3) then
         npr = 3
      else if (npro==4) then
         npr = 5
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,vx,sfms,sg)               &
!$OMP& SCHEDULE(dynamic)
      do 90 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 20 j = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! deposit fluid moments
      vx = ppart(2,j,k)
      sg(1) = 1.0
      sg(2) = vx
      if (npr > 2) then
         x = vx*vx
         sg(3) = x
         w = 0.5*x
         sg(4) = w
         sg(5) = w*vx
      endif
      do 30 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   30 continue
   40 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 60 j = 2, nn
      do 50 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   50 continue
   60 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 70 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   70 continue
      if (nn > mx) then
         do 80 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   80    continue
      endif
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine RPROFX1L(ppart,fms,kpic,ci,nppmx,idimp,npro,mx,nprd,nxv&
     &,mx1)
! for 1d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux,
! for relativistic particles
! assumes particle positions and velocities are at the same time level
! using first-order linear interpolation
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 32 flops/particle, 2 divides, 1 sqrt, 12 loads, 10 stores
! input: all, output: ppart, fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, mci = vi, where i = x
! where for i = 3, mci = vj*vk, where jk = xx
! where for i = 4, mci = gami*p2/(1.0 + gami), where p2 = px*px
! where for i = 5, mci = (gami*p2/(1.0 + gami))*vi, where i = x
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = momentum px of particle n in tile m
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic = number of particles per tile
! ci = reciprocal of velocity of light
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 2
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, npro, mx, nprd, nxv, mx1
      real ci
      real ppart, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real dxp, amx, vx, x, w
      real ci2, p2, gami
      real sfms, sg
!     dimension sfms(nprd,MXV), sg(5)
      dimension sfms(nprd,mx+1), sg(5)
      ci2 = ci*ci
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 2
      else if (npro==3) then
         npr = 3
      else if (npro==4) then
         npr = 5
      endif
      if (npr > nprd) npr = nprd
! error if local array is too small
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,vx,p2,gami,sfms,sg)       &
!$OMP& SCHEDULE(dynamic)
      do 90 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! zero out local accumulator
      do 20 j = 1, mx+1
      do 10 ii = 1, npr
      sfms(ii,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 40 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find inverse gamma
      vx = ppart(2,j,k)
      p2 = vx*vx
      gami = 1.0/sqrt(1.0 + p2*ci2)
! deposit fluid moments
      vx = vx*gami
      sg(1) = 1.0
      sg(2) = vx
      if (npr > 2) then
         sg(3) = vx*vx
         w = gami*p2/(1.0 + gami)
         sg(4) = w
         sg(5) = w*vx
      endif
      do 30 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   30 continue
   40 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 60 j = 2, nn
      do 50 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   50 continue
   60 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 70 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   70 continue
      if (nn > mx) then
         do 80 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   80    continue
      endif
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GPROFX1L(ppart,fx,fms,kpic,qbm,dt,idimp,nppmx,npro,nx, &
     &mx,nprd,nxv,mx1)
! for 1d code, this subroutine calculates fluid moments from particle
! quantities: density, momentum, momentum flux, energy, energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 32 flops/particle, 14 loads, 10 stores, if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, mci = vi, where i = x
! where for i = 3, mci = vj*vk, where jk = xx
! where for i = 4, mci = (vx**2)
! where for i = 5, mci = (vx**2+)*vi, where i = x
! velocity equations used are:
! v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
! and x(t+dt) = x(t) + v(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = velocity vx of particle n in partition in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic(k) = number of particles in tile k
! qbm = particle charge/mass
! dt = time interval between successive calculations
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 5
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, npro, nx, mx, nprd, nxv, mx1
      real qbm, dt
      real ppart, fx, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fx(nxv), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, npp, npr
      integer j, k, ii, nn
      real qtmh, dxp, amx, x, w, dx, vx
      real sfx, sfms, sg
!     dimension sfx(MXV), sfms(nprd,MXV), sg(5)
      dimension sfx(mx+1), sfms(nprd,mx+1), sg(5)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 2
      else if (npro==3) then
         npr = 3
      else if (npro==4) then
         npr = 5
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,dx,vx,sfx,sfms,sg)        &
!$OMP& SCHEDULE(dynamic)
      do 100 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      do 10 j = 1, min(mx,nx-noff)+1
      sfx(j) = fx(j+noff)
   10 continue
! zero out local accumulator
      do 30 j = 1, mx+1
      do 20 ii = 1, npr
      sfms(ii,j) = 0.0
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find acceleration
      dx = amx*sfx(nn) + dxp*sfx(nn+1)
! calculate half impulse
      x = qtmh*dx
! half acceleration
      vx = ppart(2,j,k) + x
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      if (npr > 2) then
         x = vx*vx
         sg(3) = x
         w = 0.5*x
         sg(4) = w
         sg(5) = w*vx
      endif
      do 40 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 70 j = 2, nn
      do 60 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   60 continue
   70 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 80 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   80 continue
      if (nn > mx) then
         do 90 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   90    continue
      endif
  100 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRPROFX1L(ppart,fx,fms,kpic,qbm,dt,ci,idimp,nppmx,npro,&
     &nx,mx,nprd,nxv,mx1)
! for 1d code, this subroutine calculates fluid moments from particle
! quantities: density, velocity, velocity flux, energy, energy flux,
! for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 37 flops/particle, 14 loads, 10 stores, 2 divides, 1 sqrt, 
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, mci = vi, where i = x
! where for i = 3, mci = vj*vk, where jk = xx
! where for i = 4, mci = gami*p2/(1.0 + gami), where p2 = px*px
! where for i = 5, mci = (gami*p2/(1.0 + gami))*vi, where i = x
! momentum equations used are:
! px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t))*dt,
! where q/m is charge/mass, and
! x(t+dt) = x(t) + px(t+dt/2)*dtg, where
! dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2))*ci*ci)
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = momentum px of particle n in partition in tile m
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic(k) = number of particles in tile k
! qbm = particle charge/mass
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 5
! nxv = first dimension of field array, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, npro, nx, mx, nprd, nxv, mx1
      real qbm, dt, ci
      real ppart, fx, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fx(nxv), fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, npp, npr
      integer j, k, ii, nn
      real qtmh, ci2, gami, dxp, amx, x, w, dx, vx, p2
      real sfx, sfms, sg
!     dimension sfx(MXV), sfms(nprd,MXV), sg(5)
      dimension sfx(mx+1), sfms(nprd,mx+1), sg(5)
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 2
      else if (npro==3) then
         npr = 3
      else if (npro==4) then
         npr = 5
      endif
      if (npr > nprd) npr = nprd
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,w,dxp,amx,dx,vx,p2,gami,sfx,sfms,sg)&
!$OMP& SCHEDULE(dynamic)
      do 100 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      do 10 j = 1, min(mx,nx-noff)+1
      sfx(j) = fx(j+noff)
   10 continue
! zero out local accumulator
      do 30 j = 1, mx+1
      do 20 ii = 1, npr
      sfms(ii,j) = 0.0
   20 continue
   30 continue
! loop over particles in tile
      do 50 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find acceleration
      dx = amx*sfx(nn) + dxp*sfx(nn+1)
! calculate half impulse
      x = qtmh*dx
! half acceleration
      vx = ppart(2,j,k) + x
! update inverse gamma
      p2 = vx*vx
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate average velocity
      vx = gami*vx
! deposit fluid moments
      sg(1) = 1.0
      sg(2) = vx
      if (npr > 2) then
         sg(3) = vx*vx
         w = gami*p2/(1.0 + gami)
         sg(4) = w
         sg(5) = w*vx
      endif
      do 40 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   40 continue
   50 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 70 j = 2, nn
      do 60 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   60 continue
   70 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 80 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   80 continue
      if (nn > mx) then
         do 90 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
   90    continue
      endif
  100 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,idimp,   &
     &nppmx,npro,nx,mx,nprd,nxv,mx1)
! for 1-2/2d code, this subroutine calculates fluid moments from
! particle quantities: density, momentum, momentum flux, energy,
! energy flux,
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 142 flops/particle, 1 divide, 44 loads, 28 stores
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
! where for i = 1, mci = 1.0
! where for i = 2, 4, mci = vi, where i = x,y,z
! where for i = 5, 10, mci = vj*vk, where jk = xx,xy,xz,yy,yz,zz
! where for i = 11, mci = (vx**2+vy**2+vz**2)
! where for i = 12, 14, mci = (vx**2+vy**2+vz**2)*vi, where i = x,y,z
! velocity equations at t=t+dt/2 are calculated from:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
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
! omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = velocity vx of particle n in tile m at t - dt/2
! ppart(3,n,m) = velocity vy of particle n in tile m at t - dt/2
! ppart(4,n,m) = velocity vz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic(k) = number of particles in tile k
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, npro, nx, mx, nprd, nxv, mx1
      real omx, qbm, dt
      real ppart, fxyz, byz, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real qtmh, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, w, vx, vy, vz
      real sfxyz, sbyz, sfms, sg
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfms(nprd,MXV),, sg(14)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      dimension sfms(nprd,mx+1), sg(14)
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
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,y,z,w,dxp,amx,dx,dy,dz,ox,oy,oz,vx, &
!$OMP& vy,vz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,  &
!$OMP& rot5,rot6,rot7,rot8,rot9,sfxyz,sbyz,sfms,sg) SCHEDULE(dynamic)
      do 110 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
! zero out local accumulators
      do 40 j = 1, mx+1
      do 30 ii = 1, npr
      sfms(ii,j) = 0.0
   30 continue
   40 continue
! loop over particles in tile
      do 60 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
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
         dx = vz*vz
         sg(10) = dx
         w = 0.5*(w + dx)
         sg(11) = w
         sg(12) = w*vx
         sg(13) = w*vy
         sg(14) = w*vz
      endif
      do 50 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   50 continue
   60 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 80 j = 2, nn
      do 70 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   90 continue
      if (nn > mx) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
  100    continue
      endif
  110 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine GRBPROFX13L(ppart,fxyz,byz,fms,kpic,omx,qbm,dt,ci,idimp&
     &,nppmx,npro,nx,mx,nprd,nxv,mx1)
! for 1-2/2d code, this subroutine calculates fluid moments from
! particle quantities: density, velocity, velocity flux, energy,
! energy flux, for relativistic particles
! assumes particle positions and velocities not at same time levels
! using first-order linear interpolation.
! OpenMP version using guard cells
! data deposited in tiles
! particles stored in segmented array
! 160 flops/particle, 4 divides, 2 sqrt, 44 loads, 28 stores
! if all profiles calculated
! input: all, output: fms
! fluid moments are approximated by values at the nearest grid points
! fms(i,n)=mci*(1.-dx) and fms(i,n+1)=mci*dx
! where n = nearest grid point and dx = x-n
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
! omy = (q/m)*by(x(t),y(t))*gami, and omz = (q/m)*bz(x(t),y(t))*gami.
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! similarly for fy(x), fz(x), by(x), bz(x)
! ppart(1,n,m) = position x of particle n in tile m at t
! ppart(2,n,m) = momentum px of particle n in tile m at t - dt/2
! ppart(3,n,m) = momentum py of particle n in tile m at t - dt/2
! ppart(4,n,m) = momentum pz of particle n in tile m at t - dt/2
! fxyz(1,j) = x component of force/charge at grid (j)
! fxyz(2,j) = y component of force/charge at grid (j)
! fxyz(3,j) = z component of force/charge at grid (j)
! that is, convolution of electric field over particle shape
! byz(1,j) = y component of magnetic field at grid (j)
! byz(2,j) = z component of magnetic field at grid (j)
! that is, the convolution of magnetic field over particle shape
! fms(i,j) = ith component of fluid moments at grid (j)
! kpic(k) = number of particles in tile k
! omx = magnetic field electron cyclotron frequency in x
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! npro = (1,2,3,4) = (density,velocity,velocity flux,energy+energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! mx = number of grids in sorting cell in x
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of field arrays, must be >= nx+1
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer idimp, nppmx, npro, nx, mx, nprd, nxv, mx1
      real omx, qbm, dt, ci
      real ppart, fxyz, byz, fms
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension fms(nprd,nxv)
      dimension kpic(mx1)
! local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp, npr
      integer j, k, ii, nn
      real qtmh, ci2, gami, qtmg, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, w, vx, vy, vz, p2
      real sfxyz, sbyz, sfms, sg
!     dimension sfxyz(3,MXV), sbyz(2,MXV)
!     dimension sfms(nprd,MXV),, sg(14)
      dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      dimension sfms(nprd,mx+1), sg(14)
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
!     if (mx.ge.MXV) return
! loop over tiles
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,ii,noff,npp,nn,x,y,z,w,dxp,amx,dx,dy,dz,ox,oy,oz,vx, &
!$OMP& vy,vz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,  &
!$OMP& rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,sfxyz,sbyz,sfms,sg)        &
!$OMP& SCHEDULE(dynamic)
      do 110 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
! load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
! zero out local accumulators
      do 40 j = 1, mx+1
      do 30 ii = 1, npr
      sfms(ii,j) = 0.0
   30 continue
   40 continue
! loop over particles in tile
      do 60 j = 1, npp
! find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
! find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
! find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
! calculate half impulse
      x = qtmh*dx
      y = qtmh*dy
      z = qtmh*dz
! half acceleration
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
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
! new velocity
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
      do 50 ii = 1, npr
      sfms(ii,nn) = sfms(ii,nn) + sg(ii)*amx
      sfms(ii,nn+1) = sfms(ii,nn+1) + sg(ii)*dxp
   50 continue
   60 continue
! deposit fluid moments to interior points in global array
      nn = min(mx,nxv-noff)
      do 80 j = 2, nn
      do 70 ii = 1, npr
      fms(ii,j+noff) = fms(ii,j+noff) + sfms(ii,j)
   70 continue
   80 continue
! deposit fluid moments to edge points in global array
      nn = min(mx+1,nxv-noff)
      do 90 ii = 1, npr
!$OMP ATOMIC
      fms(ii,1+noff) = fms(ii,1+noff) + sfms(ii,1)
   90 continue
      if (nn > mx) then
         do 100 ii = 1, npr
!$OMP ATOMIC
         fms(ii,nn+noff) = fms(ii,nn+noff) + sfms(ii,nn)
  100    continue
      endif
  110 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FLUIDQS13(fms,npro,nx,nprd,nxv)
! for 1-2/2d code, this subroutine calculates fluid quantities from
! fluid moments: density, velocity field, pressure tensor, energy,
! heat flux
! assumes guard cells have been added
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
! fms(i,j) = ith component of fluid moments at grid point j
! on exit
! fms(i,j) = ith component of fluid quantities at grid point j
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! nprd = maximum number of fluid components, nprd >= 14
! nxv = second dimension of current array, must be >= nx+1
      implicit none
      integer npro, nx, nprd, nxv
      real fms
      dimension fms(nprd,nxv)
! local data
      integer j, npr
      real at1, at2
      double precision dt1, dt2, dtx, dty, dtz
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
      do 10 j = 1, nx
      at1 = fms(1,j)
! calculate velocity field
      if (at1.gt.0.0) then
         at2 = 1.0/at1
         fms(2,j) = fms(2,j)*at2
         fms(3,j) = fms(3,j)*at2
         fms(4,j) = fms(4,j)*at2
      else
         fms(2,j) = 0.0
         fms(3,j) = 0.0
         fms(4,j) = 0.0
      endif
      if (npr < 10) go to 10
! calculate pressure tensor
      dt1 = dble(at1)
      dtx = dble(fms(2,j))
      dty = dble(fms(3,j))
      dtz = dble(fms(4,j))
      dt2 = dble(fms(5,j))
      fms(5,j) = dt2 - dt1*dtx*dtx
      dt2 = dble(fms(6,j))
      fms(6,j) = dt2 - dt1*dtx*dty
      dt2 = dble(fms(7,j))
      fms(7,j) = dt2 - dt1*dtx*dtz
      dt2 = dble(fms(8,j))
      fms(8,j) = dt2 - dt1*dty*dty
      dt2 = dble(fms(9,j))
      fms(9,j) = dt2 - dt1*dty*dtz
      dt2 = dble(fms(10,j))
      fms(10,j) = dt2 - dt1*dtz*dtz
! calculate heat flux
      if (npr < 14) go to 10
      dt1 = fms(11,j)
      dt2 = dtx*fms(5,j) + dty*fms(6,j) + dtz*fms(7,j)
      dt2 = fms(12,j) - dt2 - dtx*dt1
      fms(12,j) = dt2
      dt2 = dtx*fms(6,j) + dty*fms(8,j) + dtz*fms(9,j)
      dt2 = fms(13,j) - dt2 - dty*dt1
      fms(13,j) = dt2
      dt2 = dtx*fms(7,j) + dty*fms(9,j) + dtz*fms(10,j)
      dt2 = fms(14,j) - dt2 - dtz*dt1
      fms(14,j) = dt2
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FLUIDQS1(fms,npro,nx,nprd,nxv)
! for 1d code, this subroutine calculates fluid quantities from fluid
! moments: density, velocity field, pressure tensor, energy, heat flux
! assumes guard cells have been added
! OpenMP version
! 21 flops/grid, 1 divide, 15 loads, 7 stores
! if all profiles calculated
! input: all, output: fms
! fluid quantities are calculated as follows:
! fms(1,:) = density n is unchanged
! fms(2,:) = velocity field u, calculated from momentum and density n:
! u = fms(2,:)/n where n /= 0.0
! fms(3) = pressure tensor P, calculated from momentum flux, velocity
! field u and density n: P = fms(3,:) - n*[u,u]
! fms(4,:) = energy density U is unchanged
! fms(5,:) = heat flux Q, calculated from energy U, energy flux,
! velocity field u and pressure P: Q = fms(5,:) - u.P - u*U
! on entry:
! fms(i,j) = ith component of fluid moments at grid point j
! on exit
! fms(i,j) = ith component of fluid quantities at grid point j
! npro = (1,2,3,4) = (density,momentum,momentum flux,energy flux)
! if npro = n is selected, all profiles less than n are also calculated
! nx = system length in x direction
! nprd = maximum number of fluid components, nprd >= 5
! nxv = second dimension of current array, must be >= nx+1
      implicit none
      integer npro, nx, nprd, nxv
      real fms
      dimension fms(nprd,nxv)
! local data
      integer j, npr
      real at1, at2
      double precision dt1, dt2, dtx
! calculate number of components required
      npr = 0
      if (npro==1) then
         npr = 1
      else if (npro==2) then
         npr = 2
      else if (npro==3) then
         npr = 3
      else if (npro==4) then
         npr = 5
      endif
      if (npr > nprd) npr = nprd
! exit if error
      if (npr < 2) return
      do 10 j = 1, nx
      at1 = fms(1,j)
! calculate velocity field
      if (at1.gt.0.0) then
         at2 = 1.0/at1
         fms(2,j) = fms(2,j)*at2
      else
         fms(2,j) = 0.0
      endif
      if (npr < 6) go to 10
! calculate pressure tensor
      dt1 = dble(at1)
      dtx = dble(fms(2,j))
      dt2 = dble(fms(3,j))
      fms(3,j) = dt2 - dt1*dtx*dtx
! calculate heat flux
      if (npr < 5) go to 10
      dt1 = fms(4,j)
      dt2 = dtx*fms(3,j)
      dt2 = fms(5,j) - dt2 - dtx*dt1
      fms(5,j) = dt2
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine STPTRAJ1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,    &
     &nppmx,mx1,np,nprobt)
! for 1d code, this procedure sets test charge distribution by setting
! a particle id in particle location 3, whose values are between 1 and
! nprobt
! particles stored in segmented array
! input: all, output: iprobt, nprobt
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! iprobt = scratch array of size nprobt, used for nst = 2
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! idimp = size of phase space = 3
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = total number of particles in part
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer nst, idimp, nppmx, mx1, np, nprobt
      real vtx, vtsx, dvtx
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic, iprobt
      dimension kpic(mx1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, joff, it, nt, itt
      real st, at
      if (idimp < 3) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         it = np/nprobt
         itt = 1
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         st = 4.0*vtx
         it = nprobt/2
         at = real(it)
         st = at/st
         at = at + 1.0 + 0.5*(nprobt - 2*it)
         do 10 j = 1, nprobt
         iprobt(j) = 0
   10    continue
! velocity slice in x direction
      else if (nst.eq.3) then
         st = 1.0/dvtx
         itt = vtsx*st + 0.5
      endif
      joff = 0
      nt = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
! clear tag
      ppart(3,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(3,j,k) = real(nt)
            itt = itt + it
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         it = ppart(2,j,k)*st + at
         if ((it.gt.0).and.(it.le.nprobt)) then
            if (iprobt(it).eq.0) then
               nt = nt + 1
               iprobt(it) = j + joff
               ppart(3,j,k) = real(nt)
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(2,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(3,j,k) = real(nt)
         endif
      endif
   20 continue
      joff = joff + npp
   30 continue
      nprobt = nt
      return
      end
!-----------------------------------------------------------------------
      subroutine STPTRAJ13(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,   &
     &nppmx,mx1,np,nprobt)
! for 1-2/2d code, this procedure sets test charge distribution by
! setting a particle id in particle location 5, whose values are between
! 1 and nprobt
! particles stored in segmented array
! input: all, output: iprobt, nprobt
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! iprobt = scratch array of size nprobt, used for nst = 2
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = total number of particles in part
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer nst, idimp, nppmx, mx1, np, nprobt
      real vtx, vtsx, dvtx
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic, iprobt
      dimension kpic(mx1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, joff, it, nt, itt
      real st, at
      if (idimp < 5) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         it = np/nprobt
         itt = 1
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         st = 4.0*vtx
         it = nprobt/2
         at = real(it)
         st = at/st
         at = at + 1.0 + 0.5*(nprobt - 2*it)
         do 10 j = 1, nprobt
         iprobt(j) = 0
   10    continue
! velocity slice in x direction
      else if (nst.eq.3) then
         st = 1.0/dvtx
         itt = vtsx*st + 0.5
      endif
      joff = 0
      nt = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
! clear tag
      ppart(5,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
            itt = itt + it
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         it = ppart(2,j,k)*st + at
         if ((it.gt.0).and.(it.le.nprobt)) then
            if (iprobt(it).eq.0) then
               nt = nt + 1
               iprobt(it) = j + joff
               ppart(5,j,k) = real(nt)
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(2,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
         endif
      endif
   20 continue
      joff = joff + npp
   30 continue
      nprobt = nt
      return
      end
!-----------------------------------------------------------------------
      subroutine FNPTRAJ1(ppart,kpic,idimp,nppmx,mx1,nprobt)
! this procedure finds how many tagged particles are in ppart
! input: all except nprobt, output: nprobt
! ppart(3,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! idimp = size of phase space = 3
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nprobt = number of test charges whose trajectories are stored.
      implicit none
      integer idimp, nppmx, mx1, nprobt
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, npp, nt
      nprobt = 0
      if (idimp < 3) return
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,npp,nt) REDUCTION(+:nprobt)
      do 20 k = 1, mx1
      npp = kpic(k)
      nt = 0
! loop over particles in tile
      do 10 j = 1, npp
      if (ppart(3,j,k).gt.0.0) then
         nt = nt + 1
      endif
   10 continue
      nprobt = nprobt + nt
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine FNPTRAJ13(ppart,kpic,idimp,nppmx,mx1,nprobt)
! this procedure finds how many tagged particles are in ppart
! input: all except nprobt, output: nprobt
! ppart(5,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nprobt = number of test charges whose trajectories are stored.
      implicit none
      integer idimp, nppmx, mx1, nprobt
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, npp, nt
      nprobt = 0
      if (idimp < 5) return
! loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,npp,nt) REDUCTION(+:nprobt)
      do 20 k = 1, mx1
      npp = kpic(k)
      nt = 0
! loop over particles in tile
      do 10 j = 1, npp
      if (ppart(5,j,k).gt.0.0) then
         nt = nt + 1
      endif
   10 continue
      nprobt = nprobt + nt
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PTRAJ1(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, output: partt
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! idimp = size of phase space = 3
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer idimp, nppmx, mx1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mx1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mx1)
! local data
      integer i, j, k, npp, nt
      real tn
      if (idimp < 3) return
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,npp,tn)
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
      tn = ppart(3,j,k)
      if (tn.gt.0.0) then
         nt = tn
         do 10 i = 1, idimp
         partt(i,nt) = ppart(i,j,k)
   10    continue
      endif
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PTRAJ13(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, output: partt
! ppart(1,n,m) = position x of particle n in tile m
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = velocity vy of particle n in tile m
! ppart(4,n,m) = velocity vz of particle n in tile m
! ppart(5,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! partt = tagged particle coordinates
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer idimp, nppmx, mx1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mx1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mx1)
! local data
      integer i, j, k, npp, nt
      real tn
      if (idimp < 5) return
! loop over tiles
!$OMP PARALLEL DO PRIVATE(i,j,k,npp,tn)
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
      tn = ppart(5,j,k)
      if (tn.gt.0.0) then
         nt = tn
         do 10 i = 1, idimp
         partt(i,nt) = ppart(i,j,k)
   10    continue
      endif
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine STPBEAM1(part,npx,idimp,nop)
! for 1d code, this procedure marks initial beam distribution by setting
! a particle id in particle location 3, whose values are negative for
! particle n if n > npx, zero otherwise
! part(3,n) = particle id of tagged particle n
! npx = number of background particles distributed in x direction
! idimp = size of phase space = 3
! nop = number of particles
      implicit none
      integer npx, idimp, nop
      real part
      dimension part(idimp,nop)
! local data
      integer j, js
      if (idimp < 3) return
! zero tag for initial background particles
      do 10 j = 1, npx
      part(3,j) = 0.0
   10 continue
! negative tag for initial beam particles
      js = npx + 1
      do 20 j = js, nop
      part(3,j) = -1.0
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine STPBEAM13(part,npx,idimp,nop)
! for 1-2/2d code, this procedure marks initial beam distribution by
! setting a particle id in particle location 5, whose values are
! negative for particle n if n > npx, zero otherwise
! part(5,n) = particle id of tagged particle n
! npx = number of background particles distributed in x direction
! idimp = size of phase space = 5
! nop = number of particles
      implicit none
      integer npx, idimp, nop
      real part
      dimension part(idimp,nop)
! local data
      integer j, js
      if (idimp < 5) return
! zero tag for initial background particles
      do 10 j = 1, npx
      part(5,j) = 0.0
   10 continue
! negative tag for initial beam particles
      js = npx + 1
      do 20 j = js, nop
      part(5,j) = -1.0
   20 continue
      return
      end
