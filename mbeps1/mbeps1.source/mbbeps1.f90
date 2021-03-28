!-----------------------------------------------------------------------
! 1-2/2D Electromagnetic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
      program mbbeps1
      use f1
      use fb1
      use graf1
      use omplib
      implicit none
      integer :: nn, ierr
!
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      idimp = 4; ipbc = 1
!
! override default input data
      emf = 1
      relativity = 1
! read namelist
      call readnml1(iuin)
! override input data
      idcode = 2
      ndim = 3
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
!
! text output file
      fname = 'output1.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
!
      irc = 0
! nvp = number of shared memory nodes (0=default)
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call INIT_OMP(nvp)
!
! open graphics device
      irc = open_graphs(nplot)
!
! initialize scalars for standard code
! increase number of coordinates for particle tag
      if ((ntt > 0).or.((nts > 0).and.(ntsc > 0))) then
         idimp = idimp + 1
      endif
! np = total number of electrons in simulation
! nx = number of grid points in x direction
      np = npx + npxb; nx = 2**indx; nxh = nx/2
! npi = total number of ions in simulation
      if (movion > 0) npi = npxi + npxbi
      nxe = nx + 2; nxeh = nxe/2
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! nstart = initial time loop index
! ntime = current time step
      nloop = tend/dt + .0001; nstart = 1; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
         vtdzi = vtdz/sqrt(rmass*rtempdzi)
      endif
      dth = 0.0
!
! check for unimplemented features
      if (ipbc.ne.1) plist = .false.
!
! allocate field data for standard code:
! qe/qi = electron/ion charge density with guard cells
! fxe = smoothed electric field with guard cells
! ffc = form factor array for poisson solver
! mixup = bit reverse table for FFT
! sct = sine/cosine table for FFT
! cue = electron current density with guard cells
! fxyze/byze = smoothed electric/magnetic field with guard cells
! eyz/byz = transverse electric/magnetic field in fourier space
      call init_fields13()
!
! prepare fft tables: updates mixup, sct
      call mfft1_init(mixup,sct,indx)
! calculate form factors: updates ffc
      call mpois1_init(ffc,ax,affp,nx)
! initialize different ensemble of random numbers
      if (nextrand > 0) call mnextran1(nextrand,ndim,np+npi)
!
! open reset and restart files
      call open_restart1()
!
! new start
      if (nustrt==1) then
! initialize electrons: updates ppart, kpic
! ppart = tiled electron particle arrays
! kpic = number of electrons in each tile
         call init_electrons13()
!
! initialize background charge density: updates qi
         if (movion==0) then
            qi = 0.0
            call mpost1(ppart,qi,kpic,-qme,tdpost,mx)
            call maguard1(qi,tguard,nx)
         endif
!
! initialize ions: updates pparti, kipic
! pparti = tiled on particle arrays
! kipic = number of ions in each tile
! cui = ion current density with guard cells
         if (movion==1) then
            call init_ions13()
         endif
!
! initialize transverse electromagnetic fields
         eyz = cmplx(0.0,0.0)
         byz = cmplx(0.0,0.0)
!
! restart to continue a run which was interrupted
      else if (nustrt==2) then
         call bread_restart13(iur)
         if ((ntime+ntime0)> 0) dth = 0.5*dt
         nstart = ntime + 1
! start a new run with data from a previous run
      else if (nustrt==0) then
         call bread_restart13(iur0)
         if ((ntime+ntime0)> 0) dth = 0.5*dt
      endif
!
! initialize current fields
      cue = 0.0
!
! set magnitude of external transverse magnetic field
      omt = sqrt(omy*omy + omz*omz)
!
! reverse simulation at end back to start
      if (treverse==1) nloop = 2*nloop
!
! initialize all diagnostics from namelist input parameters
! wt = energy time history array=
! pkw = power spectrum for potential
! pkwdi = power spectrum for ion density
! wk = maximum frequency as a function of k for potential
! wkdi = maximum frequency as a function of k for ion density
! fmse/fmsi = electron/ion fluid moments
! fv/fvi = global electron/ion velocity distribution functions
! fvm/fvmi = electron/ion vdrift, vth, entropy for global distribution
! fvtm/fvtmi = time history of electron/ion vdrift, vth, and entropy
! fvtp = velocity distribution function for test particles
! fvmtp = vdrift, vth, and entropy for test particles
! partd = trajectory time history array
! vpkwji = power spectrum for ion current density
! vwkji = maximum frequency as a function of k for ion current
! vpkwr = power spectrum for radiative vector potential
! vwkr = maximum frequency as a function of k for radiative vector
!        potential
! vpkw = power spectrum for vector potential
! vwk = maximum frequency as a function of k for vector potential
! vpkwet = power spectrum for transverse efield
! vwket = maximum frequency as a function of k for transverse efield
! oldcu = previous current density with guard cells
      call initialize_diagnostics13(ntime)
!
! read in restart diagnostic file to continue interrupted run
      if (nustrt==2) call dread_restart13(iur)
!
! write reset file
!     call bwrite_restart13(iurr,ntime)
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (dt > 0.63*ci) then
         write (*,*) 'Info: Courant condition may be exceeded!'
      endif
!
      write (iuot,*) 'program mbbeps1'
!
! debug reset
!  10 if (irc==6) then
!        irc = 0
!        call bread_restart13(iurr)
!        call reset_diags13()
!        dth = 0.0
!        cue = 0.0
!     endif
!
! * * * start main iteration loop * * *
!
      do n = nstart, nloop 
      ntime = n - 1
      write (iuot,*) 'ntime = ', ntime
!
! debug reset
!     if (ntime==nloop/2) then
!        irc = 6
!        go to 10
!     endif
!
! save previous current in fourier space for radiative vector potential
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) oldcu = cue
      endif
!
! fluid moments diagnostic
      if (ntfm > 0) then
         it = ntime/ntfm
         if (ntime==ntfm*it) then
! updates fmse
            call efluidms_diag13(fmse)
            if (movion==1) then
! updates fmsi
               call ifluidms_diag13(fmsi)
            endif
         endif
      endif
!
! deposit electron current with OpenMP: updates ppart, kpic, and cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call deposit_ecurrent13(ppart,kpic)
!
! electron current density diagnostic:
! updates vfield=electron current
      if (ntje > 0) then
         it = ntime/ntje
         if (ntime==ntje*it) then
            call ecurrent_diag13(vfield)
! display smoothed electron current
            call dvector1(vfield,' ELECTRON CURRENT',ntime,999,0,2,nx,  &
     &irc)
            if (irc==1) exit; irc = 0
         endif
      endif
!
! deposit ion current with OpenMP: updates pparti, kipic, and cuie
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         cui = 0.0
         call dtimer(dtime,itime,1)
         tdjpost = tdjpost + real(dtime)
         call deposit_icurrent13(pparti,kipic)
      endif
!
! ion current density diagnostic:
! updates vfield=ion current, vpkwji, vwkji
      if (movion==1) then
         if (ntji > 0) then
            it = ntime/ntji
            if (ntime==ntji*it) then
               call icurrent_diag13(vfield,vpkwji,vwkji,ntime)
               if ((ndji==1).or.(ndji==3)) then
! display smoothed ion current
                  call dvector1(vfield,' ION CURRENT',ntime,999,0,2,nx, &
     &irc)
                  if (irc==1) exit; irc = 0
               endif
! ion spectral analysis
               if ((ndji==2).or.(ndji==3)) then
! display frequency spectrum
                  call dmvector1(vwkji,'ION CURRENT OMEGA VS MODE',ntime&
     &,999,2,2,modesxji,cwk,irc)
                  if (irc==1) exit; irc = 0
               endif
            endif
         endif
      endif
!
! deposit electron charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
! add guard cells: updates qe
      call maguard1(qe,tguard,nx)
!
! electron density diagnostic: updates sfield=electron density
      if (ntde > 0) then
         it = ntime/ntde
         if (ntime==ntde*it) then
            call edensity_diag1(sfield)
! display smoothed electron density
            call dscaler1(sfield,' EDENSITY',ntime,999,1,nx,irc)
            if (irc==1) exit; irc = 0
         endif
      endif
!
! deposit ion charge with OpenMP: updates qi
      if (movion==1) then
         call dtimer(dtime,itime,-1)
         qi = 0.0
         call dtimer(dtime,itime,1)
         tdpost = tdpost + real(dtime)
         call mpost1(pparti,qi,kipic,qmi,tdpost,mx)
! add guard cells: updates qi
         call maguard1(qi,tguard,nx)
      endif
!
! ion density diagnostic: updates sfield=ion density, pkwdi, wkdi
      if (movion==1) then
         if (ntdi > 0) then
            it = ntime/ntdi
            if (ntime==ntdi*it) then
               call idensity_diag1(sfield,pkwdi,wkdi,ntime)
               if ((nddi==1).or.(nddi==3)) then
! display smoothed ion density
                  call dscaler1(sfield,' ION DENSITY',ntime,999,1,nx,irc&
     &)
                  if (irc==1) exit; irc = 0
               endif
! ion spectral analysis
               if ((nddi==2).or.(nddi==3)) then
! display frequency spectrum
                  call dmscaler1(wkdi,'ION DENSITY OMEGA VS MODE',ntime,&
     &999,1,modesxdi,cwk,irc)
                  if (irc==1) exit; irc = 0
               endif
            endif
         endif
      endif
!
! add electron and ion densities: updates qe
      call maddqei1(qe,qi,tfield,nx)
!
! add electron and ion current densities: updates cue
      if (movion==1) call maddcuei1(cue,cui,tfield,nx)
!
! transform charge to fourier space: updates qe
      isign = -1
      call mfft1r(qe,isign,mixup,sct,tfft,indx)
!
! transform current to fourier space: updates cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! radiative vector potential diagnostic:
! updates vfield=radiative vector potential, vpkwr, vwkr
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) then
            call vrpotential_diag13(vfield,vpkwr,vwkr,ntime)
            if ((ndar==1).or.(ndar==3)) then
! display radiative vector potential
               call dvector1(vfield,' RADIATIVE VPOTENTIAL',ntime,999,0,&
     &2,nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndar==2).or.(ndar==3)) then
! display frequency spectrum
               call dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',&
     &ntime,999,2,2,modesxar,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! calculate electromagnetic fields in fourier space: updates eyz, byz
      if ((ntime+ntime0)==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mibpois1(cue,byz,ffc,ci,wb,tfield,nx)
         wf = 0.0
! calculate initial darwin electric field with approximate value
         call init_detfield13()
         dth = 0.5*dt
      else
! finite-difference solver
!        call mmaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wb,tfield,nx)
! analytic solver
         call mamaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wb,tfield,nx)
      endif
!
! calculate longitudinal force/charge in fourier space:
! updates fxe, we
      call mpois1(qe,fxe,ffc,we,tfield,nx)
!
! add longitudinal and transverse electric fields: updates fxyze
      call memfield1(fxyze,fxe,eyz,ffc,tfield,nx)
! copy magnetic field: updates byze
      call mbmfield1(byze,byz,ffc,tfield,nx)
!
! transform electric force to real space: updates fxyze
      isign = 1
      call mfft1rn(fxyze,isign,mixup,sct,tfft,indx)
!
! transform magnetic force to real space: updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mbaddext1(byze,tfield,omy,omz,nx)
!
! add external traveling wave field
      ts = dt*real(ntime)
      call meaddext13(fxyze,tfield,amodex,freq,ts,trmp,toff,el0,er0,ey0,&
     &ez0,nx)
!
! copy guard cells: updates fxyze, byze
      call mcguard1(fxyze,tguard,nx)
      call mcguard1(byze,tguard,nx)
!
! potential diagnostic: updates sfield=potential, pkw, wk
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
            call potential_diag1(sfield,pkw,wk,ntime)
            if ((ndp==1).or.(ndp==3)) then
! display potential
               call dscaler1(sfield,' POTENTIAL',ntime,999,0,nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndp==2).or.(ndp==3)) then
! display frequency spectrum
               call dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2, &
     &modesxp,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! longitudinal efield diagnostic: updates sfield=longitudinal efield
      if (ntel > 0) then
         it = ntime/ntel
         if (ntime==ntel*it) then
            call elfield_diag1(sfield)
! display longitudinal efield 
            call dscaler1(sfield,' ELFIELD',ntime,999,0,nx,irc)
            if (irc==1) exit; irc = 0
         endif
      endif
!
! vector potential diagnostic:updates vfield=vector potential, vpkw, vwk
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
            call vpotential_diag13(vfield,vpkw,vwk,ntime)
            if ((nda==1).or.(nda==3)) then
! display vector potential
               call dvector1(vfield,' VECTOR POTENTIAL',ntime,999,0,2,nx&
     &,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((nda==2).or.(nda==3)) then
! display frequency spectrum
               call dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime&
     &,999,2,2,modesxa,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! transverse efield diagnostic:
! updates vfield=transverse efield, vpkwet, vwket
      if (ntet > 0) then
         it = ntime/ntet
         if (ntime==ntet*it) then
            call etfield_diag13(vfield,vpkwet,vwket,ntime)
            if ((ndet==1).or.(ndet==3)) then
! display transverse efield
               call dvector1(vfield,' TRANSVERSE EFIELD',ntime,999,0,2, &
     &nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndet==2).or.(ndet==3)) then
! display frequency spectrum
               call dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',  &
     &ntime,999,2,2,modesxet,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! magnetic field diagnostic: updates vfield=bfield
      if (ntb > 0) then
         it = ntime/ntb
         if (ntime==ntb*it) then
            call bfield_diag13(vfield)
! display magnetic field
            call dvector1(vfield,' MAGNETIC FIELD',ntime,999,0,2,nx,irc)
            if (irc==1) exit; irc = 0
         endif
      endif
!
! velocity diagnostic
      if (ntv > 0) then
         it = ntime/ntv
         if (ntime==ntv*it) then
! updates fv, fe, fvm, fvtm, wkt
            call evelocity_diag13(fv,fe,fvm,fvtm,wkt)
! display electron velocity distributions
            if ((ndv==1).or.(ndv==3)) then
               if ((nvft==1).or.(nvft==3)) then
                  call displayfv1(fv,fvm,' ELECTRON',ntime,nmv,2,irc)
                  if (irc==1) exit; irc = 0
               endif
! display electron velocity distributions in cylindrical co-ordinates
               if ((nvft==4).or.(nvft==5)) then   
                  call displayfvb1(fv,fvm,' ELECTRON',ntime,nmv,2,irc)  
                  if (irc==1) exit; irc = 0    
               endif
! display electron energy distribution
               if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                  call displayfe1(fe,wkt,' ELECTRON',ntime,nmv,irc)
                  if (irc==1) exit; irc = 0
               endif
            endif
! ion distribution function
            if (movion==1) then
! updates fvi, fei, fvmi, fvtmi, wkt
               call ivelocity_diag13(fvi,fei,fvmi,fvtmi,wkt)
! display ion velocity distributions
               if ((ndv==2).or.(ndv==3)) then
                  if ((nvft==1).or.(nvft==3)) then
                     call displayfv1(fvi,fvmi,' ION',ntime,nmv,2,irc)
                     if (irc==1) exit; irc = 0
                  endif
! display ion velocity distributions in cylindrical co-ordinates
                  if ((nvft==4).or.(nvft==5)) then         
                     call displayfvb1(fvi,fvmi,' ION',ntime,nmv,2,irc) 
                     if (irc==1) exit; irc = 0 
                  endif
! display ion energy distribution
                  if ((nvft==2).or.(nvft==3).or.(nvft==5)) then
                     ts = fei(2*nmv+2,1)
                     fei(2*nmv+2,1) = rmass*ts
                     call displayfe1(fei,wkt,' ION',ntime,nmv,irc)
                     if (irc==1) exit; irc = 0
                     fei(2*nmv+2,1) = ts
                  endif
               endif
            endif
         endif
      endif
!
! trajectory diagnostic: updates partd, fvtp, fvmtp
      if (ntt > 0) then
         it = ntime/ntt
         if (ntime==ntt*it) then
            call traj_diag13(partd,fvtp,fvmtp)
            if (nst==3) then
! display test particle velocity distributions
               if (ndt==1) then
                  call displayfv1(fvtp,fvmtp,' ELECTRON',ntime,nmv,2,irc)
               else if (ndt==2) then
                  call displayfv1(fvtp,fvmtp,' ION',ntime,nmv,2,irc)
               endif
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! phase space diagnostic
      if (nts > 0) then
         it = ntime/nts
         if (ntime==nts*it) then
! calculate electron phase space distribution: updates fvs
            call ephasesp_diag1(fvs)
! plot electrons
            if ((nds==1).or.(nds==3)) then
! vx, vy, or vz versus x
               nn = nsxv; ierr = 0
               do i = 1, 3
                  if (mod(nn,2)==1) then
                     call dpmbgrasp1(ppart,kpic,' ELECTRON',ntime,999,  &
     &omx,omy,omz,nx,i+1,1,ntsc,irc)
                     if (irc==1) then
                        ierr = 1
                        exit
                     endif
                     irc = 0
                  endif
                  nn = nn/2
               enddo
               if (ierr==1) exit
! vx-vy, vx-vz or vy-vz
               nn = nsvv; ierr = 0
               do i = 1, 3
                  if (mod(nn,2)==1) then
                     call dpmbgrasp1(ppart,kpic,' ELECTRON',ntime,999,  &
     &omx,omy,omz,nx,min(i+2,4),max(i,2),ntsc,irc)
                     if (irc==1) then
                        ierr = 1
                        exit
                     endif
                     irc = 0
                  endif
                  nn = nn/2
               enddo
               if (ierr==1) exit
            endif
! ion phase space
            if (movion==1) then
! calculate ion phase space distribution: updates fvsi
               call iphasesp_diag1(fvsi)
! plot ions
               if ((nds==2).or.(nds==3)) then
! vx, vy, or vz versus x
                  nn = nsxv; ierr = 0
                  do i = 1, 3
                     if (mod(nn,2)==1) then
                        call dpmbgrasp1(pparti,kipic,' ION',ntime,999,  &
     &omx,omy,omz,nx,i+1,1,ntsc,irc)
                        if (irc==1) then
                           ierr = 1
                           exit
                        endif
                        irc = 0
                     endif
                     nn = nn/2
                  enddo
                  if (ierr==1) exit
! vx-vy, vx-vz or vy-vz
                  nn = nsvv; ierr = 0
                  do i = 1, 3
                     if (mod(nn,2)==1) then
                        call dpmbgrasp1(pparti,kipic,' ION',ntime,999,  &
     &omx,omy,omz,nx,min(i+2,4),max(i,2),ntsc,irc)
                        if (irc==1) then
                           ierr = 1
                           exit
                        endif
                        irc = 0
                     endif
                     nn = nn/2
                  enddo
                  if (ierr==1) exit
               endif
            endif
         endif
      endif
!
! push electrons with OpenMP: updates ppart, wke, kpic
      call push_electrons13(ppart,kpic)
!
! push ions with OpenMP: updates pparti, wki, kipic
      if (movion==1) then
         call push_ions13(pparti,kipic)
      endif
!
! start running simulation backwards:
! need to advance maxwell field solver one step ahead
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
            call em_time_reverse1()
         endif
      endif
!
! energy diagnostic: updates wt
      if (ntw > 0) then
         it = ntime/ntw
         if (ntime==ntw*it) then
            call energy_diag13(wt,ntime,iuot)
         endif
      endif
!
! restart file
      if (ntr > 0) then
         it = n/ntr
         if (n==ntr*it) then
            call dtimer(dtime,itime,-1)
            call bwrite_restart13(iur,n)
            call dwrite_restart13(iur)
            call dtimer(dtime,itime,1)
            tfield = tfield + real(dtime)
         endif
      endif
!
      enddo
      ntime = ntime + 1
!
! loop time
      call dtimer(dtime,ltime,1)
      tloop = tloop + real(dtime)
!
! * * * end main iteration loop * * *
!
!
      write (iuot,*)
      write (iuot,*) 'ntime, relativity = ', ntime, relativity
      if (treverse==1) write (iuot,*) 'treverse = ', treverse
!
! print timing summaries
      call print_timings13(tinit,tloop,iuot)
!
      if ((ntw > 0).or.(ntt > 0)) then
         if (nplot > 0) call reset_graphs
      endif
!
! trajectory diagnostic
      if (ntt > 0) then
         if ((nst==1).or.(nst==2)) then
            if (nplot > 0) irc = open_graphs(1)
            call displaytr1(partd,t0,dt*real(ntt),itt,2,999,irc)
            if (irc==1) stop
            call reset_nplot(nplot,irc)
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         call displayw1(wt,t0,dt*real(ntw),itw,irc)
         if (irc==1) stop
! print energy summaries
         call print_energy13(wt,iuot)
      endif
!
! velocity diagnostic
      if (ntv > 0) then
! display electron distribution time histories and entropy
         if ((ndv==1).or.(ndv==3)) then
            call displayfvt1(fvtm,' ELECT',t0,dt*real(ntv),itv,irc)
            if (irc==1) stop
         endif
! display ion distribution time histories and entropy
         if (movion==1) then
            if ((ndv==2).or.(ndv==3)) then
               call displayfvt1(fvtmi,' ION',t0,dt*real(ntv),itv,irc)
               if (irc==1) stop
            endif
         endif
      endif
!
! display final spectral analysis for ion density
      if (movion==1) then
         if (ntdi > 0) then
            if ((nddi==2).or.(nddi==3)) then
! display frequency spectrum
               call dmscaler1(wkdi,'ION DENSITY OMEGA VS MODE',ntime,999&
     &,1,modesxdi,cwk,irc)
               if (irc==1) stop
            endif
         endif
      endif
!
! display final spectral analysis for potential
      if (ntp > 0) then
         if ((ndp==2).or.(ndp==3)) then
! display frequency spectrum
            call dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,    &
     &modesxp,cwk,irc)
            if (irc==1) stop
         endif
      endif
!
! display final spectral analysis for ion current density
      if (movion==1) then
         if (ntji > 0) then
            if ((ndji==2).or.(ndji==3)) then
! display frequency spectrum
               call dmvector1(vwkji,'ION CURRENT OMEGA VS MODE',ntime,  &
     &999,2,2,modesxji,cwk,irc)
               if (irc==1) stop
            endif
         endif
      endif
!
! display final spectral analysis for vector potential
      if (nta > 0) then
         if ((nda==2).or.(nda==3)) then
! display frequency spectrum
            call dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime,  &
     &999,2,2,modesxa,cwk,irc)
            if (irc==1) stop
         endif
      endif
!
! display final spectral analysis for transverse efield
      if (ntet > 0) then
         if ((ndet==2).or.(ndet==3)) then
! display frequency spectrum
            call dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',ntime&
     &,999,2,2,modesxet,cwk,irc)
            if (irc==1) stop
         endif
      endif
!
! display final spectral analysis for radiative vector potential
      if (ntar > 0) then
         if ((ndar==2).or.(ndar==3)) then
! display frequency spectrum
            call dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',   &
     &ntime,999,2,2,modesxar,cwk,irc)
            if (irc==1) stop
         endif
      endif
!
! close diagnostics
      call close_diags13(iudm)
! close reset and restart files
      call close_restart1()
! close output file
      write (iuot,*) ' * * * q.e.d. * * *'
      close(unit=iuot)
! close graphics device
      call close_graphs
!
      stop
      end program
