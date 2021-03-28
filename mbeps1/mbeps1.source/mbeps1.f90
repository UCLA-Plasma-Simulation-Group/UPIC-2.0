!-----------------------------------------------------------------------
! 1D Electrostatic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
      program mbeps1
      use f1
      use graf1
      use omplib
      implicit none
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
! read namelist
      call readnml1(iuin)
! override input data
      idcode = 1
      ndim = 1
      if (nts > 0) then
         nsxv = min(nsxv,1); nsvv = 0
      endif
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
! increase number of coordinates for particle tag
      if ((ntt > 0).or.((nts > 0).and.(ntsc > 0))) then
         idimp = idimp + 1
      endif
!
! initialize scalars for standard code
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
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
      endif
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
      call init_fields1()
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
         call init_electrons1()
!
! initialize background charge density: updates qi
         if (movion==0) then
            qi = 0.0
            if (ndprof > 0) then
               call mpost1(ppart,qi,kpic,-qme,tdpost,mx)
               call maguard1(qi,tguard,nx)
            endif
         endif
!
! initialize ions: updates pparti, kipic
! pparti = tiled on particle arrays
! kipic = number of ions in each tile
! cui = ion current density with guard cells
         if (movion==1) then
            call init_ions1()
         endif
!
! restart to continue a run which was interrupted
      else if (nustrt==2) then
         call bread_restart1(iur)
         nstart = ntime + 1
! start a new run with data from a previous run
      else if (nustrt==0) then
         call bread_restart1(iur0)
      endif
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
      call initialize_diagnostics1(ntime)
!
! read in restart diagnostic file to continue interrupted run
      if (nustrt==2) call dread_restart1(iur)
!
! write reset file
!     call bwrite_restart1(iurr,ntime)
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      write (iuot,*) 'program mbeps1'
!
! debug reset
!  10 if (irc==6) then
!        irc = 0
!        call bread_restart1(iurr)
!        call reset_diags1()
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
! display frequency spectrum
               if ((nddi==2).or.(nddi==3)) then
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
! transform charge to fourier space: updates qe
      isign = -1
      call mfft1r(qe,isign,mixup,sct,tfft,indx)
!
! calculate force/charge in fourier space: updates fxe, we
      isign = -1
      call mpois1(qe,fxe,ffc,we,tfield,nx)
!
! transform force to real space: updates fxe
      isign = 1
      call mfft1r(fxe,isign,mixup,sct,tfft,indx)
!
! add external traveling wave field: updates fxe
      ts = dt*real(ntime)
      call meaddext1(fxe,tfield,amodex,freq,ts,trmp,toff,el0,er0,nx)
!
! copy guard cells: updates fxe
      call mdguard1(fxe,tguard,nx)
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
! display frequency spectrum
            if ((ndp==2).or.(ndp==3)) then
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
! fluid moments diagnostic
      if (ntfm > 0) then
         it = ntime/ntfm
         if (ntime==ntfm*it) then
! updates fmse
            call efluidms_diag1(fmse)
            if (movion==1) then
! updates fmsi
               call ifluidms_diag1(fmsi)
            endif
         endif
      endif
!
! velocity diagnostic
      if (ntv > 0) then
         it = ntime/ntv
         if (ntime==ntv*it) then
! updates fv, fe, fvm, fvtm, wkt
            call evelocity_diag1(fv,fe,fvm,fvtm,wkt)
! display electron velocity distributions
            if ((ndv==1).or.(ndv==3)) then
               if ((nvft==1).or.(nvft==3)) then
                  call displayfv1(fv,fvm,' ELECTRON',ntime,nmv,1,irc)
                  if (irc==1) exit; irc = 0
               endif
! display electron energy distribution
               if ((nvft==2).or.(nvft==3)) then
                  call displayfe1(fe,wkt,' ELECTRON',ntime,nmv,irc)
                  if (irc==1) exit; irc = 0
               endif
            endif
            if (movion==1) then
! updates fvi, fei, fvmi, fvtmi, wkt
               call ivelocity_diag1(fvi,fei,fvmi,fvtmi,wkt)
! display ion velocity distributions
               if ((ndv==2).or.(ndv==3)) then
                  if ((nvft==1).or.(nvft==3)) then
                     call displayfv1(fvi,fvmi,' ION',ntime,nmv,1,irc)
                     if (irc==1) exit; irc = 0
                  endif
! display ion energy distribution
                  if ((nvft==2).or.(nvft==3)) then
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
            call traj_diag1(partd,fvtp,fvmtp)
            if (nst==3) then
! display test particle velocity distributions
               if (ndt==1) then
                  call displayfv1(fvtp,fvmtp,' ELECTRON',ntime,nmv,1,irc)
               else if (ndt==2) then
                  call displayfv1(fvtp,fvmtp,' ION',ntime,nmv,1,irc)
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
! plot electrons vx versus x
            if ((nds==1).or.(nds==3)) then
               call dpmgrasp1(ppart,kpic,' ELECTRON',ntime,999,nx,2,1,  &
     &ntsc,irc)
               if (irc==1) exit; irc = 0
            endif
! ion phase space
            if (movion==1) then
! calculate ion phase space distribution: updates fvsi
               call iphasesp_diag1(fvsi)
! plot ions vx versus x
               if ((nds==2).or.(nds==3)) then
                  call dpmgrasp1(pparti,kipic,' ION',ntime,999,nx,2,1,  &
     &ntsc,irc)
                  if (irc==1) exit; irc = 0
               endif
            endif
         endif
      endif
!
! push electrons with OpenMP: updates ppart, wke, kpic
      call push_electrons1(ppart,kpic)
!
! push ions with OpenMP: updates pparti, wki, kipic
      if (movion==1) then
         call push_ions1(pparti,kipic)
      endif
!
! start running simulation backwards: updates ppart, pparti
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
            call es_time_reverse1()
         endif
      endif
!
! energy diagnostic: updates wt
      if (ntw > 0) then
         it = ntime/ntw
         if (ntime==ntw*it) then
            call energy_diag1(wt,ntime,iuot)
         endif
      endif
!
! restart file
      if (ntr > 0) then
         it = n/ntr
         if (n==ntr*it) then
            call dtimer(dtime,itime,-1)
            call bwrite_restart1(iur,n)
            call dwrite_restart1(iur)
            call writnml1(iudm)
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
      write (iuot,*)
      write (iuot,*) 'ntime, relativity = ', ntime, relativity
      if (treverse==1) write (iuot,*) 'treverse = ', treverse
!
! print timing summaries
      call print_timings1(tinit,tloop,iuot)
!
      if ((ntw>0).or.(ntt>0).or.(ntv>0).or.(ntdi>0).or.(ntp>0)) then
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
! display energy histories
         call displayw1(wt,t0,dt*real(ntw),itw,irc)
         if (irc==1) stop
! print energy summaries
         call print_energy1(wt,iuot)
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
! close diagnostics
      call close_diags1(iudm)
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
