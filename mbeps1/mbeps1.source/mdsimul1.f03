!-----------------------------------------------------------------------
! High Level library for 1-2/2D Darwin OpenMP PIC code
!
! subroutines defined:
!
! init_dfields13: allocate darwin field data for standard code
! del_dfields13: delete darwin field data for standard code
!
! calc_shift13: calculate shift constant for iteration
!
! dpush_electrons13: push darwin electrons with OpenMP
!
! dpush_ions13: push darwin ions with OpenMP
!
! d_time_reverse1: start running simulation backwards
!
! denergy_diag13: darwin energy diagnostic
!
! init_dspectrum13: allocate scratch arrays for darwin vector fields
! del_dspectrum13: delete scratch arrays for darwin vector fields
!
! init_edcurrent_diag13: initialize darwin electron current density
!                        diagnostic
! edcurrent_diag13: darwin electron current density diagnostic
! del_edcurrent_diag13: delete darwin electron current density
!                       diagnostic
!
! init_vdpotential_diag13: initialize darwin vector potential diagnostic
! vdpotential_diag13: darwin vector potential diagnostic
! del_vdpotential_diag13: delete darwin vector potential diagnostic
!
! init_detfield_diag13: initialize darwin transverse efield diagnostic
! detfield_diag13: darwin transverse efield diagnostic
! del_detfield_diag13: delete darwin transverse efield diagnostic
!
! init_dbfield_diag13: initialize darwin magnetic field diagnostic
! dbfield_diag13: darwin magnetic field diagnostic
! del_dbfield_diag13: delete darwin magnetic field diagnostic
!
! edfluidms_diag13: darwin electron fluid moments diagnostic
!
! idfluidms_diag13: darwin ion fluid moments diagnostic
!
! print_dtimings13: print darwin timing summaries
!
! reset_ddiags13: reset electrostatic/darwin diagnostics
! close_ddiags13: close darwin diagnostics
!
! initialize_ddiagnostics13: initialize all diagnostics from namelist
!                            input parameters
!
! darwin_predictor13: predictor for darwin iteration
! darwin_iteration13: darwin iteration loop
!
! bwrite_drestart13: write out basic restart file for darwin code
! bread_drestart13: read in basic restart file for darwin code
! dwrite_drestart13: write out restart diagnostic file for darwin code
! dread_drestart13: read in restart diagnostic file for darwin code
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
! update: january 13, 2021
      module fd1
      use f1
      use fb1
      use mdpush1
      implicit none
! declare scalars for standard code
      integer :: k
      real :: q2m0, wpm, wpmax, wpmin
!
! declare and initialize timing data
      real :: tdcjpost = 0.0
!
! declare arrays for standard code:
! dcu/dcui = electron/ion acceleration density with guard cells
      real, dimension(:,:), allocatable :: dcu, dcui
! amu/amui = electron/ion momentum flux with guard cells
      real, dimension(:,:), allocatable :: amu, amui
! cus = transverse electric field
      real, dimension(:,:), allocatable :: cus
! ffe = form factor arrays for poisson solvers
      complex, dimension(:), allocatable :: ffe
!
! private diagnostic arrays
! scratch array for vector field
      complex, dimension(:,:), allocatable :: vfieldc
! oldcue = previous current density
      real, dimension(:,:), allocatable :: oldcue
! vpott = store selected fourier modes for vector potential
      complex, dimension(:,:), allocatable :: vpott
! vpks = accumulated complex spectrum for vector potential
      double precision, dimension(:,:,:,:), allocatable :: vpks
! ett = store selected fourier modes for transverse efield
      complex, dimension(:,:), allocatable :: ett
! vpkset = accumulated complex spectrum for transverse efield
      double precision, dimension(:,:,:,:), allocatable :: vpkset
! bt = store selected fourier modes for magnetic field
      complex, dimension(:,:), allocatable :: bt
!
      save
!
      public :: dcu, dcui, amu, amui, cus, ffe
      private :: vfieldc, vpott, ett, bt, vpks, vpkset
!
      contains
!
!-----------------------------------------------------------------------
      subroutine init_dfields13()
! allocate darwin field data for standard code
      implicit none
! allocate electromagnetic field data: fxyze, cue, byze, eyz, byz
! allocate electrostatic field data: qe, qi, fxe, ffc, mixup, sct
      call init_fields13()
      allocate(dcu(2,nxe),cus(2,nxe),amu(2,nxe))
      allocate(ffe(nxh))
      if (movion==1) then
         allocate(dcui(2,nxe),amui(2,nxe))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_dfields13()
! delete darwin field data for standard code
      implicit none
      if (movion==1) then
         deallocate(dcui,amui)
      endif
      deallocate(dcu,cus,amu,ffe)
      call del_fields13()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine calc_shift13(iuot)
! calculate shift constant for iteration: updates q2m0
      implicit none
! iuot = output file descriptor
      integer, intent(in) :: iuot
! find maximum and minimum initial electron density
      qe = 0.0
      call mpost1(ppart,qe,kpic,qme,tdpost,mx)
      call maguard1(qe,tguard,nx)
      if (movion==1) then
         qi = 0.0
         call mpost1(pparti,qi,kipic,qmi,tdpost,mx)
         call maguard1(qi,tguard,nx)
         call mfwptminx1(qe,qi,qbme,qbmi,wpmax,wpmin,nx)
      else
         call mfwpminx1(qe,qbme,wpmax,wpmin,nx)
      endif
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (iuot,*) 'wpm = ', wpm
      q2m0 = wpm/affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dpush_electrons13(ppart,kpic)
! push darwin electrons with OpenMP
      implicit none
! ppart = tiled electron particle array
      real, dimension(:,:,:), intent(inout) :: ppart
! kpic = number of electrons in each tile
      integer, dimension(:), intent(inout) :: kpic
      wke = 0.0
! updates ppart and wke, and possibly ncl, ihole, irc
! Boris pusher
      if (mzf==0) then
         call wmbpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt,ci&
     &,wke,tpush,nx,mx,ipbc,relativity,plist,irc)
! Analytic Boris pusher
      else if (mzf==2) then
         call wmabpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt, &
     &ci,wke,tpush,nx,mx,ipbc,relativity,plist,irc)
! Exact Analytic pusher
      else if (mzf==3) then
         call wmeabpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt,&
     &ci,wke,tpush,nx,mx,ipbc,relativity,plist,irc)
! zero force: updates ppart, wke and possibly ncl, ihole, and irc
      else if (mzf==1) then
         call wmpush1zf(ppart,kpic,ncl,ihole,dt,ci,wke,tpush,nx,mx,ipbc,&
     &relativity,plist,irc)
      endif
!
! reorder electrons by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      if (irc==0) then
         call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,plist,  &
     &irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
! sanity check for electrons
      if (irc2(1)==0) then
         if (monitor > 0) then
            call mcheck1(ppart,kpic,nx,mx,irc)
            if (irc /= 0) stop
         endif
! recover from wmporder1 errors: updates ppart
      else if (irc2(1) /= 0) then
         call reorder_electrons1(irc2)
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dpush_ions13(pparti,kipic)
! push darwin ions with OpenMP
      implicit none
! pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), intent(inout) :: pparti
! kipic = number of electrons/ions in each tile
      integer, dimension(:), intent(inout) :: kipic
      wki = 0.0
! updates pparti and wki, and possibly ncl, ihole, irc
! Boris pusher
      if (mzf==0) then
         call wmbpush1(pparti,fxyze,byze,kipic,ncl,ihole,omx,qbmi,dt,dt,&
     &ci,wki,tpush,nx,mx,ipbc,relativity,plist,irc)
! Analytic Boris pusher
      else if (mzf==2) then
         call wmabpush1(pparti,fxyze,byze,kipic,ncl,ihole,omx,qbmi,dt,dt&
     &,ci,wki,tpush,nx,mx,ipbc,relativity,plist,irc)
! Exact Analytic pusher
      else if (mzf==3) then
         call wmeabpush1(pparti,fxyze,byze,kipic,ncl,ihole,omx,qbmi,dt, &
     &dt,ci,wki,tpush,nx,mx,ipbc,relativity,plist,irc)
! zero force: updates pparti, wki and possibly ncl, ihole, and irc
      else if (mzf==1) then
         call wmpush1zf(pparti,kipic,ncl,ihole,dt,ci,wki,tpush,nx,mx,   &
     &ipbc,relativity,plist,irc)
      endif
      wki = wki*rmass
!
! reorder ions by tile with OpenMP:
! updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
      if (irc==0) then
         call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,plist,&
     &irc2)
      else
         irc2(1) = 1; irc2(2) = irc; irc = 0
      endif
!
! sanity check for ions
      if (irc2(1)==0) then
         if (monitor > 0) then
            call mcheck1(pparti,kipic,nx,mx,irc)
            if (irc /= 0) stop
         endif
! recover from wmporder1 errors: updates pparti
      else if (irc2(1) /= 0) then
         call reorder_ions1(irc2)
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine d_time_reverse1()
! start running simulation backwards
! need to reverse time lag in leap-frog integration scheme
      implicit none
! still need to extrapolate cus for next iteration?
      dt = -dt
      ws = 0.0
      call wmpush1zf(ppart,kpic,ncl,ihole,dt,ci,ws,tpush,nx,mx,ipbc,    &
     &relativity,plist,irc)
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,plist,irc2)
      if (movion==1) then
         call wmpush1zf(pparti,kipic,ncl,ihole,dt,ci,ws,tpush,nx,mx,ipbc&
     &,relativity,plist,irc)
         call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,plist,&
     &irc2)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine denergy_diag13(wt,ntime,iuot)
! darwin energy diagnostic
      implicit none
! wt = energy time history array
! ntime = current time step
! iuot = output file descriptor
      real, dimension(:,:), intent(inout) :: wt
      integer, intent(in) :: ntime, iuot
      wef = we + wb
      ws = wef + wke + wki
      if (ntime==0) s(6) = ws
      if (ndw > 0) then
         write (iuot,*) 'Total Field, Kinetic and Total Energies:'
         if (movion==0) then
            write (iuot,'(3e14.7)') wef, wke, ws
         else
            write (iuot,'(4e14.7)') wef, wke, wki, ws
         endif
         write (iuot,*) 'Electric(l,t) and Magnetic Energies = '
         write (iuot,'(3e14.7)') we, wf, wb
      endif
      itw = itw + 1
! store energies in time history array
      wt(itw,:) = (/wef,wke,wki,ws,we,wf,wb/)
      s(1) = s(1) + we
      s(2) = s(2) + wke
      s(3) = s(3) + wf
      s(4) = s(4) + wb
      s(5) = s(5) + wki
      s(6) = min(s(6),dble(ws))
      s(7) = max(s(7),dble(ws))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_dspectrum13()
! allocate scratch arrays for darwin vector fields
      implicit none
      call init_spectrum13()
      allocate(vfieldc(2,nxh))
! allocate and initialize frequency array for spectral analysis
      if ((nta > 0).or.(ntet > 0)) then
         if (.not.allocated(wm)) then
            iw = (wmax - wmin)/dw + 1.5
            allocate(wm(iw))
            do it = 1, iw
               wm(it) = wmin + dw*real(it-1)
            enddo
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_dspectrum13()
! delete scratch arrays for darwin vector fields
      implicit none
      deallocate(vfieldc)
      call del_spectrum13()
      if ((nta > 0).or.(ntet > 0)) then
         if (allocated(wm)) then
            deallocate(wm)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_edcurrent_diag13()
! initialize darwin electron current density diagnostic
      implicit none
      fjename = 'curek1.'//cdrun
      modesxje = min(modesxje,nxh+1)
      allocate(oldcue(2,nxe))
      allocate(curet(2,modesxje))
! open file: updates njerec and possibly iuje
      if (njerec==0) then
         call dafopenvc1(curet,iuje,njerec,trim(fjename))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine edcurrent_diag13(vfield)
! darwin electron current density diagnostic
      implicit none
! vfield = scratch array for vector field
      real, dimension(:,:), intent(inout) :: vfield
      vfield = oldcue
! transform electron current density to fourier space: updates vfield
      isign = -1
      call mfft1rn(vfield,isign,mixup,sct,tfft,indx)
! calculate smoothed electron current in fourier space: updates vfieldc
      call msmooth13(vfield,vfieldc,ffc,tfield,nx)
! store selected fourier modes: updates curet
      call mrdvmodes1(vfieldc,curet,tfield,nx,modesxje)
! write diagnostic output: updates nderec
      call dafwritevc1(curet,tdiag,iuje,njerec,modesxje)
! transform smoothed electron current to real space: updates vfield
      call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
      call mcguard1(vfield,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_edcurrent_diag13()
! delete darwin electron current density diagnostic
      implicit none
      if (njerec > 0) then
         close(unit=iuje)
         njerec = njerec - 1
      endif
      deallocate(curet,oldcue)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_vdpotential_diag13()
! initialize darwin vector potential diagnostic
      implicit none
      faname = 'vpotk1.'//cdrun
      modesxa = min(modesxa,nxh+1)
      allocate(vpott(2,modesxa))
! open file: updates narec and possibly iua
      if (narec==0) call dafopenvc1(vpott,iua,narec,trim(faname))
! spectral analysis
      if ((nda==2).or.(nda==3)) then
         mta = (nloop - 1)/nta + 1; ita = 0
         allocate(vpkw(2,modesxa,iw,2),vpks(2,4,modesxa,iw))
         allocate(vwk(2,modesxa,2))
         vpks = 0.0d0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine vdpotential_diag13(vfield,vpkw,vwk,ntime)
! darwin vector potential diagnostic
      implicit none
! vfield = scratch array for vector field
      real, dimension(:,:), intent(inout) :: vfield
! vpkw = power spectrum for vector potential
      real, dimension(:,:,:,:), intent(inout) :: vpkw
! vwk = maximum frequency as a function of k for vector potential
      real, dimension(:,:,:), intent(inout) :: vwk
! ntime = current time step
      integer, intent(in) :: ntime
! calculate vector potential in fourier space: updates vfieldc
      call mapot1(cue,vfieldc,ffc,ci,ws,tfield,nx)
! store selected fourier modes: updates vpott
      call mrdvmodes1(vfieldc,vpott,tfield,nx,modesxa)
! write diagnostic output: updates narec
      call dafwritevc1(vpott,tdiag,iua,narec,modesxa)
! transform vector potential to real space: updates vfield
      if ((nda==1).or.(nda==3)) then
         call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
         call mcguard1(vfield,tguard,nx)
      endif
! spectral analysis
      if ((nda==2).or.(nda==3)) then
         ita = ita + 1
         ts = dt*real(ntime)
! performs frequency analysis of accumulated complex vector time series
         call mivcspect1(vpott,wm,vpkw,vpks,ts,t0,tdiag,mta,iw,modesxa, &
     &nx,1)
! find frequency with maximum power for each mode
         vwk(1,:,1) = wm(maxloc(vpkw(1,:,:,1),dim=2))
         vwk(2,:,1) = wm(maxloc(vpkw(2,:,:,1),dim=2))
         vwk(1,:,2) = wm(maxloc(vpkw(1,:,:,2),dim=2))
         vwk(2,:,2) = wm(maxloc(vpkw(2,:,:,2),dim=2))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_vdpotential_diag13()
! delete darwin vector potential diagnostic
      implicit none
      if (narec > 0) then
         close(unit=iua)
         narec = narec - 1
      endif
      deallocate(vpott)
! spectral analysis
      if ((nda==2).or.(nda==3)) then
         deallocate(vpkw,vwk,vpks)
      endif
      ceng = affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_detfield_diag13()
! initialize darwin transverse efield diagnostic
      implicit none
      fetname = 'etk1.'//cdrun
      modesxet = min(modesxet,nxh+1)
      allocate(ett(2,modesxet))
! open file: updates netrec and possibly iuet
      if (netrec==0) call dafopenvc1(ett,iuet,netrec,trim(fetname))
! spectral analysis
      if ((ndet==2).or.(ndet==3)) then
         mtet = (nloop - 1)/ntet + 1; itet = 0
         allocate(vpkwet(2,modesxet,iw,2),vpkset(2,4,modesxet,iw))
         allocate(vwket(2,modesxet,2))
         vpkset = 0.0d0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine detfield_diag13(vfield,vpkwet,vwket,ntime)
! darwin transverse efield diagnostic
      implicit none
! vfield = scratch array for vector field
      real, dimension(:,:), intent(inout) :: vfield
! vpkwet = power spectrum for transverse efield
      real, dimension(:,:,:,:), intent(inout) :: vpkwet
! vwket = maximum frequency as a function of k for transverse efield
      real, dimension(:,:,:), intent(inout) :: vwket
! ntime = current time step
      integer, intent(in) :: ntime
! calculate transverse efield in fourier space: updates vfieldc
      call metfield1(dcu,vfieldc,ffe,ci,ws,tfield,nx)
! store selected fourier modes: updates ett
      call mrdvmodes1(vfieldc,ett,tfield,nx,modesxet)
! write diagnostic output: updates netrec
      call dafwritevc1(ett,tdiag,iuet,netrec,modesxet)
! transform transverse efield to real space: updates vfield
      if ((ndet==1).or.(ndet==3)) then
         call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
         call mcguard1(vfield,tguard,nx)
      endif
! spectral analysis
      if ((ndet==2).or.(ndet==3)) then
         itet = itet + 1
         ts = dt*real(ntime)
! performs frequency analysis of accumulated complex vector time series
         call mivcspect1(ett,wm,vpkwet,vpkset,ts,t0,tdiag,mtet,iw,      &
     &modesxet,nx,0)
! find frequency with maximum power for each mode
         vwket(1,:,1) = wm(maxloc(vpkwet(1,:,:,1),dim=2))
         vwket(2,:,1) = wm(maxloc(vpkwet(2,:,:,1),dim=2))
         vwket(1,:,2) = wm(maxloc(vpkwet(1,:,:,2),dim=2))
         vwket(2,:,2) = wm(maxloc(vpkwet(2,:,:,2),dim=2))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_detfield_diag13()
! delete darwin transverse efield diagnostic
      implicit none
      if (netrec > 0) then
         close(unit=iuet)
         netrec = netrec - 1
      endif
      deallocate(ett)
! spectral analysis
      if ((ndet==2).or.(ndet==3)) then
         deallocate(vpkwet,vwket,vpkset)
      endif
      ceng = affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_dbfield_diag13()
! initialize darwin magnetic field diagnostic
      implicit none
      fbname = 'bk1.'//cdrun
      modesxb = min(modesxb,nxh+1)
      allocate(bt(2,modesxb))
! open file: updates nbrec and possibly iub
      if (nbrec==0) call dafopenvc1(bt,iub,nbrec,trim(fbname))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dbfield_diag13(vfield)
! darwin magnetic field diagnostic
      implicit none
! vfield = scratch array for vector field
      real, dimension(:,:), intent(inout) :: vfield
! calculate magnetic field in fourier space: updates vfieldc
      call mibpois1(cue,vfieldc,ffc,ci,ws,tfield,nx)
! store selected fourier modes: updates bt
      call mrdvmodes1(vfieldc,bt,tfield,nx,modesxb)
! write diagnostic output: updates nbrec
      call dafwritevc1(bt,tdiag,iub,nbrec,modesxb)
! transform magnetic field to real space: updates vfield
      call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
      call mcguard1(vfield,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_dbfield_diag13()
! delete darwin magnetic field diagnostic
      implicit none
      if (nbrec > 0) then
         close(unit=iub)
         nbrec = nbrec - 1
      endif
      deallocate(bt)
      ceng = affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine edfluidms_diag13(fmse)
! darwin electron fluid moments diagnostic
      implicit none
! fmse = electron fluid moments
      real, dimension(:,:), intent(inout) :: fmse
! calculate electron fluid moments
      if ((ndfm==1).or.(ndfm==3)) then
         call dtimer(dtime,itime,-1)
         fmse = 0.0
         call dtimer(dtime,itime,1)
         tdiag = tdiag + real(dtime)
         call wmgbprofx1(ppart,fxyze,byze,fmse,kpic,omx,qbme,dt,ci,tdiag&
     &,npro,nx,mx,relativity)
! add guard cells with OpenMP: updates fmse
         call mamcguard1(fmse,tdiag,nx)
! calculates fluid quantities from fluid moments: updates fmse
         call mfluidqs13(fmse,tdiag,npro,nx)
! write real space diagnostic output: updates nferec
         call dafwritev1(fmse,tdiag,iufe,nferec,nx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine idfluidms_diag13(fmsi)
! darwin ion fluid moments diagnostic
      implicit none
! fmsi = ion fluid moments
      real, dimension(:,:), intent(inout) :: fmsi
! calculate ion fluid moments
      if ((ndfm==2).or.(ndfm==3)) then
         call dtimer(dtime,itime,-1)
         fmsi = 0.0
         call dtimer(dtime,itime,1)
         tdiag = tdiag + real(dtime)
         call wmgbprofx1(pparti,fxyze,byze,fmsi,kipic,omx,qbmi,dt,ci,   &
     &tdiag,npro,nx,mx,relativity)
! add guard cells with OpenMP: updates fmsi
         call mamcguard1(fmsi,tdiag,nx)
! calculates fluid quantities from fluid moments: updates fmsi
         call mfluidqs13(fmsi,tdiag,npro,nx)
         fmsi = rmass*fmsi
! write real space diagnostic output: updates nfirec
         call dafwritev1(fmsi,tdiag,iufi,nfirec,nx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine print_dtimings13(tinit,tloop,iuot)
! print darwin timing summaries
      implicit none
! iuot = output file descriptor
      integer, intent(in) :: iuot
      real, intent(in) :: tinit
      real, intent(inout) :: tloop
! local data
      real :: time
      write (iuot,*)
      write (iuot,*) 'initialization time = ', tinit
      write (iuot,*) 'deposit time = ', tdpost
      write (iuot,*) 'current deposit time = ', tdjpost
      write (iuot,*) 'current derivative deposit time = ', tdcjpost
      tdpost = tdpost + tdjpost + tdcjpost
      write (iuot,*) 'total deposit time = ', tdpost
      write (iuot,*) 'guard time = ', tguard
      write (iuot,*) 'solver time = ', tfield
      write (iuot,*) 'fft time = ', tfft
      write (iuot,*) 'push time = ', tpush
      write (iuot,*) 'sort time = ', tsort
      tfield = tfield + tguard + tfft
      write (iuot,*) 'total solver time = ', tfield
      time = tdpost + tpush + tsort
      write (iuot,*) 'total particle time = ', time
      write (iuot,*) 'total diagnostic time = ', tdiag
      ws = time + tfield + tdiag
      tloop = tloop - ws
      write (iuot,*) 'total and additional time = ', ws, tloop
      write (iuot,*)
! summarize particle timings
      ws = 1.0e+09/(real(nloop)*real(np+npi))
      write (iuot,*) 'Push Time (nsec) = ', tpush*ws
      write (iuot,*) 'Deposit Time (nsec) = ', tdpost*ws
      write (iuot,*) 'Sort Time (nsec) = ', tsort*ws
      write (iuot,*) 'Total Particle Time (nsec) = ', time*ws
      write (iuot,*)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine reset_ddiags13()
! reset electrostatic/darwin diagnostics
      implicit none
! reset electrostatic diagnostics
      call reset_diags1()
      if (ntw > 0) s = 0.0d0
! reset darwin diagnositcs
      if (ntje > 0) then
         if (njerec > 1) njerec = 1
      endif
      if (nta > 0) then
         if (narec > 1) narec = 1
         if ((nda==2).or.(nda==3)) then
            ita = 0; vpks = 0.0d0
         endif
      endif
      if (ntet > 0) then
         if (netrec > 1) netrec = 1
         if ((ndet==2).or.(ndet==3)) then
            itet = 0; vpkset = 0.0d0
         endif
      endif
      if (ntb > 0) then
         if (nbrec > 1) nbrec = 1
      endif
      if (movion==1) then
         if (ntji > 0) then
            if (njirec > 1) njirec = 1
            if ((ndji==2).or.(ndji==3)) then
               itji = 0; vpksji = 0.0d0
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine close_ddiags13(iudm)
! close darwin diagnostics
! delete data, close fortran files, and write out diagnostic metafile
      implicit none
! iudm = diagnostic metafile file descriptor
      integer, intent(in) :: iudm
! electron density diagnostic
      if (ntde > 0) call del_edensity_diag1()
! potential diagnostic
      if (ntp > 0) call del_potential_diag1()
! longitudinal efield diagnostic
      if (ntel > 0) call del_elfield_diag1()
! fluid moments diagnostic
      if (ntfm > 0) then
! electrons
         call del_efluidms_diag1()
! ions
         if (movion==1) call del_ifluidms_diag1()
      endif
! darwin electron current diagnostic
      if (ntje > 0) call del_edcurrent_diag13()
! vector potential diagnostic
      if (nta > 0) then
         call del_vdpotential_diag13()
      endif
! transverse efield diagnostic
      if (ntet > 0) call del_detfield_diag13()
! magnetic field diagnostic
      if (ntb > 0) call del_dbfield_diag13()
! velocity diagnostic
      if (ntv > 0) then
         call del_evelocity_diag13()
         if (movion==1) call del_ivelocity_diag13()
      endif
! trajectory diagnostic
      if (ntt > 0) call del_traj_diag13()
! phase space diagnostic
      if (nts > 0) then
         call del_ephasesp_diag1()
         if (movion==1) call del_iphasesp_diag1()
      endif
! ion diagnostics
      if (movion==1) then
! ion density diagnostic
         if (ntdi > 0) call del_idensity_diag1()
! ion current diagnostic
         if (ntji > 0) call del_icurrent_diag13()
      endif
! write final diagnostic metafile
      call writnml1(iudm)
      close(unit=iudm)
! deallocate arrays
      call del_dfields13()
      call del_electrons1()
      if (movion==1) call del_ions1()
      deallocate(part,ppbuff,ncl,ihole)
      if ((ntde > 0).or.(ntp > 0).or.(ntel > 0).or.(ntdi > 0)) then
         call del_spectrum1()
      endif
      if ((nta>0).or.(ntet>0).or.(ntb>0).or.(ntji>0)) then
         call del_dspectrum13()
      endif
      if (ntw > 0) call del_energy_diag13()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine initialize_ddiagnostics13(ntime)
! initialize all diagnostics from namelist input parameters
      implicit none
! ntime = current time step
      integer, intent(in) :: ntime
! initialize energy diagnostic: updates wt
      if (ntw > 0) then
         call init_energy_diag13()
      endif
!
! allocate and initialize scratch arrays for scalar fields:
! allocates sfield
      if ((ntde > 0).or.(ntp > 0).or.(ntel > 0).or.(ntdi > 0)) then
         call init_spectrum1()
      endif
!
! allocate and initialize scratch arrays for vector fields:
! allocates vfield
      if ((ntje>0).or.(nta>0).or.(ntet>0).or.(ntb>0).or.(ntji>0)) then
         call init_dspectrum13()
      endif
!
! initialize electron density diagnostic
      if (ntde > 0) then
         call init_edensity_diag1()
      endif
!
! initialize ion density diagnostic: allocates pkwdi, wkdi
      if (movion==1) then
         if (ntdi > 0) then
            call init_idensity_diag1()
         endif
      endif
!
! initialize potential diagnostic: allocates pkw, wk
      if (ntp > 0) then
         call init_potential_diag1()
      endif
!
! initialize longitudinal efield diagnostic
      if (ntel > 0) then
         call init_elfield_diag1()
      endif
!
! initialize darwin electron current density diagnostic
      if (ntje > 0) then
         call init_edcurrent_diag13()
      endif
!
! initialize ion current density diagnostic: allocates vpkwji, vwkji
      if (movion==1) then
         if (ntji > 0) then
            call init_icurrent_diag13()
         endif
      endif
!
! initialize darwin vector potential diagnostic: allocates vpkw, vwk
      if (nta > 0) then
         call init_vdpotential_diag13()
      endif
!
! initialize darwin transverse efield diagnostic:
! allocates vpkwet, vwket
      if (ntet > 0) then
         call init_detfield_diag13()
      endif
!
! initialize darwin magnetic field diagnostic
      if (ntb > 0) then
         call init_dbfield_diag13()
      endif
!
! initialize fluid moments diagnostic
      if (ntfm > 0) then
! electrons: allocates fmse
         call init_efluidms_diag13()
! ions: allocates fmsi
         if (movion==1) call init_ifluidms_diag13()
      endif
!
! initialize velocity diagnostic
      if (ntv > 0) then
! electrons: allocates fv, fvm, fvtm
         call init_evelocity_diag13()
! ions: allocates fvi, fvmi, fvtmi
         if (movion==1) then
            call init_ivelocity_diag13()
         endif
      endif
!
! initialize trajectory diagnostic: allocates partd, fvtp, fvmtp
      if (ntt > 0) then
         call init_traj_diag13(ntime)
      endif
!
! initialize phase space diagnostic:
      if (nts > 0) then
! electrons: allocates fvs
         call init_ephasesp_diag13()
! ions: allocates fvsi
         if (movion==1) call init_iphasesp_diag13()
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine darwin_predictor13(q2m0)
! predictor for darwin iteration
      implicit none
! q2m0 = shift constant in darwin iteration
      real, intent(in) :: q2m0
! transform current to fourier space: updates cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! calculate magnetic field in fourier space: updates byze, wb
      call mbbpois1(cue,byze,ffc,ci,wb,tfield,nx)
!
! transform magnetic force to real space: updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field: updates byze
      call mbaddext1(byze,tfield,omy,omz,nx)
!
! copy guard cells: updates byze
      call mcguard1(byze,tguard,nx)
!
! add longitudinal and old transverse electric fields: updates fxyze
      call maddvrfield1(fxyze,cus,fxe,tfield)
!
! deposit electron acceleration density and momentum flux with OpenMP:
! updates dcu, amu
      call wmgdjpost1(ppart,fxyze,byze,dcu,amu,kpic,omx,qme,qbme,dt,ci, &
     &tdcjpost,nx,mx,relativity)
!
! deposit ion acceleration density and momentum flux with OpenMP:
! updates dcui, amui
      if (movion==1) then
         call wmgdjpost1(pparti,fxyze,byze,dcui,amui,kipic,omx,qmi,qbmi,&
     &dt,ci,tdcjpost,nx,mx,relativity)
! add electron and ion densities: updates dcu, amu
         call maddcuei1(dcu,dcui,tfield,nxe)
         call maddcuei1(amu,amui,tfield,nxe)
      endif
!
! add old scaled electric field: updates dcu
      call mascfguard1(dcu,cus,q2m0,tdcjpost,nx)
!
! add guard cells: updates dcu, amu
      call macguard1(dcu,tguard,nx)
      call macguard1(amu,tguard,nx)
!
! transform acceleration density and momentum flux to fourier space:
! updates dcu, amu
      isign = -1
      call mfft1rn(dcu,isign,mixup,sct,tfft,indx)
      call mfft1rn(amu,isign,mixup,sct,tfft,indx)
!
! take transverse part of time derivative of current: updates dcu
      call madcuperp1(dcu,amu,tfield,nx)
!
! calculate transverse electric field: updates cus, wf
      call mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)
!
! transform transverse electric field to real space: updates cus
      isign = 1
      call mfft1rn(cus,isign,mixup,sct,tfft,indx)
!
! copy guard cells: updates cus
      call mcguard1(cus,tguard,nx)
!
! add longitudinal and transverse electric fields:
! fxyze = cus + fxe, updates fxyze
! cus needs to be retained for next time step
      call maddvrfield1(fxyze,cus,fxe,tfield)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine darwin_iteration13(q2m0)
! darwin iteration loop
      implicit none
! q2m0 = shift constant in darwin iteration
      real, intent(in) :: q2m0
! deposit electron current and acceleration density and momentum flux
! with OpenMP: updates cue, dcu, amu
      call wmgdcjpost1(ppart,fxyze,byze,cue,dcu,amu,kpic,omx,qme,qbme,dt&
     &,ci,tdcjpost,nx,mx,relativity)
! add guard cells for current, acceleration density, and momentum flux:
! updates cue, dcu, amu
      call macguard1(cue,tguard,nx)
      call macguard1(dcu,tguard,nx)
      call macguard1(amu,tguard,nx)
!
! save electron current for electron current diagnostic later
      if (k==ndc) then
         if (ntje > 0) then
            it = ntime/ntje
            if (ntime==ntje*it) oldcue = cue
         endif
      endif
!
! deposit ion current and acceleration density and momentum flux
! with OpenMP: updates cui, dcui, amui
      if (movion==1) then
         call wmgdcjpost1(pparti,fxyze,byze,cui,dcui,amui,kipic,omx,qmi,&
     &qbmi,dt,ci,tdcjpost,nx,mx,relativity)
! add guard cells for current, acceleration density, and momentum flux:
! updates cui, dcui, amui
         call macguard1(cui,tguard,nx)
         call macguard1(dcui,tguard,nx)
         call macguard1(amui,tguard,nx)
! add electron and ion densities: updates cue, dcu, amu
         call maddcuei1(cue,cui,tfield,nxe)
         call maddcuei1(dcu,dcui,tfield,nxe)
         call maddcuei1(amu,amui,tfield,nxe)
      endif
!
! add scaled electric field: updates dcu
      call mascfguard1(dcu,cus,q2m0,tdcjpost,nx)
!
! transform current to fourier space: update cue
      isign = -1
      call mfft1rn(cue,isign,mixup,sct,tfft,indx)
!
! calculate magnetic field in fourier space: updates byze, wb
      call mbbpois1(cue,byze,ffc,ci,wb,tfield,nx)
!
! transform magnetic force to real space: updates byze
      isign = 1
      call mfft1rn(byze,isign,mixup,sct,tfft,indx)
!
! add constant to magnetic field: updates bzye
      call mbaddext1(byze,tfield,omy,omz,nx)
!
! transform acceleration density and momentum flux to fourier space:
! updates dcu, amu
      isign = -1
      call mfft1rn(dcu,isign,mixup,sct,tfft,indx)
      call mfft1rn(amu,isign,mixup,sct,tfft,indx)
!
! take transverse part of time derivative of current: updates dcu
      call madcuperp1(dcu,amu,tfield,nx)
!
! calculate transverse electric field: updates cus, wf
      call mepois1(dcu,cus,ffe,affp,ci,wf,tfield,nx)
!
! transform transverse electric field to real space: updates cus
      isign = 1
      call mfft1rn(cus,isign,mixup,sct,tfft,indx)
!
! copy guard cells: updates byze, cus
      call mcguard1(byze,tguard,nx)
      call mcguard1(cus,tguard,nx)
!
! add longitudinal and transverse electric fields:
! fxyze = cus + fxe, updates fxyze
! cus needs to be retained for next time step
      call maddvrfield1(fxyze,cus,fxe,tfield)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bwrite_drestart13(iur,ntime)
! write out basic restart file for darwin code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur, ntime
! local data
      integer :: i, j, irc
      integer :: nxv
! write out particles and electrostatic fields
      call bwrite_restart1(iur,ntime)
! write out shift constant for iteration
      write (iur) wpm, q2m0
! write out darwin electric field field
      nxv = size(cus,2)
      write (iur) nxv
      if (nxv > 0) then
         write (iur) ((cus(i,j),i=1,2),j=1,nxv)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_drestart13(iur)
! read in basic restart file for darwin code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j
      integer :: ios, it
! read in particles and electrostatic fields
      call bread_restart1(iur)
! allocate ion current
      if (.not.allocated(cui)) allocate(cui(2,nxe))
      cui = 0.0
! read in shift constant for iteration
      read (iur,iostat=ios) wpm, q2m0
      if (ios /= 0) then
         write (*,*) 'wpm, q2m0 restart error, ios = ', ios
         stop
      endif
! read in darwin electric field field
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'cus size restart error, ios = ', ios
         stop
      endif
      if (it > size(cus,2)) then
         write (*,*) 'cus restart error, size(cus)=',it,size(cus,2)
         stop
      endif
      if (it > 0) then
         read (iur,iostat=ios) ((cus(i,j),i=1,2),j=1,it)
         if (ios /= 0) then
            write (*,*) 'cus read error, ios = ', ios
            stop
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dwrite_drestart13(iur)
! write out restart diagnostic file for darwin code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j, k, l, it, is, ios
      character(len=32) :: fname
! write out restart diagnostic file for electrostatic code
      call dwrite_restart1(iur)
!
! write out electron current density diagnostic parameter
      write (iur) ntje
! write out record location
      if (ntje > 0) then
         write (iur) njerec
! write out record length (zero if error) and file name (if no error)
         if (njerec > 0) then
            inquire(file=fjename,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
            fname = fjename
            write (iur) fname
            endif
         endif
      endif
!
! write out vector potential diagnostic parameter
      write (iur) nta
! write out record location
      if (nta > 0) then
         write (iur) narec
! write out record length (zero if error) and file name (if no error)
         if (narec > 0) then
            inquire(file=faname,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = faname
               write (iur) fname
            endif
         endif
! write out spectrum flag
         if ((nda==2).or.(nda==3)) then
            write (iur) ita
! write out spectrum sizes and data
            if (ita > 0) then
               it = size(vpks,3); is = size(vpks,4)
               write (iur) size(vpks,1), size(vpks,2), it, is
               write (iur) ((((vpks(i,j,k,l),i=1,2),j=1,4),k=1,it),     &
     &l=1,is)
            endif
         else
            it = 0
            write (iur) it
         endif
      endif
!
! write out transverse efield diagnostic parameter
      write (iur) ntet
! write out record location
      if (ntet > 0) then
         write (iur) netrec
! write out record length (zero if error) and file name (if no error)
         if (netrec > 0) then
            inquire(file=fetname,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = fetname
               write (iur) fname
            endif
         endif
! write out spectrum flag
         if ((ndet==2).or.(ndet==3)) then
            write (iur) itet
! write out spectrum sizes and data
            if (itet > 0) then
               it = size(vpkset,3); is = size(vpkset,4)
               write (iur) size(vpkset,1), size(vpkset,2), it, is
               write (iur) ((((vpkset(i,j,k,l),i=1,2),j=1,4),k=1,it),   &
     &l=1,is)
            endif
         else
            it = 0
            write (iur) it
         endif
      endif
!
! write out magnetic field diagnostic diagnostic parameter
      write (iur) ntb
! write out record location
      if (ntb > 0) then
         write (iur) nbrec
! write out record length (zero if error) and file name (if no error)
         if (nbrec > 0) then
            inquire(file=fbname,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = fbname
               write (iur) fname
            endif
         endif
      endif
!
! write out ion current density diagnostic parameter
      if (movion==1) then
         write (iur) ntji
! write out record location
         if (ntji > 0) then
            write (iur) njirec
! write out record length (zero if error) and file name (if no error)
            if (njirec > 0) then
               inquire(file=fjiname,recl=it,iostat=ios)
               if (ios /= 0) it = 0
               write (iur) it
               if (it > 0) then
                  fname = fjiname
                  write (iur) fname
               endif
            endif
! write out spectrum flag
            if ((ndji==2).or.(ndji==3)) then
               write (iur) itji
! write out spectrum sizes and data
               if (itji > 0) then
                  it = size(vpksji,3); is = size(vpksji,4)
                  write (iur) size(vpksji,1), size(vpksji,2), it, is
                  write (iur) ((((vpksji(i,j,k,l),i=1,2),j=1,4),k=1,it),&
     &l=1,is)
               endif
            else
               it = 0
               write (iur) it
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dread_drestart13(iur)
! read in restart diagnostic file for darwin code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j, k, l, it, is, ir, iq, ios
      character(len=32) :: fname
! read in restart diagnostic file for electrostatic code
      call dread_restart1(iur)
! restore energy accumulations
      if (ntw > 0) then
         if (itw > 0) then
            s(1) = wt(1,5)
            s(2) = wt(1,2)
            s(3) = wt(1,6)
            s(4) = wt(1,7)
            s(5) = wt(1,3)
            s(6) = dble(wt(1,4))
            s(7) = s(6)
            do it = 2, itw
               s(1) = s(1) + wt(it,5)
               s(2) = s(2) + wt(it,2)
               s(3) = s(3) + wt(it,6)
               s(4) = s(4) + wt(it,7)
               s(5) = s(5) + wt(it,3)
               s(6) = min(s(6),dble(wt(it,4)))
               s(7) = max(s(7),dble(wt(it,4)))
            enddo
         endif
      endif
!
! read in electron current density diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntje restart error, ios = ', ios
         stop
      endif
      if (it /= ntje) then
         write (*,*) 'restart error: read/expected ntje=', it, ntje
         stop
      endif
! read in record location
      if (ntje > 0) then
         read (iur,iostat=ios) njerec
         if (ios /= 0) then
            write (*,*) 'njerec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (njerec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
                  write (*,*) 'ntje record length error, ios = ', ios
                  stop
            endif
            if (it==0) then
                  write (*,*) 'ntje zero length record error'
                  stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
                  write (*,*) 'ntje file name error, ios = ', ios
                  stop
            endif
            fjename = fname
         endif
      endif
!
! read in vector potential diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'nta restart error, ios = ', ios
         stop
      endif
      if (it /= nta) then
         write (*,*) 'restart error: read/expected nta=', it, nta
         stop
      endif
! read in record location
      if (nta > 0) then
         read (iur,iostat=ios) narec
         if (ios /= 0) then
            write (*,*) 'narec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (narec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'nta record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'nta zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'nta file name error, ios = ', ios
               stop
            endif
            faname = fname
         endif
! read in spectrum flag
         read (iur,iostat=ios) ita
         if (ios /= 0) then
            write (*,*) 'ita restart error, ios = ', ios
            stop
         endif
! read in spectrum sizes and data
         if (ita > 0) then
            read (iur,iostat=ios) iq, ir, it, is
            if (ios /= 0) then
               write (*,*) 'nta restart array size error, ios=', ios
               stop
            endif
            if (iq /= 2) then
               write (*,*) 'vpks size error: read/expected 2 =', iq
               stop
            endif
            if (ir /= 4) then
               write (*,*) 'vpks size error: read/expected 4 =', ir
               stop
            endif
            if (it /= modesxa) then
               write (*,*) 'vpks size error: read/expected modesxa=',it,&
     &modesxa
                  stop
            endif
            if (is /= iw) then
               write (*,*) 'vpks size error: read/expected iw=', is, iw
               stop
            endif
            read (iur,iostat=ios) ((((vpks(i,j,k,l),i=1,2),j=1,4),      &
     &k=1,it),l=1,is)
            if (ios /= 0) then
               write (*,*) 'vpks array read error, ios = ', ios
               stop
            endif
         endif
      endif
!
! read in transverse efield diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntet restart error, ios = ', ios
         stop
      endif
      if (it /= ntet) then
         write (*,*) 'restart error: read/expected ntet=', it, ntet
         stop
      endif
! read in record location
      if (ntet > 0) then
         read (iur,iostat=ios) netrec
         if (ios /= 0) then
            write (*,*) 'netrec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (netrec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntet record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntet zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntet file name error, ios = ', ios
               stop
            endif
            fetname = fname
         endif
! read in spectrum flag
         read (iur,iostat=ios) itet
         if (ios /= 0) then
            write (*,*) 'itet restart error, ios = ', ios
            stop
         endif
! read in spectrum sizes and data
         if (itet > 0) then
            read (iur,iostat=ios) iq, ir, it, is
            if (ios /= 0) then
               write (*,*) 'ntet restart array size error, ios = ', ios
               stop
            endif
            if (iq /= 2) then
               write (*,*) 'vpkset size error: read/expected 2 =', iq
               stop
            endif
            if (ir /= 4) then
               write (*,*) 'vpkset size error: read/expected 4 =', ir
               stop
            endif
            if (it /= modesxet) then
               write (*,*) 'vpkset size error: read/expected modesxet=',&
     &it, modesxet
               stop
            endif
            if (is /= iw) then
               write (*,*) 'vpkset size error: read/expected iw=',is, iw
               stop
            endif
            read (iur,iostat=ios) ((((vpkset(i,j,k,l),i=1,2),j=1,4),    &
     &k=1,it),l=1,is)
            if (ios /= 0) then
               write (*,*) 'vpkset array read error, ios = ', ios
               stop
            endif
         endif
      endif
!
! read in magnetic field diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntb restart error, ios = ', ios
         stop
      endif
      if (it /= ntb) then
         write (*,*) 'restart error: read/expected nteb=', it, ntb
         stop
      endif
! read in record location
      if (ntb > 0) then
         read (iur,iostat=ios) nbrec
         if (ios /= 0) then
            write (*,*) 'nbrec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nbrec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntb record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntb zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntb file name error, ios = ', ios
               stop
            endif
            fbname = fname
         endif
      endif
!
! read in ion current density diagnostic parameter
      if (movion==1) then
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'ntji restart error, ios = ', ios
            stop
         endif
         if (it /= ntji) then
            write (*,*) 'restart error: read/expected ntji=', it, ntji
            stop
         endif
! read in record location
         if (ntji > 0) then
            read (iur,iostat=ios) njirec
            if (ios /= 0) then
               write (*,*) 'njirec restart error, ios = ', ios
               stop
            endif
! read in record length (zero if error) and file name (if no error)
            if (njirec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'ntji record length error, ios = ', ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'ntji zero length record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'ntji file name error, ios = ', ios
                  stop
               endif
               fjiname = fname
            endif
! read in spectrum flag
            read (iur,iostat=ios) itji
            if (ios /= 0) then
               write (*,*) 'itji restart error, ios = ', ios
               stop
            endif
! read in spectrum sizes and data
            if (itji > 0) then
               read (iur,iostat=ios) iq, ir, it, is
               if (ios /= 0) then
                  write (*,*) 'ntji restart array size error, ios=', ios
                  stop
               endif
               if (iq /= 2) then
                  write (*,*) 'vpksji size error: read/expected 2 =', iq
                  stop
               endif
               if (ir /= 4) then
                  write (*,*) 'vpksji size error: read/expected 4 =', ir
                  stop
               endif
               if (it /= modesxji) then
                  write (*,*) 'vpksji size error: read/expected modesxji&
     &=',it,modesxji
                  stop
               endif
               if (is /= iwi) then
                   write (*,*) 'vpksji size error: read/expected iwi=', &
     &is,iwi
                  stop
               endif
               read (iur,iostat=ios) ((((vpksji(i,j,k,l),i=1,2),j=1,4), &
     &k=1,it),l=1,is)
               if (ios /= 0) then
                   write (*,*) 'vpksji array read error, ios = ', ios
                  stop
               endif
            endif
         endif
      endif
      end subroutine
!
      end module
