!-----------------------------------------------------------------------
! 1-2/2D Electromagnetic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
      program mbbeps1
      use in1
      use minit1
      use mpush1
      use msort1
      use mgard1
      use mfft1
      use mfield1
      use mbpush1
      use mcurd1
      use mdiag1
      use graf1
      use omplib
      implicit none
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! wke/wki/we = electron/ion kinetic energies and electrostatic energy
! wf/wb = magnetic field/transverse electric field energies
      real :: wke = 0.0, wki = 0.0, we = 0.0, wf = 0.0, wb = 0.0
      real :: zero = 0.0
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: np, nx, nxh, nxe, nxeh
      integer :: mx1, ntime, nloop, isign
      integer :: npi = 0
      real :: qbme, affp, dth, omt, ws
      real :: qbmi, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, ntmax, npbmx, irc
!
! declare scalars for diagnostics
      integer :: it, iw, mtw, mtp, mta, mtet, mtar, mtv, mtt
      integer :: itw, itp, ita, itet, itar, itv, itt
      integer :: iwi, mtdi, mtji, itdi, itji
! default Fortran unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iua = 13, iuet = 14, iub = 15, iuar = 16
      integer :: iudi = 20, iuji = 21
      real :: wef, ts
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:), allocatable :: qe, qi
! fxe = smoothed longitudinal electric field with guard cells
      real, dimension(:), allocatable :: fxe
! cue/cui = electron/ion current density with guard cells
      real, dimension(:,:), allocatable :: cue, cui
! fxyze/byze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:), allocatable :: fxyze, byze
! eyz/byz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:), allocatable :: eyz, byz
! ffc = form factor array for poisson solver
      complex, dimension(:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
! arrays required for darwin initial fields
      real, dimension(:,:), allocatable :: amu, dcu
!
! declare arrays for OpenMP (tiled) code:
! ppart/pparti = tiled electron/ion particle arrays
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppart, pparti, ppbuff
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), allocatable :: kpic, kipic
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), allocatable :: ncl
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), allocatable :: ihole
!
! diagnostic arrays
! wt = energy time history array
      real, dimension(:,:), allocatable :: wt
! s = scratch array for energies
      double precision, dimension(:), allocatable :: s
! scratch arrays for scalar field
      complex, dimension(:), allocatable :: sfieldc
      real, dimension(:), allocatable :: sfield
      complex, dimension(:,:), allocatable :: vfieldc
      real, dimension(:,:), allocatable :: vfield
! scratch arrays for spectral analysis
      real, dimension(:), allocatable :: wm, wmi
! denet/denit = store selected fourier modes for electron/ion density
      complex, dimension(:), allocatable :: denet, denit
! pkwdi = power spectrum for ion density
      real, dimension(:,:,:), allocatable :: pkwdi
! pksdi = accumulated complex spectrum for ion density
      double precision, dimension(:,:,:), allocatable :: pksdi
! wkdi = maximum frequency as a function of k for ion density
      real, dimension(:,:), allocatable :: wkdi
! pott = store selected fourier modes for potential
      complex, dimension(:), allocatable :: pott
! pkw = power spectrum for potential
      real, dimension(:,:,:), allocatable :: pkw
! pks = accumulated complex spectrum for potential
      double precision, dimension(:,:,:), allocatable :: pks
! wk = maximum frequency as a function of k for potential
      real, dimension(:,:), allocatable :: wk
! elt = store selected fourier modes for longitudinal efield
      complex, dimension(:), allocatable :: elt
! curit = store selected fourier modes for ion current density
      complex, dimension(:,:), allocatable :: curit
! vpkwji = power spectrum for ion current density
      real, dimension(:,:,:,:), allocatable :: vpkwji
! vpksji = accumulated complex spectrum for ion current density
      double precision, dimension(:,:,:,:), allocatable :: vpksji
! vwkji = maximum frequency as a function of k for ion current
      real, dimension(:,:,:), allocatable :: vwkji
! vpotr = store selected fourier modes for radiative vector potential
      complex, dimension(:,:), allocatable :: vpotr
! vpkwr = power spectrum for radiative vector potential
      real, dimension(:,:,:,:), allocatable :: vpkwr
! vpksr = accumulated complex spectrum for radiative vector potential
      double precision, dimension(:,:,:,:), allocatable :: vpksr
! vwkr = maximum frequency as a function of k for radiative vector
!        potential
      real, dimension(:,:,:), allocatable :: vwkr
! oldcu = previous current density with guard cells
      real, dimension(:,:), allocatable :: oldcu
! vpott = store selected fourier modes for vector potential
      complex, dimension(:,:), allocatable :: vpott
! vpkw = power spectrum for vector potential
      real, dimension(:,:,:,:), allocatable :: vpkw
! vpks = accumulated complex spectrum for vector potential
      double precision, dimension(:,:,:,:), allocatable :: vpks
! vwk = maximum frequency as a function of k for vector potential
      real, dimension(:,:,:), allocatable :: vwk
! ett = store selected fourier modes for transverse efield
      complex, dimension(:,:), allocatable :: ett
! vpkwet = power spectrum for transverse efield
      real, dimension(:,:,:,:), allocatable :: vpkwet
! vpkset = accumulated complex spectrum for transverse efield
      double precision, dimension(:,:,:,:), allocatable :: vpkset
! vwket = maximum frequency as a function of k for transverse efield
      real, dimension(:,:,:), allocatable :: vwket
! bt = store selected fourier modes for magnetic field
      complex, dimension(:,:), allocatable :: bt
! sfv/sfvi = electron/ion velocity distribution functions in tile
      real, dimension(:,:,:), allocatable :: sfv, sfvi
! fvm/fvmi = electron/ion vdrift, vth, entropy for global distribution
      real, dimension(:,:), allocatable :: fvm, fvmi
! fvtm/fvtmi = time history of electron/ion vdrift, vth, and entropy
!              for global distribution
      real, dimension(:,:,:), allocatable :: fvtm, fvtmi
! iprobt = scratch array 
      integer, dimension(:), allocatable :: iprobt
! partt = particle trajectories tracked
      real, dimension(:,:), allocatable :: partt
! fvtp = velocity distribution function for test particles
! fvmtp = vdrift, vth, and entropy for test particles
      real, dimension(:,:), allocatable :: fvtp, fvmtp
! partd = trajectory time history array
      real, dimension(:,:,:), allocatable :: partd
! cwk = labels for power spectrum display
      character(len=10), dimension(2) :: cwk = (/'   W > 0  ',          &
     &                                           '   W < 0  ' /)
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tdiag = 0.0
      double precision :: dtime
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
! override default input data
      emf = 1
      relativity = 1
! read namelist
      call readnml1(iuin)
! override input data
      idcode = 2
      ndim = 3
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
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
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
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
      if (list) then
         if (ipbc.ne.1) then
            write (*,*) 'ipbc /= 1 and list = .true. not yet supported'
            list = .false.
            write (*,*) 'list reset to .false.'
         endif
      endif
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxe),qi(nxe),fxe(nxe))
      allocate(fxyze(3,nxe),cue(2,nxe),byze(2,nxe))
      allocate(eyz(2,nxeh),byz(2,nxeh))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
!
! prepare fft tables
      call mfft1_init(mixup,sct,indx)
! calculate form factors
      call mpois1_init(ffc,ax,affp,nx)
! initialize different ensemble of random numbers
      if (nextrand > 0) call mnextran1(nextrand,ndim,np+npi)
!
! initialize electrons
! background electrons
      if (npx > 0) then
!        call mudistr1(part,1,npx,nx,ipbc)
         call mfdistr1(part,ampdx,scaledx,shiftdx,1,npx,nx,ipbc,ndprof)
         call wmvdistr1h(part,1,vtx,vty,vtz,vx0,vy0,vz0,npx,nvdist)
      endif
! beam electrons
      if (npxb > 0) then
         it = npx + 1
!        call mudistr1(part,it,npxb,nx,ipbc)
         call mfdistr1(part,ampdx,scaledx,shiftdx,it,npxb,nx,ipbc,ndprof&
     &)
         call wmvdistr1h(part,it,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,nvdist)
      endif
!
! mark electron beam particles
      if ((nts > 0).and.(ntsc > 0)) then
         call setmbeam1(part,npx)
      endif
!
! find number of electrons in each of mx, tiles: updates kpic, nppmx
      call mdblkp2(part,kpic,nppmx,mx,irc)
!
! allocate vector electron data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      allocate(ppart(idimp,nppmx0,mx1))
      allocate(ppbuff(idimp,npbmx,mx1))
      allocate(ncl(2,mx1))
      allocate(ihole(2,ntmax+1,mx1))
! copy ordered electron data for OpenMP: updates ppart and kpic
      call mpmovin1(part,ppart,kpic,mx,irc)
!
! sanity check for electrons
      call mcheck1(ppart,kpic,nx,mx,irc)
      deallocate(part)
!
! initialize background charge density: updates qi
      if (movion==0) then
         qi = 0.0
         call mpost1(ppart,qi,kpic,-qme,tdpost,mx)
         call maguard1(qi,tguard,nx)
      endif
!
! initialize ions
      if (movion==1) then
         allocate(part(idimp,npi),kipic(mx1),cui(2,nxe))
         cui = 0.0
! background ions
         if (npxi > 0) then
!           call mudistr1(part,1,npxi,nx,ipbc)
            call mfdistr1(part,ampdxi,scaledxi,shiftdxi,1,npxi,nx,ipbc, &
     &ndprofi)
            call wmvdistr1h(part,1,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,npxi,  &
     &nvdist)
         endif
! beam ions
         if (npxbi > 0) then
            it = npxi + 1
!           call mudistr1(part,it,npxbi,nx,ipbc)
            call mfdistr1(part,ampdxi,scaledxi,shiftdxi,it,npxbi,nx,ipbc&
     &,ndprofi)
            call wmvdistr1h(part,it,vtdxi,vtdyi,vtdzi,vdxi,vdyi,vdzi,   &
     &npxbi,nvdist)
         endif
!
! mark ion beam particles
         if ((nts > 0).and.(ntsc > 0)) then
            call setmbeam1(part,npxi)
         endif
!
! find number of ions in each of mx, tiles: updates kipic, nppmx
         call mdblkp2(part,kipic,nppmx,mx,irc)
!
! allocate vector ion data
         nppmx1 = (1.0 + xtras)*nppmx
         allocate(pparti(idimp,nppmx1,mx1))
! copy ordered ion data for OpenMP: updates pparti and kipic
         call mpmovin1(part,pparti,kipic,mx,irc)
!
! sanity check for ions
         call mcheck1(pparti,kipic,nx,mx,irc)
         deallocate(part)
      endif
!
! initialize transverse electromagnetic fields
      eyz = cmplx(0.0,0.0)
      byz = cmplx(0.0,0.0)
      cue = 0.0
! set magnitude of external transverse magnetic field
      omt = sqrt(omy*omy + omz*omz)
!
! allocate diagnostic arrays
! reverse simulation at end back to start
      if (treverse==1) nloop = 2*nloop
!
! energy time history
      if (ntw > 0) then
         mtw = (nloop - 1)/ntw + 1; itw = 0
         allocate(wt(mtw,7),s(7))
         wt = 0.0; s = 0.0d0
      endif
!
! allocate scratch arrays for scalar fields
      if ((ntde > 0).or.(ntp > 0).or.(ntel > 0).or.(ntdi > 0)) then
         allocate(sfieldc(nxh),sfield(nxe))
      endif
!
! allocate and initialize frequency array for spectral analysis
      if ((ntp > 0).or.(nta > 0).or.(ntet > 0).or.(ntar > 0)) then
         iw = (wmax - wmin)/dw + 1.5
         allocate(wm(iw))
         do it = 1, iw
         wm(it) = wmin + dw*real(it-1)
         enddo
      endif
!
! allocate and initialize frequency array for ion spectral analysis
      if (movion==1) then
         if ((ntdi > 0).or.(ntji > 0)) then
            iwi = (wimax - wimin)/dwi + 1.5
            allocate(wmi(iwi))
            do it = 1, iwi
            wmi(it) = wimin + dwi*real(it-1)
            enddo
         endif
      endif
!
! allocate scratch arrays for vector fields
      if ((nta>0).or.(ntet>0).or.(ntb>0).or.(ntar>0).or.(ntji>0)) then
         allocate(vfieldc(2,nxh),vfield(2,nxe))
      endif
!
! initialize electron density diagnostic
      if (ntde > 0) then
         fdename = 'denek1.'//cdrun
         modesxde = min(modesxde,nxh+1)
         allocate(denet(modesxde))
! open file: updates nderec and possibly iude
         if (nderec==0) call dafopenc1(denet,iude,nderec,trim(fdename))
      endif
!
! initialize ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) then
            fdiname = 'denik1.'//cdrun
            modesxdi = min(modesxdi,nxh+1)
            allocate(denit(modesxdi))
! open file: updates ndirec and possibly iudi
            if (ndirec==0) then
               call dafopenc1(denit,iudi,ndirec,trim(fdiname))
            endif
! ion spectral analysis
            if ((nddi==2).or.(nddi==3)) then
               mtdi = (nloop - 1)/ntdi + 1; itdi = 0
               allocate(pkwdi(modesxdi,iwi,2),pksdi(4,modesxdi,iwi))
               allocate(wkdi(modesxdi,2))
               pksdi = 0.0d0
            endif
         endif
      endif
!
! initialize potential diagnostic
      if (ntp > 0) then
         fpname = 'potk1.'//cdrun
         modesxp = min(modesxp,nxh+1)
         allocate(pott(modesxp))
! open file: updates nprec and possibly iup
         if (nprec==0) call dafopenc1(pott,iup,nprec,trim(fpname))
! spectral analysis
         if ((ndp==2).or.(ndp==3)) then
            mtp = (nloop - 1)/ntp + 1; itp = 0
            allocate(pkw(modesxp,iw,2),pks(4,modesxp,iw))
            allocate(wk(modesxp,2))
            pks = 0.0d0
         endif
      endif
!
! initialize longitudinal efield diagnostic
      if (ntel > 0) then
         felname = 'elk1.'//cdrun
         modesxel = min(modesxel,nxh+1)
         allocate(elt(modesxel))
! open file: updates nelrec and possibly iuel
         if (nelrec==0) call dafopenc1(elt,iuel,nelrec,trim(felname))
      endif
!
! initialize ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            fjiname = 'curik1.'//cdrun
            modesxji = min(modesxji,nxh+1)
            allocate(curit(2,modesxji))
! open file: updates njirec and possibly iuji
            if (njirec==0) then
               call dafopenvc1(curit,iuji,njirec,trim(fjiname))
            endif
! ion spectral analysis
            if ((ndji==2).or.(ndji==3)) then
               mtji = (nloop - 1)/ntji + 1; itji = 0
               allocate(vpkwji(2,modesxji,iwi,2))
               allocate(vpksji(2,4,modesxji,iwi),vwkji(2,modesxji,2))
               vpksji = 0.0d0
            endif
         endif
      endif
!
! initialize radiative vector potential diagnostic
      if (ntar > 0) then
         farname = 'vpotrk1.'//cdrun
         modesxar = min(modesxar,nxh+1)
         allocate(vpotr(2,modesxar))
! open file: updates narec and possibly iua
         if (narec==0) call dafopenvc1(vpotr,iuar,narrec,trim(farname))
! spectral analysis
         if ((ndar==2).or.(ndar==3)) then
            mtar = (nloop - 1)/ntar + 1; itar = 0
            allocate(vpkwr(2,modesxar,iw,2),vpksr(2,4,modesxar,iw))
            allocate(vwkr(2,modesxar,2),oldcu(2,nxe))
            vpksr = 0.0d0; oldcu = 0.0
         endif
      endif
!
! initialize vector potential diagnostic
      if (nta > 0) then
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
      endif
!
! initialize transverse efield diagnostic
      if (ntet > 0) then
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
      endif
!
! initialize magnetic field diagnostic
      if (ntb > 0) then
         fbname = 'bk1.'//cdrun
         modesxb = min(modesxb,nxh+1)
         allocate(bt(2,modesxb))
! open file: updates nbrec and possibly iub
         if (netrec==0) call dafopenvc1(bt,iub,nbrec,trim(fbname))
      endif
!
! initialize velocity diagnostic
      if (ntv > 0) then
         allocate(sfv(2*nmv+2,ndim,mx1+1),fvm(ndim,3))
         mtv = (nloop - 1)/ntv + 1; itv = 0
         allocate(fvtm(mtv,ndim,3))
         ws = 2.0*max(4.0*vtx+abs(vx0),4.0*vtdx+abs(vdx))
         ws = max(ws,2.0*max(4.0*vty+abs(vy0),4.0*vtdy+abs(vdy)))
         ws = max(ws,2.0*max(4.0*vtz+abs(vz0),4.0*vtdx+abs(vdz)))
         sfv(1,1,:) = ws
         sfv(1,2,:) = ws
         sfv(1,3,:) = ws
         fvtm = 0.0
! ions
         if (movion==1) then
            allocate(sfvi(2*nmv+2,ndim,mx1+1),fvmi(ndim,3))
            allocate(fvtmi(mtv,ndim,3))
            ws = 2.0*max(4.0*vtxi+abs(vxi0),4.0*vtdxi+abs(vdxi))
            ws = max(ws,2.0*max(4.0*vtyi+abs(vyi0),4.0*vtdyi+abs(vdyi)))
            ws = max(ws,2.0*max(4.0*vtzi+abs(vzi0),4.0*vtdzi+abs(vdzi)))
            sfvi(1,1,:) = ws
            sfvi(1,2,:) = ws
            sfvi(1,3,:) = ws
            fvtmi = 0.0
         endif
      endif
!
! initialize trajectory diagnostic
      if (ntt > 0) then
         allocate(iprobt(nprobt))
         call setptraj1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,np,nprobt)
         if (nprobt.gt.16777215) then
            write(*,*) 'nprobt overflow = ', nprobt
            stop
         endif
         allocate(partt(idimp,nprobt))
         if ((nst==1).or.(nst==2)) then
            mtt = (nloop - 1)/ntt + 1; itt = 0
            allocate(partd(mtt,idimp,nprobt))
         else if (nst==3) then
            allocate(fvtp(2*nmv+2,ndim),fvmtp(ndim,3))
            fvtp(1,:) = 2.0*max(4.0*vtx+abs(vx0),4.0*vtdx+abs(vdx))
         endif
      endif
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
      if (dt > 0.64*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
      write (iuot,*) 'program mbbeps1'
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
      ntime = n - 1
      write (iuot,*) 'ntime = ', ntime
!
! save previous current in fourier space for radiative vector potential
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) oldcu = cue
      endif
!
! deposit electron current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      call wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,dth,ci,tdjpost,nx,mx, &
     &ipbc,relativity,list,irc)
! add guard cells: updates cue
      call macguard1(cue,tguard,nx)
!
! reorder particles by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
! sanity check for electrons
      call mcheck1(ppart,kpic,nx,mx,irc)
!
! deposit ion current with OpenMP:
      if (movion==1) then
! updates pparti and cui, and possibly ncl, ihole, irc
         call dtimer(dtime,itime,-1)
         cui = 0.0
         call dtimer(dtime,itime,1)
         tdjpost = tdjpost + real(dtime)
         call wmdjpost1(pparti,cui,kipic,ncl,ihole,qmi,dth,ci,tdjpost,nx&
     &,mx,ipbc,relativity,list,irc)
! add guard cells: updates cui
         call macguard1(cui,tguard,nx)
!
! reorder particles by tile with OpenMP:
! updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
         call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,list, &
     &irc)
! sanity check for ions
         call mcheck1(pparti,kipic,nx,mx,irc)
      endif
!
! ion current density diagnostic
      if (movion==1) then
         if (ntji > 0) then
            it = ntime/ntji
            if (ntime==ntji*it) then
               vfield = cui
! transform ion current density to fourier space: updates vfield
               isign = -1
               call mfft1rn(vfield,isign,mixup,sct,tfft,indx)
! calculate smoothed ion current in fourier space: updates vfieldc
               call msmooth13(vfield,vfieldc,ffc,tfield,nx)
! store selected fourier modes: updates curit
               call mrdvmodes1(vfieldc,curit,tfield,nx,modesxji)
! write diagnostic output: updates ndirec
               call dafwritevc1(curit,tdiag,iuji,njirec,modesxji)
! transform smoothed ion current to real space: updates vfield
               if ((ndji==1).or.(ndji==3)) then
                  call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
                  call mcguard1(vfield,tguard,nx)
! display smoothed ion current
                  call dvector1(vfield,' ION CURRENT',ntime,999,0,2,nx, &
     &irc)
                  if (irc==1) exit; irc = 0
               endif
! ion spectral analysis
               if ((ndji==2).or.(ndji==3)) then
                  itji = itji + 1
                  ts = dt*real(ntime)
                  call mivcspect1(curit,wmi,vpkwji,vpksji,ts,t0,tdiag,  &
     &mtji,iwi,modesxji,nx,-1)
! performs frequency analysis of accumulated complex time series
                  vwkji(1,:,1) = wmi(maxloc(vpkwji(1,:,:,1),dim=2))
                  vwkji(2,:,1) = wmi(maxloc(vpkwji(2,:,:,1),dim=2))
                  vwkji(1,:,2) = wmi(maxloc(vpkwji(1,:,:,2),dim=2))
                  vwkji(2,:,2) = wmi(maxloc(vpkwji(2,:,:,2),dim=2))
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
! electron density diagnostic
      if (ntde > 0) then
         it = ntime/ntde
         if (ntime==ntde*it) then
            sfield = -qe
! transform electron density to fourier space: updates sfield
            isign = -1
            call mfft1r(sfield,isign,mixup,sct,tfft,indx)
! calculate smoothed density in fourier space: updates sfieldc
            call msmooth1(sfield,sfieldc,ffc,tfield,nx)
! store selected fourier modes: updates denet
            call mrdmodes1(sfieldc,denet,tfield,nx,modesxde)
! write diagnostic output: updates nderec
            call dafwritec1(denet,tdiag,iude,nderec,modesxde)
! transform smoothed electron density to real space: updates sfield
            call mfft1cr(sfieldc,sfield,mixup,sct,tfft,indx)
            call mdguard1(sfield,tguard,nx)
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
! ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) then
            it = ntime/ntdi
            if (ntime==ntdi*it) then
               sfield = qi
! transform ion density to fourier space: updates sfield
               isign = -1
               call mfft1r(sfield,isign,mixup,sct,tfft,indx)
! calculate smoothed density in fourier space: updates sfieldc
               call msmooth1(sfield,sfieldc,ffc,tfield,nx)
! store selected fourier modes: updates denit
               call mrdmodes1(sfieldc,denit,tfield,nx,modesxdi)
! write diagnostic output: updates ndirec
               call dafwritec1(denit,tdiag,iudi,ndirec,modesxdi)
! transform smoothed ion density to real space: updates sfield
               if ((nddi==1).or.(nddi==3)) then
                  call mfft1cr(sfieldc,sfield,mixup,sct,tfft,indx)
                  call mdguard1(sfield,tguard,nx)
! display smoothed ion density
                  call dscaler1(sfield,' ION DENSITY',ntime,999,1,nx,irc&
     &)
                  if (irc==1) exit; irc = 0
               endif
! ion spectral analysis
               if ((nddi==2).or.(nddi==3)) then
                  itdi = itdi + 1
                  ts = dt*real(ntime)
                  call micspect1(denit,wmi,pkwdi,pksdi,ts,t0,tdiag,mtdi,&
     &iwi,modesxdi,nx,-1)
! performs frequency analysis of accumulated complex time series
                  wkdi(:,1) = wmi(maxloc(pkwdi(:,:,1),dim=2))
                  wkdi(:,2) = wmi(maxloc(pkwdi(:,:,2),dim=2))
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
! radiative vector potential diagnostic
      if (ntar > 0) then
         it = ntime/ntar
         if (ntime==ntar*it) then
! average current: updates vfieldc = 0.5*(cue + oldcu)
            call mcuave1(vfieldc,cue,oldcu,tfield,nx)
! calculate radiative vector potential in fourier space: updates vfieldc
! vfieldc should contain averaged current on entry
            call mavrpot1(vfieldc,byz,ffc,ci,tfield,nx)
! store selected fourier modes: updates vpotr
            call mrdvmodes1(vfieldc,vpotr,tfield,nx,modesxar)
! write diagnostic output: updates narrec
            call dafwritevc1(vpotr,tdiag,iuar,narrec,modesxar)
! transform radiative vector potential to real space: updates vfield
            if ((ndp==1).or.(ndp==3)) then
               call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
               call mcguard1(vfield,tguard,nx)
! display radiative vector potential
               call dvector1(vfield,' RADIATIVE VPOTENTIAL',ntime,999,0,&
     &2,nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndar==2).or.(ndar==3)) then
               itar = itar + 1
               ts = dt*real(ntime)
               call mivcspect1(vpotr,wm,vpkwr,vpksr,ts,t0,tdiag,mtar,iw,&
     &modesxar,nx,1)
! performs frequency analysis of accumulated complex vector time series
               vwkr(1,:,1) = wm(maxloc(vpkwr(1,:,:,1),dim=2))
               vwkr(2,:,1) = wm(maxloc(vpkwr(2,:,:,1),dim=2))
               vwkr(1,:,2) = wm(maxloc(vpkwr(1,:,:,2),dim=2))
               vwkr(2,:,2) = wm(maxloc(vpkwr(2,:,:,2),dim=2))
! display frequency spectrum
               call dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',&
     &ntime,999,2,2,modesxar,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! calculate electromagnetic fields in fourier space:
! updates eyz, byz
      if (ntime==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mibpois1(cue,byz,ffc,ci,wb,tfield,nx)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(2,nxe),dcu(2,nxe))
         amu = 0.0
         call wmgmjpost1(ppart,amu,kpic,qme,ci,tdjpost,mx,relativity)
         call macguard1(amu,tguard,nx)
         isign = -1
         call mfft1rn(amu,isign,mixup,sct,tfft,indx)
         call mdcuperp1(dcu,amu,tfield,nx)
         deallocate(amu)
         call metfield1(dcu,eyz,ffc,ci,wf,tfield,nx)
         deallocate(dcu)
         dth = 0.5*dt
      else
         call mmaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wb,tfield,nx)
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
! add external traveling wave field
      ts = dt*real(ntime)
      call meaddext13(fxyze,tfield,amodex,freq,ts,trmp,toff,el0,er0,nx)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mbaddext1(byze,tfield,omy,omz,nx)
!
! copy guard cells: updates fxyze, byze
      call mcguard1(fxyze,tguard,nx)
      call mcguard1(byze,tguard,nx)
!
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space: updates sfieldc
            call mpot1(qe,sfieldc,ffc,ws,tfield,nx)
! store selected fourier modes: updates pott
            call mrdmodes1(sfieldc,pott,tfield,nx,modesxp)
! write diagnostic output: updates nprec
            call dafwritec1(pott,tdiag,iup,nprec,modesxp)
! transform potential to real space: updates sfield
            if ((ndp==1).or.(ndp==3)) then
               call mfft1cr(sfieldc,sfield,mixup,sct,tfft,indx)
               call mdguard1(sfield,tguard,nx)
! display potential
               call dscaler1(sfield,' POTENTIAL',ntime,999,0,nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndp==2).or.(ndp==3)) then
               itp = itp + 1
               ts = dt*real(ntime)
               call micspect1(pott,wm,pkw,pks,ts,t0,tdiag,mtp,iw,modesxp&
     &,nx,1)
! performs frequency analysis of accumulated complex time series
               wk(:,1) = wm(maxloc(pkw(:,:,1),dim=2))
               wk(:,2) = wm(maxloc(pkw(:,:,2),dim=2))
! display frequency spectrum
               call dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2, &
     &modesxp,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! longitudinal efield diagnostic
      if (ntel > 0) then
         it = ntime/ntel
         if (ntime==ntel*it) then
! calculate longitudinal efield in fourier space: updates sfieldc
            call melfield1(qe,sfieldc,ffc,ws,tfield,nx)
! store selected fourier modes: updates elt
            call mrdmodes1(sfieldc,elt,tfield,nx,modesxel)
! write diagnostic output: updates nelrec
            call dafwritec1(elt,tdiag,iuel,nelrec,modesxel)
! transform longitudinal efield to real space: updates sfield
            call mfft1cr(sfieldc,sfield,mixup,sct,tfft,indx)
            call mdguard1(sfield,tguard,nx)
! display longitudinal efield 
            call dscaler1(sfield,' ELFIELD',ntime,999,0,nx,irc)
            if (irc==1) exit; irc = 0
         endif
      endif
!
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space: updates vfieldc
            call mavpot1(byz,vfieldc,tfield,nx)
! store selected fourier modes: updates vpott
            call mrdvmodes1(vfieldc,vpott,tfield,nx,modesxa)
! write diagnostic output: updates narec
            call dafwritevc1(vpott,tdiag,iua,narec,modesxa)
! transform vector potential to real space: updates vfield
            if ((nda==1).or.(nda==3)) then
               call mfft1crn(vfieldc,vfield,mixup,sct,tfft,indx)
               call mcguard1(vfield,tguard,nx)
! display vector potential
               call dvector1(vfield,' VECTOR POTENTIAL',ntime,999,0,2,nx&
     &,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((nda==2).or.(nda==3)) then
               ita = ita + 1
               ts = dt*real(ntime)
               call mivcspect1(vpott,wm,vpkw,vpks,ts,t0,tdiag,mta,iw,   &
     &modesxa,nx,1)
! performs frequency analysis of accumulated complex vector time series
               vwk(1,:,1) = wm(maxloc(vpkw(1,:,:,1),dim=2))
               vwk(2,:,1) = wm(maxloc(vpkw(2,:,:,1),dim=2))
               vwk(1,:,2) = wm(maxloc(vpkw(1,:,:,2),dim=2))
               vwk(2,:,2) = wm(maxloc(vpkw(2,:,:,2),dim=2))
! display frequency spectrum
               call dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime&
     &,999,2,2,modesxa,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! transverse efield diagnostic
      if (ntet > 0) then
         it = ntime/ntet
         if (ntime==ntet*it) then
! store selected fourier modes: updates ett
            call mrdvmodes1(eyz,ett,tfield,nx,modesxet)
! write diagnostic output: updates netrec
            call dafwritevc1(ett,tdiag,iuet,netrec,modesxet)
! transform transverse efield to real space: updates vfield
            if ((ndet==1).or.(ndet==3)) then
               call mfft1crn(eyz,vfield,mixup,sct,tfft,indx)
               call mcguard1(vfield,tguard,nx)
! display transverse efield
               call dvector1(vfield,' TRANSVERSE EFIELD',ntime,999,0,2, &
     &nx,irc)
               if (irc==1) exit; irc = 0
            endif
! spectral analysis
            if ((ndet==2).or.(ndet==3)) then
               itet = itet + 1
               ts = dt*real(ntime)
               call mivcspect1(ett,wm,vpkwet,vpkset,ts,t0,tdiag,mtet,iw,&
     &modesxet,nx,0)
! performs frequency analysis of accumulated complex vector time series
               vwket(1,:,1) = wm(maxloc(vpkwet(1,:,:,1),dim=2))
               vwket(2,:,1) = wm(maxloc(vpkwet(2,:,:,1),dim=2))
               vwket(1,:,2) = wm(maxloc(vpkwet(1,:,:,2),dim=2))
               vwket(2,:,2) = wm(maxloc(vpkwet(2,:,:,2),dim=2))
! display frequency spectrum
               call dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',  &
     &ntime,999,2,2,modesxet,cwk,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! magnetic field diagnostic
      if (ntb > 0) then
         it = ntime/ntb
         if (ntime==ntb*it) then
! store selected fourier modes: updates bt
            call mrdvmodes1(byz,bt,tfield,nx,modesxb)
! write diagnostic output: updates nbrec
            call dafwritevc1(bt,tdiag,iub,nbrec,modesxb)
! transform magnetic field to real space: updates vfield
            call mfft1crn(byz,vfield,mixup,sct,tfft,indx)
            call mcguard1(vfield,tguard,nx)
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
! calculate electron distribution function and moments
            call mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,nmv)
! store time history electron vdrift, vth, and entropy
            itv = itv + 1
            fvtm(itv,:,:) = fvm
! display electron velocity distributions
            call displayfv1(sfv(:,:,mx1+1),fvm,' ELECTRON',ntime,nmv,2, &
     &irc)
            if (irc==1) exit; irc = 0
! ion distribution function
            if (movion==1) then
               call mvpdist1(pparti,kipic,sfvi,fvmi,tdiag,npi,nmv)
! store time history of ion vdrift, vth, and entropy
               fvtmi(itv,:,:) = fvmi
! display ion velocity distributions
               call displayfv1(sfvi(:,:,mx1+1),fvmi,' ION',ntime,nmv,2, &
     &irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! trajectory diagnostic
      if (ntt > 0) then
         it = ntime/ntt
         if (ntime==ntt*it) then
! copies trajectories to array partt
            call mptraj1(ppart,kpic,partt,tdiag)
            itt = itt + 1
            if ((nst==1).or.(nst==2)) then
               partd(itt,:,:) = partt
            else if (nst==3) then
! calculate test particle distribution function and moments
               call mvdist1(partt,fvtp,fvmtp,tdiag,nprobt,nmv)
! display test particle velocity distributions
               call displayfv1(fvtp,fvmtp,' ELECTRON',ntime,nmv,1,irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! phase space diagnostic
      if (nts > 0) then
         it = ntime/nts
         if (ntime==nts*it) then
! plot electrons vx versus x
            call dpmgrasp1(ppart,kpic,' ELECTRON',ntime,999,nx,2,1,ntsc,&
     &irc)
            if (irc==1) exit; irc = 0
! ion phase space
            if (movion==1) then
! plot ions vx versus x
               call dpmgrasp1(pparti,kipic,' ION',ntime,999,nx,2,1,ntsc,&
     &irc)
               if (irc==1) exit; irc = 0
            endif
         endif
      endif
!
! push electrons with OpenMP:
      wke = 0.0
! updates ppart and wke, and possibly ncl, ihole, irc
      if (mzf==0) then
         call wmbpush1(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dth, &
     &ci,wke,tpush,nx,mx,ipbc,relativity,list,irc)
! zero force: updates ppart, wke and possibly ncl, ihole, and irc
      else
         call wmpush1zf(ppart,kpic,ncl,ihole,dth,ci,wke,tpush,nx,mx,ipbc&
     &,relativity,list,irc)
      endif
!
! reorder electrons by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc, and possibly ihole
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,irc)
!
! sanity check for electrons
      call mcheck1(ppart,kpic,nx,mx,irc)
!
! push ions with OpenMP:
      if (movion==1) then
! updates pparti and wki, and possibly ncl, ihole, irc
         wki = 0.0
         if (mzf==0) then
            call wmbpush1(pparti,fxyze,byze,kipic,ncl,ihole,omx,qbmi,dt,&
     &dth,ci,wki,tpush,nx,mx,ipbc,relativity,list,irc)
         else
            call wmpush1zf(pparti,kipic,ncl,ihole,dth,ci,wki,tpush,nx,mx&
     &,ipbc,relativity,list,irc)
         endif
         wki = wki*rmass
!
! reorder ions by tile with OpenMP:
! updates pparti, ppbuff, kipic, ncl, irc, and possibly ihole
         call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,list, &
     &irc)
!
! sanity check for ions
         call mcheck1(pparti,kipic,nx,mx,irc)
      endif
!
! start running simulation backwards:
! need to advance maxwell field solver one step ahead
      if (treverse==1) then
         if (((ntime+1)==(nloop/2)).or.((ntime+1)==nloop)) then
! deposit electron current: updates cue
            cue = 0.0
            call wmdjpost1(ppart,cue,kpic,ncl,ihole,qme,zero,ci,tdjpost,&
     &nx,mx,ipbc,relativity,list,irc)
            call macguard1(cue,tguard,nx)
! deposit ion current: updates cui
            if (movion==1) then
               cui = 0.0
               call wmdjpost1(pparti,cui,kipic,ncl,ihole,qmi,zero,ci,   &
     &tdjpost,nx,mx,ipbc,relativity,list,irc)
               call macguard1(cui,tguard,nx)
               call maddcuei1(cue,cui,tfield,nx)
            endif
            isign = -1
            call mfft1rn(cue,isign,mixup,sct,tfft,indx)
! updates eyz, byz, wf, wb
            call mmaxwel1(eyz,byz,cue,ffc,ci,dt,wf,wb,tfield,nx)
! reverse time
            dt = -dt; dth = -dth
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         it = ntime/ntw
         if (ntime==ntw*it) then
            wef = we + wf + wb
            ws = wef + wke + wki
            if (ntime==0) s(6) = ws
            write (iuot,*) 'Total Field, Kinetic and Total Energies:'
            if (movion==0) then
               write (iuot,'(3e14.7)') wef, wke, ws
            else
               write (iuot,'(4e14.7)') wef, wke, wki, ws
            endif
            write (iuot,*) 'Electric(l,t) and Magnetic Energies = '
            write (iuot,'(3e14.7)') we, wf, wb
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
      if ((ntw > 0).or.(ntt > 0)) then
         if (nplot > 0) call reset_graphs
      endif
!
! trajectory diagnostic
      if (ntt > 0) then
         if ((nst==1).or.(nst==2)) then
            if (nplot > 0) irc = open_graphs(1)
            ts = t0 + dt*real(ntt)
            call displaytr1(partd,ts,dt*real(ntt),itt,2,3,irc)
            if (irc==1) stop
            call reset_nplot(nplot,irc)
         endif
      endif
!
! energy diagnostic
      if (ntw > 0) then
         ts = t0 + dt*real(ntw)
         call displayw1(wt,ts,dt*real(ntw),itw,irc)
         s(6) = (s(7) - s(6))/wt(1,4)
         write (iuot,*) 'Energy Conservation = ', real(s(6))
         s(1) = s(1)/real(itw)
         write (iuot,*) 'Average Field Energy <WE> = ', real(s(1))
         s(2) = s(2)/real(itw)
         write (iuot,*) 'Average Electron Kinetic Energy <WKE> = ',     &
     &real(s(2))
         write (iuot,*) 'Ratio <WE>/<WKE>= ', real(s(1)/s(2))
         s(3) = s(3)/real(itw)
         write (iuot,*) 'Average Transverse EField Energy <WF> = ',     &
     &real(s(3))
         write (iuot,*) 'Ratio <WF>/<WKE>= ', real(s(3)/s(2))
         s(4) = s(4)/real(itw)
         write (iuot,*) 'Average Magnetic Field Energy <WB> = ',        &
     &real(s(4))
         write (iuot,*) 'Ratio <WB>/<WKE>= ', real(s(4)/s(2))
      endif
!
! velocity diagnostic
      if (ntv > 0) then
         ts = t0 + dt*real(ntv)
         call displayfvt1(fvtm,' ELECT',ts,dt*real(ntv),itv,irc)
! ions
         if (movion==1) then
            call displayfvt1(fvtmi,' ION',ts,dt*real(ntv),itv,irc)
         endif
      endif
!
! display final spectral analysis for ion density
      if (movion==1) then
         if (ntdi > 0) then
            if ((nddi==2).or.(nddi==3)) then
! performs frequency analysis of accumulated complex time series
               wkdi(:,1) = wmi(maxloc(pkwdi(:,:,1),dim=2))
               wkdi(:,2) = wmi(maxloc(pkwdi(:,:,2),dim=2))
! display frequency spectrum
               call dmscaler1(wkdi,'ION DENSITY OMEGA VS MODE',ntime,   &
     &999,1,modesxdi,cwk,irc)
            endif
         endif
      endif
!
! display final spectral analysis for potential
      if (ntp > 0) then
         if ((ndp==2).or.(ndp==3)) then
! performs frequency analysis of accumulated complex time series
            wk(:,1) = wm(maxloc(pkw(:,:,1),dim=2))
            wk(:,2) = wm(maxloc(pkw(:,:,2),dim=2))
! display frequency spectrum
            call dmscaler1(wk,'POTENTIAL OMEGA VS MODE',ntime,999,2,    &
     &modesxp,cwk,irc)
         endif
      endif
!
! display final spectral analysis for ion current density
      if (movion==1) then
         if (ntji > 0) then
            if ((ndji==2).or.(ndji==3)) then
! performs frequency analysis of accumulated complex time series
               vwkji(1,:,1) = wmi(maxloc(vpkwji(1,:,:,1),dim=2))
               vwkji(2,:,1) = wmi(maxloc(vpkwji(2,:,:,1),dim=2))
               vwkji(1,:,2) = wmi(maxloc(vpkwji(1,:,:,2),dim=2))
               vwkji(2,:,2) = wmi(maxloc(vpkwji(2,:,:,2),dim=2))
! display frequency spectrum
               call dmvector1(vwkji,'ION CURRENT OMEGA VS MODE',ntime,  &
     &999,2,2,modesxji,cwk,irc)
            endif
         endif
      endif
!
! display final spectral analysis for vector potential
      if (nta > 0) then
         if ((nda==2).or.(nda==3)) then
! performs frequency analysis of accumulated complex time series
            vwk(1,:,1) = wm(maxloc(vpkw(1,:,:,1),dim=2))
            vwk(2,:,1) = wm(maxloc(vpkw(2,:,:,1),dim=2))
            vwk(1,:,2) = wm(maxloc(vpkw(1,:,:,2),dim=2))
            vwk(2,:,2) = wm(maxloc(vpkw(2,:,:,2),dim=2))
! display frequency spectrum
            call dmvector1(vwk,'VECTOR POTENTIAL OMEGA VS MODE',ntime,  &
     &999,2,2,modesxa,cwk,irc)
         endif
      endif
!
! display final spectral analysis for transverse efield
      if (ntet > 0) then
         if ((ndet==2).or.(ndet==3)) then
! performs frequency analysis of accumulated complex time series
            vwket(1,:,1) = wm(maxloc(vpkwet(1,:,:,1),dim=2))
            vwket(2,:,1) = wm(maxloc(vpkwet(2,:,:,1),dim=2))
            vwket(1,:,2) = wm(maxloc(vpkwet(1,:,:,2),dim=2))
            vwket(2,:,2) = wm(maxloc(vpkwet(2,:,:,2),dim=2))
! display frequency spectrum
            call dmvector1(vwket,'TRANSVERSE EFIELD OMEGA VS MODE',ntime&
     &,999,2,2,modesxet,cwk,irc)
         endif
      endif
!
! display final spectral analysis for radiative vector potential
      if (ntar > 0) then
         if ((ndar==2).or.(ndar==3)) then
! performs frequency analysis of accumulated complex time series
            vwkr(1,:,1) = wm(maxloc(vpkwr(1,:,:,1),dim=2))
            vwkr(2,:,1) = wm(maxloc(vpkwr(2,:,:,1),dim=2))
            vwkr(1,:,2) = wm(maxloc(vpkwr(1,:,:,2),dim=2))
            vwkr(2,:,2) = wm(maxloc(vpkwr(2,:,:,2),dim=2))
! display frequency spectrum
            call dmvector1(vwkr,'RADIATIVE VPOTENTIAL OMEGA VS MODE',   &
     &ntime,999,2,2,modesxar,cwk,irc)
         endif
      endif
!
      write (iuot,*)
      write (iuot,*) 'initialization time = ', tinit
      write (iuot,*) 'deposit time = ', tdpost
      write (iuot,*) 'current deposit time = ', tdjpost
      tdpost = tdpost + tdjpost
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
!
      ws = 1.0e+09/(real(nloop)*real(np))
      write (iuot,*) 'Push Time (nsec) = ', tpush*ws
      write (iuot,*) 'Deposit Time (nsec) = ', tdpost*ws
      write (iuot,*) 'Sort Time (nsec) = ', tsort*ws
      write (iuot,*) 'Total Particle Time (nsec) = ', time*ws
      write (iuot,*)
!
! reset parameters for final diagnostic metafile
! electron density diagnostic
      if (ntde > 0) then
         nderec = nderec - 1
      endif
! potential diagnostic
      if (ntp > 0) then
         nprec = nprec - 1; ceng = affp
      endif
! longitudinal efield diagnostic
      if (ntel > 0) then
         nelrec = nelrec - 1; ceng = affp
      endif
! radiative vector potential diagnostic
      if (ntar > 0) then
         narrec = narrec - 1; ceng = affp
      endif
! vector potential diagnostic
      if (nta > 0) then
         narec = narec - 1; ceng = affp
      endif
! transverse efield diagnostic
      if (ntet > 0) then
         netrec = netrec - 1; ceng = affp
      endif
! magnetic field diagnostic
      if (ntb > 0) then
         nbrec = nbrec - 1; ceng = affp
      endif
! ion diagnostics
      if (movion==1) then
! ion density diagnostic
         if (ntdi > 0) then
            ndirec = ndirec - 1
         endif
! ion current diagnostic
         if (ntdi > 0) then
            njirec = njirec - 1
         endif
      endif
! write final diagnostic metafile
      call writnml1(iudm)
      write (iuot,*) ' * * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
! close graphics device
      call close_graphs
!
      stop
      end program
