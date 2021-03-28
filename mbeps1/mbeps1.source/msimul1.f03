!-----------------------------------------------------------------------
! High Level library for 1D Electrostatic OpenMP PIC code
!
! subroutines defined:
!
! init_fields1: allocate field data for standard code
! del_fields1: delete field data for standard code
!
! init_electrons1: initialize electrons
! reorder_electrons1: recover from electron buffer overflow errors
! push_electrons1: push electrons with OpenMP:
! del_electrons1: delete electrons
!
! init_ions1: initialize ions
! reorder_ions1: recover from ion buffer overflow errors
! push_ions1: push ions with OpenMP:
! del_ions1: delete ions
!
! es_time_reverse1: start running simulation backwards
!
! init_energy_diag1: initialize energy diagnostic
! energy_diag1: energy diagnostic
! print_energy1: print energy summaries
! del_energy_diag1: delete energy diagnostic
!
! init_spectrum1: allocate scratch arrays for scalar fields
! del_spectrum1: delete scratch arrays for scalar fields
!
! init_edensity_diag1: initialize electron density diagnostic
! edensity_diag1: electron density diagnostic
! del_edensity_diag1: delete electron density diagnostic data
!
! init_idensity_diag1: initialize ion density diagnostic
! idensity_diag1: ion density diagnostic
! del_idensity_diag1: delete ion density diagnostic
!
! init_potential_diag1: initialize potential diagnostic
! potential_diag1: potential diagnostic
! del_potential_diag1: delete potential diagnostic
!
! init_elfield_diag1: initialize longitudinal efield diagnostic
! elfield_diag1: longitudinal efield diagnostic
! del_elfield_diag1: delete longitudinal efield diagnostic
!
! init_efluidms_diag1: initialize electron fluid moments diagnostic
! efluidms_diag1: electron fluid moments diagnostic
! del_efluidms_diag1: delete electron fluid moments diagnostic
!
! init_ifluidms_diag1: initialize ion fluid moments diagnostic
! ifluidms_diag1: ion fluid moments diagnostic
! del_ifluidms_diag1: delete ion fluid moments diagnostic
!
! init_evelocity_diag1: initialize electron velocity diagnostic
! evelocity_diag1: electron velocity diagnostic
! del_evelocity_diag1: delete electron velocity diagnostic
!
! init_ivelocity_diag1: initialize ion velocity diagnostic
! ivelocity_diag1: ion velocity diagnostic
! del_ivelocity_diag1: delete ion velocity diagnostic
!
! init_traj_diag1: initialize trajectory diagnostic
! traj_diag1: trajectory diagnostic
! del_traj_diag1: delete trajectory diagnostic
!
! init_ephasesp_diag1: initialize electron phase space diagnostic
! ephasesp_diag1: electron phase space diagnostic
! del_ephasesp_diag1: delete electron phase space diagnostic
!
! init_iphasesp_diag1: initialize ion phase space diagnostic
! iphasesp_diag1: ion phase space diagnostic
! del_iphasesp_diag1: delete ion electron phase space diagnostic
!
! print_timings1: print timing summaries
!
! reset_diags1: reset electrostatic diagnostics
! close_diags1: close diagnostics
!
! initialize_diagnostics1: initialize all diagnostics from namelist
!                          input parameters
!
! open_restart1: open reset and restart files
! bwrite_restart1: write out basic restart file for electrostatic code
! bread_restart1: read in basic restart file for electrostatic code
! dwrite_restart1: write out restart diagnostic file for electrostatic
!                  code
! dread_restart1: read in restart diagnostic file for electrostatic code
! close_restart1: close reset and restart files
!
! written by Viktor K. Decyk, UCLA
! copyright 1999-2016, regents of the university of california
! update: february 1, 2021
      module f1
      use in1
      use minit1
      use mpush1
      use msort1
      use mgard1
      use mfft1
      use mfield1
      use mdiag1
      implicit none
! idimp = number of particle coordinates = 2
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 2, ipbc = 1
! wke/wki/we = electron/ion kinetic energies and electric field energy
      real :: wke = 0.0, wki = 0.0, we = 0.0
! plist = (true,false) = list of particles leaving tiles found in push
      logical :: plist = .true.
!
! declare scalars for standard code
      integer :: n
      integer :: np, nx, nxh, nxe, nxeh
      integer :: mx1, nstart, nloop, ntime, isign
      integer :: ntime0 = 0, npi = 0
      real :: qbme, affp, ws
      real :: qbmi, vtxi, vtdxi
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nppmx1, ntmax, npbmx
      integer :: irc = 0
      integer, dimension(2) :: irc2 = 0
!
! declare scalars for diagnostics
      integer :: nmv21, it, iw, mtw, mtp, mtv, mtt
      integer :: itw = 0, itp = 0, itv = 0, itt = 0
      integer :: iwi, mtdi
      integer :: itdi = 0
! default Fortran unit numbers
      integer :: iuin = 8, iurr = 9, iuot = 18, iudm = 19
      integer :: iur = 17, iur0 = 27
      integer :: iude = 10, iup = 11, iuel = 12
      integer :: iufe = 23, iuve = 25, iut = 28, iuse = 29
      integer :: iudi = 20, iufi = 24, iuvi = 26, iusi = 30
      real :: ts, eci, wkt
      character(len=10) :: cdrun
      character(len=32) :: fname
!
! declare and initialize timing data
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0, tdiag = 0.0
      double precision :: dtime
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe/qi = electron/ion charge density with guard cells
      real, dimension(:), allocatable :: qe, qi
! fxe = smoothed electric field with guard cells
      real, dimension(:), allocatable :: fxe
! ffc = form factor array for poisson solver
      complex, dimension(:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
!
! declare arrays for OpenMP (tiled) code:
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppbuff
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), allocatable :: ncl
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), allocatable :: ihole
!
! ppart/pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), allocatable :: ppart, pparti
! kpic/kipic = number of electrons/ions in each tile
      integer, dimension(:), allocatable :: kpic, kipic
!
! public diagnostic arrays
! wt = energy time history array
      real, dimension(:,:), allocatable :: wt
! sfield = scratch array for scalar field
      real, dimension(:), allocatable :: sfield
! pkw = power spectrum for potential
      real, dimension(:,:,:), allocatable :: pkw
! pkwdi = power spectrum for ion density
      real, dimension(:,:,:), allocatable :: pkwdi
! wk = maximum frequency as a function of k for potential
      real, dimension(:,:), allocatable :: wk
! wkdi = maximum frequency as a function of k for ion density
      real, dimension(:,:), allocatable :: wkdi
! fmse/fmsi = electron/ion fluid moments
      real, dimension(:,:), allocatable :: fmse, fmsi
! fv/fvi = global electron/ion velocity distribution functions
      real, dimension(:,:), allocatable :: fv, fvi
! fe/fei = global electron/ion energy distribution functions
      real, dimension(:,:), allocatable :: fe, fei
! fvm/fvmi = electron/ion vdrift, vth, entropy for global distribution
      real, dimension(:,:), allocatable :: fvm, fvmi
! fvtm/fvtmi = time history of electron/ion vdrift, vth, and entropy
      real, dimension(:,:,:), allocatable :: fvtm, fvtmi
! fvtp = velocity distribution function for test particles
! fvmtp = vdrift, vth, and entropy for test particles
! fetp = energy distribution function for test particles
      real, dimension(:,:), allocatable :: fvtp, fvmtp, fetp
! partd = trajectory time history array
      real, dimension(:,:,:), allocatable :: partd
! fvs/fvsi = global electron/ion phase space distribution functions
      real, dimension(:,:,:), allocatable :: fvs, fvsi
! scratch arrays for spectral analysis
      real, dimension(:), allocatable :: wm, wmi
! cwk = labels for power spectrum display
      character(len=10), dimension(2) :: cwk = (/'   W > 0  ',          &
     &                                           '   W < 0  ' /)
!
! private diagnostic arrays
! s = scratch array for energies
      double precision, dimension(:), allocatable :: s
! scratch array for scalar field
      complex, dimension(:), allocatable :: sfieldc
! denet/denit = store selected fourier modes for electron/ion density
      complex, dimension(:), allocatable :: denet, denit
! pksdi = accumulated complex spectrum for ion density
      double precision, dimension(:,:,:), allocatable :: pksdi
! pott = store selected fourier modes for potential
      complex, dimension(:), allocatable :: pott
! pks = accumulated complex spectrum for potential
      double precision, dimension(:,:,:), allocatable :: pks
! elt = store selected fourier modes for longitudinal efield
      complex, dimension(:), allocatable :: elt
! sfv/sfvi = electron/ion velocity distribution functions in tile
      real, dimension(:,:,:), allocatable :: sfv, sfvi
! iprobt = scratch array 
      integer, dimension(:), allocatable :: iprobt
! partt = particle trajectories tracked
      real, dimension(:,:), allocatable :: partt
!
      save
!
      public :: qe, qi, fxe, ffc, mixup, sct
      public :: part, ppbuff, ncl, ihole
      public :: ppart, pparti, kpic, kipic
      public :: wt, sfield, pkw, pkwdi, wk, wkdi, fmse, fmsi
      public :: fv, fvi, fvm, fvmi, fvtm, fvtmi, fvtp, fvmtp, partd
      public :: fvs, fvsi
      public :: wm, wmi, cwk
      private :: s, sfieldc
      private :: denet, denit, pksdi, pott, pks, elt
      private :: sfv, sfvi, iprobt, partt
!
      contains
!
!-----------------------------------------------------------------------
      subroutine init_fields1()
! allocate field data for standard code
      implicit none
      allocate(qe(nxe),qi(nxe),fxe(nxe))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_fields1()
! delete field data for standard code
      implicit none
      deallocate(qe,qi,fxe,ffc,mixup,sct)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_electrons1()
! initialize electrons
      implicit none
! part = particle array
      allocate(part(idimp,max(np,npi)))
! background electrons
      if (npx > 0) then
! calculates initial electron co-ordinates with various density profiles
         call mfdistr1(part,ampdx,scaledx,shiftdx,1,npx,nx,ipbc,ndprof, &
     &irc)
         if (irc /= 0) stop
! initialize electron velocities
         call wmvdistr1(part,1,vtx,vx0,ci,npx,nvdist,relativity,irc)
         if (irc /= 0) stop
      endif
! beam electrons
      if (npxb > 0) then
         it = npx + 1
! calculates initial electron co-ordinates with various density profiles
         call mfdistr1(part,ampdx,scaledx,shiftdx,it,npxb,nx,ipbc,ndprof&
     &,irc)
         if (irc /= 0) stop
! initialize electron velocities
         call wmvdistr1(part,it,vtdx,vdx,ci,npxb,nvdist,relativity,irc)
         if (irc /= 0) stop
      endif
!
! mark electron beam particles
      if ((nts > 0).and.(ntsc > 0)) then
         call setmbeam1(part,npx,irc)
         if (irc /= 0) stop
      endif
!
! kpic = number of electrons in each tile
      allocate(kpic(mx1))
!
! find number of electrons in each of mx, tiles: updates kpic, nppmx
      call mdblkp2(part,kpic,nppmx,np,mx,irc)
      if (irc /= 0) stop
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
      if (irc /= 0) stop
!
! sanity check for electrons
      call mcheck1(ppart,kpic,nx,mx,irc)
      if (irc /= 0) stop
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine reorder_electrons1(irc2)
! recover from electron wmporder1 errors
! reallocates ihole, ppbuff, and ppart as necessary
      implicit none
! irc2 = error codes, returned only if error occurs, when irc2(1) /= 0
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: nter = 10
      integer :: iter
      iter = 0
      do while (irc2(1) /= 0)
         iter = iter + 1
         if (iter > nter) then
            write (*,*) 'reorder_electrons1: iteration exceeded'
            stop
         endif
! ihole overflow
         if (irc2(1)==1) then
            ntmax = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmax+1,mx1))
            irc2 = 0
            call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,     &
     &.false.,irc2)
! ppbuff overflow
         else if (irc2(1)==2) then
            npbmx = (1.0 + xtras)*irc2(2)
            deallocate(ppbuff)
            allocate(ppbuff(idimp,npbmx,mx1))
! restores initial values of ncl
            call mprsncl1(ncl,tsort)
            irc2 = 0
            call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,     &
     &.true.,irc2)
! ppart overflow
         else if (irc2(1)==3) then
! restores particle coordinates from ppbuff: updates ppart, ncl
            call mprstor1(ppart,ppbuff,ncl,ihole,tsort)
! copy ordered particles to linear array: updates part
            call mpcopyout1(part,ppart,kpic,it,irc)
            deallocate(ppart)
            nppmx0 = (1.0 + xtras)*irc2(2)
            allocate(ppart(idimp,nppmx0,mx1))
! copies unordered particles to ordered array: updates ppart
            call mpcopyin1(part,ppart,kpic,irc)
            irc2 = 0
            call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,     &
     &plist,irc2)
! electrons not in correct tile, try again
         else if (irc2(1)==(-1)) then
            irc2 = 0
            call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,     &
     &.false.,irc2)
         endif
      enddo
!
! sanity check for electrons
      if (monitor > 0) then
         call mcheck1(ppart,kpic,nx,mx,irc)
         if (irc /= 0) stop
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine push_electrons1(ppart,kpic)
! push electrons with OpenMP:
      implicit none
! ppart = tiled electron particle array
      real, dimension(:,:,:), intent(inout) :: ppart
! kpic = number of electrons in each tile
      integer, dimension(:), intent(inout) :: kpic
      wke = 0.0
! updates ppart, wke and possibly ncl, ihole, and irc
      if (mzf==0) then
         call wmpush1(ppart,fxe,kpic,ncl,ihole,qbme,dt,ci,wke,tpush,nx, &
     &mx,ipbc,relativity,plist,irc)
! zero force: updates ppart, wke and possibly ncl, ihole, and irc
      else if (mzf==1) then
         call wmpush1zf(ppart,kpic,ncl,ihole,dt,ci,wke,tpush,nx,mx,ipbc,&
     &relativity,plist,irc)
      endif
!
! reorder electrons by tile with OpenMP:
! updates ppart, ppbuff, kpic, ncl, irc2, and possibly ihole
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
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_electrons1()
! delete electrons
      implicit none
      deallocate(ppart,kpic)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_ions1()
! initialize ions
      implicit none
! part = particle array
! background ions
      if (npxi > 0) then
         call mfdistr1(part,ampdxi,scaledxi,shiftdxi,1,npxi,nx,ipbc,    &
     &ndprofi,irc)
         if (irc /= 0) stop
         call wmvdistr1(part,1,vtxi,vxi0,ci,npxi,nvdist,relativity,irc)
         if (irc /= 0) stop
      endif
! beam ions
      if (npxbi > 0) then
         it = npxi + 1
         call mfdistr1(part,ampdxi,scaledxi,shiftdxi,it,npxbi,nx,ipbc,  &
     &ndprofi,irc)
         if (irc /= 0) stop
         call wmvdistr1(part,it,vtdxi,vdxi,ci,npxbi,nvdist,relativity,  &
     &irc)
         if (irc /= 0) stop
      endif
!
! mark ion beam particles
      if ((nts > 0).and.(ntsc > 0)) then
         call setmbeam1(part,npxi,irc)
         if (irc /= 0) stop
      endif
!
! kipic = number of ions in each tile
      allocate(kipic(mx1))
!
! find number of ions in each of mx, tiles: updates kipic, nppmx
      call mdblkp2(part,kipic,nppmx,npi,mx,irc)
      if (irc /= 0) stop
!
! allocate vector ion data
      nppmx1 = (1.0 + xtras)*nppmx
      allocate(pparti(idimp,nppmx1,mx1))
! copy ordered ion data for OpenMP: updates pparti and kipic
      call mpmovin1(part,pparti,kipic,mx,irc)
      if (irc /= 0) stop
!
! sanity check for ions
      call mcheck1(pparti,kipic,nx,mx,irc)
      if (irc /= 0) stop
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine reorder_ions1(irc2)
! recover from ion wmporder1 errors
! reallocates ihole, ppbuff, and pparti as necessary
      implicit none
! irc2 = error codes, returned only if error occurs, when irc2(1) /= 0
      integer, dimension(2), intent(inout) :: irc2
! local data
      integer :: nter = 10
      integer :: iter
      iter = 0
      do while (irc2(1) /= 0)
         iter = iter + 1
         if (iter > nter) then
            write (*,*) 'reorder_ions1: iteration exceeded'
            stop
         endif
! ihole overflow
         if (irc2(1)==1) then
            ntmax = (1.0 + xtras)*irc2(2)
            deallocate(ihole)
            allocate(ihole(2,ntmax+1,mx1))
            irc2 = 0
            call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,   &
     &.false.,irc2)
! ppbuff overflow
         else if (irc2(1)==2) then
            npbmx = (1.0 + xtras)*irc2(2)
            deallocate(ppbuff)
            allocate(ppbuff(idimp,npbmx,mx1))
! restores initial values of ncl
            call mprsncl1(ncl,tsort)
            irc2 = 0
            call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,   &
     &.true.,irc2)
! pparti overflow
         else if (irc2(1)==3) then
! restores particle coordinates from ppbuff: updates ppart, ncl
            call mprstor1(pparti,ppbuff,ncl,ihole,tsort)
! copy ordered particles to linear array: updates part
            call mpcopyout1(part,pparti,kipic,it,irc)
            deallocate(pparti)
            nppmx1 = (1.0 + xtras)*irc2(2)
            allocate(pparti(idimp,nppmx1,mx1))
! copies unordered particles to ordered array: updates ppart
            call mpcopyin1(part,pparti,kipic,irc)
            irc2 = 0
            call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,   &
     &plist,irc2)
! ions not in correct tile, try again
         else if (irc2(1)==(-1)) then
            irc2 = 0
            call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,   &
     &.false.,irc2)
         endif
      enddo
!
! sanity check for ions
      if (monitor > 0) then
         call mcheck1(pparti,kipic,nx,mx,irc)
         if (irc /= 0) stop
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine push_ions1(pparti,kipic)
! push ions with OpenMP:
      implicit none
! pparti = tiled electron/ion particle arrays
      real, dimension(:,:,:), intent(inout) :: pparti
! kipic = number of electrons/ions in each tile
      integer, dimension(:), intent(inout) :: kipic
      wki = 0.0
! updates pparti, wki and possibly ncl, ihole, and irc
      if (mzf==0) then
         call wmpush1(pparti,fxe,kipic,ncl,ihole,qbmi,dt,ci,wki,tpush,nx&
     &,mx,ipbc,relativity,plist,irc)
! zero force: updates pparti, wki and possibly ncl, ihole, and irc
      else if (mzf==1) then
         call wmpush1zf(pparti,kipic,ncl,ihole,dt,ci,wki,tpush,nx,mx,   &
     &ipbc,relativity,plist,irc)
      endif
      wki = wki*rmass
!
! reorder ions by tile with OpenMP:
! updates pparti, ppbuff, kipic, ncl, irc2, and possibly ihole
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
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_ions1()
! delete ions
      implicit none
      deallocate(pparti,kipic)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine es_time_reverse1()
! start running simulation backwards:
      implicit none
! need to reverse time lag in leap-frog integration scheme
      dt = -dt
      ws = 0.0
      call wmpush1zf(ppart,kpic,ncl,ihole,dt,ci,ws,tpush,nx,mx,ipbc,    &
     &relativity,plist,irc2(1))
      call wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,plist,irc2)
      if (irc2(1) /= 0) stop
      if (movion==1) then
         call wmpush1zf(pparti,kipic,ncl,ihole,dt,ci,ws,tpush,nx,mx,ipbc&
     &,relativity,plist,irc2(1))
         call wmporder1(pparti,ppbuff,kipic,ncl,ihole,tsort,nx,mx,plist,&
     &irc2)
         if (irc2(1) /= 0) stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_energy_diag1()
! initialize energy diagnostic
      implicit none
! wt = energy time history array
      mtw = (nloop - 1)/ntw + 1; itw = 0
      allocate(wt(mtw,4),s(4))
      wt = 0.0; s = 0.0d0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine energy_diag1(wt,ntime,iuot)
! energy diagnostic
      implicit none
! wt = energy time history array
! ntime = current time step
! iuot = output file descriptor
      real, dimension(:,:), intent(inout) :: wt
      integer, intent(in) :: ntime, iuot
      ws = we + wke + wki
      if (ntime==0) s(3) = ws
      if (ndw > 0) then
         write (iuot,*) 'Field, Kinetic and Total Energies:'
         if (movion==0) then
            write (iuot,'(3e14.7)') we, wke, ws
         else
            write (iuot,'(4e14.7)') we, wke, wki, ws
         endif
      endif
      itw = itw + 1
! store energies in time history array
      wt(itw,:) = (/we,wke,wki,ws/)
      s(1) = s(1) + we
      s(2) = s(2) + wke
      s(3) = min(s(3),dble(ws))
      s(4) = max(s(4),dble(ws))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine print_energy1(wt,iuot)
! print energy summaries
      implicit none
! wt = energy time history array
! iuot = output file descriptor
      real, dimension(:,:), intent(in) :: wt
! iuot = output file descriptor
      integer, intent(in) :: iuot
      s(3) = (s(4) - s(3))/wt(1,4)
      write (iuot,*) 'Energy Conservation = ', real(s(3))
      s(1) = s(1)/real(itw)
      write (iuot,*) 'Average Field Energy <WE> = ', real(s(1))
      s(2) = s(2)/real(itw)
      write (iuot,*) 'Average Electron Kinetic Energy <WKE> = ',        &
     &real(s(2))
      write (iuot,*) 'Ratio <WE>/<WKE>= ', real(s(1)/s(2))
      write (iuot,*)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_energy_diag1()
! delete energy diagnostic
      implicit none
! wt = energy time history array
      deallocate(wt,s)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_spectrum1()
! allocate scratch arrays for scalar fields
      implicit none
      allocate(sfield(nxe),sfieldc(nxh))
! allocate and initialize frequency array for spectral analysis
      if (ntp > 0) then
         iw = (wmax - wmin)/dw + 1.5
         allocate(wm(iw))
         do it = 1, iw
            wm(it) = wmin + dw*real(it-1)
         enddo
      endif
! allocate and initialize frequency array for ion spectral analysis
      if (movion==1) then
         if (ntdi > 0) then
            iwi = (wimax - wimin)/dwi + 1.5
            allocate(wmi(iwi))
            do it = 1, iwi
               wmi(it) = wimin + dwi*real(it-1)
            enddo
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_spectrum1()
! delete scratch arrays for scalar fields
      implicit none
      deallocate(sfield,sfieldc)
! deallocate frequency array for spectral analysis
      if (ntp > 0) deallocate(wm)
! deallocate frequency array for ion spectral analysis
      if (movion==1) then
         if (ntdi > 0) deallocate(wmi)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_edensity_diag1()
! initialize electron density diagnostic
      implicit none
      fdename = 'denek1.'//cdrun
      modesxde = min(modesxde,nxh+1)
      allocate(denet(modesxde))
! open file: updates nderec and possibly iude
      if (nderec==0) call dafopenc1(denet,iude,nderec,trim(fdename))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine edensity_diag1(sfield)
! electron density diagnostic
      implicit none
! sfield = scratch array for scalar field
      real, dimension(:), intent(inout) :: sfield
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
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_edensity_diag1()
! delete electron density diagnostic data
      implicit none
      if (nderec > 0) then
         close(unit=iude)
         nderec = nderec - 1
      endif
      deallocate(denet)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_idensity_diag1()
! initialize ion density diagnostic
      implicit none
      fdiname = 'denik1.'//cdrun
      modesxdi = min(modesxdi,nxh+1)
      allocate(denit(modesxdi))
! open file: updates ndirec and possibly iudi
      if (ndirec==0) call dafopenc1(denit,iudi,ndirec,trim(fdiname))
! ion spectral analysis
      if ((nddi==2).or.(nddi==3)) then
         mtdi = (nloop - 1)/ntdi + 1; itdi = 0
         allocate(pkwdi(modesxdi,iwi,2),wkdi(modesxdi,2))
         allocate(pksdi(4,modesxdi,iwi))
         pksdi = 0.0d0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine idensity_diag1(sfield,pkwdi,wkdi,ntime)
! ion density diagnostic
      implicit none
! sfield = scratch array for scalar field
      real, dimension(:), intent(inout) :: sfield
! pkwdi = power spectrum for ion density
      real, dimension(:,:,:), intent(inout) :: pkwdi
! wkdi = maximum frequency as a function of k for ion density
      real, dimension(:,:), intent(inout) :: wkdi
! ntime = current time step
      integer, intent(in) :: ntime
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
      endif
! ion spectral analysis
      if ((nddi==2).or.(nddi==3)) then
         itdi = itdi + 1
         ts = dt*real(ntime)
! performs frequency analysis of accumulated complex time series
! zero out mode 0
         denit(1) = cmplx(0.0,0.0)
         call micspect1(denit,wmi,pkwdi,pksdi,ts,t0,tdiag,mtdi,iwi,     &
     &modesxdi,nx,-1)
! find frequency with maximum power for each mode
         wkdi(:,1) = wmi(maxloc(pkwdi(:,:,1),dim=2))
         wkdi(:,2) = wmi(maxloc(pkwdi(:,:,2),dim=2))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_idensity_diag1()
! delete ion density diagnostic
      implicit none
      if (ndirec > 0) then
         close(unit=iudi)
         ndirec = ndirec - 1
      endif
      deallocate(denit)
! spectral analysis
      if ((nddi==2).or.(nddi==3)) then
         deallocate(pkwdi,wkdi,pksdi)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_potential_diag1()
! initialize potential diagnostic
      implicit none
      fpname = 'potk1.'//cdrun
      modesxp = min(modesxp,nxh+1)
      allocate(pott(modesxp))
! open file: updates nprec and possibly iup
      if (nprec==0) call dafopenc1(pott,iup,nprec,trim(fpname))
! spectral analysis
      if ((ndp==2).or.(ndp==3)) then
         mtp = (nloop - 1)/ntp + 1; itp = 0
         allocate(pkw(modesxp,iw,2),wk(modesxp,2))
         allocate(pks(4,modesxp,iw))
         pks = 0.0d0
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine potential_diag1(sfield,pkw,wk,ntime)
! potential diagnostic
      implicit none
! sfield = scratch array for scalar field
      real, dimension(:), intent(inout) :: sfield
! pkw = power spectrum for potential
      real, dimension(:,:,:), intent(inout) :: pkw
! wk = maximum frequency as a function of k for potential
      real, dimension(:,:), intent(inout) :: wk
! ntime = current time step
      integer, intent(in) :: ntime
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
      endif
! spectral analysis
      if ((ndp==2).or.(ndp==3)) then
         itp = itp + 1
         ts = dt*real(ntime)
! performs frequency analysis of accumulated complex time series
         call micspect1(pott,wm,pkw,pks,ts,t0,tdiag,mtp,iw,modesxp,nx,1)
! find frequency with maximum power for each mode
         wk(:,1) = wm(maxloc(pkw(:,:,1),dim=2))
         wk(:,2) = wm(maxloc(pkw(:,:,2),dim=2))
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_potential_diag1()
! delete potential diagnostic
      implicit none
      if (nprec > 0) then
         close(unit=iup)
         nprec = nprec - 1
      endif
      deallocate(pott)
! spectral analysis
      if ((ndp==2).or.(ndp==3)) then
         deallocate(pkw,wk,pks)
      endif
      ceng = affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_elfield_diag1()
! initialize longitudinal efield diagnostic
      implicit none
      felname = 'elk1.'//cdrun
      modesxel = min(modesxel,nxh+1)
      allocate(elt(modesxel))
! open file: updates nelrec and possibly iuel
      if (nelrec==0) call dafopenc1(elt,iuel,nelrec,trim(felname))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine elfield_diag1(sfield)
! longitudinal efield diagnostic
      implicit none
! sfield = scratch array for scalar field
      real, dimension(:), intent(inout) :: sfield
! calculate longitudinal efield in fourier space: updates sfieldc
      call melfield1(qe,sfieldc,ffc,ws,tfield,nx)
! store selected fourier modes: updates elt
      call mrdmodes1(sfieldc,elt,tfield,nx,modesxel)
! write diagnostic output: updates nelrec
      call dafwritec1(elt,tdiag,iuel,nelrec,modesxel)
! transform longitudinal efield to real space: updates sfield
      call mfft1cr(sfieldc,sfield,mixup,sct,tfft,indx)
      call mdguard1(sfield,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_elfield_diag1()
! delete longitudinal efield diagnostic
      implicit none
      if (nelrec > 0) then
         close(unit=iuel)
         nelrec = nelrec - 1
      endif
      deallocate(elt)
      ceng = affp
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_efluidms_diag1()
! initialize electron fluid moments diagnostic
      implicit none
! calculate first dimension of fluid arrays
      if (npro==1) then
         nprd = 1
      else if (npro==2) then
         nprd = 2
      else if (npro==3) then
         nprd = 3
      else if (npro==4) then
         nprd = 5
      endif
      if ((ndfm==1).or.(ndfm==3)) then
         allocate(fmse(nprd,nxe))
! open file for real data: updates nferec and possibly iufe
         ffename = 'fmer1.'//cdrun
         if (nferec==0) then
            call dafopenv1(fmse,nx,iufe,nferec,trim(ffename))
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine efluidms_diag1(fmse)
! electron fluid moments diagnostic
      implicit none
! fmse = electron fluid moments
      real, dimension(:,:), intent(inout) :: fmse
! calculate electron fluid moments
      if ((ndfm==1).or.(ndfm==3)) then
         call dtimer(dtime,itime,-1)
         fmse = 0.0
         call dtimer(dtime,itime,1)
         tdiag = tdiag + real(dtime)
         call wmgprofx1(ppart,fxe,fmse,kpic,qbme,dt,ci,tdiag,npro,nx,mx,&
     &relativity)
! add guard cells with OpenMP: updates fmse
         call mamcguard1(fmse,tdiag,nx)
! calculates fluid quantities from fluid moments: updates fmse
         call mfluidqs1(fmse,tdiag,npro,nx)
! write real space diagnostic output: updates nferec
          call dafwritev1(fmse,tdiag,iufe,nferec,nx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_efluidms_diag1()
! delete electron fluid moments diagnostic
      implicit none
      if (nferec > 0) then
         close(unit=iufe)
         nferec = nferec - 1
      endif
      if (allocated(fmse)) deallocate(fmse)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_ifluidms_diag1()
! initialize ion fluid moments diagnostic
      implicit none
! calculate first dimension of fluid arrays
      if (npro==1) then
         nprd = 1
      else if (npro==2) then
         nprd = 2
      else if (npro==3) then
         nprd = 3
      else if (npro==4) then
         nprd = 5
      endif
      if ((ndfm==2).or.(ndfm==3)) then
         allocate(fmsi(nprd,nxe))
! open file for real data: updates nfirec and possibly iufi
         ffiname = 'fmir1.'//cdrun
         if (nfirec==0) then
            call dafopenv1(fmsi,nx,iufi,nfirec,trim(ffiname))
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ifluidms_diag1(fmsi)
! ion fluid moments diagnostic
      implicit none
! fmsi = ion fluid moments
      real, dimension(:,:), intent(inout) :: fmsi
! calculate ion fluid moments
      if ((ndfm==2).or.(ndfm==3)) then
         call dtimer(dtime,itime,-1)
         fmsi = 0.0
         call dtimer(dtime,itime,1)
         tdiag = tdiag + real(dtime)
         call wmgprofx1(pparti,fxe,fmsi,kipic,qbmi,dt,ci,tdiag,npro,nx, &
     &mx,relativity)
! add guard cells with OpenMP: updates fmsi
         call mamcguard1(fmsi,tdiag,nx)
! calculates fluid quantities from fluid moments: updates fmsi
         call mfluidqs1(fmsi,tdiag,npro,nx)
         fmsi = rmass*fmsi
! write real space diagnostic output: updates nfirec
         call dafwritev1(fmsi,tdiag,iufi,nfirec,nx)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_ifluidms_diag1()
! delete ion fluid moments diagnostic
      implicit none
      if (nfirec > 0) then
         close(unit=iufi)
         nfirec = nfirec - 1
      endif
      if (allocated(fmsi)) deallocate(fmsi)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_evelocity_diag1()
! initialize electron velocity diagnostic
      implicit none
      nfvd = 0; nfed = 0
      if ((nvft==1).or.(nvft==3)) then
         nfvd = ndim
      endif
      if ((nvft==2).or.(nvft==3)) then
         nfed = 1
      endif
      mtv = (nloop - 1)/ntv + 1; itv = 0
      eci = ci; if (relativity==0) eci = 0.0
      wkt = 0.0
      if ((ndv==1).or.(ndv==3)) then
! estimate maximum electron velocity or momentum
         ws = 0.0
         if (npx > 0) then
            ws = 4.0*vtx+abs(vx0)
         endif
         if (npxb > 0) then
            ws = max(ws,4.0*vtdx+abs(vdx))
         endif
! cylindrical not supported in electrostatic code
         if (nvft < 4) then
            allocate(fv(2*nmv+2,nfvd),fvm(ndim,3),fe(2*nmv+2,nfed))
            allocate(sfv(2*nmv+2,ndim,mx1+1))
            fvm = 0.0
! open file for electron velocity data: updates nverec and possibly iuve
            fvename = 'fve1.'//cdrun
            if (nverec==0) then
               call dafopenfv1(fvm,fv,fe,wkt,iuve,nverec,trim(fvename))
            endif
         endif
! cartesian distribution
         if ((nvft==1).or.(nvft==3)) then
            allocate(fvtm(mtv,ndim,3))
            fvtm = 0.0
! set velocity or momentum scale
            fv(2*nmv+2,1) = 2.0*ws
         endif
! energy distribution
         if ((nvft==2).or.(nvft==3)) then
! set energy scale for electrons
            ws = ws*ws
            fe(2*nmv+2,1) = ws/(1.0 + sqrt(1.0 + ws*eci*eci))
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine evelocity_diag1(fv,fe,fvm,fvtm,wkt)
! electron velocity diagnostic
      implicit none
! fv = global electron velocity distribution function
      real, dimension(:,:), intent(inout) :: fv
! fe = global electron energy distribution function
      real, dimension(:,:), intent(inout) :: fe
! fvm = electron vdrift, vth, entropy for global distribution
      real, dimension(:,:), intent(inout)  :: fvm
! fvtm = time history of electron vdrift, vth, and entropy
      real, dimension(:,:,:), intent(inout) :: fvtm
! wkt = total energy contained in distribution
      real, intent(inout) :: wkt
      if ((ndv==1).or.(ndv==3)) then
! calculate electron cartesian distribution function and moments
         if ((nvft==1).or.(nvft==3)) then
            sfv(2*nmv+2,:,mx1+1) = fv(2*nmv+2,:)
            call mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,nmv)
            fv = sfv(:,:,mx1+1)
! store time history electron vdrift, vth, and entropy
            itv = itv + 1
            fvtm(itv,:,:) = fvm
         endif
! electron energy distribution
         if ((nvft==2).or.(nvft==3)) then
            sfv(2*nmv+2,1,mx1+1) = fe(2*nmv+2,1)
            call merpdist1(ppart,kpic,sfv,eci,wkt,tdiag,nmv)
            fe(:,1) = sfv(:,1,mx1+1)
         endif
! cylindrical not supported in electrostatic code
         if (nvft < 4) then
! write electron velocity-space diagnostic output: updates nverec
            call dafwritefv1(fvm,fv,fe,wkt,tdiag,iuve,nverec)
         endif
      endif   
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_evelocity_diag1()
! delete electron velocity diagnostic
      implicit none
      if (nverec > 0) then
         close(unit=iuve)
         nverec = nverec - 1
      endif
      if ((ndv==1).or.(ndv==3)) deallocate(fv,fvm,fe)
      if (allocated(fvtm)) deallocate(fvtm)
      if (allocated(sfv)) deallocate(sfv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_ivelocity_diag1()
! initialize ion velocity diagnostic
      implicit none
      nfvd = 0; nfed = 0
      if ((nvft==1).or.(nvft==3)) then
         nfvd = ndim
      endif
      if ((nvft==2).or.(nvft==3)) then
         nfed = 1
      endif
      mtv = (nloop - 1)/ntv + 1; itv = 0
      eci = ci; if (relativity==0) eci = 0.0
      wkt = 0.0
      if ((ndv==2).or.(ndv==3)) then
! estimate maximum ion velocity or momentum
         ws = 0.0
         if (npxi > 0) then
            ws = 4.0*vtxi+abs(vxi0)
         endif
         if (npxbi > 0) then
            ws = max(ws,4.0*vtdxi+abs(vdxi))
         endif
! cylindrical not supported in electrostatic code
         if (nvft < 4) then
            allocate(fvi(2*nmv+2,nfvd),fvmi(ndim,3),fei(2*nmv+2,nfed))
            allocate(sfvi(2*nmv+2,ndim,mx1+1))
            fvmi = 0.0
! open file for ion velocity data: updates nvirec and possibly iuvi
            fviname = 'fvi1.'//cdrun
            if (nvirec==0) then
               call dafopenfv1(fvmi,fvi,fei,wkt,iuvi,nvirec,            &
     &trim(fviname))
            endif
         endif
! cartesian distribution
         if ((nvft==1).or.(nvft==3)) then
            allocate(fvtmi(mtv,ndim,3))
            fvtmi = 0.0
! set velocity or momentum scale
            fvi(2*nmv+2,1) = 2.0*ws
         endif
! energy distribution
         if ((nvft==2).or.(nvft==3)) then
! set energy scale for ions
            ws = ws*ws
            fei(2*nmv+2,1) = ws/(1.0 + sqrt(1.0 + ws*eci*eci))
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ivelocity_diag1(fvi,fei,fvmi,fvtmi,wkt)
! ion velocity diagnostic
      implicit none
! fvi = global ion velocity distribution function
      real, dimension(:,:), intent(inout) :: fvi
! fei = global ion energy distribution function
      real, dimension(:,:), intent(inout) :: fei
! fvmi = ion vdrift, vth, entropy for global distribution
      real, dimension(:,:), intent(inout)  :: fvmi
! fvtmi = time history of ion vdrift, vth, and entropy
      real, dimension(:,:,:), intent(inout) :: fvtmi
! wkt = total energy contained in distribution
      real, intent(inout) :: wkt
      if ((ndv==2).or.(ndv==3)) then
! calculate ion cartesian distribution function and moments
         if ((nvft==1).or.(nvft==3)) then
            sfvi(2*nmv+2,:,mx1+1) = fvi(2*nmv+2,:)
            call mvpdist1(pparti,kipic,sfvi,fvmi,tdiag,npi,nmv)
            fvi = sfvi(:,:,mx1+1)
! update time step if electrons have not been calculated
            if (ndv==2) itv = itv + 1
! store time history of ion vdrift, vth, and entropy
            fvtmi(itv,:,:) = fvmi
         endif
! ion energy distribution
         if ((nvft==2).or.(nvft==3)) then
            sfvi(2*nmv+2,1,mx1+1) = fei(2*nmv+2,1)
            call merpdist1(pparti,kipic,sfvi,eci,wkt,tdiag,nmv)
            fei(:,1) = sfvi(:,1,mx1+1)
            wkt = rmass*wkt
         endif
! cylindrical not supported in electrostatic code
         if (nvft < 4) then
! write ion velocity-space diagnostic output: updates nvirec
            call dafwritefv1(fvmi,fvi,fei,wkt,tdiag,iuvi,nvirec)
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_ivelocity_diag1()
! delete ion velocity diagnostic
      implicit none
      if (nvirec > 0) then
         close(unit=iuvi)
         nvirec = nvirec - 1
      endif
      if ((ndv==2).or.(ndv==3)) deallocate(fvi,fvmi,fei)
      if (allocated(fvtmi)) deallocate(fvtmi)
      if (allocated(sfvi)) deallocate(sfvi)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_traj_diag1(ntime)
! initialize trajectory diagnostic
      implicit none
! ntime = current time step
      integer, intent(in) :: ntime
! set initial test trajectories
      if ((ntime+ntime0)==0) then
         if ((ndt==2).and.(movion==0)) ndt = 0
         if ((ndt==1).or.(ndt==2)) then
            allocate(iprobt(nprobt))
         endif
! electron trajectories
         if (ndt==1) then
! sets electron test charge distribution: updates ppart, iprobt, nprobt
            call setptraj1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,np,nprobt&
     &,irc)
            if (irc /= 0) then
               write (*,*) 'esetptraj1 error: irc=', irc
               stop
            endif
! estimate maximum electron velocity or momentum
            if (nst==3) then
               ws = 0.0
               if (npx > 0) then
                  ws = 4.0*vtx+abs(vx0)
               endif
               if (npxb > 0) then
                  ws = max(ws,4.0*vtdx+abs(vdx))
               endif
            endif
! ion trajectories
         else if (ndt==2) then
            if (movion==1) then
! sets ion test charge distribution: updates ppart, iprobt, nprobt
               call setptraj1(pparti,kipic,iprobt,nst,vtxi,vtsx,dvtx,npi&
     &,nprobt,irc)
               if (irc /= 0) then
                  write (*,*) 'isetptraj1 error: irc=', irc
                  stop
               endif
! estimate maximum ion velocity or momentum
               if (nst==3) then
                  ws = 0.0
                  if (npxi > 0) then
                     ws = 4.0*vtxi+abs(vxi0)
                  endif
                  if (npxbi > 0) then
                     ws = max(ws,4.0*vtdxi+abs(vdxi))
                  endif
               endif
            endif
         endif
         if (allocated(iprobt)) deallocate(iprobt)
! find number of existing test trajectories: updates nprobt
      else
         if (ndt==1) then
            call mfnptraj1(ppart,kpic,nprobt,irc)
            if (nst==3) then
               ws = 0.0
               if (npx > 0) then
                  ws = 4.0*vtx+abs(vx0)
               endif
               if (npxb > 0) then
                  ws = max(ws,4.0*vtdx+abs(vdx))
               endif
            endif
         else if (ndt==2) then
            call mfnptraj1(pparti,kipic,nprobt,irc)
            if (nst==3) then
               ws = 0.0
               if (npxi > 0) then
                   ws = 4.0*vtxi+abs(vxi0)
               endif
               if (npxbi > 0) then
                  ws = max(ws,4.0*vtdxi+abs(vdxi))
               endif
            endif
         endif
         if (irc /= 0) then
            write (*,*) 'mfnptraj1 error: irc=', irc
            stop
         endif
      endif
! electron or ion trajectories
      if ((ndt==1).or.(ndt==2)) then
         if (nprobt > 16777215) then
               write(*,*) 'nprobt overflow = ', nprobt
               stop
         endif
         ndimp = idimp
         allocate(partt(idimp,nprobt))
         ftname = 'tr1.'//cdrun
! track particle trajectories
         if ((nst==1).or.(nst==2)) then
            mtt = (nloop - 1)/ntt + 1
            itt = 0
            allocate(partd(mtt,idimp,nprobt))
            partd = 0.0
! open file for trajectory data: updates ntrec and possibly iut
            if (ntrec==0) then
               call dafopentr1(partt,iut,ntrec,trim(ftname))
            endif
! calculate test particle distribution function and moments
         else if (nst==3) then
            allocate(fvtp(2*nmv+2,ndim),fvmtp(ndim,3))
            allocate(fetp(2*nmv+2,0))
            fvtp(2*nmv+2,:) = 2.0*ws
! open file for test particle diagnostic: updates ntrec and possibly iut
            if (ntrec==0) then
               ws = 0.0
               call dafopenfv1(fvmtp,fvtp,fetp,ws,iut,ntrec,            &
     &trim(ftname))
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine traj_diag1(partd,fvtp,fvmtp)
! trajectory diagnostic
      implicit none
! partd = trajectory time history array
      real, dimension(:,:,:), intent(inout) :: partd
! fvtp = velocity distribution function for test particles
! fvmtp = vdrift, vth, and entropy for test particles
      real, dimension(:,:), intent(inout) :: fvtp, fvmtp
! copies tagged electron coordinates to array partt
      if (ndt==1) then
         call mptraj1(ppart,kpic,partt,tdiag,irc)
! copies tagged ion coordinatess to array partti
      else if (ndt==2) then
         if (movion==1) then
            call mptraj1(pparti,kipic,partt,tdiag,irc)
         endif
      endif
      if (irc /= 0) stop
! electron or ion trajectories
      if ((ndt==1).or.(ndt==2)) then
! store particle trajectories
         if ((nst==1).or.(nst==2)) then
! write trajectory diagnostic output: updates ntrec
            call dafwritetr1(partt,tdiag,iut,ntrec)
            itt = itt + 1
            partd(itt,:,:) = partt
! calculate test particle distribution function and moments
         else if (nst==3) then
! calculate test particle distribution function and moments
            call mvdist1(partt,fvtp,fvmtp,tdiag,nprobt,nmv)
! write test particle diagnostic output: updates ntrec
             ws = 0.0
             call dafwritefv1(fvmtp,fvtp,fetp,ws,tdiag,iut,ntrec)
         endif
      endif 
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_traj_diag1()
! delete trajectory diagnostic
      implicit none
      if (ntrec > 0) then
         close(unit=iut)
         ntrec = ntrec - 1
      endif
      if (allocated(partt)) deallocate(partt)
      if ((nst==1).or.(nst==2)) then
         deallocate(partd)
      else if (nst==3) then
         deallocate(fvtp,fvmtp,fetp)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_ephasesp_diag1()
! initialize electron phase space diagnostic
      implicit none
      nmv21 = 2*nmv + 1
      mvx = min(mvx,nx)
      nsxb = (nx - 1)/mvx + 1
! electron phase space diagnostic
      if ((nds==1).or.(nds==3)) then
! estimate maximum electron velocity or momentum
         ws = 0.0
         if (npx > 0) then
            ws = 4.0*vtx+abs(vx0)
         endif
         if (npxb > 0) then
            ws = max(ws,4.0*vtdx+abs(vdx))
         endif
         allocate(fvs(nmv21+1,ndim,nsxb))
         fvs = 0.0
         fvs(nmv21+1,:,1) = 1.25*ws
! open file for electron phase space data:
! updates nserec and possibly iuse
! opens a new fortran unformatted stream file
         if (nserec==0) then
            fsename = 'pse1.'//cdrun
            iuse =  get_funit(iuse)
            call fnopens1(iuse,trim(fsename))
            nserec = 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ephasesp_diag1(fvs)
! electron phase space diagnostic
! fvs = global electron phase space distribution function
      real, dimension(:,:,:), intent(inout) :: fvs
      if ((nds==1).or.(nds==3)) then
! calculates velocity distribution in different regions of space:
! updates fvs
         call mvspdist1(ppart,kpic,fvs,tdiag,nmv,mvx)
! write phase space diagnostic output: updates nserec
         if (nserec > 0) then
            call mwrfvsdata1(fvs,tdiag,iuse)
            nserec = nserec + 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_ephasesp_diag1()
! delete electron phase space diagnostic
      implicit none
      if (nserec > 0) then
         close(unit=iuse)
         nserec = nserec - 1
      endif
      if (allocated(fvs)) deallocate(fvs)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine init_iphasesp_diag1()
! initialize ion phase space diagnostic
      implicit none
      nmv21 = 2*nmv + 1
      mvx = min(mvx,nx)
      nsxb = (nx - 1)/mvx + 1
      if ((nds==2).or.(nds==3)) then
! estimate maximum ion velocity or momentum
         ws = 0.0
         if (npxi > 0) then
            ws = 4.0*vtxi+abs(vxi0)
         endif
         if (npxbi > 0) then
            ws = max(ws,4.0*vtdxi+abs(vdxi))
         endif
         allocate(fvsi(nmv21+1,ndim,nsxb))
         fvsi = 0.0
         fvsi(nmv21+1,:,1) = 1.25*ws
! open file for ion phase space data:
! updates nsirec and possibly iusi
! opens a new fortran unformatted stream file
         if (nsirec==0) then
            fsiname = 'psi1.'//cdrun
            iusi =  get_funit(iusi)
            call fnopens1(iusi,trim(fsiname))
            nsirec = 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine iphasesp_diag1(fvsi)
! ion phase space diagnostic
! fvsi = global ion phase space distribution function
      real, dimension(:,:,:), intent(inout) :: fvsi
      if ((nds==2).or.(nds==3)) then
! calculates velocity distribution in different regions of space:
! updates fvsi
         call mvspdist1(pparti,kipic,fvsi,tdiag,nmv,mvx)
! write phase space diagnostic output: updates nsirec
         if (nsirec > 0) then
            call mwrfvsdata1(fvsi,tdiag,iusi)
            nsirec = nsirec + 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine del_iphasesp_diag1()
! delete ion phase space diagnostic
      implicit none
      if (nsirec > 0) then
         close(unit=iusi)
         nsirec = nsirec - 1
      endif
      if (allocated(fvsi)) deallocate(fvsi)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine print_timings1(tinit,tloop,iuot)
! print timing summaries
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
      subroutine reset_diags1()
! reset electrostatic diagnostics
      implicit none
      if (ntw > 0) then
         itw = 0; wt = 0.0
         if (allocated(s)) s = 0.0d0
      endif
      if (ntde > 0) then
         if (nderec > 1) nderec = 1
      endif
      if (ntp > 0) then
         if (nprec > 1) nprec = 1
         if ((ndp==2).or.(ndp==3)) then
            itp = 0; pks = 0.0d0
         endif
      endif
      if (movion==1) then
         if (ntdi > 0) then
            if (ndirec > 1) ndirec = 1
            if ((nddi==2).or.(nddi==3)) then
               itdi = 0; pksdi = 0.0d0
            endif
         endif
      endif
      if (ntel > 0) then
         if (nelrec > 1) nelrec = 1
      endif
      if (ntfm > 0) then
         if (nferec > 1) nferec = 1
         if (movion==1) then
            if (nfirec > 1) nfirec = 1
         endif
      endif
      if (ntv > 0) then
         if (nverec > 1) nverec = 1
         itv = 0; fvtm = 0.0
         if (movion==1) then
            if (nvirec > 1) nvirec = 1
            fvtmi = 0.0
         endif
      endif
      if (ntt > 0) then
         if (ntrec > 1) ntrec = 1
         itt = 0
         if (allocated(partd)) partd = 0.0
      endif
      if (nts > 0) then
         if (nserec > 1) nserec = 1
         if (movion==1) then
            if (nsirec > 1) nsirec = 1
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine close_diags1(iudm)
! close diagnostics
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
! velocity diagnostic
      if (ntv > 0) then
         call del_evelocity_diag1()
         if (movion==1) call del_ivelocity_diag1()
      endif
! trajectory diagnostic
      if (ntt > 0) call del_traj_diag1()
! phase space diagnostic
      if (nts > 0) then
         call del_ephasesp_diag1()
         if (movion==1) call del_iphasesp_diag1()
      endif
! ion density diagnostic
      if (movion==1) then
         if (ntdi > 0) call del_idensity_diag1()
      endif
! write final diagnostic metafile
      call writnml1(iudm)
! deallocate arrays
      call del_fields1()
      call del_electrons1()
      if (movion==1) call del_ions1()
      deallocate(part,ppbuff,ncl,ihole)
      if ((ntde > 0).or.(ntp > 0).or.(ntel > 0).or.(ntdi > 0)) then
         call del_spectrum1()
      endif
      if (ntw > 0) call del_energy_diag1()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine initialize_diagnostics1(ntime)
! initialize all diagnostics from namelist input parameters
      implicit none
! ntime = current time step
      integer, intent(in) :: ntime
! initialize energy diagnostic: updates wt
      if (ntw > 0) then
         call init_energy_diag1()
      endif
!
! allocate and initialize scratch arrays for scalar fields:
! allocates sfield
      if ((ntde > 0).or.(ntp > 0).or.(ntel > 0).or.(ntdi > 0)) then
         call init_spectrum1()
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
! initialize fluid moments diagnostic
      if (ntfm > 0) then
! electrons: allocates fmse
         call init_efluidms_diag1()
! ions: allocates fmsi
         if (movion==1) call init_ifluidms_diag1()
      endif
!
! initialize velocity diagnostic:
      if (ntv > 0) then
! electrons: allocates fv, fvm, fvtm
         call init_evelocity_diag1()
! ions: allocates fvi, fvmi, fvtmi
         if (movion==1) call init_ivelocity_diag1()
      endif
!
! initialize trajectory diagnostic: allocates partd, fvtp, fvmtp, fetp
      if (ntt > 0) then
         call init_traj_diag1(ntime)
      endif
!
! initialize phase space diagnostic:
      if (nts > 0) then
! electrons: allocates fvs
         call init_ephasesp_diag1()
! ions: allocates fvsi
         if (movion==1) call init_iphasesp_diag1()
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine open_restart1()
! open reset and restart files
      implicit none
! iur, iurr = restart, reset, old restart file descriptors
! local data
      character(len=10) :: cdrun0
      character(len=32) :: fname
!     character(len=10) :: a = 'sequential'
      character(len=6) :: a = 'stream'
      character(len=11) :: f = 'unformatted'
! reset file
!     fname = 'reset1'
!     iurr = get_funit(iurr)
!     open(iurr,file=trim(fname),access=a,form=f,status='unknown')
! restart file
      fname = 'rstrt1.'//cdrun
      iur = get_funit(iur)
! start a new run from random numbers
      if (nustrt==1) then
         if (ntr > 0) then
            open(iur,file=trim(fname),access=a,form=f,status= 'unknown')
         endif
! continue a run which was interrupted
      else if (nustrt==2) then
         open(iur,file=trim(fname),access=a,form=f,status='old')
! start a new run with data from a previous run
      else if (nustrt==0) then
         if (ntr > 0) then
            open(iur,file=trim(fname),access=a,form=f,status= 'unknown')
         endif
         if (idrun /= idrun0) then
            write (cdrun0,'(i10)') idrun0
            cdrun0 = adjustl(cdrun0)
            fname = 'rstrt1.'//cdrun0
            open(iur0,file=trim(fname),access=a,form=f,status='old')
         else
            write (*,*) 'restart warning: old, new idruns identical'
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bwrite_restart1(iur,ntime)
! write out basic restart file for electrostatic code
      implicit none
! iur = restart file descriptor
! ntime = current time step
      integer, intent(in) :: iur, ntime
! local data
      integer :: i, j
      integer :: idimp, npp, nxv
      idimp = size(part,1)
      rewind iur
! write out current and initial time
      write (iur) ntime, ntime0
! copy ordered particles to linear array: updates part, npp
      call mpcopyout1(part,ppart,kpic,npp,irc)
! write out size of electron array
      write (iur) npp, idimp
! write out electrons, if non-zero
      if (npp > 0) then
         write (iur) ((part(i,j),i=1,idimp),j=1,npp)
      endif
! write out if ions are moving
      write (iur) movion
      if (movion > 0) then
! copy ordered particles to linear array: updates part, npp
         call mpcopyout1(part,pparti,kipic,npp,irc)
! write out size of ion array
         write (iur) npp, idimp
! write out ions, if non-zero
         if (npp > 0) then
            write (iur) ((part(i,j),i=1,idimp),j=1,npp)
         endif
! write out ion density, if ions are not moving
      else
         nxv = size(qi)
         write (iur) nxv
         if (nxv > 0) then
            write (iur) (qi(j),j=1,nxv)
         endif
      endif
! write out electric field parameter
      write (iur) emf
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine bread_restart1(iur)
! read in basic restart file for electrostatic code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j
      integer :: ndimp, npp, ios, it
      integer :: nppmx, nppmx0, ntmax, npbmx, nppmx1
! part = linear particle array
      if (.not.allocated(part)) allocate(part(idimp,max(np,npi)))
!
! read in current and initial time
      rewind iur
      read (iur,iostat=ios) ntime, ntime0
      if (ios /= 0) then
         write (*,*) 'ntime restart error, ios = ', ios
         stop
      endif
      if (nustrt==0) then
         write (*,*) 'restarting from ntime, idrun0 = ', ntime, idrun0
      else
         write (*,*) 'restarting from ntime = ', ntime
      endif
! read in size of electron array
      read (iur,iostat=ios) npp, ndimp
      if (ios /= 0) then
         write (*,*) 'np restart error, ios = ', ios
         stop
      endif
      if (ndimp /= size(part,1)) then
         write (*,*) 'restart error, idimp=', ndimp, size(part,1)
         stop
      endif
      if (npp /= np) then
         write (*,*) 'restart warning: new np/old np=', npp, np
         np = npp
      endif
! read in electrons, if non-zero
      if (npp > 0) then
         read (iur,iostat=ios) ((part(i,j),i=1,idimp),j=1,npp)
         if (ios /= 0) then
            write (*,*) 'electron array read error, ios = ', ios
            stop
         endif
      endif
! kpic = number of electrons in each tile
      if (.not.allocated(kpic)) allocate(kpic(mx1))
! find number of electrons in each of mx, tiles: updates kpic, nppmx
      call mdblkp2(part,kpic,nppmx,np,mx,irc)
! allocate vector electron data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      if (.not.allocated(ppart)) allocate(ppart(idimp,nppmx0,mx1))
      if (.not.allocated(ppbuff)) allocate(ppbuff(idimp,npbmx,mx1))
      if (.not.allocated(ncl)) allocate(ncl(2,mx1))
      if (.not.allocated(ihole)) allocate(ihole(2,ntmax+1,mx1))
! copy ordered electron data for OpenMP: updates ppart and kpic
      call mpmovin1(part,ppart,kpic,mx,irc)
! sanity check for electrons
      call mcheck1(ppart,kpic,nx,mx,irc)
      if (irc /= 0) stop
! read in movion to determine if ions are moving
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'movion restart error, ios = ', ios
         stop
      endif
      if (it /= movion) then
         write (*,*) 'movion restart error, movion = ', it, movion
         stop
      endif
! ions are moving
      if (movion > 0) then
! read in size of ion array
         read (iur,iostat=ios) npp, ndimp
         if (ios /= 0) then
            write (*,*) 'npi restart error, ios = ', ios
            stop
         endif
         if (ndimp /= size(part,1)) then
            write (*,*) 'ion restart error, idimp=',ndimp,size(part,1)
            stop
         endif
         if (npp /= npi) then
            write (*,*) 'restart warning: new npi/old npi=', npp, npi
            npi = npp
         endif
! read in ions, if non-zero
         if (npi > 0) then
            read (iur,iostat=ios) ((part(i,j),i=1,idimp),j=1,npp)
            if (ios /= 0) then
               write (*,*) 'ion array read error, ios = ', ios
               stop
            endif
         endif
! kipic = number of ions in each tile
         if (.not.allocated(kipic)) allocate(kipic(mx1))
! find number of ions in each of mx, tiles: updates kipic, nppmx
         call mdblkp2(part,kipic,nppmx,npi,mx,irc)
! allocate vector ion data
         nppmx1 = (1.0 + xtras)*nppmx
         if (.not.allocated(pparti)) allocate(pparti(idimp,nppmx1,mx1))
! copy ordered ion data for OpenMP: updates pparti and kipic
         call mpmovin1(part,pparti,kipic,mx,irc)
! sanity check for ions
         call mcheck1(pparti,kipic,nx,mx,irc)
         if (irc /= 0) stop
! ions are not moving, read in ion density
      else
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'qi size restart error, ios = ', ios
            stop
         endif
         if (it > size(qi)) then
            write (*,*) 'qi restart error, size(qi)=',it,size(qi)
            stop
         endif
         if (it > 0) then
            read (iur,iostat=ios) (qi(j),j=1,it)
            if (ios /= 0) then
               write (*,*) 'qi read error, ios = ', ios
               stop
            endif
         endif
      endif
! read in electric field parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'emf restart error, ios = ', ios
         stop
      endif
      if (it /= emf) then
         write (*,*) 'warning: emf values differ, emf=',it,emf
      endif
      ntime0 = ntime0 + ntime
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dwrite_restart1(iur)
! write out restart diagnostic file for electrostatic code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j, k, it, is, ios
      character(len=32) :: fname
! write out energy diagnostic parameter
      write (iur) ntw
      if (ntw > 0) then
         write (iur) itw
! write out time history array sizes and data
         if (itw > 0) then
            it = size(wt,2)
            write (iur) size(wt,1), it
            write (iur) ((wt(i,j),i=1,itw),j=1,it)
         endif
      endif
!
! write out electron density diagnostic parameter
      write (iur) ntde
! write out record location
      if (ntde > 0) then
         write (iur) nderec
! write out record length (zero if error) and file name (if no error)
         if (nderec > 0) then
            inquire(file=fdename,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = fdename
               write (iur) fname
            endif
         endif
      endif
!
! write out potential diagnostic parameter
      write (iur) ntp
! write out record location
      if (ntp > 0) then
         write (iur) nprec
! write out record length (zero if error) and file name (if no error)
         if (nprec > 0) then
            inquire(file=fpname,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = fpname
               write (iur) fname
            endif
         endif
! write out spectrum flag
         if ((ndp==2).or.(ndp==3)) then
            write (iur) itp
! write out spectrum sizes and data
            if (itp > 0) then
               it = size(pks,2); is = size(pks,3)
               write (iur) size(pks,1), it, is
               write (iur) (((pks(i,j,k),i=1,4),j=1,it),k=1,is)
            endif
         else
            it = 0
            write (iur) it
         endif
      endif
!
! write out longitudinal efield diagnostic parameter
      write (iur) ntel
! write out record location
      if (ntel > 0) then
         write (iur) nelrec
! write out record length (zero if error) and file name (if no error)
         if (nelrec > 0) then
            inquire(file=felname,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = felname
               write (iur) fname
            endif
         endif
      endif
!
! write out ion density diagnostic parameter
      if (movion > 0) then
         write (iur) ntdi
! write out record location
         if (ntdi > 0) then
            write (iur) ndirec
! write out record length (zero if error) and file name (if no error)
            if (ndirec > 0) then
               inquire(file=fdiname,recl=it,iostat=ios)
               if (ios /= 0) it = 0
               write (iur) it
               if (it > 0) then
                  fname = fdiname
                  write (iur) fname
               endif
            endif
! write out spectrum flag
            if ((nddi==2).or.(nddi==3)) then
               write (iur) itdi
! write out spectrum sizes and data
               if (itdi > 0) then
                  it = size(pksdi,2); is = size(pksdi,3)
                  write (iur) size(pksdi,1), it, is
                  write (iur) (((pksdi(i,j,k),i=1,4),j=1,it),k=1,is)
               endif
            else
               it = 0
               write (iur) it
            endif
         endif
      endif
!
! write out fluid moments diagnostic parameter
      write (iur) ntfm
      if (ntfm > 0) then
! write out electron record location
         write (iur) nferec
! write out record length (zero if error) and file name (if no error)
         if (nferec > 0) then
            inquire(file=ffename,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = ffename
               write (iur) fname
            endif
         endif
         if (movion > 0) then
! write out ion record location
            write (iur) nfirec
! write out record length (zero if error) and file name (if no error)
            if (nfirec > 0) then
               inquire(file=ffiname,recl=it,iostat=ios)
               if (ios /= 0) it = 0
               write (iur) it
               if (it > 0) then
                  fname = ffiname
                  write (iur) fname
               endif
            endif
         endif
      endif
!
! write out velocity diagnostic parameter
      write (iur) ntv
      if (ntv > 0) then
! write out electron record location
         write (iur) nverec
! write out record length (zero if error) and file name (if no error)
         if (nverec > 0) then
            inquire(file=fvename,recl=it,iostat=ios)
            if (ios /= 0) it = 0
            write (iur) it
            if (it > 0) then
               fname = fvename
               write (iur) fname
            endif
         endif
         write (iur) itv, ndv
! write out electron time history array sizes and data
         if (itv > 0) then
            if ((ndv==1).or.(ndv==3)) then
               it = size(fvtm,2)
               write (iur) size(fvtm,1), it, size(fvtm,3)
               write (iur) (((fvtm(i,j,k),i=1,itv),j=1,it),k=1,3)
            endif
         endif
         if (movion > 0) then
! write out ion record location
            write (iur) nvirec
! write out record length (zero if error) and file name (if no error)
            if (nvirec > 0) then
               inquire(file=fviname,recl=it,iostat=ios)
               if (ios /= 0) it = 0
               write (iur) it
               if (it > 0) then
                  fname = fviname
                  write (iur) fname
               endif
            endif
! write out ion time history array sizes and data
            if (itv > 0) then
               if ((ndv==2).or.(ndv==3)) then
                  it = size(fvtmi,2)
                  write (iur) size(fvtmi,1), it, size(fvtmi,3)
                  write (iur) (((fvtmi(i,j,k),i=1,itv),j=1,it),k=1,3)
               endif
            endif
         endif
      endif
!
! write out trajectory diagnostic parameter
      write (iur) ntt
      if (ntt > 0) then
! electron or ion trajectories
         if ((ndt==1).or.(ndt==2)) then
! write out trajectory record location
            write (iur) ntrec
! write out record length (zero if error) and file name (if no error)
            if (ntrec > 0) then
               inquire(file=ftname,recl=it,iostat=ios)
               if (ios /= 0) it = 0
               write (iur) it
               if (it > 0) then
                  fname = ftname
                  write (iur) fname
               endif
            endif
! store particle trajectories
            if ((nst==1).or.(nst==2)) then
               write (iur) itt
! write out time history sizes and data
               if (itt > 0) then
                  it = size(partd,2); is = size(partd,3)
                  write (iur) size(partd,1), it, is
                  write (iur) (((partd(i,j,k),i=1,itt),j=1,it),k=1,is)
               endif
            endif
         endif
      endif
!
! write out phase space diagnostic parameter
      write (iur) nts
      if (nts > 0) then
! write out electron phase space record location
         write (iur) nserec
! write out record length (zero if error) and file name (if no error)
         if (nserec > 0) then
            it = size(fvs)
            write (iur) it
            if (it > 0) then
               fname = fsename
               write (iur) fname
            endif
         endif
         if (movion > 0) then
! write out ion phase space record location
            write (iur) nsirec
! write out record length (zero if error) and file name (if no error)
            if (nsirec > 0) then
               it = size(fvsi)
               write (iur) it
               if (it > 0) then
                  fname = fsiname
                  write (iur) fname
               endif
            endif
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dread_restart1(iur)
! read in restart diagnostic file for electrostatic code
      implicit none
! iur = restart file descriptor
      integer, intent(in) :: iur
! local data
      integer :: i, j, k, nd, it, is, ir, ios
      character(len=32) :: fname
! read in energy diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntw restart error, ios = ', ios
         stop
      endif
      if (it /= ntw) then
         write (*,*) 'restart error: read/expected ntw=', it, ntw
         stop
      endif
      if (ntw > 0) then
         read (iur,iostat=ios) itw
         if (ios /= 0) then
            write (*,*) 'itw restart error, ios = ', ios
            stop
         endif
! read in time history array sizes and data
         if (itw > 0) then
            read (iur,iostat=ios) is, it
            if (ios /= 0) then
               write (*,*) 'ntw restart array size error, ios = ', ios
               stop
            endif
            if (is /= mtw) then
               write (*,*) 'restart error: read/expected mtw=', is, mtw
               stop
            endif
            if (it > size(wt,2)) then
               write (*,*) 'wt size error read/expected=', it,size(wt,2)
               stop
            endif
            read (iur,iostat=ios) ((wt(i,j),i=1,itw),j=1,it)
            if (ios /= 0) then
               write (*,*) 'wt array read error, ios = ', ios
               stop
            endif
! restore energy accumulations
            if (allocated(s)) then
               s = 0.0
               s(1) = wt(1,1)
               s(2) = wt(1,2)
               s(3) = dble(wt(1,4))
               s(4) = s(3)
               do it = 2, itw
                  s(1) = s(1) + wt(it,1)
                  s(2) = s(2) + wt(it,2)
                  s(3) = min(s(3),dble(wt(it,4)))
                  s(4) = max(s(4),dble(wt(it,4)))
               enddo
            endif
         endif
      endif
!
! read in electron density diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntde restart error, ios = ', ios
         stop
      endif
      if (it /= ntde) then
         write (*,*) 'restart error: read/expected ntde=', it, ntde
         stop
      endif
! read in record location
      if (ntde > 0) then
         read (iur,iostat=ios) nderec
         if (ios /= 0) then
            write (*,*) 'nderec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nderec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntde record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntde zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntde file name error, ios = ', ios
               stop
            endif
            fdename = fname
         endif
      endif
!
! read in potential diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntp restart error, ios = ', ios
         stop
      endif
      if (it /= ntp) then
         write (*,*) 'restart error: read/expected ntp=', it, ntp
         stop
      endif
! read in record location
      if (ntp > 0) then
         read (iur,iostat=ios) nprec
         if (ios /= 0) then
            write (*,*) 'nprec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nprec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntp record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntp zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntp file name error, ios = ', ios
               stop
            endif
            fpname = fname
         endif
! read in spectrum flag
         read (iur,iostat=ios) itp
         if (ios /= 0) then
            write (*,*) 'itp restart error, ios = ', ios
            stop
         endif
! read in spectrum sizes and data
         if (itp > 0) then
            read (iur,iostat=ios) ir, it, is
            if (ios /= 0) then
               write (*,*) 'ntp restart array size error, ios = ', ios
               stop
            endif
            if (ir /= 4) then
               write (*,*) 'pks size error: read/expected 4 =', ir
               stop
            endif
            if (it /= modesxp) then
               write (*,*) 'pks size error: read/expected modesxp=', it,&
     &modesxp
               stop
            endif
            if (is /= iw) then
               write (*,*) 'pks size error: read/expected iw=', is, iw
               stop
            endif
            read (iur,iostat=ios) (((pks(i,j,k),i=1,4),j=1,it),k=1,is)
            if (ios /= 0) then
                  write (*,*) 'pks array read error, ios = ', ios
                  stop
            endif
         endif
      endif
!
! read in longitudinal efield diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntel restart error, ios = ', ios
         stop
      endif
      if (it /= ntel) then
         write (*,*) 'restart error: read/expected ntel=', it, ntel
         stop
      endif
! read in record location
      if (ntel > 0) then
         read (iur,iostat=ios) nelrec
         if (ios /= 0) then
            write (*,*) 'nelrec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nelrec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntel record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntel zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntel file name error, ios = ', ios
               stop
            endif
            felname = fname
         endif
      endif
!
! read in ion density diagnostic parameter
      if (movion > 0) then
         read (iur,iostat=ios) it
         if (ios /= 0) then
            write (*,*) 'ntdi restart error, ios = ', ios
            stop
         endif
         if (it /= ntdi) then
            write (*,*) 'restart error: read/expected ntdi=', it, ntdi
            stop
         endif
! read in record location
         if (ntdi > 0) then
            read (iur,iostat=ios) ndirec
            if (ios /= 0) then
               write (*,*) 'ndirec restart error, ios = ', ios
               stop
            endif
! read in record length (zero if error) and file name (if no error)
            if (ndirec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'ntdi record length error, ios = ', ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'ntdi zero length record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'ntdi file name error, ios = ', ios
               stop
               endif
               fdiname = fname
            endif
! read in spectrum flag
            read (iur,iostat=ios) itdi
            if (ios /= 0) then
               write (*,*) 'itdi restart error, ios = ', ios
               stop
            endif
! read in spectrum sizes and data
            if (itdi > 0) then
               read (iur,iostat=ios) ir, it, is
               if (ios /= 0) then
                  write (*,*) 'ntdi restart array size error, ios=', ios
                  stop
               endif
               if (ir /= 4) then
                  write (*,*) 'pksdi size error: read/expected 4 =', ir
                  stop
               endif
               if (it /= modesxdi) then
                  write (*,*) 'pksdi size error: read/expected modesxdi=&
     &', it, modesxdi
                  stop
               endif
               if (is /= iwi) then
                  write (*,*) 'pksdi size error: read/expected iwi=',is,&
     &iwi
                  stop
               endif
               read (iur,iostat=ios) (((pksdi(i,j,k),i=1,4),j=1,it),    &
     &k=1,is)
               if (ios /= 0) then
                  write (*,*) 'pksdi array read error, ios = ', ios
                  stop
               endif
            endif
         endif
      endif
!
! read in fluid moments diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntfm restart error, ios = ', ios
         stop
      endif
      if (it /= ntfm) then
         write (*,*) 'restart error: read/expected ntfm=', it, ntfm
         stop
      endif
      if (ntfm > 0) then
! read in electron record location
         read (iur,iostat=ios) nferec
         if (ios /= 0) then
            write (*,*) 'nferec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nferec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntfm record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntfm zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntfm file name error, ios = ', ios
               stop
            endif
            ffename = fname
         endif
! read in ion data
         if (movion > 0) then
! read in ion record location
            read (iur,iostat=ios) nfirec
            if (ios /= 0) then
               write (*,*) 'nfirec restart error, ios = ', ios
               stop
            endif
! read in ion record length (zero if error) and file name (if no error)
            if (nfirec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'ntfm ion record length error, ios = ',ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'ntfm zero length ion record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'ntfm ion file name error, ios = ', ios
                  stop
               endif
               ffiname = fname
            endif
         endif
      endif
!
! read in velocity diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntv restart error, ios = ', ios
         stop
      endif
      if (it /= ntv) then
         write (*,*) 'restart error: read/expected ntv=', it, ntv
         stop
      endif
      if (ntv > 0) then
! read in electron record location
         read (iur,iostat=ios) nverec
         if (ios /= 0) then
            write (*,*) 'nverec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nverec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'ntv record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'ntv zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'ntv file name error, ios = ', ios
               stop
            endif
            fvename = fname
         endif
         read (iur,iostat=ios) itv, nd
         if (ios /= 0) then
            write (*,*) 'itv restart error, ios = ', ios
            stop
         endif
! read in electron time history array sizes and data
         if (itv > 0) then
            if ((nd==1).or.(nd==3)) then
               read (iur,iostat=ios) is, it, ir
               if (ios /= 0) then
                  write (*,*) 'ntv restart array size error, ios = ',ios
                  stop
               endif
               if (is /= mtv) then
                  write (*,*) 'restart error: read/expected mtv=',is,mtv
                  stop
               endif
               if (it /= ndim) then
                  write (*,*) 'fvtm size error read/expected ndim=', it,&
     &ndim
                  stop
               endif
               if (ir /= 3) then
                  write (*,*) 'fvtm size error read/expected 3=', ir
                  stop
               endif
               read (iur,iostat=ios) (((fvtm(i,j,k),i=1,itv),j=1,it),k=1,3)
               if (ios /= 0) then
                  write (*,*) 'fvtm array read error, ios = ', ios
                  stop
               endif
            endif
         endif
         if (movion > 0) then
! read in ion record location
            read (iur,iostat=ios) nvirec
            if (ios /= 0) then
               write (*,*) 'nvirec restart error, ios = ', ios
               stop
            endif
! read in record length (zero if error) and file name (if no error)
            if (nvirec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'ntv ion record length error, ios = ', ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'ntv ion zero length record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'ntv ion file name error, ios = ', ios
                  stop
               endif
               fviname = fname
            endif
! read in ion time history array sizes and data
            if (itv > 0) then
               if ((nd==2).or.(nd==3)) then
                  read (iur,iostat=ios) is, it, ir
                  if (ios /= 0) then
                     write (*,*) 'ntvi restart array size error,ios= ', &
     &ios
                     stop 
                  endif
                  if (is /= mtv) then
                     write (*,*) 'ion restart error: read/expected mtv='&
     &, is, mtv
                     stop
                  endif
                  if (it /= ndim) then
                     write (*,*) 'fvtmi size error read/expected ndim=',&
     & it, ndim
                     stop
                  endif
                  if (ir /= 3) then
                     write (*,*) 'fvtmi size error read/expected 3=', ir
                     stop
                  endif
                  read (iur,iostat=ios) (((fvtmi(i,j,k),i=1,itv),j=1,it)&
     &,k=1,3)
                  if (ios /= 0) then
                     write (*,*) 'fvtmi array read error, ios = ', ios
                     stop
                  endif
               endif
            endif
         endif
      endif
!
! read in trajectory diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'ntt restart error, ios = ', ios
         stop
      endif
      if (it /= ntt) then
         write (*,*) 'restart error: read/expected ntt=', it, ntt
         stop
      endif
      if (ntt > 0) then
! electron or ion trajectories
         if ((ndt==1).or.(ndt==2)) then
! read in trajectory record location
            read (iur,iostat=ios) ntrec
            if (ios /= 0) then
               write (*,*) 'ntrec restart error, ios = ', ios
               stop
            endif
! read in record length (zero if error) and file name (if no error)
            if (ntrec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'ntt record length error, ios = ', ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'ntt zero length record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'ntt file name error, ios = ', ios
                  stop
               endif
               ftname = fname
            endif
! restore particle trajectories
            if ((nst==1).or.(nst==2)) then
               read (iur,iostat=ios) itt
               if (ios /= 0) then
                  write (*,*) 'itt restart error, ios = ', ios
                  stop
               endif
! read in time history sizes and data
               if (itt > 0) then
                  read (iur,iostat=ios) ir, it, is
                  if (ios /= 0) then
                     write (*,*) 'ntt restart array size error, ios = ',&
     &ios
                     stop
                  endif
                  if (ir /= mtt) then
                     write (*,*) 'restart error: read/expected mtt=', ir&
     &,mtt
                     stop
                  endif
                  if (it /= idimp) then
                     write (*,*) 'partd size error read/expected idimp='&
     &,it, size(partd,2)
                     stop
                  endif
                  if (is /= nprobt) then
                     write (*,*) 'partd size err read/expected nprobt=',&
     &is, nprobt
                     stop
                  endif
                  read (iur,iostat=ios) (((partd(i,j,k),i=1,itt),j=1,it)&
     &,k=1,is)
                  if (ios /= 0) then
                     write (*,*) 'partd array read error, ios = ', ios
                     stop
                  endif
               endif
            endif
         endif
      endif
!
! read in phase space diagnostic parameter
      read (iur,iostat=ios) it
      if (ios /= 0) then
         write (*,*) 'nts restart error, ios = ', ios
         stop
      endif
      if (it /= nts) then
         write (*,*) 'restart error: read/expected nts=', it, nts
         stop
      endif
      if (nts > 0) then
! read in electron record location
         read (iur,iostat=ios) nserec
         if (ios /= 0) then
            write (*,*) 'nserec restart error, ios = ', ios
            stop
         endif
! read in record length (zero if error) and file name (if no error)
         if (nserec > 0) then
            read (iur,iostat=ios) it
            if (ios /= 0) then
               write (*,*) 'nts record length error, ios = ', ios
               stop
            endif
            if (it==0) then
               write (*,*) 'nts zero length record error'
               stop
            endif
            read (iur,iostat=ios) fname
            if (ios /= 0) then
               write (*,*) 'nts file name error, ios = ', ios
               stop
            endif
            fsename = fname
! reposition stream file
            call fnsets1(nserec-1,it,fsename)
         endif
! read in ion data
         if (movion > 0) then
! read in ion record location
            read (iur,iostat=ios) nsirec
            if (ios /= 0) then
               write (*,*) 'nsirec restart error, ios = ', ios
               stop
            endif
! read in ion record length (zero if error) and file name (if no error)
            if (nsirec > 0) then
               read (iur,iostat=ios) it
               if (ios /= 0) then
                  write (*,*) 'nts ion record length error, ios = ',ios
                  stop
               endif
               if (it==0) then
                  write (*,*) 'nts zero length ion record error'
                  stop
               endif
               read (iur,iostat=ios) fname
               if (ios /= 0) then
                  write (*,*) 'nts ion file name error, ios = ', ios
                  stop
               endif
               fsiname = fname
! reposition stream file
               call fnsets1(nsirec-1,it,fsiname)
            endif
         endif
      endif
!
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine close_restart1()
! close reset and restart files
      implicit none
! iur, iurr = restart, reset, old restart file descriptors
!     close(unit=iurr)
      if (nustrt==1) then
         if (ntr > 0) close(unit=iur)
      else if (nustrt==2) then
         close(unit=iur)
      else if (nustrt==0) then
         if (ntr > 0) close(unit=iur)
         if (idrun /= idrun0) close(unit=iur0)
      endif
      end subroutine
!
      end module
