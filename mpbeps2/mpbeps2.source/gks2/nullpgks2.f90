!-----------------------------------------------------------------------
! Null general 2d parallel gks graphics library
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: february 5, 2017
      module plibgks2
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine IPLTCOMM(nplt)
! initialize required elements of common block plotcm,
! to allow only one MPI node to display graphics.
      implicit none
      integer, intent(in) :: nplt
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PGRCLOSE
! this subroutine deactivates workstation and closes gks
      implicit none
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PDSYNC(irc)
! this subroutine synchronizes irc and iplot element in plotcm
      implicit none
      integer, intent(inout) :: irc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PCARPET(f,g,nyp,nvp,label,isc,ist,nx,ny,nxv,nypmx,chr, &
     &ntc,irc)
! this subroutine displays an array f as a color raster image, for
! distributed data
! a 256 color palette must have been defined prior to this call. 
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = distributed field array to be plotted
! g = scratch array for receiving messages
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plots, determined by the absolute value of f
! nx/ny = length of field f in x/y direction
! nxv = first dimension of field array f, must be >= nx
! nypmx = maximum size of particle partition, including guard cells
! chr = additional long character string comment for plot
! ntc = number of valid colors, should be power of 2, <= 256
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, isc, ist, nx, ny, nxv, nypmx, ntc
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv,nypmx), intent(in) :: f
      real, dimension(nxv,nypmx), intent(inout) ::  g
      irc = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PCONTUR(f,g,lf,nyp,nvp,label,isc,ist,nx,ny,nxv,nypmx,  &
     &chr,nc,irc)
! this subroutine displays an array f as a contour plot.
! a maximum of ncols colors are used, used in order from lowest to
! highest contour: blue, green, cyan, foreground, yellow, magenta, red
! multiple plots per page can be displayed by dividing the screen into
! n x n subregions, where n*n is the next largest integer >= nplot
! the location (ix,iy) of a plot in the subregions is determined by
! the parameter iplot = ix + iy*n
! f = field array to be plotted
! lf = scratch field array
! g = scratch array for receiving messages
! nvp = number of real or virtual processors requested
! label = long character string label for plot
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! the range of values of f are given by fmax and fmin.
! if ist = 0, then fmax = 2**isc and fmin = -2**isc.
! if ist = 1, then fmax = 2**isc and fmin = 0.
! if ist = -1, then fmax = 0 and fmin = -2**isc.
! if ist = 2, then fmax = fmin + 2**ir,
! where fmin/fmax are the function minimum/maximum, 
! and ir = power of 2 scale for (fmax - fmin)
! if abs(isc) < 116, then the isc value passed is used for scale.
! if abs(isc) > 116, then the program finds the minimum value of isc
! which will contain the plots, determined by the absolute value of f
! nx/ny = length of field f in x/y direction
! nxv = first dimension of field array f, must be >= nx
! nypmx = maximum size of particle partition, including guard cells
! chr = additional long character string comment for plot
! nc = number of contour lines
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: nyp, nvp, isc, ist, nx, ny, nxv, nypmx, nc
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label, chr
      real, dimension(nxv,nypmx), intent(in) :: f
      real, dimension(nxv,nypmx), intent(inout) ::  g
      integer, dimension(nxv,nypmx), intent(inout) ::  lf
      irc = 0
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PGRASP23(part,f,npp,label,itime,isc,nx,ny,iyp,ixp,idimp&
     &,npmax,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! f = scratch array for receiving messages
! npp = number of particles in partition
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: npp
      integer, intent(in) :: itime, isc, nx, ny, idimp, npmax, iyp, ixp
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(idimp,npmax), intent(in) :: part
      real, dimension(2*npmax), intent(inout) :: f
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PBGRASP23(part,f,npp,label,itime,isc,omx,omy,omz,nx,ny,&
     &iyp,ixp,idimp,npmax,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! for magnetized plasma, rotating cartesian co-ordinates so that B
! points in the z direction.
! if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
! if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! part(5,n) = velocity vz of particle n in partition
! f = scratch array for receiving messages
! npp = number of particles in partition
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! npmax = maximum number of particles in each partition
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: npp
      integer, intent(in) :: itime, isc, nx, ny, idimp, npmax, iyp, ixp
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(idimp,npmax), intent(in) :: part
      real, dimension(2*npmax), intent(inout) :: f
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPGRASP23(ppart,f,kpic,label,itime,isc,nx,ny,iyp,ixp,  &
     &idimp,nppmx,mxyp1,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! ppart(1,n,m) = position x of particle n in tile m in partition
! ppart(2,n,m) = position y of particle n in tile m in partition
! ppart(3,n,m) = velocity vx of particle n in tile m in partition
! ppart(4,n,m) = velocity vy of particle n in tile m in partition
! ppart(5,n,m) = velocity vz of particle n in tile m in partition
! kpic = number of particles per tile
! f = scratch array for receiving messages
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! and where mx1 = (system length in x direction - 1)/mx + 1
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp, idimp, nppmx
      integer, intent(in) :: mxyp1
      integer, intent(inout) :: irc
      character(len=*), intent(in) :: label
      real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
      real, dimension(2*nppmx*mxyp1), intent(inout) :: f
      integer, dimension(mxyp1) :: kpic
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBGRASP23(ppart,f,kpic,label,itime,isc,omx,omy,omz,nx,&
     &ny,iyp,ixp,idimp,nppmx,mxyp1,irc)
! for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
! for magnetized plasma, rotating cartesian co-ordinates so that B
! points in the z direction.
! if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
! if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
! ppart(1,n,m) = position x of particle n in tile m in partition
! ppart(2,n,m) = position y of particle n in tile m in partition
! ppart(3,n,m) = velocity vx of particle n in tile m in partition
! ppart(4,n,m) = velocity vy of particle n in tile m in partition
! ppart(5,n,m) = velocity vz of particle n in tile m in partition
! kpic = number of particles per tile
! f = scratch array for receiving messages
! label = species label
! itime = current time step
! isc = power of 2 scale of range of values of velocity
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
! nx/ny = system length in x/y direction
! iyp/ixp = phase space coordinates to be displayed
! idimp = size of phase space = 4 or 5
! nppmx = maximum number of particles in tile
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! and where mx1 = (system length in x direction - 1)/mx + 1
! irc = return code (0 = normal return)
      implicit none
      integer, intent(in) :: itime, isc, nx, ny, iyp, ixp, idimp, nppmx
      integer, intent(in) :: mxyp1
      integer, intent(inout) :: irc
      real, intent(in) :: omx, omy, omz
      character(len=*), intent(in) :: label
      real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
      real, dimension(2*nppmx*mxyp1), intent(inout) :: f
      integer, dimension(mxyp1) :: kpic
      end subroutine
!
      end module


