c Null general 1d gks graphics library
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: december 12, 2017
c-----------------------------------------------------------------------
      subroutine GROPEN
c this subroutine opens gks and activates standard workstation
c colors and maximum size of display surface are also set
c if a string (keyboard) device is available, it is initialized
c if a locator (mouse) device is available, it is initialized
      return
      end
c-----------------------------------------------------------------------
      subroutine GROPEN0(iecho)
c this subroutine opens gks and activates workstation
c colors and maximum size of display surface are also set
c if a string (keyboard) device is available, it is initialized
c if a locator (mouse) device is available, it is initialized
      return
      end
c-----------------------------------------------------------------------
      subroutine GRCLOSE
c this subroutine deactivates workstation and closes gks
      return
      end
c-----------------------------------------------------------------------
      subroutine SETNPLT(nplt,irc)
c this subroutine resets the maximum number of plots per page
c if requested number is negative, it is set to default (=1)
c the current plot location is also reset to the next available location
c if next available location is iplot = 0 and the old location was not
c 0, then the workstation is updated and a return code can be generated.
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,chrs
     1,irc)
c this subroutine plots ngs subarrays of the array f, on a common graph,
c each plot with nx points, versus a linear function in x,
c where xmin < x < xmax.
c depending on the number of colors in the display device, each subarray
c is plotted with a different color, given in order by:
c blue, red, yellow, cyan, magenta, green, and foreground.
c after all the colors are cycled through, then different line styles
c are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
c or different marker types if mks>0: dot, plus, star, circle, cross.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = array containing subarrays to be plotted
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plots have a common scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nx = number of points plotted in each subarray
c nxv = first dimension of array f, nxv >= nx
c ngs = second dimension of array f, number of subarrays to be plotted
c chr = additional long character string comment for plot
c chrs = array of ngs short character labels used by subroutine tickd
c to label individual line or marker samples
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      character*(*) label, chr
      character*(*) chrs(ngs)
      dimension f(nxv,ngs)
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine DISPC(f,g,label,zsc,zst,mks,nx,nxv,ngs,chr,chrs,irc)
c this subroutine plots ngs subarrays of the array f, on a common graph,
c each plot with nx points, versus the corresponding subarray of the
c array g.
c depending on the number of colors in the display device, each subarray
c is plotted with a different color, given in order by:
c blue, red, yellow, cyan, magenta, green, and foreground
c after all the colors are cycled through, then different lines styles
c are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
c or different marker types if mks>0: dot, plus, star, circle, cross.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f, g = arrays containing subarrays to be plotted
c label = long character string label for plot
c real(zsc)/aimag(zsc) = power of 2 scale of x/y coordinate for plot
c real(zst)/aimag(zst) = flag for positive and/or negative x/y values.
c the plots have a common scale in y given by ymax and ymin,
c where isc = int(aimag(zsc)) and ist = int(aimag(zst)), as follows:
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c the plots have a common scale in x given by xmax and xmin,
c where isc = int(real(zsc)) and ist = int(real(zst)), as follows:
c if ist = 0, then xmax = 2**isc and xmin = -2**isc.
c if ist = 1, then xmax = 2**isc and xmin = 0.
c if ist = -1, then xmax = 0 and xmin = -2**isc.
c if ist = 2, then xmin = gmin, xmax = gmin + 2**ir,
c where gmin/gmax are the function minimum/maximum, 
c and ir = power of 2 scale for (gmax - gmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of g
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nx = number of points plotted in each subarray
c nxv = first dimension of array f, nxv >= nx
c ngs = second dimension of array f, number of subarrays to be plotted
c chr = additional character string comment for plot
c chrs = array of ngs short character labels used by subroutine tickd
c to label individual line or marker samples
c irc = return code (0 = normal return)
      character*(*) label, chr
      character*(*) chrs(ngs)
      complex zsc, zst
      dimension f(nxv,ngs), g(nxv,ngs)
      irc = 0
      return
      end
      subroutine DISPS(f,label,xmin,xmax,isc,ist,nx,chr,irc)
c this subroutine plots an array f versus a linear function in x,
c where xmin < x < xmax.  It is plotted in solid line style, in blue
c if color is available.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = array to be plotted
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plot has a scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of f
c nx = dimension of array f, the number of points to be plotted.
c chr = additional long character string comment for plot
c irc = return code (0 = normal return)
      dimension f(nx)
      character*(*) label, chr
      irc = 0
      return
      end
      subroutine DISPD(f,g,label,zsc,zst,nx,chr,irc)
c this subroutine plots an array f versus an array g, in a solid line
c style, in blue if color is available.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f, g = arrays to be plotted
c label = long character string label for plot
c real(zsc)/aimag(zsc) = power of 2 scale of x/y coordinate for plot
c real(zst)/aimag(zst) = flag for positive and/or negative x/y values.
c the plot has a scale in y given by ymax and ymin,
c where isc = int(aimag(zsc)) and ist = int(aimag(zst)), as follows:
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist > 0, then ymax = 2**isc and ymin = 0.
c if ist < 0, then ymax = 0 and ymin = -2**isc.
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of f
c the plot has a scale in x given by xmax and xmin,
c where isc = int(real(zsc)) and ist = int(real(zst)), as follows:
c if ist = 0, then xmax = 2**isc and xmin = -2**isc.
c if ist > 0, then xmax = 2**isc and xmin = 0.
c if ist < 0, then xmax = 0 and xmin = -2**isc.
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of g
c nx = dimension of arrays f and g, the number of points to be plotted.
c chr = additional character string comment for plot
c irc = return code (0 = normal return)
      dimension f(nx), g(nx)
      character*(*) label, chr
      complex zsc, zst
      irc = 0
      return
      end
      subroutine CLRSCRN
c clear screen
      return
      end
      subroutine RSTSCRN
c clear screen and force input device into request mode for next time
      return
      end
      subroutine GTINPUT(label,prompt,input,irc)
c display label and prompt and request input
c irc = return code (0 = normal return)
      integer irc
      character*(*) label, prompt, input
      irc = 0
      end
      subroutine PTOTPUT(label,prompt)
c display label and prompt
      character*(*) label, prompt
      end
      subroutine GTMINPUT(label,prompt,input,ndim,irc)
c display labels and prompt and request input
c ndim = number of labels
c irc = return code (0 = normal return)
      integer ndim, irc
      character*(*) label, prompt, input
      dimension label(ndim)
      irc = 0
      end
      subroutine idtran
c this subroutine performs the identify transformation number 1
c rx, ry = ndc coordinates of upper-right corner of workstation window
      return
      end
      subroutine dfplps (chh,nfont,iprec)
c this subroutine sets default plotting parameters
c chh = character height, in world coordinates
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
c ifrg = index of foreground color
      return
      end
      subroutine readrc(irc)
c this subroutine reads return code from input device
c if mouse = 0, string device is used if present, otherwise locator
c if mouse = 1, locator device is used if present, otherwise string
c if neither device is present, the subroutine exits with irc = 0
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np,ir
     1c)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c plots background particles in blue and beam particles in red
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 5
c npx = number of background particles
c np = number of particles
c irc = return code (0 = normal return)
      dimension part(idimp,np)
      character*(*) label
      return
      end
c-----------------------------------------------------------------------
      subroutine BGRASP13(part,label,itime,isc,omx,omy,omz,nx,iyp,ixp,  
     1idimp,npx,np,irc)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c for magnetized plasma, rotating cartesian co-ordinates so that B
c points in the z direction.
c if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
c if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
c plots background particles in blue and beam particles in red
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4
c npx = number of background particles
c np = number of particles
c irc = return code (0 = normal return)
      dimension part(idimp,np)
      character*(*) label
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp,  
     1nppmx,mx1,irc)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
      return
      end
c-----------------------------------------------------------------------
      subroutine PBGRASP13(ppart,kpic,label,itime,isc,omx,omy,omz,nx,iyp
     1,ixp,idimp,nppmx,mx1,irc)
c for 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c for magnetized plasma, rotating cartesian co-ordinates so that B
c points in the z direction.
c if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
c if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c irc = return code (0 = normal return)
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
      return
      end
c-----------------------------------------------------------------------
      subroutine PMGRASP13(ppart,kpic,label,itime,isc,nx,iyp,ixp,idimp, 
     1nppmx,mx1,ltag,irc)
c for 1d or 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c with beam  particles marked for color
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c ppart(ltag,n,m) = particle id of tagged particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 3
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c ltag = location in particle array of tags (3 or 5)
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
      return
      end
c-----------------------------------------------------------------------
      subroutine PMBGRASP13(ppart,kpic,label,itime,isc,omx,omy,omz,nx,  
     1iyp,ixp,idimp,nppmx,mx1,ltag,irc)
c for 1d or 1-2/2d code, this subroutine displays (iyp-ixp) phase space
c with beam  particles marked for color
c for magnetized plasma, rotating cartesian co-ordinates so that B
c points in the z direction.
c if iyp=2, plot vperp1, if iyp=3, plot vperp2, if iyp=4, plot vparallel
c if ixp=2, plot vperp1, if ixp=3, plot vperp2, if ixp=4, plot vparallel
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c ppart(ltag,n,m) = particle id of tagged particle n in tile m
c kpic = number of particles per tile
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx = system length in x direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 3
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c ltag = location in particle array of tags
c irc = return code (0 = normal return
      dimension ppart(idimp,nppmx,mx1), kpic(mx1)
      character*(*) label
      return
      end
