c Null general 2d gks graphics library
c written by viktor k. decyk, ucla
c copyright 1996, regents of the university of california
c update: february 8, 2017
c-----------------------------------------------------------------------
      subroutine CARPET(f,label,isc,ist,nx,ny,nxv,chr,ntc,irc)
c this subroutine displays an array f as a color raster image.
c a 256 color palette must have been defined prior to this call. 
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c ntc = number of valid colors, should be power of 2, <= 256
c irc = return code (0 = normal return)
      character*(*) label, chr
      dimension f(nxv,ny)
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine CARPETL(f,label,xmin,xmax,ymin,ymax,isc,ist,nx,ny,nxv,c
     1hr,ntc,irc)
c this subroutine displays an array f as a color raster image.
c a 256 color palette must have been defined prior to this call. 
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c ntc = number of valid colors, should be power of 2, <= 256
c irc = return code (0 = normal return)
c npald = number of palette entries
      character*(*) label, chr
      dimension f(nxv,ny)
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine CONTUR(f,lf,label,isc,ist,nx,ny,nxv,chr,nc,irc)
c this subroutine displays an array f as a contour plot.
c a maximum of ncols colors are used, used in order from lowest to
c highest contour: blue, green, cyan, foreground, yellow, magenta, red
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c lf = scratch field array
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c chr = additional long character string comment for plot
c nc = number of contour lines
c irc = return code (0 = normal return)
      character*(*) label, chr
      dimension f(nxv,ny), lf(nxv,ny)
      irc = 0
      return
      end
      subroutine GRASP23(part,label,itime,isc,nx,ny,iyp,ixp,idimp,npxy,n
     1p,irc)
c for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
c plots background particles in blue and beam particles in red
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx/ny = system length in x/y direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 5
c npxy = number of background particles
c np = number of particles
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      dimension part(idimp,np)
      character*(*) label
      irc = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine STPALIT(idpal)
c this subroutine selects one from three available palettes
c idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
      return
      end

