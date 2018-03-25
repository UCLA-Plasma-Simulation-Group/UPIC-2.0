c     
c     smile_f77.f by Fred Hucht (C) 1993-2007.
c     Short and simple example program for FORTRAN bindings.
c     Compile with
c     make smile_f77 (uses f77 compiler)
c     or
c     make smile_f2c (uses f2c compiler)
c     
c     NOTE: This example uses the FORTRAN bindings from Ygl,
c           so most routine- and parameter-names are truncated
c           to six characters.
c     
c     $Id: smile_f77.f,v 3.3 2007-05-09 12:27:49+02 fred Exp $
c     
      program smile
      integer*4 win

      include "../irisgl/fgl.h"


      
      call prefsi(100, 100)
      win = winope('Smile!', 6)

c     write(*,*) win

c     The background
      call color(BLACK)
      call clear()

c     The face
      call color(RED)
      call circfi(50, 50, 40)
  
c     The eyes
      call color(WHITE)
      call circfi(30, 60, 10)
      call circfi(70, 60, 10)
      call color(BLACK)
      call circfi(30, 60,  5)
      call circfi(70, 60,  5)
      
c     The smile
      call arci(50, 50, 25, 2000, 3400)

c     The sleep
      call sleep(2)

c     The twinkle
      call color(RED)
      call circfi(30, 65, 10)
      call sleep(1)
      call color(WHITE)
      call circfi(30, 60, 10)
      call color(BLACK)
      call circfi(30, 60,  5)

c     The end
      call sleep(2)
      call gexit()
      end
