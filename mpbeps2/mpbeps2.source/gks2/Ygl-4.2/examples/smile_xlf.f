c     
c     smile.f by Fred Hucht (C) %VAL(1993-2007).
c     Short and simple example program. This program
c     works under AIX. 
c     Compile with
c     
c     make smile_xlf
c     
c     NOTE: This example uses the C routines from libYgl.a.
c     
c     
c     $Id: smile_xlf.f,v 3.3 2007-05-09 12:27:55+02 fred Exp $
c     
c     
      program smile
      integer*4 win
      
      include "../irisgl/fgl.h"
c     fgl.h only declares winope, so...
      integer*4 winopen
      
      call prefsize(%VAL(100),%VAL(100))
      win = winopen('Smile!')
      
c     write(*,*) win
      
c     The background
      call color(%VAL(BLACK))
      call clear()
      
c     The face
      call color(%VAL(RED))
      call circfi(%VAL(50), %VAL(50), %VAL(40))
      
c     The eyes
      call color(%VAL(WHITE))
      call circfi(%VAL(30), %VAL(60), %VAL(10))
      call circfi(%VAL(70), %VAL(60), %VAL(10))
      call color(%VAL(BLACK))
      call circfi(%VAL(30), %VAL(60),  %VAL(5))
      call circfi(%VAL(70), %VAL(60),  %VAL(5))
      
c     The smile
      call arci(%VAL(50), %VAL(50), %VAL(25), %VAL(2000), %VAL(3400))
      
c     The sleep
      call sleep(%VAL(2))
      
c     The twinkle
      call color(%VAL(RED))
      call circfi(%VAL(30), %VAL(65), %VAL(10))
      call sleep(%VAL(1))
      call color(%VAL(WHITE))
      call circfi(%VAL(30), %VAL(60), %VAL(10))
      call color(%VAL(BLACK))
      call circfi(%VAL(30), %VAL(60),  %VAL(5))
      
c     The end
      call sleep(%VAL(2))
      call gexit()
      end
