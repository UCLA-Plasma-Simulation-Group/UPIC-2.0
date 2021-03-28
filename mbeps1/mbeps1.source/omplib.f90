! OpenMP utility library
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 12, 2016
      module omplib
      use omp_lib
      implicit none
!
! common data for parallel processing
! nthreads = number of threads used
      integer :: nthreads
      save
!
      private
      public :: INIT_OMP, SETNTHSIZE, GETNTHSIZE
!
      contains
!
!-----------------------------------------------------------------------
      subroutine INIT_OMP(nth)
! initialize openmp library
! use nth threads if nth > 0; otherwise, use the number found
      implicit none
      integer, intent(in) :: nth
! local data
      integer :: ncpus
! determine how many processors are available
      ncpus = omp_get_num_procs()
      write (*,*) 'number of cpus found = ', ncpus
      nthreads = omp_get_max_threads()
      write (*,*) 'maximum number of threads = ', nthreads
      if (nth > 0) nthreads = nth
      call omp_set_num_threads(nthreads)
      write (*,*) 'using ',  nthreads, ' thread(s)'
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine SETNTHSIZE(nth)
! set number of threads
      implicit none
      integer, intent(in) :: nth
      if (nth > 0) nthreads = nth
      call omp_set_num_threads(nthreads)
      end subroutine
!
!-----------------------------------------------------------------------
      function GETNTHSIZE()
! get number of threads
      implicit none
      integer :: GETNTHSIZE
      GETNTHSIZE = nthreads
      end function
!
      end module omplib
