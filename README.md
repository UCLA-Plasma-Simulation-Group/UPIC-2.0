# UPIC 2.0
2D MPI/OpenMP and 3D OpenMP Particle-in-Cell (PIC) codes  
BEPS2 and BEPS3  
by Viktor K. Decyk, UCLA  
Copyright 1994-2017, Regents of the University of California

These codes are part of the UPIC 2.0.2 Framework.  The primary purpose
of this framework is to provide trusted components for building and
customizing parallel Particle-in-Cell codes for plasma simulation.
The framework provides libraries as well as reference
applications illustrating their use.  The libraries are designed to
work with up to 3 levels of parallelization (using MPI for distributed
memory, OpenMP for shared memory, local vectorization).  Currently only
the first two are supported.  The reference applications are designed to
provide basic functionality common to many PIC codes, and are intended
to be customized and specialized by users for specific applications.
Currently, only spectral field solvers supported, but the particle
methods are designed to be gneral.  In addition, other test codes will
also be provided, illustrating the use of some specific component.  This
is primarily intended to users who wish to incorporate some component
into an already existing code. Currently, only one test code is
provided.

# Upon cloning the repository

If you clone this repository, we ask that you __please contact__ Ben Winjum (bwinjum@ucla.edu). The development of UPIC relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies.

Please also feel free to send an email to this address if you would like assistance using the code, if you have any questions related to the code base, or if you have suggestions for improvements. We are in the process of establishing a user forum for discussions related to the use and development of UPIC and will be more than happy to include you in such a forum if you desire.

# Further Details

Directions for compiling and executing the codes are available in the README files of each subfolder.

A Software License for use of these codes are in the documents:  
mpbeps2/Documents/License(CommercialReservation)  
mpbeps3/Documents/License(CommercialReservation)

Descriptions of the design of these codes and their components are in the
documents:  
mpbeps2/Documents/BEPS2Design.pdf  
mpbeps3/Documents/BEPS3Design.pdf

Description of the current capabilities of these codes and their
components are in the documents:  
mpbeps2/Documents/BEPS2Overview.pdf  
mpbeps3/Documents/BEPS3Overview.pdf

Details about the mathematical equations and units used in these codes are
given in the documents:  
mpbeps2/Documents/UPICModels.pdf  
mpbeps3/Documents/UPICModels.pdf

Please see the READMEs in each folder for further details.
