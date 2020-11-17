!=====================================================================
!     ****** LBE/check_case
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple check of the configuration...
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
        subroutine check_case
!
        use storage
        use timing
!
#ifdef _OPENMP
        use omp_lib
#endif
!
#ifdef CUDAFOR
        use cudafor
#endif
!
        implicit none
!
        integer:: ierr
        real(mykind) :: a
!
#ifdef _OPENMP
        INTEGER:: nthreads, threadid
#endif
!
#ifdef CUDAFOR
        INTEGER:: ierr, nDevices
#endif
!
! Error section  (compiling will fail)
#if defined DRIVEN && defined COUETTE
        write(6,*) " ERROR: defined both DRIVEN and COUETTE	! the string is not closed: so error compiling!!!
        stop
#endif

#if defined ORIGINAL && defined TWO_FIELDS
        write(6,*) " ERROR: defined both ORIGINAL and TWO_FIELDS ! the string is not closed: so error compiling!!!
        stop
#endif

#if defined DRIVEN && defined INLET
        write(6,*) " ERROR: defined both DRIVEN and INLET     ! the string is not closed: so error compiling!!!
        stop
#endif

#if defined COUETTE && defined INLET
        write(6,*) " ERROR: defined both COUETTE and INLET     ! the string is not closed: so error compiling!!!
        stop
#endif
!
! Warning section  (run will proceed...) 
        if((u0.GT.0).AND.(u_inflow.GT.0)) then
           write(16,*) "WARNING: inflow ad volume force at the same & 
                        time", myrank, u0, u_inflow 
        endif
!
! Info section  
        if(myrank == 0) then 
#ifdef REGOLARIZED
           write(16,*) "INFO: Eegolarized collision is used" 
#endif

#ifdef MB
           write(16,*) "INFO: Multi-block approach is used" 
#endif
!
#ifdef OBSTACLES
           write(16,*) "INFO: The test case has obstacles"
#endif
!
#ifdef SPONGES
           write(16,*) "INFO: The test case has noslip-freeslip, & 
                              inflow-outflow bc "
           write(16,*) "WARNING: Hand made fix for serial free-slip"
#endif
!
#ifdef DRIVEN
           write(16,*) "INFO: The test case is driven cavity" 
#else
# ifdef INLET
           write(16,*) "INFO: The test case is channel with inlet" 
# else
#  ifdef COUETTE
           write(16,*) "INFO: The test case is couette"
#  else
#   ifdef PERIODIC
           write(16,*) "INFO: The test case is triperiodic"
#   else
           write(16,*) "INFO: The test case is channel"
#   endif
#  endif
# endif
#endif
!
#ifdef SLAB
           write(16,*) "INFO: using slab "
#endif
!
#ifdef WALL
           write(16,*) "INFO: using wall "
#endif
!
#ifdef FREESLIP
           write(16,*) "INFO: free slip simulation (still to fix)"
#endif
!
#ifdef DUCT
           write(16,*) "INFO: The test case is a duct (inflow/outflow)"
#endif
!
#ifdef CHECK2D
           write(16,*) "INFO:  2D-channel simulation (inflow/outflow)"
#endif
!
#ifdef FORCING_Y
           write(16,*) "INFO: volume forcing along y direction"
#else
# ifdef FORCING_Z
           write(16,*) "INFO: volume forcing along z direction"
# else
           write(16,*) "INFO: volume forcing along x direction"
# endif
#endif

#ifdef NO_OUTPUT
           write(16,*) "INFO: no output mode enabled (few I/O) "
#endif

#ifdef HPC
           write(16,*) "INFO: HPC mode enabled (all binary dump) "
#endif
 
#ifdef CILINDER
           write(16,*) "INFO: using cilinder "
#endif
!
#ifdef TEST1
           write(16,*) "INFO: TEST1 "
#endif
!
#ifdef QUAD_P
           write(16,*) "INFO: using quad precision"
#else
# ifdef DOUBLE_P
           write(16,*) "INFO: using double precision"
# else
#  ifdef HALF_P
           write(16,*) "INFO: using half precision"
#  else
           write(16,*) "INFO: using single precision"
#  endif
# endif
#endif
!
#ifdef MIXEDPRECISION
           write(16,*) "INFO: using mixed precision"
#endif
!
           write(16,*) "INFO: mykind=", mykind, "range  =", range(a)
           write(16,*) "INFO: mykind=", mykind, "huge   =", huge(a)
           write(16,*) "INFO: mykind=", mykind, "epsilon=", epsilon(a)
!
#ifdef HDF5
           write(16,*) "INFO: using I/O with HDF5"
#else
# ifdef MPIIO
           write(16,*) "INFO: MPI-IO version"
# else
           write(16,*) "INFO: using RAW I/O"
# endif
#endif
!
#ifdef PGI
           write(16,*) "INFO: using PGI compiler"
           write(16,*) "      quad precision not supported"
#endif
!
#ifdef ARM
           write(16,*) "INFO: using ARM version"
           write(16,*) "      quad precision not supported"
#endif
!
#ifdef OPENACC
           write(16,*) "INFO: using PGI openacc"
#endif

#ifdef CUDAFOR
           write(16,*) "INFO: using PGI cudafortran"
           ierr = cudaGetDeviceCount (nDevices)
           if (ierr /= cudaSuccess) then
              write (* ,*) cudaGetErrorString (ierr)
           else
              write(16,*) "INFO: GPU found"
           endif
 #ifdef DRIVEN
! nothing to do
 #else
  #ifdef PERIODIC
! still working on...
  #else
           write(16,*) "CUDAFOR version works only with Driven Cavity or Peridic BC  ! the string is not closed: so error compiling!!!  
  #endif
 #endif
#endif
!
#ifdef SENDRECV
# ifdef SENDRECV_TRY
           write(16,*) "INFO: using sendrecv comm. (not packed)"
# else
           write(16,*) "INFO: using sendrecv comm. (packed)"
# endif
#endif
!
#ifdef DUPLEX
           write(16,*) "INFO: using duplex comm. "
# endif
!
#ifdef NOBLOCK
# ifdef NOBLOCK_TRY
           write(16,*) "INFO: using nonblocking comm. (not packed)"
# else
           write(16,*) "INFO: using nonblocking comm. (packed)"
# endif
#endif

#ifdef ORIGINAL
           write(16,*) "INFO: using Move & collide (Original) version"
#else
# ifdef TWO_FIELDS 
           write(16,*) "INFO: using Move & collide (Two field) version"
# else
           write(16,*) "INFO: using Move+Collide (FUSED) version"
# endif
#endif

#ifdef SERIAL
           write(16,*) "INFO: serial version"
#else
           write(16,*) "INFO: MPI version"
#endif
!
#ifdef CYLINDER_INF
           if(myrank == 0) then
              write(6,*) "INFO: Infinite Cylinder test case"
              write(16,*) "INFO: Infinite Cylinder test case"
           endif
# ifdef DRAG
           if(myrank == 0) then
              write(6,*) "INFO: drag check activated"
              write(6,*) "WARNING: radius value handmade fix", radius
              write(16,*) "INFO: drag check activated"
              write(16,*) "WARNING: radius value handmade fix", radius
           endif
# endif
#endif
!
#ifdef _OPENMP
!$omp parallel
           nthreads = OMP_GET_NUM_THREADS()
           threadid = OMP_GET_THREAD_NUM()
           if(threadid.eq.0) then
              write(16,*) "INFO: using OpenMP version with threads = ", nthreads
           endif
!$omp end parallel
#endif
        endif

#ifdef TRYBOX
        write( 6,*) "WARNING: TRYBOX option activated border -->", border
        write(16,*) "WARNING: TRYBOX option activated border -->", border 
#endif
!
#ifdef SERIAL
! do nothing        
#else
        call MPI_barrier(MPI_COMM_WORLD,ierr)
#endif
!
#ifdef DEBUG_3
        if(myrank == 0) then
           write(6,*) "INFO: DEBUG3 mode enabled"
        endif
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "INFO: DEBUG2 mode enabled"
        endif
#endif
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "INFO: DEBUG1 mode enabled"
           write(6,*) "DEBUG1: Exiting from sub. check_case"
        endif
#endif
!
#ifdef MEM_CHECK
        if(myrank == 0) then
           mem_stop = get_mem();
           write(6,*) "MEM_CHECK: after sub check_case mem =", mem_stop
        endif
#endif
!
! error section, run will be killed
        if(rear(1) == 4) then
           if(u0.NE.0) then
              write(6,*) "ERROR: both inflow and volume forcing" 
              write(16,*) "ERROR: both inflow and volume forcing" 
              stop
           endif
        endif
!
        if(up(1) == 2) then
           if(u0.NE.0) then
              write(6,*) "ERROR: both moving wall and volume forcing" 
              write(16,*) "ERROR: both moving wall and volume forcing" 
              stop
           endif
        endif
        return
!
        end subroutine check_case
