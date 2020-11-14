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
        implicit none
!
        integer:: ierr
        real(mykind) :: a
!
#ifdef _OPENMP
        INTEGER:: nthreads, threadid
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
#ifdef MB
           write(16,*) "INFO: Multi-block approach is used" 
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
#ifdef MPIIO
           write(16,*) "INFO: MPI-IO version"
#else
           write(16,*) "INFO: using RAW I/O"
#endif
!
#ifdef PGI
           write(16,*) "INFO: using PGI compiler"
           write(16,*) "      quad precision not supported"
#endif
!
#ifdef OPENACC
           write(16,*) "INFO: using openacc"
#endif
!
#ifdef SENDRECV
           write(16,*) "INFO: using sendrecv comm. (packed)"
#endif
!
#ifdef NOBLOCK
           write(16,*) "INFO: using nonblocking comm. (packed)"
#endif
!
#ifdef NOMANAGED
           write(16,*) "INFO: using no-maneged data movement (openacc)"
#endif
!
#ifdef ORIGINAL
           write(16,*) "INFO: using Move & collide (Original) version"
#else
           write(16,*) "INFO: using Move+Collide (FUSED) version"
#endif
!
#ifdef SERIAL
           write(16,*) "INFO: serial version"
#else
           write(16,*) "INFO: MPI version"
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
