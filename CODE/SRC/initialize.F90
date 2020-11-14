!=====================================================================
!     ****** LBE/initialize
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different initializations routines..
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
       subroutine initialize(itrestart,init_v,itfin,itstart,ivtim,isignal,itsave,icheck,tstep)
!
       use storage
       use timing
!
       implicit none
!
       integer:: i
       integer:: ierr
       integer:: itstart,itrestart,init_v
       integer:: itfin,ivtim,isignal,itsave,icheck,tstep
!
! some other initializations
!
       time_coll  = 0.0
       time_move  = 0.0
       time_obs   = 0.0
       time_bc    = 0.0
       time_mp    = 0.0
       time_dg    = 0.0
       time_io    = 0.0
       time_dev   = 0.0
       time_coll1 = 0.0
       time_move1 = 0.0
       time_obs1  = 0.0
       time_bc1   = 0.0
       time_mp1   = 0.0
       time_dg1   = 0.0
       time_io1   = 0.0
       time_dev1  = 0.0
       timeZ      = 0.0
       timeY      = 0.0
       timeX      = 0.0
       old1       = 0.0
       old2       = 0.0
       old3       = 0.0
       i_shift = 1
!
! set boundary condition flags..
!
       call build_bcond
!
! check
#ifdef SERIAL
!do nothing
#else
#ifdef CHECK
       if(myrank == 0) then
             open(37,file='bgk.topo'); close(37)
       endif
!
       do i = 0, nprocs
          if(i == myrank) then
             write(*,*) "writing topo file for task", myrank
             open(37,file='bgk.topo', access='append')
             write(37,*) myrank, ":my x coord is -->", mpicoords(1)
             write(37,*) myrank, ":my z coord is -->", mpicoords(2)
             write(37,*) myrank, ":my front   is -->", front(2), front(1)
             write(37,*) myrank, ":my rear    is -->", rear(2),  rear(1)
             write(37,*) myrank, ":my up      is -->", up(2),    up(1)
             write(37,*) myrank, ":my down    is -->", down(2),  down(1)
             close(37)
          endif
          call mpi_barrier(lbecomm,ierr)
       enddo
#endif
#endif
!
! build obstacles
#ifdef OBSTACLES
       call build_obs
#endif
!
!
! build starting configuration
!
       if(itrestart.eq.1) then
!
          call time(tcountF0)
          call SYSTEM_CLOCK(countF0, count_rate, count_max)
!
#ifdef HDF5
# ifdef SERIAL
          call restore_hdf5_serial(itstart)
# else
          call restore_hdf5_parallel(itstart)
# endif
#else
# ifdef MPIIO
          call restore_mpiio(itstart)
# else
          call restore_raw(itstart)
# endif
#endif
!
          if(myrank == 0) then
             write(*,*) "restoring at timestep ", itstart
          endif
!
          call SYSTEM_CLOCK(countF1, count_rate, count_max)
          call time(tcountF1)
          time_io = time_io + real(countF1-countF0)/(count_rate)
          time_io1 = time_io1 + tcountF1-tcountF0
!
#ifdef DEBUG_1
          if(myrank == 0) then
             write(6,*) "DEBUG1: I/O time (1)", real(countF1-countF0)/(count_rate)
             write(6,*) "DEBUG1: I/O time (2)", tcountF1-tcountF0
          endif
#endif
!
       else
          if(itrestart.ge.2) then
             itstart = 0
             write(6,*) "still not implemented"
             stop
          else
             itstart = 0
             call init(init_v)
             call diagno(itstart)
#ifdef NO_OUTPUT
! do nothing
#else
             call varm(itstart)
             call prof_i(itstart,m/2,n/2)
             call prof_j(itstart,l/2,n/2)
             call prof_k(itstart,l/2,m/2)
#endif
          endif
       endif
!
! compute collision parameters
       call hencol
!
       if(myrank==0) then
          call outdat(itfin,itstart,ivtim,isignal,itsave,icheck,tstep)
       endif
!
#ifdef NO_OUTPUT
! do nothing
#else
       call prof_i(0,m/2,n/2)
       call prof_j(0,l/2,n/2)
       call prof_k(0,l/2,m/2)
#endif
!
#ifdef DEBUG_1
       if(myrank == 0) then
          write(6,*) "DEBUG1: Exiting from sub. initialize"
       endif
#endif
!
# ifdef MEM_CHECK
       if(myrank == 0) then
          mem_stop = get_mem();
          write(6,*) "MEM_CHECK: after intialize. mem =", mem_stop
       endif
# endif

       return
       end subroutine initialize
