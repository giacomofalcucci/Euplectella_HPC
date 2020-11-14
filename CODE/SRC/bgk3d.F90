! =====================================================================
!     ****** LBE/bgk3D
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bgk2d
!     DESCRIPTION
!       main program for LBM 3D
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables used: itfin, itstart, ivtim
!                               itime, itsave, icheck, itrestart
!                               isignal
!       real variables used: tempo1, tempo2     
!       open the following unit: 16 (bgk.log)
!                                60 (prof_k.dat)
!                                61 (prof_i.dat)
!                                62 (u_med.dat)
!                                63 (diagno.dat)
!                                68 (probe.dat)
!                                69 (bgk.perf)
!       velocity directions:
!        direction  1    unit vector = ( 1,-1, 0)   tipo 2
!        direction  2    unit vector = ( 1, 0,-1)   tipo 2
!        direction  3    unit vector = ( 1, 1, 0)   tipo 2
!        direction  4    unit vector = ( 1, 0, 1)   tipo 2
!        direction  5    unit vector = ( 1, 0, 0)   tipo 1
!        direction  6    unit vector = ( 0, 0, 1)   tipo 1
!        direction  7    unit vector = ( 0, 1, 1)   tipo 2
!        direction  8    unit vector = ( 0, 1, 0)   tipo 1
!        direction  9    unit vector = ( 0, 1,-1)   tipo 2
!        direction 10    unit vector = (-1,-1, 0)   tipo 2
!        direction 11    unit vector = (-1, 0,-1)   tipo 2
!        direction 12    unit vector = (-1, 1, 0)   tipo 2
!        direction 13    unit vector = (-1, 0, 1)   tipo 2
!        direction 14    unit vector = (-1, 0, 0)   tipo 1
!        direction 15    unit vector = ( 0, 0,-1)   tipo 1
!        direction 16    unit vector = ( 0,-1,-1)   tipo 2
!        direction 17    unit vector = ( 0,-1, 0)   tipo 1
!        direction 18    unit vector = ( 0,-1, 1)   tipo 2
!        direction 19    unit vector = ( 0, 0, 0)   tipo 0
!                              
!     *****
! =====================================================================
!
      program bgk3d
!
      use storage
      use timing
      use real_kinds
#ifdef OPENACCOLD
      use openacc
#endif
!
!
#ifdef _OPENMP
      use omp_lib
#endif
!
      implicit none
!
      INTEGER:: itfin, itstart, ivtim, ierr
      INTEGER:: itime, itsave, icheck, itrestart, init_v
      INTEGER:: isignal
      integer:: tstep, tt
!
#ifdef LIKWID
      call likwid_markerInit()
!
! set up the simulation...
      call likwid_markerStartRegion("init")
#endif
!
! set up the simulation...
      call setup(itfin,ivtim,isignal,itsave,icheck,itrestart,init_v,tstep)
!
! initialize the flow...
!
      call initialize(itrestart,init_v,itfin,itstart,ivtim,isignal,itsave,icheck,tstep)
!
#ifdef LIKWID
      call likwid_markerStopRegion("init")
#endif
!
      call SYSTEM_CLOCK(countE0, count_rate, count_max)
      call time(tcountE0)
!
      call SYSTEM_CLOCK(countD0, count_rate, count_max)
      call time(tcountD0)
!
!      call check_isend
! main loop starts here.....
!
      if(myrank==0) then 
         call system("date       >> time.log")
      endif
!
#ifdef LIKWID
      call likwid_markerStartRegion("loop")
#endif
!
#ifdef NOMANAGED
!$acc data copy(a01,a02,a03,a04 &
!$acc          ,a05,a06,a07,a08 &
!$acc          ,a09,a10,a11,a12 &
!$acc          ,a13,a14,a15,a16 &
!$acc          ,a17,a18,a19     &
!$acc          ,b01,b02,b03,b04 &
!$acc          ,b05,b06,b07,b08 &
!$acc          ,b09,b10,b11,b12 &
!$acc          ,b13,b14,b15,b16 &
!$acc          ,b17,b18)
#endif
!
      do itime=itstart+1,itfin
!
#ifdef DEBUG_2
         if(myrank == 0) then
            write(6,*) "DEBUG2: starting time step =", itime
         endif
#endif
!
!         call SYSTEM_CLOCK(countD0, count_rate, count_max)
!         call time(tcountD0)
!
#ifdef LIKWID
         call likwid_markerStartRegion("boundaries")
#endif
!
         call boundaries         ! boundary conditions
!
#ifdef LIKWID
         call likwid_markerStopRegion("boundaries")
#endif
!
         call propagation        ! propagation step
!
#ifdef LIKWID
         call likwid_markerStartRegion("collision")
#endif
!
         call collision(itime)   ! collision step
!
#ifdef LIKWID
         call likwid_markerStopRegion("collision")
#endif
!
! get macroscopic values
         call diagnostic(itime,ivtim,icheck,itsave)
!
! get timing/profiling values
         if (mod(itime,isignal).eq.0) then
            if (myrank == 0 ) then 
               call profile(itime,itfin,isignal) 
            endif
         endif
      enddo
!#endif
!
#ifdef NOMANAGED
!$acc end data 
#endif
!
#ifdef LIKWID
      call likwid_markerStopRegion("loop")
#endif
!
!     some global timings
      call SYSTEM_CLOCK(countE1, count_rate, count_max)
      call time(tcountE1)
      time_loop = real(countE1-countE0)/(count_rate)
      time_loop1 = tcountE1-tcountE0
!
      call finalize(itstart,itfin)     ! finalize all
!
#ifdef LIKWID
      call likwid_markerClose()
#endif
!
      if(myrank==0) then
         call system("date       >> time.log")
         write(6,*) "That's all folks!!!!"
      endif
!
      end program bgk3d
