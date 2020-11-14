!=====================================================================
!     ****** LBE/profile
!
!     COPYRIGHT
!       (c) 2018-20?? by CINECA/G.Amati
!     NAME
!       profile
!     DESCRIPTION
!       prifle subroutine:
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine profile(itime,itfin,tstep)
!
      use storage
      use timing
      implicit none
!
      integer:: itime, itfin, tstep
! 
! here I am signal
      call SYSTEM_CLOCK(countD1, count_rate, count_max)
      call time(tcountD1)
      time_inn_loop = real(countD1-countD0)/(count_rate)
      time_inn_loop1 = time_inn_loop1 + (tcountD1-tcountD0)
      write(6,1001)(time_inn_loop)/tstep,itime,itfin
!
#ifdef MEM_CHECK
      mem_stop = get_mem();
      write(6,*) "MEM_CHECK: iteration", itime, " mem =", mem_stop
#endif
!
! some info about time...
       write(100,2001) itime, (time_mp-old1),   &
                              (time_coll-old2), & 
                              (time_move-old3)
!
       old1 = time_mp
       old2 = time_coll
       old3 = time_move
!
       call SYSTEM_CLOCK(countD0, count_rate, count_max)
       call time(tcountD0)
!
! formats...
1001  format(" Mean   time",1(e14.6,1x),i8,"/",i8)
2001  format(i8, 3(e14.6,1x))
!
      return
      end subroutine profile
