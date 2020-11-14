!=====================================================================
!     ****** LBE/collision
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for movef routine...
!       - movef (if -DORIGINAL)
!       - none  else
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
        subroutine collision(itime)
!
        use timing
        use storage
!
        implicit none
!
        integer:: itime, ierr, istat
!
! start timing
        call SYSTEM_CLOCK(countC0, count_rate, count_max)
        call time(tcountC0)
!
#ifdef ORIGINAL
        call col(itime)
#else
        call col_MC(itime)
#endif
!
! stop timing
        call time(tcountC1)
        call SYSTEM_CLOCK(countC1, count_rate, count_max)
        time_coll = time_coll + real(countC1-countC0)/(count_rate)
        time_coll1 = time_coll1 + (tcountC1-tcountC0)
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. collision"
        endif
#endif
!
        return
        end subroutine collision
