!=====================================================================
!     ****** LBE/propagation
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
        subroutine propagation
!
        use timing
        use storage
!
        implicit none
!

#ifdef ORIGINAL
! start timing
        call SYSTEM_CLOCK(countB0, count_rate, count_max)
        call time(tcountB0)
        call movef
!
! stop timing
        call time(tcountB1)
        call SYSTEM_CLOCK(countB1, count_rate, count_max)
        time_move = time_move + real(countB1-countB0)/(count_rate)
        time_move1 = time_move1 + (tcountB1-tcountB0)
#endif

#ifdef TWO_FIELDS
! start timing
        call SYSTEM_CLOCK(countB0, count_rate, count_max)
        call time(tcountB0)
!
        call movef_TF
!
! stop timing
        call time(tcountB1)
        call SYSTEM_CLOCK(countB1, count_rate, count_max)
        time_move = time_move + real(countB1-countB0)/(count_rate)
        time_move1 = time_move1 + (tcountB1-tcountB0)
#endif


        return
        end subroutine propagation
