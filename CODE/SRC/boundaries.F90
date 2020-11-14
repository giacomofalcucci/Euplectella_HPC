!=====================================================================
!     ****** LBE/boundaries
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for boundaries routine...
!       - call bcond  else
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
        subroutine boundaries
!
        use timing
        use storage
!
        implicit none
        integer:: ierr, i
        real(mykind):: temp1, temp2
!
#ifdef OBSTACLES
! start timing
        call bcond_comm_noblock_packed_try
        call SYSTEM_CLOCK(countO0, count_rate, count_max)
        call time(tcountO0)
!
        call bcond_obs
!
! stop timing
        call time(tcountO1)
        call SYSTEM_CLOCK(countO1, count_rate, count_max)
        time_obs = time_obs + real(countO1-countO0)/(count_rate)
        time_obs1 = time_obs1 + (tcountO1-tcountO0)
#endif
!
#ifdef SENDRECV
# ifdef SENDRECV_TRY
! standard sendrecv
!        call bcond_comm
!
! sendrecv + section (1 for every sendrecv)
        call bcond_comm_openmp
        write(6,*) "ERROR: Wrong comm"
        stop
# else
        call bcond_comm_packed
# endif
#endif
!
#ifdef NOBLOCK
# ifdef NOBLOCK_TRY
        call bcond_comm_noblock
!        write(6,*) "ERROR: Wrong comm"
!        stop
# else
        call bcond_comm_noblock_packed
!        call bcond_comm_noblock_packed_overlap
# endif
#endif
!
!
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. boundaries"
        endif
#endif
!
        return
      end subroutine boundaries
