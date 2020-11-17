!=====================================================================
!     ****** LBE/bcond_comm_duplex
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_comm_duplex
!     DESCRIPTION
!       
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       experimental
!
!     *****
!=====================================================================
!
        subroutine bcond_comm_duplex
!
        use timing
        use storage
#ifdef SERIAL
! do nothing
#else
        use mpi
#endif
!
        implicit none
!
        integer      :: i,j,k 
        integer      :: tag, ierr
       
        
#ifdef SERIAL
! do nothing
#else
        integer      :: status_front(MPI_STATUS_SIZE)
        integer      :: reqs_front(10)
        integer      :: status_rear(MPI_STATUS_SIZE)
        integer      :: reqs_rear(10)
        integer      :: status_left(MPI_STATUS_SIZE)
        integer      :: reqs_left(10)
        integer      :: status_right(MPI_STATUS_SIZE)
        integer      :: reqs_right(10)
        integer      :: status_up(MPI_STATUS_SIZE)
        integer      :: reqs_up(10)
        integer      :: status_down(MPI_STATUS_SIZE)
        integer      :: reqs_down(10)
#endif
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
#ifdef SERIAL
! do nothing
#else
        write(6,*) "zzzzzz0"
!
! first isend/recv from front
        if(front(1) == 0) then
           tag = 01
           call mpi_isend(a01(l,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(1), ierr)
           tag = 02
           call mpi_isend(a02(l,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(2), ierr)
           tag = 03
           call mpi_isend(a03(l,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(3), ierr)
           tag = 04
           call mpi_isend(a04(l,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(4), ierr)
           tag = 05
           call mpi_isend(a05(l,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(5), ierr)
!!
           tag = 10
           call mpi_irecv(a10(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(6), ierr)
           tag = 11
           call mpi_irecv(a11(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(7), ierr)
           tag = 12
           call mpi_irecv(a12(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(8), ierr)
           tag = 13
           call mpi_irecv(a13(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(9), ierr)
           tag = 14
           call mpi_irecv(a14(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, reqs_front(10), ierr)
!        endif
!
!        if(front(1) == 0) then
           call mpi_wait(1,reqs_front( 1), ierr)
           call mpi_wait(1,reqs_front( 2), ierr)
           call mpi_wait(1,reqs_front( 3), ierr)
           call mpi_wait(1,reqs_front( 4), ierr)
           call mpi_wait(1,reqs_front( 5), ierr)
           call mpi_wait(1,reqs_front( 6), ierr)
           call mpi_wait(1,reqs_front( 7), ierr)
           call mpi_wait(1,reqs_front( 8), ierr)
           call mpi_wait(1,reqs_front( 9), ierr)
           call mpi_wait(1,reqs_front(10), ierr)
        endif
!
! second recv from rear
        if(rear(1) == 0) then
           tag = 01
           call mpi_irecv(a01(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(1), ierr)
           tag = 02
           call mpi_irecv(a02(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(2), ierr)
           tag = 03
           call mpi_irecv(a03(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(3), ierr)
           tag = 04
           call mpi_irecv(a04(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(4), ierr)
           tag = 05
           call mpi_irecv(a05(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(5), ierr)
!!
           tag = 10
           call mpi_isend(a10(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(6), ierr)
           tag = 11
           call mpi_isend(a11(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(7), ierr)
           tag = 12
           call mpi_isend(a12(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(8), ierr)
           tag = 13
           call mpi_isend(a13(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(9), ierr)
           tag = 14
           call mpi_isend(a14(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(10), ierr)
        endif
!
        if(rear(1) == 0) then
           call mpi_waitall(1,reqs_rear(1), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(2), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(3), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(4), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(5), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(6), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(7), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(8), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(9), status_rear, ierr)
           call mpi_waitall(1,reqs_rear(10), status_rear, ierr)
        endif
!
!        call mpi_barrier(lbecomm,ierr)
!
        write(6,*) "zzzzzz10"
!--------------------------------------------------
! first isend from right

        if(right(1) == 0) then 
           tag = 03
           call mpi_isend(a03(0,m,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(1), ierr)
           tag = 07
           call mpi_isend(a07(0,m,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(2), ierr)
           tag = 08
           call mpi_isend(a08(0,m,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(3), ierr)
           tag = 09
           call mpi_isend(a09(0,m,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(4), ierr)
           tag = 12
           call mpi_isend(a12(0,m,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(5), ierr)
!!
           tag = 01
           call mpi_irecv(a01(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(6), ierr)
           tag = 10
           call mpi_irecv(a10(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(7), ierr)
           tag = 16
           call mpi_irecv(a16(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(8), ierr)
           tag = 17
           call mpi_irecv(a17(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(9), ierr)
           tag = 18
           call mpi_irecv(a18(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, reqs_right(10), ierr)
        endif
!
        write(6,*) "zzzzzz11"
        if(right(1) == 0) then
           call mpi_waitall(1,reqs_right(1), status_left, ierr)
           call mpi_waitall(1,reqs_right(2), status_left, ierr)
           call mpi_waitall(1,reqs_right(3), status_left, ierr)
           call mpi_waitall(1,reqs_right(4), status_left, ierr)
           call mpi_waitall(1,reqs_right(5), status_left, ierr)
           call mpi_waitall(1,reqs_right(6), status_left, ierr)
           call mpi_waitall(1,reqs_right(7), status_left, ierr)
           call mpi_waitall(1,reqs_right(8), status_left, ierr)
           call mpi_waitall(1,reqs_right(9), status_left, ierr)
           call mpi_waitall(1,reqs_right(10), status_left, ierr)
        endif
!
        write(6,*) "zzzzzz12"
! second recv from left
        if(left(1) == 0) then 
           tag = 03
           call mpi_irecv(a03(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(1), ierr)
           tag = 07
           call mpi_irecv(a07(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(2), ierr)
           tag = 08
           call mpi_irecv(a08(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(3), ierr)
           tag = 09
           call mpi_irecv(a09(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(4), ierr)
           tag = 12
           call mpi_irecv(a12(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(5), ierr)
!!
           tag = 01
           call mpi_isend(a01(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(6), ierr)
           tag = 10
           call mpi_isend(a10(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(7), ierr)
           tag = 16
           call mpi_isend(a16(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(8), ierr)
           tag = 17
           call mpi_isend(a17(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(9), ierr)
           tag = 18
           call mpi_isend(a18(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(10), ierr)
        endif
!
        write(6,*) "zzzzzz13"
        if(left(1) == 0) then
           call mpi_waitall(1,reqs_left(1), status_left, ierr)
           call mpi_waitall(1,reqs_left(2), status_left, ierr)
           call mpi_waitall(1,reqs_left(3), status_left, ierr)
           call mpi_waitall(1,reqs_left(4), status_left, ierr)
           call mpi_waitall(1,reqs_left(5), status_left, ierr)
           call mpi_waitall(1,reqs_left(6), status_left, ierr)
           call mpi_waitall(1,reqs_left(7), status_left, ierr)
           call mpi_waitall(1,reqs_left(8), status_left, ierr)
           call mpi_waitall(1,reqs_left(9), status_left, ierr)
           call mpi_waitall(1,reqs_left(10), status_left, ierr)
        endif
!
        call mpi_barrier(lbecomm,ierr)
!
        write(6,*) "zzzzzz2"
!--------------------------------------------------
! first isend from up
        if(up(1) == 0) then
           tag = 04
           call mpi_isend(a04(0,0,n), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(1), ierr)
           tag = 06
           call mpi_isend(a06(0,0,n), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(2), ierr)
           tag = 07
           call mpi_isend(a07(0,0,n), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(3), ierr)
           tag = 13
           call mpi_isend(a13(0,0,n), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(4), ierr)
           tag = 18
           call mpi_isend(a18(0,0,n), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(5), ierr)
!!
           tag = 02
           call mpi_irecv(a02(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(6), ierr)
           tag = 09
           call mpi_irecv(a09(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(7), ierr)
           tag = 11
           call mpi_irecv(a11(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(8), ierr)
           tag = 15
           call mpi_irecv(a15(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(9), ierr)
           tag = 16
           call mpi_irecv(a16(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, reqs_up(10), ierr)
        endif
!
! second recv from down
        if(down(1) == 0) then
           tag = 04
           call mpi_irecv(a04(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(1), ierr)
           tag = 06
           call mpi_irecv(a06(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(2), ierr)
           tag = 07
           call mpi_irecv(a07(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(3), ierr)
           tag = 13
           call mpi_irecv(a13(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(4), ierr)
           tag = 18
           call mpi_irecv(a18(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(5), ierr)
!!
           tag = 02
           call mpi_isend(a02(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(6), ierr)
           tag = 09
           call mpi_isend(a09(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(7), ierr)
           tag = 11
           call mpi_isend(a11(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(8), ierr)
           tag = 15
           call mpi_isend(a15(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(9), ierr)
           tag = 16
           call mpi_isend(a16(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(10), ierr)
        endif
!
        if(up(1) == 0) then
           call mpi_waitall(1,reqs_up(1), status_up, ierr)
           call mpi_waitall(1,reqs_up(2), status_up, ierr)
           call mpi_waitall(1,reqs_up(3), status_up, ierr)
           call mpi_waitall(1,reqs_up(4), status_up, ierr)
           call mpi_waitall(1,reqs_up(5), status_up, ierr)
           call mpi_waitall(1,reqs_up(6), status_up, ierr)
           call mpi_waitall(1,reqs_up(7), status_up, ierr)
           call mpi_waitall(1,reqs_up(8), status_up, ierr)
           call mpi_waitall(1,reqs_up(9), status_up, ierr)
           call mpi_waitall(1,reqs_up(10), status_up, ierr)
        endif
        if(down(1) == 0) then
           call mpi_waitall(1,reqs_down(1), status_down, ierr)
           call mpi_waitall(1,reqs_down(2), status_down, ierr)
           call mpi_waitall(1,reqs_down(3), status_down, ierr)
           call mpi_waitall(1,reqs_down(4), status_down, ierr)
           call mpi_waitall(1,reqs_down(5), status_down, ierr)
           call mpi_waitall(1,reqs_down(6), status_down, ierr)
           call mpi_waitall(1,reqs_down(7), status_down, ierr)
           call mpi_waitall(1,reqs_down(8), status_down, ierr)
           call mpi_waitall(1,reqs_down(9), status_down, ierr)
           call mpi_waitall(1,reqs_down(10), status_down, ierr)
        endif
!
        call mpi_barrier(lbecomm,ierr)
!
#endif
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_mp = time_mp + real(countA1-countA0)/(count_rate)
        time_mp1 = time_mp1 + (tcountA1-tcountA0)
!
! now local bc...
        call bcond_bc
!
!

!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_comm_duplex"
        endif
#endif
!
        return
        end subroutine bcond_comm_duplex

