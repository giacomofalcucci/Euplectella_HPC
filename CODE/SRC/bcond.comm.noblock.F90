!=====================================================================
!     ****** LBE/bcond_comm_noblock
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       
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
        subroutine bcond_comm_noblock
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
        integer      :: reqs_front(5)
        integer      :: status_rear(MPI_STATUS_SIZE)
        integer      :: reqs_rear(5)
        integer      :: status_left(MPI_STATUS_SIZE)
        integer      :: reqs_left(5)
        integer      :: status_right(MPI_STATUS_SIZE)
        integer      :: reqs_right(5)
        integer      :: status_up(MPI_STATUS_SIZE)
        integer      :: reqs_up(5)
        integer      :: status_down(MPI_STATUS_SIZE)
        integer      :: reqs_down(5)
#endif
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
#ifdef SERIAL
! do nothing
#else
!
        call time(tcountX0)
! first isend from front
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
        endif
!
! second recv from rear
        if(rear(1) == 0) then
           tag = 01
           call mpi_recv(a01(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, status_front, ierr)
           tag = 02
           call mpi_recv(a02(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, status_front, ierr)
           tag = 03
           call mpi_recv(a03(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, status_front, ierr)
           tag = 04
           call mpi_recv(a04(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, status_front, ierr)
           tag = 05
           call mpi_recv(a05(0,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, status_front, ierr)
        endif
!
! third isend from rear
        if(rear(1) == 0) then
           tag = 10
           call mpi_isend(a10(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(1), ierr)
           tag = 11
           call mpi_isend(a11(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(2), ierr)
           tag = 12
           call mpi_isend(a12(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(3), ierr)
           tag = 13
           call mpi_isend(a13(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(4), ierr)
           tag = 14
           call mpi_isend(a14(1,0,0), 1, yzplane, rear(2), tag, &
                             lbecomm, reqs_rear(5), ierr)
        endif
!
! forth recv from front
        if(front(1) == 0) then
           tag = 10
           call mpi_recv(a10(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, status_rear, ierr)
           tag = 11
           call mpi_recv(a11(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, status_rear, ierr)
           tag = 12
           call mpi_recv(a12(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, status_rear, ierr)
           tag = 13
           call mpi_recv(a13(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, status_rear, ierr)
           tag = 14
           call mpi_recv(a14(l+1,0,0), 1, yzplane, front(2), tag, &
                             lbecomm, status_rear, ierr)
        endif
!
        if(front(1) == 0) then
           call mpi_wait(reqs_front(1), status_front, ierr)
           call mpi_wait(reqs_front(2), status_front, ierr)
           call mpi_wait(reqs_front(3), status_front, ierr)
           call mpi_wait(reqs_front(4), status_front, ierr)
           call mpi_wait(reqs_front(5), status_front, ierr)
        endif
        if(rear(1) == 0) then
           call mpi_wait(reqs_rear(1), status_rear, ierr)
           call mpi_wait(reqs_rear(2), status_rear, ierr)
           call mpi_wait(reqs_rear(3), status_rear, ierr)
           call mpi_wait(reqs_rear(4), status_rear, ierr)
           call mpi_wait(reqs_rear(5), status_rear, ierr)
        endif
!
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
!        call mpi_barrier(lbecomm,ierr)
!--------------------------------------------------
! first isend from right

        call time(tcountY0)
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
        endif
!
! second recv from left
        if(left(1) == 0) then 
           tag = 03
           call mpi_recv(a03(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status_right, ierr)
           tag = 07
           call mpi_recv(a07(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status_right, ierr)
           tag = 08
           call mpi_recv(a08(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status_right, ierr)
           tag = 09
           call mpi_recv(a09(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status_right, ierr)
           tag = 12
           call mpi_recv(a12(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status_right, ierr)
        endif
!
! third isend from left
        if(left(1) == 0) then 
           tag = 01
           call mpi_isend(a01(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(1), ierr)
           tag = 10
           call mpi_isend(a10(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(2), ierr)
           tag = 16
           call mpi_isend(a16(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(3), ierr)
           tag = 17
           call mpi_isend(a17(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(4), ierr)
           tag = 18
           call mpi_isend(a18(0,1,0), 1, xzplane, left(2), tag, &
                          lbecomm, reqs_left(5), ierr)
        endif
!
! forth recv from right
        if(right(1) == 0) then
           tag = 01
           call mpi_recv(a01(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status_left, ierr)
           tag = 10
           call mpi_recv(a10(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status_left, ierr)
           tag = 16
           call mpi_recv(a16(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status_left, ierr)
           tag = 17
           call mpi_recv(a17(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status_left, ierr)
           tag = 18
           call mpi_recv(a18(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status_left, ierr)

        endif
!
        if(right(1) == 0) then
           call mpi_wait(reqs_right(1), status_left, ierr)
           call mpi_wait(reqs_right(2), status_left, ierr)
           call mpi_wait(reqs_right(3), status_left, ierr)
           call mpi_wait(reqs_right(4), status_left, ierr)
           call mpi_wait(reqs_right(5), status_left, ierr)
        endif
        if(left(1) == 0) then
           call mpi_wait(reqs_left(1), status_left, ierr)
           call mpi_wait(reqs_left(2), status_left, ierr)
           call mpi_wait(reqs_left(3), status_left, ierr)
           call mpi_wait(reqs_left(4), status_left, ierr)
           call mpi_wait(reqs_left(5), status_left, ierr)
        endif
!
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!
!        call mpi_barrier(lbecomm,ierr)
!
!--------------------------------------------------
! first isend from up
        call time(tcountZ0)
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
        endif
!
! second recv from down
        if(down(1) == 0) then
           tag = 04
           call mpi_recv(a04(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, status_up, ierr)
           tag = 06
           call mpi_recv(a06(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, status_up, ierr)
           tag = 07
           call mpi_recv(a07(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, status_up, ierr)
           tag = 13
           call mpi_recv(a13(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, status_up, ierr)
           tag = 18
           call mpi_recv(a18(0,0,0), 1, xyplane, down(2), tag, &
                             lbecomm, status_up, ierr)
        endif
!
! third isend from down
        if(down(1) == 0) then
           tag = 02
           call mpi_isend(a02(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(1), ierr)
           tag = 09
           call mpi_isend(a09(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(2), ierr)
           tag = 11
           call mpi_isend(a11(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(3), ierr)
           tag = 15
           call mpi_isend(a15(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(4), ierr)
           tag = 16
           call mpi_isend(a16(0,0,1), 1, xyplane, down(2), tag, &
                             lbecomm, reqs_down(5), ierr)
        endif
!
! forth recv from up
        if(up(1) == 0) then
           tag = 02
           call mpi_recv(a02(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, status_down, ierr)
           tag = 09
           call mpi_recv(a09(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, status_down, ierr)
           tag = 11
           call mpi_recv(a11(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, status_down, ierr)
           tag = 15
           call mpi_recv(a15(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, status_down, ierr)
           tag = 16
           call mpi_recv(a16(0,0,n+1), 1, xyplane, up(2), tag, &
                             lbecomm, status_down, ierr)
                             
        endif
!
#ifdef QQQQQQ
        if(up(1) == 0) then
           call mpi_wait(reqs_up(1), status_up, ierr)
           call mpi_wait(reqs_up(2), status_up, ierr)
           call mpi_wait(reqs_up(3), status_up, ierr)
           call mpi_wait(reqs_up(4), status_up, ierr)
           call mpi_wait(reqs_up(5), status_up, ierr)
        endif
        if(down(1) == 0) then
           call mpi_wait(reqs_down(1), status_down, ierr)
           call mpi_wait(reqs_down(2), status_down, ierr)
           call mpi_wait(reqs_down(3), status_down, ierr)
           call mpi_wait(reqs_down(4), status_down, ierr)
           call mpi_wait(reqs_down(5), status_down, ierr)
        endif
#endif
!
!
#endif
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
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
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
#ifdef SERIAL
! do nothing
#else
        if(up(1) == 0) then
           call mpi_wait(reqs_up(1), status_up, ierr)
           call mpi_wait(reqs_up(2), status_up, ierr)
           call mpi_wait(reqs_up(3), status_up, ierr)
           call mpi_wait(reqs_up(4), status_up, ierr)
           call mpi_wait(reqs_up(5), status_up, ierr)
        endif
        if(down(1) == 0) then
           call mpi_wait(reqs_down(1), status_down, ierr)
           call mpi_wait(reqs_down(2), status_down, ierr)
           call mpi_wait(reqs_down(3), status_down, ierr)
           call mpi_wait(reqs_down(4), status_down, ierr)
           call mpi_wait(reqs_down(5), status_down, ierr)
        endif
!
!        call mpi_barrier(lbecomm,ierr)
#endif
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_mp = time_mp + real(countA1-countA0)/(count_rate)
        time_mp1 = time_mp1 + (tcountA1-tcountA0)
!
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_comm_noblock"
        endif
#endif
!
        return
        end subroutine bcond_comm_noblock

