!=====================================================================
!     ****** LBE/bcond_comm
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_comm
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
        subroutine bcond_comm
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
        integer      :: status(MPI_STATUS_SIZE)
#endif
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
#ifdef SERIAL
!do nothing
#else
!------------------------------------------------------------------------
! comms along z + 
!
        call time(tcountZ0)

        tag = 04
        call mpi_sendrecv(a04(0,0,n), 1, xyplane, up(2), tag,      &
                          a04(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 06
        call mpi_sendrecv(a06(0,0,n), 1, xyplane, up(2), tag,      &
                          a06(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 07
        call mpi_sendrecv(a07(0,0,n), 1, xyplane, up(2), tag,      &
                          a07(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 13
        call mpi_sendrecv(a13(0,0,n), 1, xyplane, up(2), tag,      &
                          a13(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 18
        call mpi_sendrecv(a18(0,0,n), 1, xyplane, up(2), tag,      &
                          a18(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
!
!
! comms along z - 
!
        tag = 02
        call mpi_sendrecv(a02(0,0,  1), 1, xyplane, down(2), tag,  &
                          a02(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 09
        call mpi_sendrecv(a09(0,0,  1), 1, xyplane, down(2), tag,  &
                          a09(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 11
        call mpi_sendrecv(a11(0,0,  1), 1, xyplane, down(2), tag,  &
                          a11(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 15
        call mpi_sendrecv(a15(0,0,  1), 1, xyplane, down(2), tag,  &
                          a15(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 16
        call mpi_sendrecv(a16(0,0,  1), 1, xyplane, down(2), tag,  &
                          a16(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
!
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
!!!        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along x + 
        call time(tcountX0)
        tag = 01
        call mpi_sendrecv(a01(l,0,0), 1, yzplane, front(2), tag,   &
                          a01(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 02
        call mpi_sendrecv(a02(l,0,0), 1, yzplane, front(2), tag,   &
                          a02(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 03
        call mpi_sendrecv(a03(l,0,0), 1, yzplane, front(2), tag,   &
                          a03(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 04
        call mpi_sendrecv(a04(l,0,0), 1, yzplane, front(2), tag,   &
                          a04(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 05
        call mpi_sendrecv(a05(l,0,0), 1, yzplane, front(2), tag,   &
                          a05(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
!!!        call mpi_barrier(lbecomm,ierr)
!
! comms along x - 
        tag = 10
        call mpi_sendrecv(a10(  1,0,0), 1, yzplane, rear(2), tag,  &
                          a10(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 11
        call mpi_sendrecv(a11(  1,0,0), 1, yzplane, rear(2), tag,  &
                          a11(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 12
        call mpi_sendrecv(a12(  1,0,0), 1, yzplane, rear(2), tag,  &
                          a12(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 13
        call mpi_sendrecv(a13(  1,0,0), 1, yzplane, rear(2), tag,  &
                          a13(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 14
        call mpi_sendrecv(a14(  1,0,0), 1, yzplane, rear(2), tag,  &
                          a14(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
!
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
!!!        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along y + 
        call time(tcountY0)
        tag = 3
        call mpi_sendrecv(a03(0,m,0), 1, xzplane, right(2), tag,  &
                          a03(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 7
        call mpi_sendrecv(a07(0,m,0), 1, xzplane, right(2), tag,  &
                          a07(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 8
        call mpi_sendrecv(a08(0,m,0), 1, xzplane, right(2), tag,  &
                          a08(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 9
        call mpi_sendrecv(a09(0,m,0), 1, xzplane, right(2), tag,  &
                          a09(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 12
        call mpi_sendrecv(a12(0,m,0), 1, xzplane, right(2), tag,  &
                          a12(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
!!!        call mpi_barrier(lbecomm,ierr)
!
! comms along y - 
        tag = 1
        call mpi_sendrecv(a01(0,  1,0), 1, xzplane, left(2), tag,  &
                          a01(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 10
        call mpi_sendrecv(a10(0,  1,0), 1, xzplane, left(2), tag,  &
                          a10(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 16
        call mpi_sendrecv(a16(0,  1,0), 1, xzplane, left(2), tag,  &
                          a16(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 17
        call mpi_sendrecv(a17(0,  1,0), 1, xzplane, left(2), tag,  &
                          a17(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 18
        call mpi_sendrecv(a18(0,  1,0), 1, xzplane, left(2), tag,  &
                          a18(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
!
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!
        call mpi_barrier(lbecomm,ierr)
!
!
#endif
!
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_mp = time_mp + real(countA1-countA0)/(count_rate)
        time_mp1 = time_mp1 + (tcountA1-tcountA0)
!
! now local bc...
        call bcond_bc
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_comm"
        endif
#endif
        return
        end subroutine bcond_comm
