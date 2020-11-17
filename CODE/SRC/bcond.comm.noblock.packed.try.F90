!=====================================================================
!     ****** LBE/bcond_comm_noblock_packed_try
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
!       It is nothing else then bcond_comm_noblock_packed without call for bc. 
!       to use only for obstacles
!       total time can reduced checking if there's an obstacle between tasks...
!
!     *****
!=====================================================================
!
        subroutine bcond_comm_noblock_packed_try
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
! 
#ifdef SERIAL
! do nothing
#else
        integer      :: status_front(MPI_STATUS_SIZE)
        integer      :: reqs_front(1)
        integer      :: status_rear(MPI_STATUS_SIZE)
        integer      :: reqs_rear(1)
        integer      :: status_left(MPI_STATUS_SIZE)
        integer      :: reqs_left(1)
        integer      :: status_right(MPI_STATUS_SIZE)
        integer      :: reqs_right(1)
        integer      :: status_up(MPI_STATUS_SIZE)
        integer      :: reqs_up(1)
        integer      :: status_down(MPI_STATUS_SIZE)
        integer      :: reqs_down(1)
!
        real(mykind), dimension(0:m+1,0:n+1,5) :: bufferXIN, bufferXIN2
        real(mykind), dimension(0:l+1,0:n+1,5) :: bufferYIN, bufferYIN2
        real(mykind), dimension(0:l+1,0:m+1,5) :: bufferZIN, bufferZIN2
        real(mykind), dimension(0:m+1,0:n+1,5) :: bufferXOUT,bufferXOUT2
        real(mykind), dimension(0:l+1,0:n+1,5) :: bufferYOUT,bufferYOUT2
        real(mykind), dimension(0:l+1,0:m+1,5) :: bufferZOUT,bufferZOUT2
!
        integer      :: msgsizeX
        integer      :: msgsizeY
        integer      :: msgsizeZ
!
        msgsizeX = (n+2)*(m+2)*5
        msgsizeY = (l+2)*(n+2)*5
        msgsizeZ = (l+2)*(m+2)*5
!
#endif
!
! start global timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
#ifdef SERIAL
! do nothing
#else
!
! --------------------------------------------------------------------------
! comms along x 
! --------------------------------------------------------------------------
!
        call time(tcountX0)
! first isend 
! isend from front
        if(front(1) == 0) then
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXIN,a01,a02,a03,a04,a05,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do j = 0,m+1
                 bufferXIN(j,k,1)=a01(l,j,k)
                 bufferXIN(j,k,2)=a02(l,j,k)
                 bufferXIN(j,k,3)=a03(l,j,k)
                 bufferXIN(j,k,4)=a04(l,j,k)
                 bufferXIN(j,k,5)=a05(l,j,k)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 01
           call mpi_isend(bufferXIN(0,0,1),msgsizeX,MYMPIREAL,front(2), tag, &
                             lbecomm, reqs_front(1), ierr)
        endif
!
! isend from rear
        if(rear(1) == 0) then
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXIN2,a10,a11,a12,a13,a14,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do j = 0,m+1
                 bufferXIN2(j,k,1)=a10(1,j,k)
                 bufferXIN2(j,k,2)=a11(1,j,k)
                 bufferXIN2(j,k,3)=a12(1,j,k)
                 bufferXIN2(j,k,4)=a13(1,j,k)
                 bufferXIN2(j,k,5)=a14(1,j,k)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 10
           call mpi_isend(bufferXIN2(0,0,1), msgsizex, MYMPIREAL, rear(2), tag, &
                             lbecomm, reqs_rear(1), ierr)
        endif
!
! second recv 
! recv from rear
        if(rear(1) == 0) then
!
           tag = 01
           call mpi_recv(bufferXOUT(0,0,1),msgsizeX,MYMPIREAL,rear(2), tag, &
                             lbecomm, status_front, ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXOUT,a01,a02,a03,a04,a05,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do j = 0,m+1
                 a01(0,j,k) = bufferXOUT(j,k,1)
                 a02(0,j,k) = bufferXOUT(j,k,2)
                 a03(0,j,k) = bufferXOUT(j,k,3)
                 a04(0,j,k) = bufferXOUT(j,k,4)
                 a05(0,j,k) = bufferXOUT(j,k,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
        endif

! forth recv from front
        if(front(1) == 0) then
!
           tag = 10
           call mpi_recv(bufferXOUT2(0,0,1), msgsizeX, MYMPIREAL, front(2), tag, &
                             lbecomm, status_rear, ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXOUT2,a10,a11,a12,a13,a14,m,n,l) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do j = 0,m+1
                 a10(l+1,j,k) = bufferXOUT2(j,k,1)
                 a11(l+1,j,k) = bufferXOUT2(j,k,2)
                 a12(l+1,j,k) = bufferXOUT2(j,k,3)
                 a13(l+1,j,k) = bufferXOUT2(j,k,4)
                 a14(l+1,j,k) = bufferXOUT2(j,k,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
        endif
!
        if(front(1) == 0) then
           call mpi_wait(reqs_front(1), status_front, ierr)
        endif
        if(rear(1) == 0) then
           call mpi_wait(reqs_rear(1), status_rear, ierr)
        endif
!
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
!        call mpi_barrier(lbecomm,ierr)
!
! --------------------------------------------------------------------------
! comms along y 
! --------------------------------------------------------------------------
! first isend 

        call time(tcountY0)
        if(right(1) == 0) then 
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYIN,a03,a07,a08,a09,a12,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do i = 0,l+1
                 bufferYIN(i,k,1)=a03(i,m,k)
                 bufferYIN(i,k,2)=a07(i,m,k)
                 bufferYIN(i,k,3)=a08(i,m,k)
                 bufferYIN(i,k,4)=a09(i,m,k)
                 bufferYIN(i,k,5)=a12(i,m,k)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 03
           call mpi_isend(bufferYIN(0,0,1), msgsizeY, MYMPIREAL, right(2), tag, &
                          lbecomm, reqs_right(1), ierr)
!
        endif
!
! isend from left
        if(left(1) == 0) then 
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYIN2,a01,a10,a16,a17,a18,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do i = 0,l+1
                 bufferYIN2(i,k,1)=a01(i,1,k)
                 bufferYIN2(i,k,2)=a10(i,1,k)
                 bufferYIN2(i,k,3)=a16(i,1,k)
                 bufferYIN2(i,k,4)=a17(i,1,k)
                 bufferYIN2(i,k,5)=a18(i,1,k)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 01
           call mpi_isend(bufferYIN2(0,0,1), msgsizey, MYMPIREAL, left(2), tag, &
                          lbecomm, reqs_left(1), ierr)
!
        endif
!

! second recv from left
        if(left(1) == 0) then
!
           tag = 03
           call mpi_recv(bufferYOUT(0,0,1), msgsizey, MYMPIREAL, left(2), tag, &
                          lbecomm, status_right, ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYOUT,a03,a07,a08,a09,a12,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do i = 0,l+1
                 a03(i,0,k)=bufferYOUT(i,k,1)
                 a07(i,0,k)=bufferYOUT(i,k,2)
                 a08(i,0,k)=bufferYOUT(i,k,3)
                 a09(i,0,k)=bufferYOUT(i,k,4)
                 a12(i,0,k)=bufferYOUT(i,k,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
        endif


! forth recv from right
        if(right(1) == 0) then
!
           tag = 01
           call mpi_recv(bufferYOUT2(0,0,1), msgsizey, MYMPIREAL, right(2), tag, &
                          lbecomm, status_left, ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYOUT2,a01,a10,a16,a17,a18,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do k = 0,n+1
!$acc loop independent
              do i = 0,l+1
                 a01(i,m+1,k)=bufferYOUT2(i,k,1)
                 a10(i,m+1,k)=bufferYOUT2(i,k,2)
                 a16(i,m+1,k)=bufferYOUT2(i,k,3)
                 a17(i,m+1,k)=bufferYOUT2(i,k,4)
                 a18(i,m+1,k)=bufferYOUT2(i,k,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
        endif
!
        if(right(1) == 0) then
           call mpi_wait(reqs_right(1), status_left, ierr)
        endif
        if(left(1) == 0) then
           call mpi_wait(reqs_left(1), status_left, ierr)
        endif
!
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!
        call mpi_barrier(lbecomm,ierr)
! --------------------------------------------------------------------------
! comms along z 
! --------------------------------------------------------------------------
! first isend 
        call time(tcountZ0)
        if(up(1) == 0) then
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZIN,a04,a06,a07,a13,a18,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do j = 0,m+1
!$acc loop independent
              do i = 0,l+1
                 bufferZIN(i,j,1)=a04(i,j,n)
                 bufferZIN(i,j,2)=a06(i,j,n)
                 bufferZIN(i,j,3)=a07(i,j,n)
                 bufferZIN(i,j,4)=a13(i,j,n)
                 bufferZIN(i,j,5)=a18(i,j,n)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 04
           call mpi_isend(bufferZIN(0,0,1), msgsizez, MYMPIREAL, up(2), tag, &
                             lbecomm, reqs_up(1), ierr)
        endif
!
! isend from down
        if(down(1) == 0) then
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZIN2,a02,a09,a11,a15,a16,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do j = 0,m+1
!$acc loop independent
              do i = 0,l+1
                 bufferZIN2(i,j,1)=a02(i,j,1)
                 bufferZIN2(i,j,2)=a09(i,j,1)
                 bufferZIN2(i,j,3)=a11(i,j,1)
                 bufferZIN2(i,j,4)=a15(i,j,1)
                 bufferZIN2(i,j,5)=a16(i,j,1)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
!
           tag = 02
           call mpi_isend(bufferZIN2(0,0,1), msgsizez, MYMPIREAL, down(2), tag, &
                             lbecomm, reqs_down(1), ierr)
        endif
!
! recv from up
        if(up(1) == 0) then
!
           tag = 02
           call mpi_recv(bufferZOUT2(0,0,1), msgsizez, MYMPIREAL, up(2), tag, &
                             lbecomm, status_down, ierr)
!!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZOUT2,a02,a09,a11,a15,a16,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do j = 0,m+1
!$acc loop independent
              do i = 0,l+1
                 a02(i,j,n+1) = bufferZOUT2(i,j,1)
                 a09(i,j,n+1) = bufferZOUT2(i,j,2)
                 a11(i,j,n+1) = bufferZOUT2(i,j,3)
                 a15(i,j,n+1) = bufferZOUT2(i,j,4)
                 a16(i,j,n+1) = bufferZOUT2(i,j,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL

        endif
!
! recv from down
        if(down(1) == 0) then
!
           tag = 04
           call mpi_recv(bufferZOUT(0,0,1), msgsizez, MYMPIREAL, down(2), tag, &
                             lbecomm, status_up, ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZOUT,a04,a06,a07,a13,a18,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
           do j = 0,m+1
!$acc loop independent
              do i = 0,l+1
                 a04(i,j,0)=bufferZOUT(i,j,1)
                 a06(i,j,0)=bufferZOUT(i,j,2)
                 a07(i,j,0)=bufferZOUT(i,j,3)
                 a13(i,j,0)=bufferZOUT(i,j,4)
                 a18(i,j,0)=bufferZOUT(i,j,5)
              enddo
           enddo
!$acc end kernels
!$OMP END PARALLEL
        endif

        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
#endif
!
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_mp = time_mp + real(countA1-countA0)/(count_rate)
        time_mp1 = time_mp1 + (tcountA1-tcountA0)
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
        endif
        if(down(1) == 0) then
           call mpi_wait(reqs_down(1), status_down, ierr)
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
           write(6,*) "DEBUG2: Exiting from sub. bcond_comm_noblock_packed_try"
        endif
#endif
!
        return
        end subroutine bcond_comm_noblock_packed_try

