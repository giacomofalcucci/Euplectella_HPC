!=====================================================================
!     ****** LBE/bcond_comm_packed
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
        subroutine bcond_comm_packed
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
        integer      :: msgsizeX
        integer      :: msgsizeY
        integer      :: msgsizeZ
!
        real(mykind), dimension(0:m+1,0:n+1,5) :: bufferXIN 
        real(mykind), dimension(0:m+1,0:n+1,5) :: bufferXOUT 
        real(mykind), dimension(0:l+1,0:n+1,5) :: bufferYIN 
        real(mykind), dimension(0:l+1,0:n+1,5) :: bufferYOUT 
        real(mykind), dimension(0:l+1,0:m+1,5) :: bufferZIN 
        real(mykind), dimension(0:l+1,0:m+1,5) :: bufferZOUT 
!
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
        msgsizeX = (n+2)*(m+2)*5
        msgsizeY = (l+2)*(n+2)*5
        msgsizeZ = (l+2)*(m+2)*5
!------------------------------------------------------------------------
! comms along z + 
!
        call time(tcountZ0)
        tag = 04
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZIN,a04,a06,a07,a13,a18,l,m,n) 
!$OMP DO
!$acc kernels
        do j = 0,m+1
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
        call mpi_sendrecv(bufferZIN(0,0,1),msgsizeZ,MYMPIREAL,up(2),tag,&
                          bufferZOUT(0,0,1),msgsizeZ,MYMPIREAL,down(2),tag,&
                          lbecomm,status,ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZOUT,a04,a06,a07,a13,a18,l,m,n) 
!$OMP DO
!$acc kernels
        do j = 0,m+1
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
!
!!!        call mpi_barrier(lbecomm,ierr)
! comms along z - 
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZIN,a02,a09,a11,a15,a16,l,m,n) 
!$OMP DO
!$acc kernels
        do j = 0,m+1
           do i = 0,l+1
              bufferZIN(i,j,1)=a02(i,j,1) 
              bufferZIN(i,j,2)=a09(i,j,1)
              bufferZIN(i,j,3)=a11(i,j,1)
              bufferZIN(i,j,4)=a15(i,j,1)
              bufferZIN(i,j,5)=a16(i,j,1)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
!
        call mpi_sendrecv(bufferZIN(0,0,1),msgsizeZ,MYMPIREAL,down(2),tag,&
                          bufferZOUT(0,0,1),msgsizeZ,MYMPIREAL,up(2),tag,&
                          lbecomm,status,ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,i)  &
!$OMP SHARED(bufferZOUT,a02,a09,a11,a15,a16,l,m,n) 
!$OMP DO
!$acc kernels
        do j = 0,m+1
           do i = 0,l+1
              a02(i,j,n+1) = bufferZOUT(i,j,1)
              a09(i,j,n+1) = bufferZOUT(i,j,2)
              a11(i,j,n+1) = bufferZOUT(i,j,3)
              a15(i,j,n+1) = bufferZOUT(i,j,4)
              a16(i,j,n+1) = bufferZOUT(i,j,5)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
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
!
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXIN,a01,a02,a03,a04,a05,l,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
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
        call mpi_sendrecv(bufferXIN(0,0,1),msgsizeX,MYMPIREAL,front(2),tag,&
                          bufferXOUT(0,0,1),msgsizeX,MYMPIREAL,rear(2),tag,&
                          lbecomm,status,ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXOUT,a01,a02,a03,a04,a05,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
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
!
!!!        call mpi_barrier(lbecomm,ierr)
!
! comms along x - 
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXIN,a10,a11,a12,a13,a14,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
           do j = 0,m+1
              bufferXIN(j,k,1)=a10(1,j,k)
              bufferXIN(j,k,2)=a11(1,j,k)
              bufferXIN(j,k,3)=a12(1,j,k)
              bufferXIN(j,k,4)=a13(1,j,k)
              bufferXIN(j,k,5)=a14(1,j,k)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
!
        call mpi_sendrecv(bufferXIN(0,0,1),msgsizeX,MYMPIREAL,rear(2),tag,&
                          bufferXOUT(0,0,1),msgsizeX,MYMPIREAL,front(2),tag,&
                          lbecomm,status,ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(j,k)  &
!$OMP SHARED(bufferXOUT,a10,a11,a12,a13,a14,m,n,l) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
           do j = 0,m+1
              a10(l+1,j,k) = bufferXOUT(j,k,1)
              a11(l+1,j,k) = bufferXOUT(j,k,2)
              a12(l+1,j,k) = bufferXOUT(j,k,3)
              a13(l+1,j,k) = bufferXOUT(j,k,4)
              a14(l+1,j,k) = bufferXOUT(j,k,5)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
!
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
!!!        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along y + 
        call time(tcountY0)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYIN,a03,a07,a08,a09,a12,l,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
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
        tag = 3
        call mpi_sendrecv(bufferYIN(0,0,1),msgsizeY,MYMPIREAL,right(2),tag,&
                          bufferYOUT(0,0,1),msgsizeY,MYMPIREAL,left(2),tag,&
                          lbecomm,status,ierr)
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYOUT,a03,a07,a08,a09,a12,l,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
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
!!!        call mpi_barrier(lbecomm,ierr)
!
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYIN,a01,a10,a16,a17,a18,l,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
           do i = 0,l+1
              bufferYIN(i,k,1)=a01(i,1,k)
              bufferYIN(i,k,2)=a10(i,1,k)
              bufferYIN(i,k,3)=a16(i,1,k)
              bufferYIN(i,k,4)=a17(i,1,k)
              bufferYIN(i,k,5)=a18(i,1,k)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
!
        tag = 3
        call mpi_sendrecv(bufferYIN(0,0,1),msgsizeY,MYMPIREAL,left(2),tag,&
                          bufferYOUT(0,0,1),msgsizeY,MYMPIREAL,right(2),tag,&
                          lbecomm,status,ierr)
!
!
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,k)  &
!$OMP SHARED(bufferYOUT,a01,a10,a16,a17,a18,l,m,n) 
!$OMP DO
!$acc kernels
        do k = 0,n+1
           do i = 0,l+1
              a01(i,m+1,k)=bufferYOUT(i,k,1)
              a10(i,m+1,k)=bufferYOUT(i,k,2)
              a16(i,m+1,k)=bufferYOUT(i,k,3)
              a17(i,m+1,k)=bufferYOUT(i,k,4)
              a18(i,m+1,k)=bufferYOUT(i,k,5)
           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL
!
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!
        call mpi_barrier(lbecomm,ierr)
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
           write(6,*) "DEBUG2: Exiting from sub. bcond_comm_packed"
        endif
#endif
        return
        end subroutine bcond_comm_packed
