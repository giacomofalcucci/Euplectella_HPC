!=====================================================================
!     ****** LBE/check_isend
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       check_isend
!     DESCRIPTION
!       Simple wrapper for different BCs
!       - channel (default)
!       - couette (if -DCOUETTE)
!       - driven  (if -DDRIVEN)
!		- first order
!		- second order (bounce-back)
!       - inlet   (if -DINLET)
!	- bi-periodic   (if -DPERIODIC)
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
        subroutine check_isend
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
#ifdef SERIAL
! do nothing
#else
        integer      :: reqs_up
        integer      :: reqs_down
        integer      :: in1, in2, out1, out2
        integer      :: tag, ierr
        integer      :: status(MPI_STATUS_SIZE)
!
! set some values..

        in1 = myrank+100
        in2 = myrank-100
        out1 = -1
        out2 = -1
#endif
!
#ifdef SERIAL
!do nothing
#else
!
        write(6,*) "--> pre-barrier", myrank, in1, in2, out1, out2
        call mpi_barrier(lbecomm,ierr)
        write(6,*) "--> post-barrier", myrank, in1, in2, out1, out2
!
! first irecv from up
        if(up(1) == 0) then
            write(6,*) "--> pre isend-up", myrank
            tag = 04
            call mpi_isend(in1, 1, MPI_INTEGER, up(2), tag, &
                             lbecomm, reqs_up, ierr)
            write(6,*) "--> post isend-up", myrank
        endif
!
! second send from down
        if(down(1) == 0) then
            write(6,*) "--> pre recv-down", myrank
            tag = 04
            call mpi_recv(out1, 1, MPI_INTEGER, down(2), tag, &
                             lbecomm, status, ierr)
            write(6,*) "--> post recv-down", myrank
        endif
!
! first irecv from down
        if(down(1) == 0) then
            tag = 05
            write(6,*) "--> pre isend-down", myrank
            call mpi_isend(in2, 1, MPI_INTEGER, down(2), tag, &
                             lbecomm, reqs_down, ierr)
            write(6,*) "--> post isend-down", myrank
        endif
!
! second send from up
        if(up(1) == 0) then
            write(6,*) "--> pre recv-up", myrank
            tag = 05
            call mpi_recv(out2, 1, MPI_INTEGER, up(2), tag, &
                             lbecomm, status, ierr)
            write(6,*) "--> post recv-up", myrank
        endif

        call mpi_barrier(lbecomm,ierr)
        if(up(1) == 0) then
           write(6,*) "--> pre waitall-up", myrank
           call mpi_wait(reqs_up, status, ierr)
           write(6,*) "--> post waitall-up", myrank
        endif
        if(down(1) == 0) then
           write(6,*) "--> pre waitall-down", myrank
           call mpi_wait(reqs_down, status, ierr)
           write(6,*) "--> post waitall-up", myrank
        endif
!
#endif
!
#ifdef SERIAL
!do nothing
#else
        call mpi_barrier(lbecomm,ierr)
!!        write(6,*) "out -->", myrank, in1, in2, out1, out2
!!        call mpi_finalize(ierr)
!!        stop
#endif

        end subroutine check_isend
