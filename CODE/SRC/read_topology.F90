! =====================================================================
!     ****** LBE/read_topology
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       input
!     DESCRIPTION
!       read topology information
!       read from unit 41 (topo.input)
!     INPUTS
!       none
!     OUTPUT
!
!     TODO
!	
!     NOTES
!       integer variables used: 
!       real variables defined: 
!
!     *****
! =====================================================================
!
      subroutine read_topology()
!
#ifdef SERIAL
!do nothing
#else
      use mpi
      use storage
      implicit none
!
      real(mykind):: dt             ! not yet used
!
      integer:: ios, i, check, ierr, tag, tag1, icode
      integer:: nblock, blockid, taskid
      integer:: status(MPI_STATUS_SIZE)
!
! temporary integers..
      integer:: temp(18)
!
      tag = 1
      if(myrank == 0) then
        open(41,FILE='topo.input',STATUS='old',IOSTAT=ios)
        if(ios == 0) then
          read(41,*) nblock
          if(nblock == nprocs) then
            read(41,*) lx, ly, lz
            do i = 0, nprocs-1
              read(41,*) taskid
              if(taskid == 0) then 
                read(41,*) l,m,n,start_idx(1),start_idx(2),start_idx(3)
                read(41,*) front(1), rear(1), right(1), left(1), up(1), down(1)
                read(41,*) front(2), rear(2), right(2), left(2), up(2), down(2)
              else
                write(6,*) "INFO: I am task", myrank, "and I send to", taskid
                read(41,*) temp( 1),temp( 2),temp( 3),temp( 4),temp( 5),temp( 6)
                read(41,*) temp( 7),temp( 8),temp( 9),temp(10),temp(11),temp(12)
                read(41,*) temp(13),temp(14),temp(15),temp(16),temp(17),temp(18)
!
                tag = taskid
                call mpi_send(temp(1), 18, MPI_INTEGER, taskid, tag, MPI_COMM_WORLD, ierr)
                write(6,*) "INFO: I am task", myrank, "and I have sent to", taskid
              endif
            enddo
            read(41,*) check
            if(check == 666) then
              write(6,*) "INFO: topology input file OK", myrank
            else
              write(6,*) "ERROR: topology input file not correct!", myrank
              call MPI_finalize(ierr)
            endif
          else
            write(6,*) "ERROR: #blocks differs from mpi tasks!!"
            call MPI_abort(MPI_COMM_WORLD,icode,ierr)
          endif
        else
          write(6,*) "ERROR: topology file doesn't exist..."
          call MPI_abort(MPI_COMM_WORLD,icode,ierr)
        endif
        close(41)
      else
        write(6,*) "INFO: I am task", myrank, "and I'm waiting"
        tag1 = myrank
        call mpi_recv(temp(1), 18, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
        write(6,*) "INFO: I am task", myrank, "and I have received"
!
        l           = temp(1) 
        m           = temp(2) 
        n           = temp(3) 
        start_idx(1)= temp(4) 
        start_idx(2)= temp(5) 
        start_idx(3)= temp(6) 
        front(1)    = temp(7)
        rear(1)     = temp(8)
        right(1)    = temp(9)
        left(1)     = temp(10)
        up(1)       = temp(11)
        down(1)     = temp(12)
        front(2)    = temp(13)
        rear(2)     = temp(14)
        right(2)    = temp(15)
        left(2)     = temp(16)
        up(2)       = temp(17)
        down(2)     = temp(18)

      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
# ifdef DEBUG_1
       if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. read_topology"
       endif
# endif
#endif
!
       return
       end subroutine read_topology
