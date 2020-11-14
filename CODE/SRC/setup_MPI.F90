!=====================================================================
!     ****** LBE/setup_mpi
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different mpi initializations 
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
      subroutine setup_mpi
!
      use storage
      use timing
!
#ifdef SERIAL
!do nothing
#else
      use mpi
#endif
!
#ifdef OPENACC
      use openacc
#endif
!
      implicit none
!
      integer:: i, uni
      integer:: ierr                    ! mpi variables
!
      real(mykind):: knorm
!
#ifdef SERIAL
! do nothing
# ifdef OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
      call acc_set_device_num(mydev,acc_device_nvidia)
      write(6,*) "INFO: using GPU",mydev, ndev
# endif
#else
!
# ifdef OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
!
! check 1 (max task 4)
!      if(nprocs.gt.4) then
!         write(6,*) "ERROR: this version works with only 4 task", ndev
!         call MPI_Abort(MPI_COMM_WORLD, 1,ierr)
!      endif
!
! check 2 (every node has 4 gpus)   --> Marconi100
! check 2 (every node has 2 gpus)   --> Galileo K80
! check 2 (every node has 1 gpus)   --> Galileo V100
!
      if(ndev.ne.4) then
         write(6,*) "ERROR: ndev.ne.4", ndev
         call MPI_Abort(MPI_COMM_WORLD, 1,ierr)
      else
!
! set the gpu to the task id
!
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, localcomm, ierr)
        call MPI_Comm_rank(localcomm, mydev, ierr)
        call acc_set_device_num(mydev,acc_device_nvidia)
        write(6,*) "INFO:  MPI rank",myrank," using GPU",mydev
!         idev = myrank
!        call acc_set_device_num(idev,acc_device_nvidia)
!       write(6,*) "INFO --> myrank, mygpu", myrank, idev
      endif
# endif
!
# ifdef DOUBLE_P
       knorm = 8.0/1024.0
# else
       knorm = 4.0/1024.0
# endif
!
! check
      if ((proc_x*proc_y*proc_z).ne.nprocs) then
        if (myrank.eq.0) then
          write(*,*) 'ERROR: decomposed for', & 
                             proc_x*proc_y*proc_z, 'procs'
          write(*,*) 'ERROR: launched on', nprocs, 'processes'
        end if
        call MPI_finalize(ierr)
        stop
      end if
!
      rreorder=.false.
!
      periodic(1) = .true.
      periodic(2) = .true.
      periodic(3) = .true.
!
      prgrid(1) = proc_x
      prgrid(2) = proc_y
      prgrid(3) = proc_z
!
!        write(*,*) proc_x, proc_z
!
! building virtual topology
      call MPI_cart_create(mpi_comm_world, mpid, prgrid, &
                              periodic,rreorder,lbecomm,ierr)

      call MPI_comm_rank(lbecomm, myrank, ierr)
      call MPI_cart_coords(lbecomm, myrank, mpid, &
                            mpicoords, ierr)
!
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
# ifdef MB
! nothing to do...
# else
! x dir  & y dir
      call MPI_cart_shift(lbecomm, 0, 1, rear(2), front(2), ierr)
      call MPI_cart_shift(lbecomm, 1, 1, left(2), right(2), ierr)
      call MPI_cart_shift(lbecomm, 2, 1, down(2), up(2), ierr)
# endif
!
! yz plane is composed by single point (stride.ne.1)
      call MPI_type_vector((n+2)*(m+2), 1, l+2, MYMPIREAL, yzplane, ierr)
      call MPI_type_commit(yzplane,ierr)
      if(myrank.eq.0) then
         write(16,*) "INFO: yzplane (KB)-->", (n+2)*(m+2)*knorm
      endif
!
! xz plane is composed by single arrays (stride.ne.1)
      call MPI_type_vector(n+2, l+2, (m+2)*(l+2), MYMPIREAL, xzplane, ierr)
      call MPI_type_commit(xzplane,ierr)
      if(myrank.eq.0) then
         write(16,*) "INFO: xzplane (KB)-->", (n+2)*(l+2)*knorm
      endif
!
! xy plane is contiguous arrays (stride.eq.1)
      call MPI_type_contiguous((l+2)*(m+2), MYMPIREAL, xyplane, ierr)
      call MPI_type_commit(xyplane,ierr)
      if(myrank.eq.0) then
         write(16,*) "INFO: xyplane (KB)-->", (m+2)*(l+2)*knorm
      endif
!
! create subarray (for save-restore)
!
! setting some offsets...
      gsizes(1) = lx   ! global size 
      gsizes(2) = ly
      gsizes(3) = lz
!
      lsizes(1) = l       ! local size 
      lsizes(2) = m
      lsizes(3) = n
!
# ifdef MB
! nothing to do (already defined)
# else
      start_idx(1) = mpicoords(1)*l   ! offset 
      start_idx(2) = mpicoords(2)*m
      start_idx(3) = mpicoords(3)*n
# endif
!
      offset(1) = start_idx(1)   ! offset 
      offset(2) = start_idx(2)
      offset(3) = start_idx(3)
!
!      file_offset = myrank*l*n*8
      file_offset = 0    !to check
!
      buffer_size = l*m*n
!
! 
# ifdef MB
! 3dump is not possible....
# else
      call MPI_Type_create_subarray(mpid, gsizes, lsizes, start_idx, &
           MPI_ORDER_FORTRAN, MYMPIREAL, dump3d,ierr)
      call MPI_Type_commit(dump3d,ierr)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
# endif
#endif

#ifdef DEBUG_1
# ifdef MPIIO
      write(6,*)  'DEBUG1: task ', myrank, file_offset
      write(16,*) 'DEBUG1: task ', myrank, start_idx(1)
      write(16,*) 'DEBUG1: task ', myrank, start_idx(2)
      write(16,*) 'DEBUG1: task ', myrank, start_idx(3)
# endif
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. setup_mpi"
      endif
#endif

# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. setup_MPI mem =", mem_stop
      endif
# endif
!
      return
      end subroutine setup_mpi
