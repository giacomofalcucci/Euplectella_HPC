! =====================================================================
!     ****** LBE/save
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save
!     DESCRIPTION
!       Simple wrapper for saving in:
!        - raw
!        - hdf5
!     INPUTS
!       itime --> timestep
!     OUTPUT
!
!     TODO
!	
!     NOTES
!
!     *****
! =====================================================================
!
      subroutine save (itime)
!
      use storage
      use timing
      use real_kinds
!
      implicit none
!
      integer:: itime, ierr
!
#ifdef SERIAL
      call system("date")
#else
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(myrank == 0) then
         call system("date")
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
      call time(tcountF0)
      call SYSTEM_CLOCK(countF0, count_rate, count_max)
!
#ifdef HDF5
# ifdef SERIAL
      call save_hdf5_serial(itime)
# else
      call save_hdf5_parallel(itime)
# endif
#else
# ifdef MPIIO
      call save_mpiio(itime)
# else
      call save_raw(itime)
# endif
#endif
!
      call SYSTEM_CLOCK(countF1, count_rate, count_max)
      call time(tcountF1)
      time_io = time_io + real(countF1-countF0)/(count_rate)
      time_io1 = time_io1 + tcountF1-tcountF0
!
#ifdef SERIAL
! do nothing
#else
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(myrank == 0) then
         call system("date")
         write(6,*) "DEBUG1: I/O time (1)", real(countF1-countF0)/(count_rate)
         write(6,*) "DEBUG1: I/O time (2)", tcountF1-tcountF0 
         write(6,*) "DEBUG1: Exiting from sub. save"
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. save mem =", mem_stop
      endif
# endif
      return
      end subroutine save
