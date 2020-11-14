!=====================================================================
!     ****** LBE/setup
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different initializations routines..
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
      subroutine setup(itfin,ivtim,isignal,itsave,icheck,itrestart,init_v,tstep)
!
      use storage
      use timing
! 
#ifdef OPENACC
      use openacc
#endif
!
      implicit none
!
      INTEGER:: itfin, itstart, ivtim
      INTEGER:: itime, itsave, icheck, itrestart, init_v
      INTEGER:: isignal
      INTEGER:: tstep
      INTEGER:: ierr
      INTEGER:: provided
      INTEGER:: idev
      character*19 file_name1, file_name2, file_name3
      character*23 file_name4
      character*15 file_name5
      character*19 file_name6
      character*21 file_name7
      character*21 file_name8		! drag
!
#ifdef SERIAL
# ifdef OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
!
!     forcing to run on a particular gpu
!
      idev = 0
      call acc_set_device_num(idev,acc_device_nvidia)
      write(6,*) "INFO --> acc_device_nvidia", ndev, idev
# endif
#endif
!
#ifdef SERIAL
! set values for serial version...
      myrank = 0
      nprocs = 0
      mpicoords(1) = 0
      mpicoords(2) = 0
      mpicoords(3) = 0
!
      call system("date       > time.log")
!
#else
!
# ifdef MEM_CHECK
      if(myrank == 0) then
        mem_stop = get_mem();
        write(6,*) "MEM_CHECK: before setup mpi stuff mem =", mem_stop
      endif
# endif
!
! setup MPI
!      call mpi_init(ierr)
!      call mpi_init_thread(MPI_THREAD_SINGLE,provided,ierr)
#ifdef _OPENMP
      call mpi_init_thread(MPI_THREAD_SINGLE,provided,ierr)
!      call mpi_init_thread(MPI_THREAD_FUNNELED,provided,ierr)
!      call mpi_init_thread(MPI_THREAD_SERIALIZED,provided,ierr)
!      call mpi_init_thread(MPI_THREAD_MULTIPLE,provided,ierr)
!
! check
!!      if (provided < MPI_THREAD_MULTIPLE) then
!!        write(6,*) "ERROR: --> MPI_THREAD_MULTIPLE"
!!        call MPI_Abort(MPI_COMM_WORLD, 1,ierr)
!!      endif
#else
      call mpi_init(ierr)
#endif
      call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
!
      if(myrank==0) then
         call system("date       > time.log")
      endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after setup mpi stuff mem =", mem_stop
      endif
# endif
!
# ifdef DEBUG_1
      write(*,*) "DEBUG1: I'm of procs=",nprocs ," myrank=",myrank
# endif
!
#endif
!
! probe
      file_name3 = 'probe.xxxx.dat'
      write(file_name3(7:10),4000) myrank
!
! probe_g
      file_name7 = 'probe_g.xxxx.dat'
      write(file_name7(9:12),4000) myrank
!
! drag
      file_name8 = 'drag.xxxx.dat'
      write(file_name8(6:9),4000) myrank
!
! prof_k
      file_name1 = 'prof_k.xxxx.dat'
      write(file_name1(8:11),4000) myrank
!
! prof_j
      file_name6 = 'prof_j.xxxx.dat'
      write(file_name6(8:11),4000) myrank
!
! prof_i
      file_name2 = 'prof_i.xxxx.dat'
      write(file_name2(8:11),4000) myrank
!
! task
      file_name5 = 'task.xxxxxx.log'
      write(file_name5(6:11),3100) myrank
!
      if(myrank == 0) then 
         open(16,file='bgk.log',  status='unknown')
#ifdef _OPENMP
        if(myrank == 0) then
           write(6,*)  "INFO --> provided", provided
           write(16,*) "INFO --> provided", provided
        endif
#endif
      endif
!
#ifdef NO_OUTPUT
! do nothing
#else
      open(60,file=file_name1, status='unknown')        ! prof_k
      open(61,file=file_name2, status='unknown')        ! prof_i
      open(64,file=file_name6, status='unknown')        ! prof_j
      open(66,file=file_name8, status='unknown')        ! drag
      open(67,file=file_name7, status='unknown')        ! probe_gb
      open(68,file=file_name3, status='unknown')        ! probe
      open(38,file=file_name5, status='unknown')	! task.XXXXXX.log
#endif
!
      if(myrank==0) then
         open(63,file='diagno.dat',status='unknown')
      endif
!
      if(myrank==0) then 
         write(6,*) "========================="
         write(6,*) " version 0: uniform LBE  "
         write(6,*) "========================="
      endif
!
! check, read input data and data allocation
!
      call check_case
!
      call input(itfin,ivtim,isignal,itsave,icheck, & 
                   itrestart,init_v,tstep)
!
#ifdef SERIAL
! do nothing...
#else
      call setup_mpi
#endif
!
#ifdef MB
      file_name4 = 'u_med.xxxx.dat'
      write(file_name4(7:10),4000) myrank
#else
      file_name4 = 'u_med.xxx.xxx.xxx.dat'
      write(file_name4(7:9),4100) mpicoords(1)
      write(file_name4(11:13),4100) mpicoords(2)
      write(file_name4(15:17),4100) mpicoords(3)
#endif

#ifdef NO_OUTPUT
! do nothing...
#else
      open(62,file=file_name4, status='unknown')
#endif
!
      call alloca

#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. setup"
      endif
#endif

# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. setup mem =", mem_stop
      endif
# endif
!
! format
4000  format(i4.4)
4100  format(i3.3)
3100  format(i6.6)
!
      return
      end subroutine setup
