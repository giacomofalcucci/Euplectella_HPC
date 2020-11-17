!=====================================================================
!     ****** LBE/build_bcond
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       set BC flags (for multiblock & mpi)
!       if flag == 0 periodic b.c.
!               == 1 solid wall
!               == 2 moving wall
!               == 3 no-slip wall
!               == 4 inflow
!               == 5 outflow
!               == other --> error
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
      subroutine build_bcond
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
      integer:: uni, len, ierr, i
!
      character*15 hname
!
#ifdef MB
! nothing to do
#else
! set default 
      up(1)   = 0
      down(1) = 0
      right(1)= 0
      left(1) = 0
      front(1)= 0
      rear(1) = 0
!
# ifdef DRIVEN
!
      u00 = u0             ! boundary condition
      u0 = 0.0             ! volume force
!
      if(mpicoords(3)==(proc_z-1)) then 
         up(1)   = 2
      endif
      if(mpicoords(3)==0) then 
         down(1) = 1
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 1
      endif
      if(mpicoords(2)==0) then
         left(1) = 1
      endif
      if(mpicoords(1)==(proc_x-1)) then 
         front(1)= 1
      endif
      if(mpicoords(1)==0) then 
         rear(1) = 1
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> DRIVEN", u0, u00
      endif
!
# else
#  ifdef COUETTE  
      u00 = u0             ! boundary condition
      u0 = 0.0             ! volume force
      if(mpicoords(3)==(proc_z-1)) then 
         up(1)   = 2
      endif
      if(mpicoords(3)==0) then 
         down(1) = 1
      endif
      right(1)= 0
      left(1) = 0
      front(1)= 0
      rear(1) = 0
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> COUETTE", u0, u00
      endif
!
#  else
#   ifdef INLET
!
      u_inflow = u0
      u0 = 0
!
      if(mpicoords(3)==(proc_z-1)) then 
         up(1)   = 0
      endif
      if(mpicoords(3)==0) then 
         down(1) = 0
      endif
      if(mpicoords(1)==(proc_x-1)) then 
         front(1)= 5
      endif
      if(mpicoords(1)==0) then 
         rear(1) = 4
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> INLET", u0, u_inflow
      endif
!
#   else
#    ifdef PERIODIC
      up(1)   = 0
      down(1) = 0
      right(1)= 0
      left(1) = 0
      front(1)= 0
      rear(1) = 0
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> PERIODIC"
      endif
!
#    else 
      if(mpicoords(3)==(proc_z-1)) then  ! i.e. Channel
         up(1)   = 1
#     ifdef DEBUG_1
         if(myrank == 0) then
            write(6,*) "DEBUG1: --> UP", myrank, up(1)
         endif
#     endif
      endif
      if(mpicoords(3)==0) then 
         down(1) = 1
#     ifdef DEBUG_1
         if(myrank == 0) then
            write(6,*) "DEBUG1: --> DOWN", myrank, down(1)
         endif
#     endif
      endif
      right(1)= 0
      left(1) = 0
      front(1)= 0
      rear(1) = 0
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> CHANNEL"
      endif
!
#    endif
#   endif
#  endif
# endif

# ifdef CHECK2D
!
      u_inflow = u0
      u0 = 0
!
      if(mpicoords(3)==(proc_z-1)) then
         up(1)   = 0
      endif
      if(mpicoords(3)==0) then
         down(1) = 0
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 1
      endif
      if(mpicoords(2)==0) then
         left(1) = 1
      endif
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         rear(1) = 4
      endif
!
# endif
!
# ifdef DUCT
!
      front(1)= 0
      rear(1) = 0
!
      u_inflow= u0
      u0 = 0
!
      if(mpicoords(3)==(proc_z-1)) then
         up(1)   = 1
      endif
      if(mpicoords(3)==0) then
         down(1) = 1
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 1
      endif
      if(mpicoords(2)==0) then
         left(1) = 1
      endif
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         rear(1) = 4
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> DUCT", u0, u_inflow
      endif
!
# endif 
!
!
# ifdef MIXING
!
      front(1)= 0
      rear(1) = 0
!
      if(mpicoords(3)==(proc_z-1)) then
         up(1)   = 1
      endif
      if(mpicoords(3)==0) then
         down(1) = 1
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 0
      endif
      if(mpicoords(2)==0) then
         left(1) = 0
      endif
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         u_inflow= u0
         rear(1) = 4
      endif
# endif 
!
# ifdef SPONGES
!
      u_inflow= u0
      u0 = 0
!
      if(mpicoords(3)==(proc_z-1)) then
         up(1)   = 3
      endif
      if(mpicoords(3)==0) then
         down(1) = 1
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 0
      endif
      if(mpicoords(2)==0) then
         left(1) = 0
      endif
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         rear(1) = 4
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> for SPONGES"
      endif
!
# endif
!
# ifdef FREESLIP
!
      u_inflow= u0
      u0 = 0
!
      if(mpicoords(3)==(proc_z-1)) then
         up(1)   = 3
      endif
      if(mpicoords(3)==0) then
         down(1) = 3
      endif
      if(mpicoords(2)==(proc_y-1)) then
         right(1)= 3
      endif
      if(mpicoords(2)==0) then
         left(1) = 3
      endif
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         u_inflow= u0
         rear(1) = 4
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> for FREESLIP"
      endif
!
# endif
! closing MB ifdef...
#endif 
!
#ifdef CYLINDER_INF
!
      u_inflow= u0
      u0 = 0
!
      up(1)   = 0
      down(1) = 0
      right(1)= 0
      left(1) = 0
      if(mpicoords(1)==(proc_x-1)) then
         front(1)= 5
      endif
      if(mpicoords(1)==0) then
         rear(1) = 4
      endif
!
      if(myrank == 0) then
         write(16,*) "INFO: build bcond --> for CYLINDER_INF"
      endif
!
#endif
!
! log
#ifdef SERIAL
! do nothing
#else
# ifdef FERMI_FIX
      hname = " --> IBM BG/Q"
# else
      call MPI_GET_PROCESSOR_NAME(hname, len,ierr)
# endif
#endif
       
#ifdef NO_OUTPUT
! do nothing
#else
# ifdef SERIAL
! do nothing
# else
      write(38,*) "#", myrank, ":my host is ------>", hname
#  ifdef OPENACC
      write(38,*) "#", myrank, ":my GPU  is ------>", mydev, ndev
#  endif
      write(38,*) "#", myrank, ":my rank is ------>", myrank
      write(38,*) "#", myrank, ":my x mpicoord --->", mpicoords(1)
      write(38,*) "#", myrank, ":my y mpicoord --->", mpicoords(2)
      write(38,*) "#", myrank, ":my z mpicoord --->", mpicoords(3)
# endif
#  ifdef OPENACC
      write(38,*) "#", myrank, ":my GPU  is ------>", mydev, ndev
#  endif
      write(38,*) "#", myrank, ":my size (x) is -->", l
      write(38,*) "#", myrank, ":my size (y) is -->", m
      write(38,*) "#", myrank, ":my size (z) is -->", n
      write(38,*) "#", myrank, ":my up BC is ----->", up(1), up(2) 
      write(38,*) "#", myrank, ":my down BC is --->", down(1), down(2)
      write(38,*) "#", myrank, ":my right BC is -->", right(1), right(2)
      write(38,*) "#", myrank, ":my left BC is --->", left(1), left(2)
      write(38,*) "#", myrank, ":my front BC is -->", front(1), front(2)
      write(38,*) "#", myrank, ":my rear BC is --->", rear(1), rear(2)
      write(38,*) "#", myrank, ":Memory (start) -->", mem_start
!
! write info for topology
      write(38,*) offset(1)  , offset(2)  , offset(3)
      write(38,*) offset(1)  , offset(2)  , offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)+l, offset(2)  , offset(3)
      write(38,*) offset(1)+l, offset(2)  , offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)+l, offset(2)+m, offset(3)
      write(38,*) offset(1)+l, offset(2)+m, offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)  , offset(2)+m, offset(3)
      write(38,*) offset(1)  , offset(2)+m, offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)  , offset(2)  , offset(3)
      write(38,*) offset(1)  , offset(2)  , offset(3)+n
      write(38,*) " "
#endif
!
! write topo.input back-up file...
!
#ifdef SERIAL
! do nothing....
#else
      if(myrank == 0) then
         open(27,file='topo.input.bak',  status='unknown')
         write(27,*) nprocs
         write(27,*) lx,ly,lz
         close(27)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
      do i = 0, nprocs-1
         if(myrank == i) then
            open(27,file='topo.input.bak',  access='append')
            write(27,*) myrank
            write(27,*) l,m,n,start_idx(1),start_idx(2),start_idx(3)
            write(27,*) front(1), rear(1), right(1), left(1), up(1), down(1)
            write(27,*) front(2), rear(2), right(2), left(2), up(2), down(2)
            close(27)
         endif
         call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo
!
      if(myrank == 0) then
         open(27,file='topo.input.bak',  access='append')
         write(27,*) 666
         close(27)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
!       
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. build_bcond"
      endif
#endif
      return
      end subroutine build_bcond
