! =====================================================================
!     ****** LBE/input
!
!     COPYRIGHT
!       (c) 2000-2008 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       input
!     DESCRIPTION
!       read input parameters
!       read from unit 15 (bgk.input)
!     INPUTS
!       none
!     OUTPUT
!       size     --> size_x, size_y
!       itfin    --> end of the run
!       ivtim    --> interval between two different visualization
!       isignal  --> interval between two timings 
!       itsave   --> time of saving
!       icheck   --> interval between two different diagnostic
!       irestart --> restart from file (1=yes,other=no)
!       init_v   --> initial velocity condition (see subroutine init)
!     TODO
!	
!     NOTES
!       integer variables used: itfin,ivtim,isignal,itsave,icheck,irestart
!                               init_v, tstep
!       real variables defined: dt (not used)
!
!     *****
! =====================================================================
!
      subroutine input (itfin,ivtim,isignal,itsave,icheck,irestart, & 
                          init_v,tstep)
!
      use storage
      use timing
! 
      implicit none
!
      real(mykind):: dt             ! not yet used
!
      integer:: itfin,ivtim,isignal,itsave,icheck
      integer:: irestart,init_v, tstep
!
      namelist /parameters/ svisc, u0, itfin, ivtim, isignal, & 
                            itsave, icheck, irestart, init_v, &
                            lx, ly, lz, proc_x, proc_y, proc_z, tstep, &
                            flag1, flag2, flag3, ipad, jpad, kpad
!
! default values
      flag1 = 2         ! creating obstacles  (1-file/2-creating)
      flag2 = 2         ! obstacles           (1-sphere/2-cilinder)
      flag3 = 1 
!
      ipad  = 0         ! no memory padding (x)
      jpad  = 0         ! no memory padding (y)
      kpad  = 0         ! no memory padding (z)
!
      open(15,FILE='bgk.input',STATUS='old')
      read(15,parameters)
      close(15)

#ifdef SERIAL
      proc_x = 1
      proc_y = 1
      proc_z = 1
      nprocs = 1
      l = lx
      m = ly
      n = lz
#else
# ifdef MB
#  ifdef DEBUG_1
      write(6,*) "INFO: reading topology file"
#  endif
      write(16,*) "INFO: reading topology file"
      call read_topology
# else
#  ifdef DEBUG_1
      write(6,*) "INFO: simple decomposition"
#  endif
      l = lx/proc_x
      m = ly/proc_y
      n = lz/proc_z
!
! some check
      if((l*proc_x).NE.lx) then
         write(6,*) "ERROR: global and local size along x not valid!!", lx, l, proc_x
         stop
      endif
!
      if((m*proc_y).NE.ly) then
         write(6,*) "ERROR: global and local size along y not valid!!", ly, m, proc_y
         stop
      endif
!
      if((n*proc_z).NE.lz) then
         write(6,*) "ERROR: global and local size along z not valid!!", lz, n, proc_z
         stop
      endif
# endif
#endif

#ifdef TRYBOX
       if((border.ge.l/2).OR.(border.ge.m/2).OR.(border.ge.n/2)) then
           write( 6,*) "ERROR: border too big", border, l/2, m/2, n/2
           stop
       endif
#endif

! set tstep to 1... (to remove)
      tstep = unoi
!
      if(myrank == 0) then 
           write(16,*) "INFO: tstep --->", tstep
      endif
!
! some fix..
      l1 = l+1
      m1 = m+1
      n1 = n+1
!
! default value: volume forcing along x
!
      u0x = 0.0
      u0y = 0.0
      u0z = 0.0
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. input"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. input mem =", mem_stop
      endif
# endif
!
      return
      end subroutine input
