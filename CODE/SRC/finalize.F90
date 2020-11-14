!=====================================================================
!     ****** LBE/finalize
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different finalizations..
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       the MLUPS value is correct for equal-size bloks, else to fix
!
!     *****
!=====================================================================
!
      subroutine finalize(itstart,itfin)
!
      use storage
      use timing
!
      implicit none
!
      integer:: itstart, itfin
      integer:: ierr
!
      real(mykind):: knorm
!
#ifdef SERIAL
!do nothing
#else
# ifdef DOUBLE_P
       knorm = 8.0/1024.0
# else
       knorm = 4.0/1024.0
# endif
#endif


#ifdef PGI
#else
# ifdef FUGAKU
# else
      mem_stop = get_mem();
# endif
#endif
!!!!      call dealloca()
!
! write on external file some performance figure...
!
      if(myrank == 0) then 
         open(69,file='bgk.perf',  status='unknown')
         write(69,9999) 
         write(69,*)  "# Run info "
         write(69,*)  lx, ly, lz
         write(69,*)  l, m, n
         write(69,9999) 
         write(69,*)  "# Time for section "
         write(69,1101) time_loop, time_loop1
         write(69,1102) time_coll, time_coll1
#ifdef ORIGINAL
         write(69,1103) time_move, time_move1
#endif
#ifdef TWO_FIELDS
         write(69,1103) time_move, time_move1
#endif
         write(69,1104) time_bc, time_bc1
         write(69,1105) time_io, time_io1
         write(69,1114) time_dg, time_dg1
#ifdef OBSTACLES
         write(69,1116) time_obs, time_obs1
#endif
#ifdef SERIAL
! do nothing
#else
         write(69,1115) time_mp, time_mp1
#endif
         write(69,1117) time_loop-(time_coll+time_move+time_bc+time_io+time_dg+time_mp+time_obs)
         write(69,9999)
         write(69,*)  "# Ratio per section "
#ifdef SERIAL
! do nothing
#else
         write(69,1201) time_mp/time_loop, time_mp1/time_loop1
#endif
         write(69,1202) time_io/time_loop, time_io1/time_loop1
         write(69,1203) time_bc/time_loop, time_bc1/time_loop1
         write(69,1204) time_coll/time_loop, time_coll1/time_loop1
         write(69,1205) time_dg/time_loop, time_dg1/time_loop1
#ifdef ORIGINAL
         write(69,1206) time_move/time_loop, time_move/time_loop1
#endif
#ifdef OBSTACLES
         write(69,1207) time_obs/time_loop, time_obs/time_loop1
#endif
#ifdef TWO_FIELDS
         write(69,1206) time_move/time_loop, time_move/time_loop1
#endif
         write(69,9999) 
         write(69,*)  "# Derived (global) metrics "
         write(69,1106) float(lx)*float(ly)*float(lz)* & 
                        (itfin-itstart)/(time_loop)/1000.0/1000.0
         write(69,1107) float(215)*float(lx)*float(ly)*float(lz)*  &
                        (itfin-itstart)/(time_coll)/1000.0/1000.0
         write(69,1108) float(19*8)*float(lx)*float(ly)*float(lz)* &
                        (itfin-itstart)/(time_coll)/1000.0/1000.0
#ifdef ORIGINAL
         write(69,1109) float(18*8)*float(lx)*float(ly)*float(lz)* &
                        (itfin-itstart)/(time_move)/1000.0/1000.0
#endif
         write(69,9999) 
         write(69,*)  " # Memory (task 0) metrics "
         write(69,1110) mem_start, mem_stop
         write(69,*) (l/1024.0)*(m/1024.0)*(n/1024.0)*19*2*4
         write(69,9999) 
         close(69)   ! bgk.perf
!
         call system("git log | grep commit | tail -n 1 >> bgk.perf")
!
      endif
!#endif
!
!     closing files...
!
      close(16)   ! bgk.log
      close(31)   ! restore_vel.bin
      close(60)   ! prof_k.XX.dat
      close(61)   ! prof_i.XX.dat
      close(64)   ! prof_j.XX.dat
      close(62)   ! u_med.dat
      close(67)   ! probe_g.xxxx.dat
      close(68)   ! probe.xxxx.dat
      if(myrank==0) then
         close(63)   ! diagno.dat
      endif
!
#ifdef NO_OUTPUT
! do nothing
#else 
# ifdef SERIAL
! do nothing
# else
      write(38,*)  "# Time for section "
      write(38,1101) time_loop, time_loop1
      write(38,1102) time_coll, time_coll1
#  ifdef ORIGINAL
      write(38,1103) time_move, time_move1
#  endif
      write(38,1104) time_bc, time_bc1
      write(38,1105) time_io, time_io1
      write(38,1114) time_dg, time_dg1
#  ifdef SERIAL
! do nothing
#  else
      write(38,1115) time_mp, time_mp1
#  endif
      write(38,9999)
      write(38,*)  "# Derived (global) metrics "
      write(38,1106) float(l*m*n)*(itfin-itstart)/(time_loop)/1000.0/1000.0
      write(38,1107) float(215*l*m*n)*(itfin-itstart)/(time_coll)/1000.0/1000.0
      write(38,1108) float(l*m*n*19*8)*(itfin-itstart)/(time_coll)/1000.0/1000.0
#  ifdef ORIGINAL
      write(38,1109) float(l*m*n*18*8)*(itfin-itstart)/(time_move)/1000.0/1000.0
#  endif
      write(38,9999)
      write(38,*) "# MPI time/BW (z) -->", timeZ,   &
                               (m+2)*(l+2)*(itfin-itstart)*knorm/timeZ
      write(38,*) "# MPI time/BW (y) -->", timeY,   &
                               (n+2)*(l+2)*(itfin-itstart)*knorm/timeY
      write(38,*) "# MPI time/BW (x) -->", timeX,   &
                               (n+2)*(m+2)*(itfin-itstart)*knorm/timeX
      write(38,*) "#", myrank, ":Memory (stop) --->", mem_stop
      close(38)   ! task.XXXXXXXXX.log
# endif
#endif
!
#ifdef SERIAL
! do nothing
#else
! free derived datatype
      call MPI_type_free(xyplane,ierr)
      call MPI_type_free(yzplane,ierr)
      call MPI_type_free(xzplane,ierr)
# ifdef MB
! do nothing
# else
      call MPI_Type_free(dump3d,ierr)
# endif
!
      call mpi_barrier(lbecomm,ierr)
!
! finalize all!!
      call MPI_finalize(ierr)
#endif
!
! formats
!
9999  format(" #--------------------------------")
1100  format(" # init   time",1(e14.6,1x))
1101  format(" # loop   time",2(e14.6,1x))
1102  format(" # coll   time",2(e14.6,1x))
1103  format(" # move   time",2(e14.6,1x))
1104  format(" # bc     time",2(e14.6,1x))
1105  format(" # I/O    time",2(e14.6,1x))
1114  format(" # diagno time",2(e14.6,1x))
1115  format(" # MPI    time",2(e14.6,1x))
1116  format(" # Obst   time",2(e14.6,1x))
1117  format(" # Check      ",1(e14.6,1x))
1201  format(" # Ratio MPI  ",2(f7.3,1x))
1202  format(" # Ratio I/O  ",2(f7.3,1x))
1203  format(" # Ratio BC   ",2(f7.3,1x))
1204  format(" # Ratio Coll ",2(f7.3,1x))
1205  format(" # Ratio Diag.",2(f7.3,1x))
1206  format(" # Ratio Move ",2(f7.3,1x))
1207  format(" # Ratio Obs  ",2(f7.3,1x))
1106  format(" # Mlups      ",1(f14.6,1x))
1107  format(" # coll   Flop",1(e14.6,1x), "MFlops")
1108  format(" # coll   BW  ",1(e14.6,1x), "GB/s")    ! double precision only
1109  format(" # move   BW  ",1(e14.6,1x), "GB/s")    ! double precision only
1110  format(" # Memory (start,stop)",2(f14.6,1x), "MB")  ! double precision only
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. finalize"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. finalize mem =", mem_stop
      endif
# endif

        return
        end subroutine finalize
