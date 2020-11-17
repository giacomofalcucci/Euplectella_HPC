!=====================================================================
!     ****** LBE/diagnostic
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different diagnostic routines..
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
        subroutine diagnostic(itime,ivtim,icheck,itsave)
!
        use storage
        use timing
#ifdef CUDAFOR
        use cudafor
#endif
        implicit none
!
        integer:: itime,ivtim,icheck,itsave
        integer:: istart,istop
        integer:: jstart,jstop
        integer:: kstart,kstop
!
        real(mykind) :: fluxX, fluxY, fluxZ
!
! get macroscopic values
!
        if (mod(itime,ivtim).eq.0) then
!
! start timing...
           call SYSTEM_CLOCK(countA0, count_rate, count_max)
           call time(tcountA0)
!
!          
#ifdef CUDAFOR
! copy back (all) to host for diagnostic
           a02 = dev_a02
           a04 = dev_a04
           a05 = dev_a05
           a06 = dev_a06
           a11 = dev_a11
           a13 = dev_a13
           a14 = dev_a14
           a15 = dev_a15
           a19 = dev_a19
#endif

#ifdef HDF5
# ifdef SERIAL
           call hdf5_2d_xz_serial(itime)
!           call vtk_3d_bin(itime)               ! to remove
           call vtk_3d_bin(itime)               ! to remove
# else
           call vtk_xy_bin(itime,n-1)
           call vtk_xz_bin(itime,l/2-1)
!           call hdf5_2d_xy_parallel(itime)
!           call hdf5_2d_xz_parallel(itime)
!           call hdf5_2d_xy_parallel(itime)   (Da fare!!!!!)
           call vtk_3d_bin(itime)               ! to remove
# endif
#else

# ifdef NO_OUTPUT
!do nothing
# else
#  ifdef HPC
           call vtk_3d_bin(itime)
           call vtk_xy_bin(itime,n/2+1)
!           call vtk_om_bin(itime,n/2+1)
#  else
!           call vtk_yz(itime)
!           call vtk_3d(itime)
!           call vtk_xz(itime,m/2)
           call vtk_xy_bin(itime,n/2-1)
           call vtk_xz_bin(itime,m/2-1)
           call vtk_yz_bin(itime,l/2-1)
           call vtk_om_bin(itime,n/2+1)
!           call vtk_3d_bin(itime)
#  endif
# endif
#endif
!
! stop timing
           call time(tcountA1)
           call SYSTEM_CLOCK(countA1, count_rate, count_max)
           time_dg = time_dg + real(countA1-countA0)/(count_rate)
           time_dg1 = time_dg1 + (tcountA1-tcountA0)
!
        end if
!
        if (mod(itime,icheck).eq.0) then
!
! start timing...
           call SYSTEM_CLOCK(countA0, count_rate, count_max)
           call time(tcountA0)
!
! copy back (all) to host for diagnostic
#ifdef CUDAFOR
           a02 = dev_a02
           a04 = dev_a04
           a05 = dev_a05
           a06 = dev_a06
           a11 = dev_a11
           a13 = dev_a13
           a14 = dev_a14
           a15 = dev_a15
           a19 = dev_a19
#endif
!
!GA           call diagno(itime)
!
#ifdef NO_OUTPUT
! do nothing
#else
!
! drag & drop
!
#ifdef CYLINDER_INF
# ifdef DRAG
           call drag(itime)
# endif 
#endif 
!
           call varm(itime)
           call prof_i(itime,m/2,n/2)
           call prof_j(itime,l/2,n/2)
           call prof_k(itime,l/2,m/2)
!
!            call probe(itime,l-2,m/2,n/2)
! global probe: 
! x ---> LX = 20  x_p = 8
! y ---> LX = 16  x_p = 9
            call probe_global(itime,(8*lx/20),(9*ly/16),(lz/2-5))
!            call probe_global(itime,(3*lx/4),(ly/2),(lz/2-5))
!
#endif
!           call flush(61)            ! flush prof_i.dat
!           call flush(68)            ! flush probe.dat
!           call flush(88)            ! flush fort.88 (convergence)
!
! stop timing
           call time(tcountA1)
           call SYSTEM_CLOCK(countA1, count_rate, count_max)
           time_dg = time_dg + real(countA1-countA0)/(count_rate)
           time_dg1 = time_dg1 + (tcountA1-tcountA0)
!
        end if
!
        if (mod(itime,itsave).eq.0) then
!
#ifdef CUDAFOR
! copy back (all) to host for diagnostic
           a02 = dev_a02
           a04 = dev_a04
           a05 = dev_a05
           a06 = dev_a06
           a11 = dev_a11
           a13 = dev_a13
           a14 = dev_a14
           a15 = dev_a15
           a19 = dev_a19
#endif
           call save(itime)
!
        end if
!
! trick, in this way a (light) diagnostic call is done every tstep 
!
#ifdef OPENACCOLD
!!!         call diagno(itime)
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. diagnostic"
        endif
#endif
!
        return
        end subroutine diagnostic
