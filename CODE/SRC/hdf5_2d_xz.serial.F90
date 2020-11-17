!=======================================================================
!     ****** LBE/hdf5_2d_xz_serial
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       hdf5_2d_xz_serial
!     DESCRIPTION
!       dump physical variables using hdf5 (using serial hdf5)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!
!     *****
!=======================================================================
!
      subroutine hdf5_2d_xz_serial(itime)
!
      use storage
      use hdf5
!
      implicit none
!
      character(len=16) :: file_name="u_xz_XXXXXXXX.h5"
      integer           :: itime,i,k,j0
      integer(hid_t)    :: file_id
      integer(hid_t)    :: group_a, group_c
      integer(hid_t)    :: dset_a,  dset_c
      integer(hid_t)    :: space_a, space_c
      integer(hsize_t)  :: dims_a(1), dims_c(2)
      integer           :: hdferr       
!
      real(sp), dimension(:,:), allocatable:: u,v,w,den  ! always single precision...
!
      dims_a(1)=1
!
      dims_c(1)=l
      dims_c(2)=n
!
      j0 = m/2
!
      allocate(den(dims_c(1),dims_c(2)))
      allocate(  u(dims_c(1),dims_c(2)))
      allocate(  w(dims_c(1),dims_c(2)))
      allocate(  v(dims_c(1),dims_c(2)))
!
      do k = 1,dims_c(2)
         do i = 1,dims_c(1)
            den(i,k)=                                                         &
                  a01(i,j0,k)+a02(i,j0,k)+a03(i,j0,k)+a04(i,j0,k)+a05(i,j0,k) &
                 +a06(i,j0,k)+a07(i,j0,k)+a08(i,j0,k)+a09(i,j0,k)+a10(i,j0,k) &
                 +a11(i,j0,k)+a12(i,j0,k)+a13(i,j0,k)+a14(i,j0,k)+a15(i,j0,k) &
                 +a16(i,j0,k)+a17(i,j0,k)+a18(i,j0,k)+a19(i,j0,k)
!
            u(i,k) =                                                          &
                 (a01(i,j0,k)+a02(i,j0,k)+a03(i,j0,k)+a04(i,j0,k)+a05(i,j0,k) &
                 -a10(i,j0,k)-a11(i,j0,k)-a12(i,j0,k)-a13(i,j0,k)-a14(i,j0,k))&
                 /den(i,k)
!
            w(i,k) =                                                          &
                 (a03(i,j0,k)+a07(i,j0,k)+a08(i,j0,k)+a09(i,j0,k)+a12(i,j0,k) &
                 -a01(i,j0,k)-a10(i,j0,k)-a16(i,j0,k)-a17(i,j0,k)-a18(i,j0,k))&
                 /den(i,k)
!
            v(i,k) =                                                          &
                 (a04(i,j0,k)+a06(i,j0,k)+a07(i,j0,k)+a13(i,j0,k)+a18(i,j0,k) &
                 -a02(i,j0,k)-a09(i,j0,k)-a11(i,j0,k)-a15(i,j0,k)-a16(i,j0,k))&
                 /den(i,k)
         enddo
      enddo
!
      write(file_name(6:13),4000) itime
!
      write(16,*) 'saving t=', itime, file_name
      write(6,*)  'saving t=', itime, file_name
      write(6,*)  'writing file_name -->', file_name
!
! open the fortran interface for hdf5
      call h5open_f(hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in h5open_f", hdferr
         stop
      endif
!
! create the file...
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, hdferr, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5fcreate_f", hdferr
         stop
      endif
!
! create the groups
! create group time, its dataspace and write...
      call h5gcreate_f(file_id, "time", group_a, hdferr, &
                     OBJECT_NAMELEN_DEFAULT_F, H5P_DEFAULT_F, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5gcreate_f", hdferr
         stop
      endif
      call h5screate_simple_f(1, dims_a, space_a, hdferr)
      call h5dcreate_f(group_a, "timestep", H5T_STD_I32LE, space_a, & 
                     dset_a, hdferr)
      call h5dwrite_f(dset_a, H5T_STD_I32LE, itime, dims_a, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (timestep)", hdferr
         stop
      endif
!
! create group field, its dataspace and write...
      call h5gcreate_f(file_id, "field", group_c, hdferr, &
                     OBJECT_NAMELEN_DEFAULT_F, H5P_DEFAULT_F, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5gcreate_f", hdferr
         stop
      endif
!
! create dataspace for file
      call h5screate_simple_f(2, dims_c, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for velocity along x) 
      call h5dcreate_f(group_c, "u", H5T_NATIVE_REAL, space_c, &
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (u)", hdferr
         stop
      endif
!
! write velocity along x
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, u, dims_c, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (u)", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for velocity along y)
      call h5dcreate_f(group_c, "w", H5T_NATIVE_REAL, space_c, &
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (w)", hdferr
         stop
      endif
!
! write velocity along y
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, w, dims_c, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (w)", hdferr
         stop
      endif
!

!-----------------------------------------------------------------------------
! create dataspace (for velocity along z)
      call h5dcreate_f(group_c, "v", H5T_NATIVE_REAL, space_c, &
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (v)", hdferr
         stop
      endif
!
! write velocity along y
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, v, dims_c, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (v)", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for density)
      call h5dcreate_f(group_c, "rho", H5T_NATIVE_REAL, space_c, &
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (rho)", hdferr
         stop
      endif

      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, den, dims_c, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (rho)", hdferr
         stop
      endif
!
! close the group
      call h5gclose_f(group_a, hdferr)
      call h5gclose_f(group_c, hdferr)
!
! close the file...
      call h5fclose_f(file_id, hdferr)
! 
! close the interface
      call h5close_f(hdferr);
!
! free memory...
      deallocate(den)
      deallocate(u)
      deallocate(w)
      deallocate(v)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. hdf5_2d_xz_serial"
        endif
#endif
!
4000  format(i8.8)

      return
      end subroutine hdf5_2d_xz_serial
