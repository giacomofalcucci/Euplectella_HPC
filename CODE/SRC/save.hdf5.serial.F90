!=======================================================================
!     ****** LBE/save_hdf5_serial
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save_hdf5_serial
!     DESCRIPTION
!       save all microscopic variables (all populations) in hdf5 format
!       using serial hdf5..
!       write on unit 20 (save.xxxxxxxx.bin, unformatted, xxxxxxxx is timestep)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!       character*18 file_name
!       integer itime,i,k
!       max time allowed  99'999'999
!
!     *****
!=======================================================================
!
      subroutine save_hdf5_serial(itime)
!
      use timing
      use storage
      use hdf5
!
      implicit none
!
      character*17 file_name
!
      integer           :: itime,i,k
      integer(hid_t)    :: file_id
      integer(hid_t)    :: group_a, group_b, group_c
      integer(hid_t)    :: dset_a,  dset_b,  dset_c
      integer(hid_t)    :: space_a, space_b, space_c
      integer    :: myrealtype
      integer(hsize_t)  :: dims_a(1), dims_b(1), dims_c(3)
      integer(hsize_t)  :: dims(3)
      integer           :: buffer_b(3)
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c01
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c02
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c03
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c04
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c05
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c06
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c07
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c08
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c09
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c10
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c11
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c12
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c13
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c14
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c15
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c16
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c17
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c18
      real(mykind), dimension(:,:,:), allocatable    :: buffer_c19
      integer           :: hdferr       
!
! this trick doesn't work....(why?)
!!#ifdef DOUBLE_P
!!      myrealtype = H5T_NATIVE_DOUBLE
!!#else
!!      myrealtype = H5T_NATIVE_REAL
!!#endif
!
      file_name = 'save.xxxxxxxx.h5'
!
      write(file_name(6:13),4000) itime
      open(20,file=file_name,form="unformatted",status='unknown')
!
      write(16,*) 'saving t=', itime, file_name
      write(6,*)  'saving t=', itime, file_name
!
      dims_a(1) = 1
      dims_b(1) = 3
      dims(1) = l+2             ! to remove????
      dims(2) = n+2             ! to remove????
      dims(3) = 9               ! to remove????
      dims_c(1) = l+2
      dims_c(2) = m+2
      dims_c(3) = n+2
!
      buffer_b(1) = l
      buffer_b(2) = m
      buffer_b(3) = n
!
      allocate(buffer_c01(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c02(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c03(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c04(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c05(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c06(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c07(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c08(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c09(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c10(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c11(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c12(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c13(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c14(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c15(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c16(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c17(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c18(0:l+1,0:m+1,0:n+1))
      allocate(buffer_c19(0:l+1,0:m+1,0:n+1))
!
      buffer_c01(:,:,:) = a01(:,:,:)
      buffer_c02(:,:,:) = a02(:,:,:)
      buffer_c03(:,:,:) = a03(:,:,:)
      buffer_c04(:,:,:) = a04(:,:,:)
      buffer_c05(:,:,:) = a05(:,:,:)
      buffer_c06(:,:,:) = a06(:,:,:)
      buffer_c07(:,:,:) = a07(:,:,:)
      buffer_c08(:,:,:) = a08(:,:,:)
      buffer_c09(:,:,:) = a09(:,:,:)
      buffer_c10(:,:,:) = a10(:,:,:)
      buffer_c11(:,:,:) = a11(:,:,:)
      buffer_c12(:,:,:) = a12(:,:,:)
      buffer_c13(:,:,:) = a13(:,:,:)
      buffer_c14(:,:,:) = a14(:,:,:)
      buffer_c15(:,:,:) = a15(:,:,:)
      buffer_c16(:,:,:) = a16(:,:,:)
      buffer_c17(:,:,:) = a17(:,:,:)
      buffer_c18(:,:,:) = a18(:,:,:)
      buffer_c19(:,:,:) = a19(:,:,:)
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
!
! create group size, its dataspace and write...
      call h5gcreate_f(file_id, "size", group_b, hdferr, &
                     OBJECT_NAMELEN_DEFAULT_F, H5P_DEFAULT_F, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5gcreate_f", hdferr
         stop
      endif
      call h5screate_simple_f(1, dims_b, space_b, hdferr)
      call h5dcreate_f(group_b, "dimension", H5T_STD_I32LE, space_b, &
                     dset_b, hdferr)
      call h5dwrite_f(dset_b, H5T_STD_I32LE, buffer_b, dims_b, &
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! create group popolation, its dataspace and write...
      call h5gcreate_f(file_id, "popolations", group_c, hdferr, &
                     OBJECT_NAMELEN_DEFAULT_F, H5P_DEFAULT_F, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5gcreate_f", hdferr
         stop
      endif
!
! -----------------------------------------------------------------
! dump popolation f01
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f01", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f01", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c01, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c01, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f02
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f02", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f02", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c02, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c02, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f03
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f03", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f03", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c03, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c03, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f04
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f04", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f04", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c04, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c04, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f05
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f05", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f05", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c05, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c05, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
! -----------------------------------------------------------------
! dump popolation f06
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f06", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f06", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c06, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c06, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f07
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f07", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f07", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c07, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c07, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f08
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f08", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f08", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c08, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c08, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f09
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f09", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f09", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c09, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c09, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f10
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f10", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f10", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c10, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c10, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f11
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f11", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f11", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c11, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c11, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f12
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f12", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f12", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c12, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c12, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f13
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f13", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f13", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c13, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c13, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f14
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f14", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f14", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c14, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c14, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f15
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f15", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f15", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c15, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c15, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f16
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f16", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f16", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c16, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c16, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f17
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f17", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f17", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c17, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c17, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f18
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f18", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f18", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c18, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c18, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f19
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f19", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f19", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c19, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c19, dims_c, &
#endif
                     hdferr, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
!
! -----------------------------------------------------------------
!
! close the dataset
      call h5dclose_f(dset_c, hdferr)
      call h5dclose_f(dset_a, hdferr)
      call h5dclose_f(dset_b, hdferr)
!
! close the group
      call h5gclose_f(group_a, hdferr)
      call h5gclose_f(group_b, hdferr)
      call h5gclose_f(group_c, hdferr)
!
! close the file...
      call h5fclose_f(file_id, hdferr)
! 
! close the interface
      call h5close_f(hdferr);
!
! free memory...
      deallocate(buffer_c01)
      deallocate(buffer_c02)
      deallocate(buffer_c03)
      deallocate(buffer_c04)
      deallocate(buffer_c05)
      deallocate(buffer_c06)
      deallocate(buffer_c07)
      deallocate(buffer_c08)
      deallocate(buffer_c09)
      deallocate(buffer_c10)
      deallocate(buffer_c11)
      deallocate(buffer_c12)
      deallocate(buffer_c13)
      deallocate(buffer_c14)
      deallocate(buffer_c15)
      deallocate(buffer_c16)
      deallocate(buffer_c17)
      deallocate(buffer_c18)
      deallocate(buffer_c19)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. save_hdf5_serial"
        endif
#endif
!
! formats
4000  format(i8.8)

      return
      end subroutine save_hdf5_serial
