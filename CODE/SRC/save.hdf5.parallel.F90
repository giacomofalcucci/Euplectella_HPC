!=======================================================================
!     ****** LBE/save_hdf5_parallel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save_hdf5_parallel
!     DESCRIPTION
!       save all microscopic variables (all populations) in hdf5 format
!       using parallel hdf5..
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
      subroutine save_hdf5_parallel(itime)
!
      use timing
      use storage
#ifdef SERIAL
! noting to include
#else
      use mpi
#endif
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
      integer(hid_t)    :: memspace        ! dataspace id (task memory)
      integer(hid_t)    :: plist1_id        ! property access file list id
      integer(hid_t)    :: plist3_id       ! property access dataset list id
      integer           :: myrealtype
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
      integer(hsize_t)  :: offset3D(3),count(3)
      integer           :: hdferr       
      integer           :: ierr
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
      if(myrank == 0) then
         write(16,*) 'saving t=', itime, file_name
         write(6,*)  'saving t=', itime, file_name
      endif
!
      dims_a(1) = 1 ! one value (timestep)
      dims_b(1) = 3 ! two values (size x,y & z )
      dims(1) = lx  ! global size
      dims(2) = ly  ! global size
      dims(3) = lz  ! global size
      dims_c(1) = l ! task size
      dims_c(2) = m ! task size
      dims_c(3) = n ! task size
!
      buffer_b(1) = lx
      buffer_b(2) = ly
      buffer_b(3) = lz
!
! offset or starting point
      offset3d(1) = mpicoords(1)*l;
      offset3d(2) = mpicoords(2)*m;
      offset3d(3) = mpicoords(3)*n;
!
      count(1)  = dims_c(1)
      count(2)  = dims_c(2)
      count(3)  = dims_c(3)
!
      allocate(buffer_c01(1:l,1:m,1:n))
      allocate(buffer_c02(1:l,1:m,1:n))
      allocate(buffer_c03(1:l,1:m,1:n))
      allocate(buffer_c04(1:l,1:m,1:n))
      allocate(buffer_c05(1:l,1:m,1:n))
      allocate(buffer_c06(1:l,1:m,1:n))
      allocate(buffer_c07(1:l,1:m,1:n))
      allocate(buffer_c08(1:l,1:m,1:n))
      allocate(buffer_c09(1:l,1:m,1:n))
      allocate(buffer_c10(1:l,1:m,1:n))
      allocate(buffer_c11(1:l,1:m,1:n))
      allocate(buffer_c12(1:l,1:m,1:n))
      allocate(buffer_c13(1:l,1:m,1:n))
      allocate(buffer_c14(1:l,1:m,1:n))
      allocate(buffer_c15(1:l,1:m,1:n))
      allocate(buffer_c16(1:l,1:m,1:n))
      allocate(buffer_c17(1:l,1:m,1:n))
      allocate(buffer_c18(1:l,1:m,1:n))
      allocate(buffer_c19(1:l,1:m,1:n))
!
      buffer_c01(1:l,1:m,1:n) = a01(1:l,1:m,1:n)
      buffer_c02(1:l,1:m,1:n) = a02(1:l,1:m,1:n)
      buffer_c03(1:l,1:m,1:n) = a03(1:l,1:m,1:n)
      buffer_c04(1:l,1:m,1:n) = a04(1:l,1:m,1:n)
      buffer_c05(1:l,1:m,1:n) = a05(1:l,1:m,1:n)
      buffer_c06(1:l,1:m,1:n) = a06(1:l,1:m,1:n)
      buffer_c07(1:l,1:m,1:n) = a07(1:l,1:m,1:n)
      buffer_c08(1:l,1:m,1:n) = a08(1:l,1:m,1:n)
      buffer_c09(1:l,1:m,1:n) = a09(1:l,1:m,1:n)
      buffer_c10(1:l,1:m,1:n) = a10(1:l,1:m,1:n)
      buffer_c11(1:l,1:m,1:n) = a11(1:l,1:m,1:n)
      buffer_c12(1:l,1:m,1:n) = a12(1:l,1:m,1:n)
      buffer_c13(1:l,1:m,1:n) = a13(1:l,1:m,1:n)
      buffer_c14(1:l,1:m,1:n) = a14(1:l,1:m,1:n)
      buffer_c15(1:l,1:m,1:n) = a15(1:l,1:m,1:n)
      buffer_c16(1:l,1:m,1:n) = a16(1:l,1:m,1:n)
      buffer_c17(1:l,1:m,1:n) = a17(1:l,1:m,1:n)
      buffer_c18(1:l,1:m,1:n) = a18(1:l,1:m,1:n)
      buffer_c19(1:l,1:m,1:n) = a19(1:l,1:m,1:n)
!
! open the fortran interface for hdf5
      call h5open_f(hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in h5open_f", hdferr
         stop
      endif
!
#ifdef SERIAL
! do noting
#else
!
! Set up file access property list for MPI-IO access 
      call H5Pcreate_f(H5P_FILE_ACCESS_F, plist1_id, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pcreate_F (1)", hdferr
         stop
      endif
!
      call H5Pset_fapl_mpio_f(plist1_id, lbecomm, MPI_INFO_NULL, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pset_fapl_mpio", hdferr
         stop
      endif
!
#endif
!
#ifdef SERIAL
! do nothing
#else
! Create property list for collective dataset write 
      call H5Pcreate_f(H5P_DATASET_XFER_F, plist3_id, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pcreate_F (2)", hdferr
         stop
      endif
!
      call H5Pset_dxpl_mpio_f(plist3_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pset_dxpl_mpio", hdferr
         stop
      endif
#endif
!
! Each task defines dataset in memory (to writes it to the hyperslab in the file)
      call h5screate_simple_f(3, dims_c, memspace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f (memspace)", hdferr
         stop
      endif
!
!
! create the file...
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, hdferr, &
                     H5P_DEFAULT_F, plist1_id)
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
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (dimension)", hdferr
         stop
      endif
!
!      mem_start = get_mem(); write(6,*) "MEM -->", mem_start
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
!-----------------------------------------------------------------------------
! dump popolation f01
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f01", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f01", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f01)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c01, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c01, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f02
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f02", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f02", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f02)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c02, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c02, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f03
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f03", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f03", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f03)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c03, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c03, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f04
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f04", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f04", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f04)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c04, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c04, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f05
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f05", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f05", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f05)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c05, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c05, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f06
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f06", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f06", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f06)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c06, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c06, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f07
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f07", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f07", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f07)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c07, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c07, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f08
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f08", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f08", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f08)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c08, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c08, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f09
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f09", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f09", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f09)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c09, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c09, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!!!!                     hdferr, memspace, space_c, H5P_DEFAULT_F)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f10
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f10", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f10", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f10)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c10, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c10, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f11
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f11", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f11", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f11)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c11, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c11, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f12
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f12", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f12", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f12)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c12, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c12, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f13
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f13", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f13", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f13)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c13, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c13, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f14
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f14", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f14", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f14)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c14, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c14, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f15
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f15", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f15", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f15)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c15, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c15, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f16
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f16", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f16", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f16)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c16, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c16, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f17
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f17", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f17", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f17)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c17, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c17, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f18
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f18", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f18", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f18)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c18, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c18, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------
! dump popolation f19
      call h5screate_simple_f(3, dims, space_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "f19", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "f19", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (f19)", hdferr
         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3d, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c19, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_c19, dims_c, &
#endif
                     hdferr, memspace, space_c, plist3_id)
!
! close the dataspace
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! -----------------------------------------------------------------






!
! close memspace
      call h5sclose_f(memspace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in close (memspace)", hdferr
         stop
      endif

! close the property list
      call h5pclose_f(plist1_id, hdferr)
      call h5pclose_f(plist3_id, hdferr)
!
! close the dataset
!      call h5dclose_f(dset_c, hdferr)
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
      if(hdferr.ne.0) then
         write(*,*) "problem in h5fclose_f", hdferr
         stop
      endif
!
! 
! close the interface
      call h5close_f(hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in h5close_f", hdferr
         stop
      endif
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
! formats
4000  format(i8.8)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. save_hdf5_parallel"
        endif
#endif
!
      return
      end subroutine save_hdf5_parallel
