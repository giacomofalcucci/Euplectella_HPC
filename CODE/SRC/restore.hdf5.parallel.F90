!=======================================================================
!     ****** LBE/restore_hdf5_parallel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       restore_hdf5_parallel
!     DESCRIPTION
!       restore all microscopic variables (all populations) in hdf5 format
!       parallel hdf5 version
!       read from unit 30 (restore.h5)
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
      subroutine restore_hdf5_parallel(itime)
!
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
      integer    :: myrealtype
      integer(hid_t)    :: file_id
      integer(hid_t)    :: group_a, group_b, group_c
      integer(hid_t)    :: dset_a,  dset_b,  dset_c
      integer(hid_t)    :: space_a, space_b, space_c
      integer(hsize_t)  :: dims_a(1), dims_b(1), dims_c(3)
      integer(hsize_t)  :: dims(3)
      integer(hid_t)    :: memspace        ! dataspace id (task memory)

      integer(hid_t)    :: plist1_id        ! property access file list id
      integer(hid_t)    :: plist3_id       ! property access dataset list id
      character(len=16) :: dsetname_a, dsetname_b
      character(len=18) :: dsetname_c01
      character(len=18) :: dsetname_c02
      character(len=18) :: dsetname_c03
      character(len=18) :: dsetname_c04
      character(len=18) :: dsetname_c05
      character(len=18) :: dsetname_c06
      character(len=18) :: dsetname_c07
      character(len=18) :: dsetname_c08
      character(len=18) :: dsetname_c09
      character(len=18) :: dsetname_c10
      character(len=18) :: dsetname_c11
      character(len=18) :: dsetname_c12
      character(len=18) :: dsetname_c13
      character(len=18) :: dsetname_c14
      character(len=18) :: dsetname_c15
      character(len=18) :: dsetname_c16
      character(len=18) :: dsetname_c17
      character(len=18) :: dsetname_c18
      character(len=18) :: dsetname_c19
      integer           :: buffer_a(1), buffer_b(3)
      integer(hsize_t)  :: offset3D(3), count(3)
      integer           :: hdferr       
!
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
!
! this trick doesn't work .... (why???)
!!#ifdef DOUBLE_P
!!      myrealtype = H5T_NATIVE_DOUBLE
!!#else
!!      myrealtype = H5T_NATIVE_REAL
!!#endif
!
      file_name = 'restore.h5'
!
      if (myrank==0) then
         write(16,*) 'restoring from ', file_name
         write(6,*)  'restoring from ', file_name
      endif
!
      dsetname_a = "/time/timestep/"
      dims_a(1) = 1
!
      dsetname_b = "/size/dimension/"
      dims_b(1) = 3
!
      dsetname_c01 = "/popolations/f01/"
      dsetname_c02 = "/popolations/f02/"
      dsetname_c03 = "/popolations/f03/"
      dsetname_c04 = "/popolations/f04/"
      dsetname_c05 = "/popolations/f05/"
      dsetname_c06 = "/popolations/f06/"
      dsetname_c07 = "/popolations/f07/"
      dsetname_c08 = "/popolations/f08/"
      dsetname_c09 = "/popolations/f09/"
      dsetname_c10 = "/popolations/f10/"
      dsetname_c11 = "/popolations/f11/"
      dsetname_c12 = "/popolations/f12/"
      dsetname_c13 = "/popolations/f13/"
      dsetname_c14 = "/popolations/f14/"
      dsetname_c15 = "/popolations/f15/"
      dsetname_c16 = "/popolations/f16/"
      dsetname_c17 = "/popolations/f17/"
      dsetname_c18 = "/popolations/f18/"
      dsetname_c19 = "/popolations/f19/"
!
      dims(1) = lx
      dims(2) = ly
      dims(3) = lz
!
      dims_c(1) = l
      dims_c(2) = m
      dims_c(3) = n
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
! open the fortran interface for hdf5
      call h5open_f(hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in h5open_f", hdferr
         stop
      endif
!
#ifdef SERIAL
! do nothing
#else
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
! Create property list for collective dataset write/read
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
! open the file...
      call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr,plist1_id)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5fopen_f", hdferr
         stop
      endif
!
! datafile for timestep
      if(myrank==0) then 
         write(6,*) "reading -->", dsetname_a
      endif
! open the datafile
      call h5dopen_f(file_id, dsetname_a, dset_a, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
         stop
      endif
!
! read the datafile...
      call h5dread_f(dset_a, H5T_STD_I32LE, buffer_a, dims_a, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f (timestep)", hdferr
      endif
!
      itime = buffer_a(1)
!
! datafile for time
#ifdef DEBUG_1
      write(*,*) "DEBUG1: time -->", itime
#endif
!
      call h5dclose_f(dset_a, hdferr)
!
! datafile for size
      if(myrank==0) then 
         write(6,*) "reading -->", dsetname_b
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_b, dset_b, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
         stop
      endif
!
! read the datafile...
      call h5dread_f(dset_b, H5T_STD_I32LE, buffer_b, dims_b, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f (size)", hdferr
         stop
      endif
!
! check 
      if(buffer_b(1).ne.lx) then
         write(*,*) "ERROR: size doesn't match", lx, buffer_b(1)
         stop
      endif
!
      if(buffer_b(2).ne.ly) then
         write(*,*) "ERROR: size doesn't match", ly, buffer_b(2)
         stop
      endif
!
      if(buffer_b(3).ne.lz) then
         write(*,*) "ERROR: size doesn't match", lz, buffer_b(3)
         stop
      endif
!
#ifdef DEBUG_1
      write(*,*) "DEBUG1: global size ->", buffer_b(1), buffer_b(2), buffer_b(3)
      write(*,*) "DEBUG1: local size -->", dims_c(1), dims_c(2), dims_c(3)
      write(*,*) "DEBUG1: offset ------>", offset3d(1), offset3d(2), offset3d(3)
#endif
!
      call h5dclose_f(dset_b, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c01
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c01, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f01)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c01, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c01, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c01", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c02
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c02, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f02)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c02, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c02, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c02", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
          write(6,*) "reading -->", dsetname_c03
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c03, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f03)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c03, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c03, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c03", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c04
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c04, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f04)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c04, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c04, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c04", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c05
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c05, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f05)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c05, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c05, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c05", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c06
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c06, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f06)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c06, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c06, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c06", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c07
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c07, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f07)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c07, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c07, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c07", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c08
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c08, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f08)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c08, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c08, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c08", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c09
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c09, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f09)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c09, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c09, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c09", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c10
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c10, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f10)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c10, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c10, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c10", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c11
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c11, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!     
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f11)", hdferr
         stop
      endif
! 
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c11, dims_c, hdferr, &
#else    
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c11, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
! 
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c11", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c12
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c12, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f12)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c12, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c12, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c12", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c13
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c13, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f13)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c13, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c13, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c13", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c14
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c14, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f14)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c14, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c14, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c14", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c15
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c15, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f15)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c15, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c15, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c15", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c16
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c16, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f16)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c16, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c16, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c16", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c17
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c17, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f17)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c17, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c17, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c17", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c18
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c18, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f18)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c18, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c18, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c18", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      if (myrank==0) then
         write(6,*) "reading -->", dsetname_c19
      endif
! open the datafile...
      call h5dopen_f(file_id, dsetname_c19, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! Select hyperslab in the file
      call H5Dget_space_F(dset_c,space_c,hdferr);
      call H5Sselect_hyperslab_F(space_c, H5S_SELECT_SET_F, offset3D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (f19)", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c19, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c19, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c19", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------









! closing dataspace
      call h5sclose_f(memspace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5sclose_f", hdferr
         stop
      endif
!
! closing property lists...
      call h5pclose_f(plist1_id, hdferr)
      call h5pclose_f(plist3_id, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5pclose_f", hdferr
         stop
      endif
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
! coping back values...
      a01(1:l,1:m,1:n) = buffer_c01(1:l,1:m,1:n)
      a02(1:l,1:m,1:n) = buffer_c02(1:l,1:m,1:n)
      a03(1:l,1:m,1:n) = buffer_c03(1:l,1:m,1:n)
      a04(1:l,1:m,1:n) = buffer_c04(1:l,1:m,1:n)
      a05(1:l,1:m,1:n) = buffer_c05(1:l,1:m,1:n)
      a06(1:l,1:m,1:n) = buffer_c06(1:l,1:m,1:n)
      a07(1:l,1:m,1:n) = buffer_c07(1:l,1:m,1:n)
      a08(1:l,1:m,1:n) = buffer_c08(1:l,1:m,1:n)
      a09(1:l,1:m,1:n) = buffer_c09(1:l,1:m,1:n)
      a10(1:l,1:m,1:n) = buffer_c10(1:l,1:m,1:n)
      a11(1:l,1:m,1:n) = buffer_c11(1:l,1:m,1:n)
      a12(1:l,1:m,1:n) = buffer_c12(1:l,1:m,1:n)
      a13(1:l,1:m,1:n) = buffer_c13(1:l,1:m,1:n)
      a14(1:l,1:m,1:n) = buffer_c14(1:l,1:m,1:n)
      a15(1:l,1:m,1:n) = buffer_c15(1:l,1:m,1:n)
      a16(1:l,1:m,1:n) = buffer_c16(1:l,1:m,1:n)
      a17(1:l,1:m,1:n) = buffer_c17(1:l,1:m,1:n)
      a18(1:l,1:m,1:n) = buffer_c18(1:l,1:m,1:n)
      a19(1:l,1:m,1:n) = buffer_c19(1:l,1:m,1:n)
!
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
! check

#ifdef DEBUG_1
      write(6,*) "DEBUG1: check 01 -->", a01(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 02 -->", a02(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 03 -->", a03(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 04 -->", a04(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 05 -->", a05(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 06 -->", a06(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 07 -->", a07(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 08 -->", a08(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 09 -->", a09(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 10 -->", a10(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 11 -->", a11(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 12 -->", a12(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 13 -->", a13(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 14 -->", a14(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 15 -->", a15(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 16 -->", a16(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 17 -->", a17(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 18 -->", a18(l/2,m/2,n/2), myrank
      write(6,*) "DEBUG1: check 19 -->", a19(l/2,m/2,n/2), myrank

      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. restore_hdf5_parallel"
      endif
#endif
!
4000  format(i8.8)
!
      return
      end subroutine restore_hdf5_parallel
