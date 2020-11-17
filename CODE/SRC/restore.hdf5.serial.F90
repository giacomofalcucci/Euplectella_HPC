!=======================================================================
!     ****** LBE/restore_hdf5_serial
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       restore_hdf5_serial
!     DESCRIPTION
!       restore all microscopic variables (all populations) in hdf5 format
!       serial hdf5 version
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
      subroutine restore_hdf5_serial(itime)
!
      use storage
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
      write(16,*) 'restoring t=', itime, file_name
      write(6,*)  'restoring t=', itime, file_name
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
      dims_c(1) = l+2
      dims_c(2) = m+2
      dims_c(3) = n+2
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
! open the fortran interface for hdf5
      call h5open_f(hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in h5open_f", hdferr
         stop
      endif
!
! open the file...
      call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr, &
                     H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5fopen_f", hdferr
         stop
      endif
!
! datafile for timestep
      write(6,*) "reading -->", dsetname_a
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
         write(*,*) "problem in h5dread_f", hdferr
         stop
      endif
!
      itime = buffer_a(1)
!
! datafile for size
      write(6,*) "reading -->", dsetname_b
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
         write(*,*) "problem in h5dread_f", hdferr
         stop
      endif
!
! check 
      if(buffer_b(1).ne.l) then
         write(*,*) "ERROR: size doesn't match", l, buffer_b(1)
         stop
      endif
!
      if(buffer_b(2).ne.m) then
         write(*,*) "ERROR: size doesn't match", m, buffer_b(2)
         stop
      endif
!
      if(buffer_b(3).ne.n) then
         write(*,*) "ERROR: size doesn't match", n, buffer_b(3)
         stop
      endif
!
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c01
! open the datafile...
      call h5dopen_f(file_id, dsetname_c01, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c01, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c01, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c01", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c02
! open the datafile...
      call h5dopen_f(file_id, dsetname_c02, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c02, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c02, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c02", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c03
! open the datafile...
      call h5dopen_f(file_id, dsetname_c03, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c03, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c03, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c03", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c04
! open the datafile...
      call h5dopen_f(file_id, dsetname_c04, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c04, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c04, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c04", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c05
! open the datafile...
      call h5dopen_f(file_id, dsetname_c05, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c05, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c05, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c05", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c06
! open the datafile...
      call h5dopen_f(file_id, dsetname_c06, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c06, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c06, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c06", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c07
! open the datafile...
      call h5dopen_f(file_id, dsetname_c07, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c07, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c07, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c07", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c08
! open the datafile...
      call h5dopen_f(file_id, dsetname_c08, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c08, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c08, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c08", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c09
! open the datafile...
      call h5dopen_f(file_id, dsetname_c09, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c09, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c09, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c09", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c10
! open the datafile...
      call h5dopen_f(file_id, dsetname_c10, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c10, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c10, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c10", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c11
! open the datafile...
      call h5dopen_f(file_id, dsetname_c11, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c11, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c11, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c11", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c12
! open the datafile...
      call h5dopen_f(file_id, dsetname_c12, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c12, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c12, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c12", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
      write(6,*) "reading -->", dsetname_c13
! open the datafile...
      call h5dopen_f(file_id, dsetname_c13, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c13, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c13, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c13", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c14
! open the datafile...
      call h5dopen_f(file_id, dsetname_c14, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c14, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c14, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c14", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c15
! open the datafile...
      call h5dopen_f(file_id, dsetname_c15, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c15, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c15, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c15", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c16
! open the datafile...
      call h5dopen_f(file_id, dsetname_c16, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c16, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c16, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c16", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c17
! open the datafile...
      call h5dopen_f(file_id, dsetname_c17, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c17, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c17, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c17", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c18
! open the datafile...
      call h5dopen_f(file_id, dsetname_c18, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c18, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c18, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c18", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
     write(6,*) "reading -->", dsetname_c19
! open the datafile...
      call h5dopen_f(file_id, dsetname_c19, dset_c, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dopen_f", hdferr
!         stop
      endif
!
! read the datafile...
#ifdef DOUBLE_P
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_c19, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_c19, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c19", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
! --------------------------------------------------------
!
      a01(:,:,:) = buffer_c01(:,:,:)
      a02(:,:,:) = buffer_c02(:,:,:)
      a03(:,:,:) = buffer_c03(:,:,:)
      a04(:,:,:) = buffer_c04(:,:,:)
      a05(:,:,:) = buffer_c05(:,:,:)
      a06(:,:,:) = buffer_c06(:,:,:)
      a07(:,:,:) = buffer_c07(:,:,:)
      a08(:,:,:) = buffer_c08(:,:,:)
      a09(:,:,:) = buffer_c09(:,:,:)
      a10(:,:,:) = buffer_c10(:,:,:)
      a11(:,:,:) = buffer_c11(:,:,:)
      a12(:,:,:) = buffer_c12(:,:,:)
      a13(:,:,:) = buffer_c13(:,:,:)
      a14(:,:,:) = buffer_c14(:,:,:)
      a15(:,:,:) = buffer_c15(:,:,:)
      a16(:,:,:) = buffer_c16(:,:,:)
      a17(:,:,:) = buffer_c17(:,:,:)
      a18(:,:,:) = buffer_c18(:,:,:)
      a19(:,:,:) = buffer_c19(:,:,:)
!
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
! check
#ifdef DEBUG_1
      write(6,*) "DEBUG1: check 01 -->", a01(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 02 -->", a02(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 03 -->", a03(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 04 -->", a04(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 05 -->", a05(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 06 -->", a06(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 07 -->", a07(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 08 -->", a08(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 09 -->", a09(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 10 -->", a10(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 11 -->", a11(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 12 -->", a12(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 13 -->", a13(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 14 -->", a14(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 15 -->", a15(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 16 -->", a16(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 17 -->", a17(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 18 -->", a18(l/2,m/2,n/2)
      write(6,*) "DEBUG1: check 19 -->", a19(l/2,m/2,n/2)

      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. restore_hdf5_serial"
      endif
#endif
!
! fix
!
!      itime = itime -1 

4000  format(i8.8)

      return
      end subroutine restore_hdf5_serial
