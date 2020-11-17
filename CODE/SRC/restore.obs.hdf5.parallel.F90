!=======================================================================
!     ****** LBE/restore_obs_hdf5_parallel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       restore_obs_hdf5_parallel
!     DESCRIPTION
!       restore obstacle in hdf5 format
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
      subroutine restore_obs_hdf5_parallel
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
      integer           :: buffer_a(1), buffer_b(3)
      integer(hsize_t)  :: offset3D(3), count(3)
      integer           :: hdferr       
!
      real(mykind), dimension(:,:,:), allocatable    :: buffer_obs
!
! this trick doesn't work .... (why???)
!!#ifdef DOUBLE_P
!!      myrealtype = H5T_NATIVE_DOUBLE
!!#else
!!      myrealtype = H5T_NATIVE_REAL
!!#endif
!
      file_name = 'restore.obs.h5'
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
      dsetname_c01 = "/obstacles/obs/"
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
      allocate(buffer_obs(1:l,1:m,1:n))
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
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_obs, dims_c, hdferr, &
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_obs, dims_c, hdferr, &
#endif
                     memspace, space_c, plist3_id)
!
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / obs", hdferr
         stop
      endif
!
      call h5sclose_f(space_c, hdferr)
      call h5dclose_f(dset_c, hdferr)
!

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
      obs(1:l,1:m,1:n) = buffer_obs(1:l,1:m,1:n)
!
! free memory...
      deallocate(buffer_obs)
!
! check

#ifdef DEBUG_1
      write(6,*) "DEBUG1: check 01 -->", obs(l/2,m/2,n/2), myrank

      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. restore_hdf5_parallel"
      endif
#endif
!
4000  format(i8.8)
!
      return
      end subroutine restore_obs_hdf5_parallel
