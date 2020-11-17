!=======================================================================
!     ****** LBE/restore_obs_hdf5_serial
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       restore_obs_hdf5_serial
!     DESCRIPTION
!       restore obstacle hdf5 format
!       serial hdf5 version
!       read from unit 300 (restore.obs.h5)
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!
!     *****
!=======================================================================
!
      subroutine restore_obs_hdf5_serial
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
      integer           :: buffer_a(1), buffer_b(3)
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
      write(16,*) 'restoring obs ', file_name
      write(6,*)  'restoring obs ', file_name
!
      dsetname_a = "/time/timestep/"
      dims_a(1) = 1
!
      dsetname_b = "/size/dimension/"
      dims_b(1) = 3
!
      dsetname_c01 = "/obstacles/obs/"
!
      dims_c(1) = l+2
      dims_c(2) = m+2
      dims_c(3) = n+2
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
      call h5dread_f(dset_c, H5T_NATIVE_DOUBLE, buffer_obs, dims_c, hdferr)
#else
      call h5dread_f(dset_c, H5T_NATIVE_REAL, buffer_obs, dims_c, hdferr)
#endif
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dread_f / c01", hdferr
         stop
      endif
!
      call h5dclose_f(dset_c, hdferr)
!
! --------------------------------------------------------
!
      obs(:,:,:) = buffer_obs(:,:,:)
!
! close the file...
      call h5fclose_f(file_id, hdferr)
! 
! close the interface
      call h5close_f(hdferr);
!
! free memory...
      deallocate(buffer_obs)
!
! check
#ifdef DEBUG_1
      write(6,*) "DEBUG1: check obs -->", obs(l/2,m/2,n/2)

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
      end subroutine restore_obs_hdf5_serial
