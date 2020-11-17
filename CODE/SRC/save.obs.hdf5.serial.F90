!=======================================================================
!     ****** LBE/save_obs_hdf5_serial
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       LBE/save_obs_hdf5_serial
!     DESCRIPTION
!       save obstacle in hdf5 format
!       using serial hdf5..
!       write on unit 200 (save.obs.h5)
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!       character*18 file_name
!       integer itime,i,k
!
!     *****
!=======================================================================
!
      subroutine save_obs_hdf5_serial
!
      use timing
      use storage
      use hdf5
!
      implicit none
!
      character*11 file_name
!
      integer           :: itime,i,k
      integer(hid_t)    :: file_id
      integer(hid_t)    :: group_a, group_b, group_c
      integer(hid_t)    :: dset_a,  dset_b,  dset_c
      integer(hid_t)    :: space_a, space_b, space_c
      integer    :: myrealtype
      integer(hsize_t)  :: dims_a(1), dims_b(1), dims_c(3)
      integer           :: buffer_b(3)
      real(mykind), dimension(:,:,:), allocatable    :: buffer_obs
      integer           :: hdferr       
!
! this trick doesn't work....(why?)
!!#ifdef DOUBLE_P
!!      myrealtype = H5T_NATIVE_DOUBLE
!!#else
!!      myrealtype = H5T_NATIVE_REAL
!!#endif
!
      file_name = 'save.obs.h5'
      itime = 0                          ! starting time...
!
      open(200,file=file_name,form="unformatted",status='unknown')
!
      write(16,*) 'dumping obstacle =', file_name
      write(6,*)  'dumping obstacle =', file_name
!
      dims_a(1) = 1
      dims_b(1) = 3
      dims_c(1) = l
      dims_c(2) = m
      dims_c(3) = n
!
      buffer_b(1) = l
      buffer_b(2) = m
      buffer_b(3) = n
!
      allocate(buffer_obs(1:l,1:m,1:n))
!
      buffer_obs(:,:,:) = obs(:,:,:)
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
      call h5gcreate_f(file_id, "obstacles", group_c, hdferr, &
                     OBJECT_NAMELEN_DEFAULT_F, H5P_DEFAULT_F, &
                     H5P_DEFAULT_F, H5P_DEFAULT_F)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5gcreate_f", hdferr
         stop
      endif
!
! -----------------------------------------------------------------
! dump obstacles 
      call h5screate_simple_f(3, dims_c, space_c, hdferr)
#ifdef DOUBLE_P
      call h5dcreate_f(group_c, "obs", H5T_NATIVE_DOUBLE, space_c, &
#else
      call h5dcreate_f(group_c, "obs", H5T_NATIVE_REAL, space_c, &
#endif
                     dset_c, hdferr)
#ifdef DOUBLE_P
      call h5dwrite_f(dset_c, H5T_NATIVE_DOUBLE, buffer_obs, dims_c, &
#else
      call h5dwrite_f(dset_c, H5T_NATIVE_REAL, buffer_obs, dims_c, &
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
      deallocate(buffer_obs)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. save_obs_hdf5_serial"
        endif
#endif
!
      return
      end subroutine save_obs_hdf5_serial
