!=======================================================================
!     ****** LBE/hdf5_2d_xy_parallel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       hdf5_2d_xy_parallel
!     DESCRIPTION
!       dump physical variables using hdf5 (parallel)
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
      subroutine hdf5_2d_xy_parallel(itime)
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
      character(len=20) :: file_name="u_xy_XXX_XXXXXXXX.h5"
      integer           :: itime,i,j,k0
      integer(hid_t)    :: file_id
      integer(hid_t)    :: plist_id        ! property access file list id
      integer(hid_t)    :: plist2_id       ! property access dataset list id
      integer(hid_t)    :: dset_a,  dset_c
      integer(hid_t)    :: dset_u, dset_v, dset_w, dset_rho  
      integer(hid_t)    :: group_a         ! group_id (timestep)
      integer(hid_t)    :: group_c         ! group_id (field)
      integer(hid_t)    :: space_a         ! dataspace_id (timestep)
      integer(hid_t)    :: memspace        ! dataspace id (task memory)
      iNteger(hid_t)    :: filespace       ! filespace id (to create dataspace)
      integer(hid_t)    :: filespace1      ! filespace id (to copy hyperslab)
      integer(hsize_t)  :: dims_a(1), dims_c(2), dims(2)
      integer(hsize_t)  :: offset2D(2), count(2)
      integer           :: hdferr
      integer           :: info
!
      real(sp), dimension(:,:), allocatable:: u,v,w,rho  ! always single precision...
!
      dims_a(1) = 1
!
! global size
      dims(1)   = lx
      dims(2)   = ly
!
! task size
      dims_c(1) = l
      dims_c(2) = m
!
! offset or starting point
      offset2D(1) = mpicoords(1)*l;  
      offset2D(2) = mpicoords(2)*m;  
!
      k0 = n/2
!
      count(1)  = dims_c(1)
      count(2)  = dims_c(2)
!
      allocate(rho(dims_c(1),dims_c(2)))
      allocate(  u(dims_c(1),dims_c(2)))
      allocate(  v(dims_c(1),dims_c(2)))
      allocate(  w(dims_c(1),dims_c(2)))
!
      do j = 1,dims_c(2)
         do i = 1,dims_c(1)
            rho(i,j) =                           &
               a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
              +a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a10(i,j,k0) &
              +a11(i,j,k0)+a12(i,j,k0)+a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0) &
              +a16(i,j,k0)+a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
            u(i,j) = &
              (a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
              -a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0)-a13(i,j,k0)-a14(i,j,k0))&
                 /rho(i,j)
!
            w(i,j) = &
              (a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a12(i,j,k0) &
              -a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0)-a17(i,j,k0)-a18(i,j,k0))&
                 /rho(i,j)
!
            v(i,j) = &
              (a04(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a13(i,j,k0)+a18(i,j,k0) &
              -a02(i,j,k0)-a09(i,j,k0)-a11(i,j,k0)-a15(i,j,k0)-a16(i,j,k0))&
                 /rho(i,j)
         enddo
      enddo
!
      write(file_name(6:8),3000) k0 ! to fix!!!!!!
      write(file_name(10:17),4000) itime
!
      if(myrank == 0) then
         write(16,*) 'saving t=', itime, file_name
         write(6,*)  'saving t=', itime, file_name
!         write(6,*)  'writing file_name -->', file_name
      endif
         write(6,*)  'writing file_name -->', file_name
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
      call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pcreate_F (1)", hdferr
         stop
      endif
!
      call H5Pset_fapl_mpio_f(plist_id, lbecomm, MPI_INFO_NULL, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pset_fapl_mpio", hdferr
         stop
      endif
!
#endif
!
#ifdef SERIAL
! do notihng
#else
! Create property list for collective dataset write 
      call H5Pcreate_f(H5P_DATASET_XFER_F, plist2_id, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pcreate_F (2)", hdferr
         stop
      endif
!
      call H5Pset_dxpl_mpio_f(plist2_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Pset_dxpl_mpio", hdferr
         stop
      endif
#endif
!
! Each task defines dataset in memory (to writes it to the hyperslab in the file)
      call h5screate_simple_f(2, dims_c, memspace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f (memspace)", hdferr
         stop
      endif
!
! create the file...
      call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, hdferr, &
                     H5P_DEFAULT_F, plist_id)
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
                     hdferr, H5S_ALL_F, H5S_ALL_F, plist2_id)
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
      call h5screate_simple_f(2, dims, filespace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5screate_simple_f", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for velocity along x) 
      call h5dcreate_f(group_c, "u", H5T_NATIVE_REAL, filespace, dset_u, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (u)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_u,filespace1,hdferr);
      call H5Sselect_hyperslab_F(filespace1, H5S_SELECT_SET_F, offset2D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (u)", hdferr
         stop
      endif
!
! write velocity along x
      call h5dwrite_f(dset_u, H5T_NATIVE_REAL, u, dims_c, hdferr, memspace, filespace1, plist2_id)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (u)", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for velocity along y) 
      call h5dcreate_f(group_c, "v", H5T_NATIVE_REAL, filespace, dset_v, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (v)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_v,filespace1,hdferr);
      call H5Sselect_hyperslab_F(filespace1, H5S_SELECT_SET_F, offset2D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (v)", hdferr
         stop
      endif
!
! write velocity along y
      call h5dwrite_f(dset_v, H5T_NATIVE_REAL, v, dims_c, hdferr, memspace, filespace1, plist2_id)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (v)", hdferr
         stop
      endif

!-----------------------------------------------------------------------------
! create dataspace (for velocity along z) 
      call h5dcreate_f(group_c, "w", H5T_NATIVE_REAL, filespace, dset_w, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (w)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_w,filespace1,hdferr);
      call H5Sselect_hyperslab_F(filespace1, H5S_SELECT_SET_F, offset2D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (w)", hdferr
         stop
      endif
!
! write velocity along z
      call h5dwrite_f(dset_w, H5T_NATIVE_REAL, w, dims_c, hdferr, memspace, filespace1, plist2_id)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (w)", hdferr
         stop
      endif
!
!-----------------------------------------------------------------------------
! create dataspace (for density) 
      call h5dcreate_f(group_c, "rho", H5T_NATIVE_REAL, filespace, dset_rho, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dcreate_f (rho)", hdferr
         stop
      endif
!
! Select hyperslab in the file 
      call H5Dget_space_F(dset_rho,filespace1,hdferr);
      call H5Sselect_hyperslab_F(filespace1, H5S_SELECT_SET_F, offset2D, count, hdferr);
      if(hdferr.ne.0) then
         write(*,*) "problem in H5Sselect_hyperslab_F (rho)", hdferr
         stop
      endif
!
! write density 
      call h5dwrite_f(dset_rho, H5T_NATIVE_REAL, rho, dims_c, hdferr, memspace, filespace1, plist2_id)
      if(hdferr.ne.0) then
         write(*,*) "problem in h5dwrite_f (rho)", hdferr
         stop
      endif
!
! close the property list
      call h5pclose_f(plist_id, hdferr)
      call h5pclose_f(plist2_id, hdferr)
!
! close the group
      call h5gclose_f(group_a, hdferr)
      call h5gclose_f(group_c, hdferr)
!
! close the dataset/dataspace
      call h5dclose_f(dset_a, hdferr)
      call h5dclose_f(dset_u, hdferr)
      call h5dclose_f(dset_v, hdferr)
      call h5dclose_f(dset_w, hdferr)
      call h5dclose_f(dset_rho, hdferr)
!
      call h5sclose_f(space_a, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in close (space_a)", hdferr
         stop
      endif
!
      call h5sclose_f(filespace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in close (filespace)", hdferr
         stop
      endif
!
      call h5sclose_f(filespace1, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in close (filespace1)", hdferr
         stop
      endif
!
      call h5sclose_f(memspace, hdferr)
      if(hdferr.ne.0) then
         write(*,*) "problem in close (memspace)", hdferr
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
!
! free memory...
      deallocate(rho)
      deallocate(u)
      deallocate(w)
      deallocate(v)
!
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. hdf5_2d_xy_parallel"
        endif
#endif
!
3000  format(i3.3)
4000  format(i8.8)

      return
      end subroutine hdf5_2d_xy_parallel
