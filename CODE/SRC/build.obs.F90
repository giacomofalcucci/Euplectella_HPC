!=====================================================================
!     ****** LBE/build_obs
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       set obstacles 
!       no obstacle (flag1=0)
!       from file   (flag1=1)
!       creating    (flag1=2)
!                    sphere   ---> flag2=1
!                    cilinder ---> flag2=2
!                    h-cilinder -> flag2=3
!       obs dump    (flag3)
!                    none    ---> flag3=0
!                    vtk     ---> flag3=1
!                    binary  ---> flag3=2 (restart)
!                    hdf5    ---> flag3=3 (restart)
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
      subroutine build_obs
!
      use timing
      use storage
#ifdef SERIAL
! do nothing
#else
      use mpi
#endif
!
      implicit none
!
      integer:: i, j, k, ierr
      integer:: icoord, jcoord, kcoord
      integer:: sizeX, sizeY, sizeZ
      integer:: startX, startY, startZ
      integer:: stopX, stopY, stopZ
      integer:: offsetX, offsetY, offsetZ
      integer:: itime
      integer:: check, index
      real(mykind):: d2, R2a, R2b, R
      character*24 file_nameVTK
      character*21 file_nameBIN
      character*14 file_nameHDF5
!
      file_nameVTK  = 'tec_ob.xxxx.xxxxxxxx.vtk'
      file_nameBIN  = 'obs.xxxx.xxxxxxxx.bin'
      file_nameHDF5 = 'restore.obs.h5'
!
      imin = l
      jmin = m
      kmin = n
      imax = 0
      jmax = 0
      kmax = 0
      nobs = -100
!
#ifdef SERIAL
      myrank = 0
#endif
!
      itime  = 0
!
! some check
      if(flag1.gt.2) then
        write(6,*) "ERROR: flag1 not supported -->", flag1
        stop 
      endif
!
      if(flag2.gt.3) then
        write(6,*) "ERROR: flag2 not supported -->", flag2
        stop 
      endif
!
      if(flag3.gt.3) then
        write(6,*) "ERROR: flag3 not supported -->", flag3
        stop 
      endif
!
!
! Section 1 (flag1): reading section 
!
! let's start...
      if(flag1 == 0) then
         write(6,*) "INFO: no obstacle" 
      else
!
! reading from external file (one for task)
         if(flag1 == 1) then
!
#ifdef HDF5
# ifdef SERIAL
            write(6,*) "INFO: reading obstacle from (serial)", file_nameHDF5
            call restore_obs_hdf5_serial
# else 
            write(6,*) "INFO: reading obstacle from (parallel)", file_nameHDF5
            call restore_obs_hdf5_parallel
# endif
!
! check
            check = 0
            do k = 1, n
               do j = 1, m
                  do i = 1, l
 
                  if(obs(i,j,k)==1) then 
                     check = check+1
!
                     imin = min(imin,i)
                     jmin = min(jmin,j)
                     kmin = min(kmin,k)
!
                     imax = max(imax,i)
                     jmax = max(jmax,j)
                     kmax = max(kmax,k)
                  endif
                  enddo
               enddo
            enddo
            nobs = check
            write(6,*) "INFO: num. obs       -->",nobs, myrank
            write(6,*) "INFO: ratio obs/size -->",nobs/float(l*n*m), myrank
            write(6,*) "INFO: obs (x)           ",imax, imin, myrank
            write(6,*) "INFO: obs (y)           ",jmax, jmin, myrank
            write(6,*) "INFO: obs (z)           ",kmax, kmin, myrank
            write(16,*) "INFO: num. obs       -->",check, myrank
            write(16,*) "INFO: ratio obs/size -->",check/float(l*n*m), myrank
            write(16,*) "INFO: obs (x)           ",imax, imin, myrank
            write(16,*) "INFO: obs (y)           ",jmax, jmin, myrank
            write(16,*) "INFO: obs (z)           ",kmax, kmin, myrank
!
!
#else
!
            write(file_nameBIN(5:8),3100) myrank
            write(file_nameBIN(10:17),4000) itime
            write(6,*) "INFO: reading obstacle from ", file_nameBIN
            rewind(53)
            open(53,file=file_nameBIN, &
                    form="unformatted",status='unknown')
!
            check = 0
            read(53) nobs
            write(6,*) "INFO: obstacles in the task",nobs, myrank
!
#ifdef SERIAL
! do nothing
#else
            call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
!
            do index = 1,nobs
               read(53) i,j,k
               obs(i,j,k)=1
               check = check+1
!
               imin = min(imin,i)
               jmin = min(jmin,j)
               kmin = min(kmin,k)
!
               imax = max(imax,i)
               jmax = max(jmax,j)
               kmax = max(kmax,k)

            end do
!
! check
            if (check.NE.nobs) then
                write(6,*) "WARNING: obstacles not correct->", check,nobs
            else
                write(6,*)  "INFO: reading obstacle from file done" 
                write(16,*) "INFO: reading obstacle from file done" 
                close(53) ! binary obs. dump file
            endif
!
            write(6,*) "INFO: num. obs       -->",nobs
            write(6,*) "INFO: ratio obs/size -->",nobs/float(l*n*m), myrank
            write(6,*) "INFO: obs (x)           ",imax, imin, myrank
            write(6,*) "INFO: obs (y)           ",jmax, jmin, myrank
            write(6,*) "INFO: obs (z)           ",kmax, kmin, myrank
!
#endif
         endif
!
!
! Section 2 (flag1): creating  obstacle 
!
! creating obstacles
         if(flag1 == 2) then
            if(myrank==0) then 
               write(6,*) "INFO: creating obstacle" 
            endif
!
! sphere
            if(flag2 == 1) then
               write(16,*) "INFO: building obstacle (sphere) "
               icoord = lx/4
               jcoord = ly/2
               kcoord = lz/2
#ifdef TEST1
               radius = 5
               radius = min(ly/4,lz/4)   
#else
               radius = min(ly/4,lz/4)   
               radius = 20                      ! GA: to remove
#endif
!
! sphere
#ifdef NO_OUTPUT
! do nothing
#else
               write(6,*) "INFO: sphere radius    -->", radius
               write(6,*) "INFO: sphere icoord    -->", icoord
               write(6,*) "INFO: sphere jcoord    -->", jcoord
               write(6,*) "INFO: sphere kcoord   --->", kcoord
#endif
!
               R = radius*1.0
               R2a = (R-2)*(R-2)   ! lower
               R2b = (R+2)*(R+2)   ! upper
!
               offsetX = mpicoords(1)*l
               offsetY = mpicoords(2)*m
               offsetZ = mpicoords(3)*n
!
               do k = 1, n
                  do j = 1, m
                     do i = 1, l
!
                        d2 = (icoord-(i+offsetX))*(icoord-(i+offsetX))  &
                            +(jcoord-(j+offsetY))*(jcoord-(j+offsetY))  &
                            +(kcoord-(k+offsetZ))*(kcoord-(k+offsetZ))
                        if((d2.gt.R2a).and.(d2.lt.R2b)) then
!
                           obs(i,j,k) = 1
                           nobs = nobs + 1
!
                           imin = min(imin,i)
                           jmin = min(jmin,j)
                           kmin = min(kmin,k)
!
                           imax = max(imax,i)
                           jmax = max(jmax,j)
                           kmax = max(kmax,k)
!
                        endif
                     end do
                  end do
               end do
!
               write(6,*) "INFO: num. obs         -->", nobs
               write(6,*) "INFO: ratio obs/size   -->", nobs/float(l*n*m)
               write(6,*) "INFO: obs (x)             ",imax, imin
               write(6,*) "INFO: obs (y)             ",jmax, jmin
               write(6,*) "INFO: obs (z)             ",kmax, kmin
            endif
!
! cilinder...
            if (flag2 == 2) then
!
               write(16,*) "INFO: building obstacle (cilinder) "
               SizeY  = ly/8
               SizeX  = SizeY
#ifdef CHECK2D
!               radius = ly/8 
!               icoord = lx/8
#else
# ifdef CYLINDER_INF
               radius = lx/40
               write(6,*) "WARNING: hand made fix, to solve"
               radius = 100
               write(6,*) "WARNING: hand made fix, to solve"
               icoord = lx/4
               jcoord = ly/2
               kcoord = 0
               SizeZ  = lz+1
# else
               radius = min(100,lx/8)
!
! to fix
               write(6,*) "WARNING: hand made fix, to solve"
               radius = min(20,lx/8)
               write(6,*) "WARNING: hand made fix, to solve"
               icoord = lx/8
               SizeZ  = lz ! full column
               icoord = lx/4
               jcoord = ly/2
               kcoord = 0
# endif
#endif
!
               startX = icoord-mpicoords(1)*l-radius-1
               startY = jcoord-mpicoords(2)*m-radius+1
               startZ = kcoord-mpicoords(3)*n
!
               stopX = startX+radius
               stopY = startY+radius
               stopZ = startZ+sizeZ
!              
               offsetX = mpicoords(1)*l
               offsetY = mpicoords(2)*m
               offsetZ = mpicoords(3)*n
!              
#ifdef NO_OUTPUT
! do nothing
#else
               write(6,*) "INFO: column radius       -->", radius
               write(6,*) "INFO: column sizeZ        -->", sizeZ
               write(6,*) "INFO: column startX,stopX -->", startX, stopX
               write(6,*) "INFO: column startY,stopX -->", startY, stopY
               write(6,*) "INFO: column startZ,stopZ -->", startZ, stopZ
#endif
!
               R = radius*1.0
!               R2a = (R-2)*(R-2)   ! lower
               R2a = 0
               R2b = (R+2)*(R+2)   ! upper
!
               do k = 1, n
                  if((k.ge.startZ).AND.(k.lt.stopZ)) then
                  do j = 1, m
                     do i = 1, l
!
                        d2 = (icoord-(i+offsetX))*(icoord-(i+offsetX))  &
                             +(jcoord-(j+offsetY))*(jcoord-(j+offsetY))  
!
                        if((d2.ge.R2a).and.(d2.lt.R2b)) then
                             obs(i,j,k) = 1
                             nobs = nobs + 1
!
                             imin = min(imin,i)
                             jmin = min(jmin,j)
                             kmin = min(kmin,k)
!
                             imax = max(imax,i)
                             jmax = max(jmax,j)
                             kmax = max(kmax,k)
!
                        endif
                     end do
                  end do
                  endif
               end do
!
#ifdef NO_OUTPUT
! do nothing
#else
               write(6,*) "INFO: num. obs            -->",nobs
               write(6,*) "INFO: ratio obs/size      -->",nobs/float(l*n*m)
               write(6,*) "INFO: obs (x)                ",imax, imin
               write(6,*) "INFO: obs (y)                ",jmax, jmin
               write(6,*) "INFO: obs (z)                ",kmax, kmin
#endif
!
            endif
!
! half cilinder
            if (flag2 == 3) then
!
               write(16,*) "INFO: building obstacle (half-cilinder) "
               SizeY  = ly/8
               SizeX  = SizeY
!              SizeZ  = lz/2 ! half column
               SizeZ  = 60
               radius = min(20,lz/8)
!               radius = min(20,lz/8)
               icoord = lx/4
               jcoord = ly/2
               kcoord = 0
!
               startX = icoord-mpicoords(1)*l-radius-1
               startY = jcoord-mpicoords(2)*m-radius+1
               startZ = kcoord-mpicoords(3)*n
!
!
               stopX = startX+radius
               stopY = startY+radius
               stopZ = startZ+sizeZ
!              
               offsetX = mpicoords(1)*l
               offsetY = mpicoords(2)*m
               offsetZ = mpicoords(3)*n
!              
               write(6,*) "INFO: column radius       -->", radius
               write(6,*) "INFO: column sizeZ        -->", sizeZ
               write(6,*) "INFO: column startX,stopX -->", startX, stopX
               write(6,*) "INFO: column startY,stopX -->", startY, stopY
               write(6,*) "INFO: column startZ,stopZ -->", startZ, stopZ
!
               R = radius*1.0
               R2a = 0
               R2b = (R+2)*(R+2)   ! upper
!
               do k = 1, n
                  if((k.ge.startZ).AND.(k.lt.stopZ)) then
                  do j = 1, m
                     do i = 1, l
                        d2 = (icoord-(i+offsetX))*(icoord-(i+offsetX)) &
                            +(jcoord-(j+offsetY))*(jcoord-(j+offsetY))  
!
                        if((d2.ge.R2a).and.(d2.lt.R2b)) then
                             obs(i,j,k) = 1
                             nobs = nobs + 1
!
                             imin = min(imin,i)
                             jmin = min(jmin,j)
                             kmin = min(kmin,k)
!
                             imax = max(imax,i)
                             jmax = max(jmax,j)
                             kmax = max(kmax,k)
                        endif
                     end do
                  end do
                  endif
               end do
!
               write(6,*) "INFO: num. obs        -->",nobs
               write(6,*) "INFO: ratio obs/size  -->",nobs/float(l*n*m)
               write(6,*) "INFO: obs (x)            ",imax, imin
               write(6,*) "INFO: obs (y)            ",jmax, jmin
               write(6,*) "INFO: obs (z)            ",kmax, kmin
!
            endif
!
         endif
      endif
!
!
! Section 3 (flag3): dumping obstacle 
!
! hdf5/binary/vtk dump section
!
#ifdef NO_OUTPUT
      if(myrank==0) then 
         write(6,*) "INFO: no output mode enabled, no dump at all"
         write(16,*) "INFO: no output mode enabled, no dump at all"
      endif
      flag3 = 0
#endif
!
      if(flag3 == 0) then
! nothing to do
      else
         if(flag3 == 1) then            ! vtk dump
!
            if(nobs.le.0) then
               write(6,*) "task", myrank, "has no obstacles"
            else 
               write(6,*) "task", myrank, "has ", nobs, " obstacles"
! 
               write(file_nameVTK(8:11),3100) myrank
               write(file_nameVTK(13:20),4000) itime
               write(6,*) "INFO: obstacle vtk dump ", file_nameVTK
!
               open(52,file=file_nameVTK,status='unknown')
               write(52,'(A26)')'# vtk DataFile Version 2.0'
               write(52,'(A5)')'Campo'
               write(52,'(A5)')'ASCII'
               write(52,'(A24)')'DATASET RECTILINEAR_GRID'
               write(52,'(A11,I10,A1,I10,A1,I10)')  & 
                                   'DIMENSIONS ',l,' ',m,' ',n
!
               write(52,'(A14,I10,A7)')'X_COORDINATES ',l,' double'
               do i = 1,l
                  write(52, *) i + offset(1)
               enddo
!
               write(52,'(A14,I10,A7)')'Y_COORDINATES ',m,' double'
               do j = 1,m
                  write(52, *) j + offset(2)
               enddo
!
               write(52,'(A14,I10,A7)')'Z_COORDINATES ',n,' double'
               do k = 1,n
                  write(52, *) k + offset(3)
               enddo
!
               write(52,'(A10,I10)')'POINT_DATA ',l*m*n
               write(52,'(A21)')'SCALARS obs double'
               write(52,'(A20)')'LOOKUP_TABLE default'
               do k = 1,n
                  do j = 1,m
                     do i = 1,l
                        write(52,*) obs(i,j,k)*1.0
                     end do
                  end do
               end do
               write(6,*) "INFO: obs vtk dump done"
               close(52) ! vtk obs. dump file
            endif 
         endif 
!
         if(flag3 == 2) then            ! bin dump
!
! by definition all obstacles are no-slip.....
!
            write(file_nameBIN(5:8),3100) myrank
            write(file_nameBIN(10:17),4000) itime
#ifdef NO_OUTPUT
! do nothing
#else
            write(6,*) "INFO: obstacle binary dump ", file_nameBIN
            write(6,*) "INFO_ writing obstacles for task", nobs, myrank
#endif
            open(53,file=file_nameBIN, & 
                    form="unformatted",status='unknown')
            check = 0
!
            write(53) nobs
            do k = 1,n
               do j = 1,m
                  do i = 1,l
                     if (obs(i,j,k)==1) then
                        write(53) i,j,k
                        check = check+1
                     endif
                  end do
               end do
            end do
!
! check
            if (check.NE.nobs) then 
                write(6,*) "ERROR: obstacles not correct->", check,nobs
            else
                write(6,*) "INFO: obs binary dump done"
                close(53) ! binary obs. dump file
            endif
         endif 
!
         if(flag3 == 3) then    ! HDF5 dump
#ifdef HDF5
# ifdef SERIAL
           write(6,*) "INFO: serial hdf5 dump"
           call save_obs_hdf5_serial
# else
           write(6,*) "INFO: parallel hdf5 dump"
           write(6,*) "ERROR: still not implemented"
           stop
# endif
#else
           write(6,*) "ERROR: HDF5 preproc flag not activated..."
           stop
#endif
         endif 
      endif
!
      write(38,*) "#", myrank, ":num. obstacles--->", nobs
      write(38,*) "#", myrank, ":ratio obs      -->", float(nobs)/(l*m*n)
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. build_obs"
      endif
#endif
!
3100    format(i4.4)
4000    format(i8.8)
!
      return
      end subroutine build_obs
