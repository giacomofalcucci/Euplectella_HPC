!=======================================================================
!     ****** LBE/save_mpiio
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save
!     DESCRIPTION
!       save all microscopic variables (all populations)
!       using mpiio
!       write on unit 20 (aYY.xxxxxxxx.bin, unformatted, xxxxxxxx is timestep)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!       now 1 file for population, in the future 1 for all the
!       populations...
!     NOTES
!       character*18 file_name
!       integer itime,i,k
!       max time allowed  99'999'999
!     *****
!=======================================================================
!
      subroutine save_mpiio(itime)
!
#ifdef SERIAL
! do nothing
#else
      use storage
      use mpi
!
      implicit none
!
      character*24 file_name01
      character*24 file_name02
      character*24 file_name03
      character*24 file_name04
      character*24 file_name05
      character*24 file_name06
      character*24 file_name07
      character*24 file_name08
      character*24 file_name09
      character*24 file_name10
      character*24 file_name11
      character*24 file_name12
      character*24 file_name13
      character*24 file_name14
      character*24 file_name15
      character*24 file_name16
      character*24 file_name17
      character*24 file_name18
      character*24 file_name19
!
      integer:: itime,i,k,ierr
!
      integer:: myfile01
      integer:: myfile02
      integer:: myfile03
      integer:: myfile04
      integer:: myfile05
      integer:: myfile06
      integer:: myfile07
      integer:: myfile08
      integer:: myfile09
      integer:: myfile10
      integer:: myfile11
      integer:: myfile12
      integer:: myfile13
      integer:: myfile14
      integer:: myfile15
      integer:: myfile16
      integer:: myfile17
      integer:: myfile18
      integer:: myfile19
!
      real(mykind), dimension(:,:,:), allocatable :: buffer
!
      file_name01 = 'a01.xxxxxxxx.bin'
      file_name02 = 'a02.xxxxxxxx.bin'
      file_name03 = 'a03.xxxxxxxx.bin'
      file_name04 = 'a04.xxxxxxxx.bin'
      file_name05 = 'a05.xxxxxxxx.bin'
      file_name06 = 'a06.xxxxxxxx.bin'
      file_name07 = 'a07.xxxxxxxx.bin'
      file_name08 = 'a08.xxxxxxxx.bin'
      file_name09 = 'a09.xxxxxxxx.bin'
      file_name10 = 'a10.xxxxxxxx.bin'
      file_name11 = 'a11.xxxxxxxx.bin'
      file_name12 = 'a12.xxxxxxxx.bin'
      file_name13 = 'a13.xxxxxxxx.bin'
      file_name14 = 'a14.xxxxxxxx.bin'
      file_name15 = 'a15.xxxxxxxx.bin'
      file_name16 = 'a16.xxxxxxxx.bin'
      file_name17 = 'a17.xxxxxxxx.bin'
      file_name18 = 'a18.xxxxxxxx.bin'
      file_name19 = 'a19.xxxxxxxx.bin'
!
      write(file_name01(5:12),4000) itime
      write(file_name02(5:12),4000) itime
      write(file_name03(5:12),4000) itime
      write(file_name04(5:12),4000) itime
      write(file_name05(5:12),4000) itime
      write(file_name06(5:12),4000) itime
      write(file_name07(5:12),4000) itime
      write(file_name08(5:12),4000) itime
      write(file_name09(5:12),4000) itime
      write(file_name10(5:12),4000) itime
      write(file_name11(5:12),4000) itime
      write(file_name12(5:12),4000) itime
      write(file_name13(5:12),4000) itime
      write(file_name14(5:12),4000) itime
      write(file_name15(5:12),4000) itime
      write(file_name16(5:12),4000) itime
      write(file_name17(5:12),4000) itime
      write(file_name18(5:12),4000) itime
      write(file_name19(5:12),4000) itime
! 
! allocate & fill buffer vector
      allocate(buffer(1:l,1:m,1:n))
!
! ----------------- 1st population -------------------------------
!
      buffer(1:l,1:m,1:n) = a01(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name01
      call MPI_FILE_OPEN(lbecomm, file_name01, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                       MPI_INFO_NULL, myfile01, ierr) 
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile01, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile01, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_FILE_CLOSE(myfile01,ierr)
!
! ----------------- 2nd population -------------------------------
!
      buffer(1:l,1:m,1:n) = a02(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name02
      call MPI_FILE_OPEN(lbecomm, file_name02, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile02, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile02, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile02, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile02,ierr)
!
! ----------------- 3rd population -------------------------------
!
      buffer(1:l,1:m,1:n) = a03(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name03
      call MPI_FILE_OPEN(lbecomm, file_name03, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile03, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile03, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile03, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile03,ierr)
!
! ----------------- 4th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a04(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name04
      call MPI_FILE_OPEN(lbecomm, file_name04, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile04, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile04, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile04, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile04,ierr)
!
! ----------------- 5th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a05(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name05
      call MPI_FILE_OPEN(lbecomm, file_name05, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile05, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile05, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile05, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile05,ierr)
!
! ----------------- 6th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a06(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name06
      call MPI_FILE_OPEN(lbecomm, file_name06, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile06, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile06, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile06, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile06,ierr)
!
! ----------------- 7th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a07(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name07
      call MPI_FILE_OPEN(lbecomm, file_name07, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile07, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile07, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile07, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile03,ierr)
!
! ----------------- 8th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a08(1:l,1:m,1:n)
!
! open file...
     if (myrank == 0)  write(6,*) 'Saving using MPI-IO to ', file_name08
      call MPI_FILE_OPEN(lbecomm, file_name08, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile08, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile08, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile08, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile08,ierr)
!
! ----------------- 9th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a09(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name09
      call MPI_FILE_OPEN(lbecomm, file_name09, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile09, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile09, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile09, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile09,ierr)
!
! ----------------- 10th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a10(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name10
      call MPI_FILE_OPEN(lbecomm, file_name10, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile10, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile10, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile10, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile10,ierr)
!
! ----------------- 11th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a11(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name11
      call MPI_FILE_OPEN(lbecomm, file_name11, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile11, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile11, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile11, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile11,ierr)
!
! ----------------- 12th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a12(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name12
      call MPI_FILE_OPEN(lbecomm, file_name12, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile12, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile12, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile12, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile12,ierr)
!
! ----------------- 13th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a13(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name13
      call MPI_FILE_OPEN(lbecomm, file_name13, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile13, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile13, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile13, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile13,ierr)
!
! ----------------- 14th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a14(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name14
      call MPI_FILE_OPEN(lbecomm, file_name14, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile14, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile14, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile14, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile14,ierr)
!
! ----------------- 15th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a15(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name15
      call MPI_FILE_OPEN(lbecomm, file_name15, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile15, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile15, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile15, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile15,ierr)
!
! ----------------- 16th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a16(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name16
      call MPI_FILE_OPEN(lbecomm, file_name16, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile16, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile16, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile16, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile16,ierr)
!
! ----------------- 17th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a17(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name17
      call MPI_FILE_OPEN(lbecomm, file_name17, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile17, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile17, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile17, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile17,ierr)
!
! ----------------- 18th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a18(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name18
      call MPI_FILE_OPEN(lbecomm, file_name18, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile18, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile18, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile18, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile18,ierr)
!
! ----------------- 19th population -------------------------------
!
      buffer(1:l,1:m,1:n) = a19(1:l,1:m,1:n)
!
! open file...
      if (myrank == 0) write(6,*) 'Saving using MPI-IO to ', file_name19
      call MPI_FILE_OPEN(lbecomm, file_name19, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, myfile19, ierr)
!
! define view & write
      call MPI_FILE_SET_VIEW(myfile19, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_WRITE_ALL(myfile19, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile19,ierr)
!
!-------------------------------------------------------------
! barrier...
      call mpi_barrier(lbecomm,ierr)
!
! free memory...
      deallocate(buffer)
#endif
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. save_mpiio"
      endif
#endif
!
4000    format(i8.8)

      return
      end subroutine save_mpiio
