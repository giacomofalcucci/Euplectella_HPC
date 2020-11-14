!=======================================================================
!     ****** LBE/restore_mpiio
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
!       3D decomposition
!	
!     NOTES
!       character*18 file_name
!       integer itime,i,k
!       max time allowed  99'999'999
!
!     *****
!=======================================================================
!
      subroutine restore_mpiio(itime)
!
#ifdef SERIAL
! do nothig
#else
      use storage
      use mpi
!
      implicit none
!
      LOGICAL :: file_exists
!
      character*24 file_name00
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
      real(mykind), dimension(:,:,:), allocatable:: buffer
!
      file_name00 = 'time.restore.txt'
      file_name01 = 'a01.XXXXXXXX.bin'
      file_name02 = 'a02.XXXXXXXX.bin'
      file_name03 = 'a03.XXXXXXXX.bin'
      file_name04 = 'a04.XXXXXXXX.bin'
      file_name05 = 'a05.XXXXXXXX.bin'
      file_name06 = 'a06.XXXXXXXX.bin'
      file_name07 = 'a07.XXXXXXXX.bin'
      file_name08 = 'a08.XXXXXXXX.bin'
      file_name09 = 'a09.XXXXXXXX.bin'
      file_name10 = 'a10.XXXXXXXX.bin'
      file_name11 = 'a11.XXXXXXXX.bin'
      file_name12 = 'a12.XXXXXXXX.bin'
      file_name13 = 'a13.XXXXXXXX.bin'
      file_name14 = 'a14.XXXXXXXX.bin'
      file_name15 = 'a15.XXXXXXXX.bin'
      file_name16 = 'a16.XXXXXXXX.bin'
      file_name17 = 'a17.XXXXXXXX.bin'
      file_name18 = 'a18.XXXXXXXX.bin'
      file_name19 = 'a19.XXXXXXXX.bin'
!
! allocate buffer vector
      allocate(buffer(1:l,1:m,1:n))
!
! ----------------- timestep ----- -------------------------------
      if(myrank == 0) write(6,*) 'Restoring timestep using MPI-IO ', file_name00
      open(72,file=file_name00, status='OLD')
      read(72,*) itime
      close(72) 
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
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: restoring at time", itime
         write(6,*) "DEBUG1: file ", file_name01
      endif
#endif
! ----------------- 1st population -------------------------------
! open file...
      INQUIRE(FILE=file_name01, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name01,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name01
      call MPI_FILE_OPEN(lbecomm, file_name01, &
                       MPI_MODE_RDONLY, & 
                       MPI_INFO_NULL, myfile01, ierr) 
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile01, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile01, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile01,ierr)
!
! fill from buffer vector
      a01(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 2nd population -------------------------------
! open file...
      INQUIRE(FILE=file_name02, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name02,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name02
      call MPI_FILE_OPEN(lbecomm, file_name02, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile02, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile02, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile02, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile02,ierr)
!
! fill from buffer vector
      a02(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 3rd population -------------------------------
! open file...
      INQUIRE(FILE=file_name03, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name03,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name03
      call MPI_FILE_OPEN(lbecomm, file_name03, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile03, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile03, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile03, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile03,ierr)
!
! fill from buffer vector
      a03(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 4th population -------------------------------
! open file...
      INQUIRE(FILE=file_name04, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name04,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name04
      call MPI_FILE_OPEN(lbecomm, file_name04, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile04, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile04, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile04, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile04,ierr)
!
! fill from buffer vector
      a04(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 5th population -------------------------------
! open file...
      INQUIRE(FILE=file_name05, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name05,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name05
      call MPI_FILE_OPEN(lbecomm, file_name05, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile05, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile05, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile05, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile05,ierr)
!
! fill from buffer vector
      a05(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 6th population -------------------------------
! open file...
      INQUIRE(FILE=file_name06, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name06,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name06
      call MPI_FILE_OPEN(lbecomm, file_name06, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile06, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile06, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile06, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile06,ierr)
!
! fill from buffer vector
      a06(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 7th population -------------------------------
! open file...
      INQUIRE(FILE=file_name07, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name07,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name07
      call MPI_FILE_OPEN(lbecomm, file_name07, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile07, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile07, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile07, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile07,ierr)
!
! fill from buffer vector
      a07(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 8th population -------------------------------
! open file...
      INQUIRE(FILE=file_name08, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name08,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name08
      call MPI_FILE_OPEN(lbecomm, file_name08, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile08, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile08, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile08, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile08,ierr)
!
! fill from buffer vector
      a08(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 9th population -------------------------------
! open file...
      INQUIRE(FILE=file_name09, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name09,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name09
      call MPI_FILE_OPEN(lbecomm, file_name09, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile09, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile09, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile09, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile09,ierr)
! 
! fill from buffer vector
      a09(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 10th population -------------------------------
! open file...
      INQUIRE(FILE=file_name10, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name10,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name10
      call MPI_FILE_OPEN(lbecomm, file_name10, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile10, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile10, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile10, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile10,ierr)
!
! fill from buffer vector
      a10(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 11th population -------------------------------
! open file...
      INQUIRE(FILE=file_name11, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name11,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name11
      call MPI_FILE_OPEN(lbecomm, file_name11, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile11, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile11, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile11, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile11,ierr)
!
! fill from buffer vector
      a11(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 12th population -------------------------------
! open file...
      INQUIRE(FILE=file_name12, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name12,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name12
      call MPI_FILE_OPEN(lbecomm, file_name12, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile12, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile12, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile12, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile12,ierr)
!
! fill from buffer vector
      a12(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 13th population -------------------------------
! open file...
      INQUIRE(FILE=file_name13, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name13,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name13
      call MPI_FILE_OPEN(lbecomm, file_name13, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile13, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile13, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile13, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile13,ierr)
!
! fill from buffer vector
      a13(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 14th population -------------------------------
! open file...
      INQUIRE(FILE=file_name14, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name14,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name14
      call MPI_FILE_OPEN(lbecomm, file_name14, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile14, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile14, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile14, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile14,ierr)
!
! fill from buffer vector
      a14(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 15th population -------------------------------
! open file...
      INQUIRE(FILE=file_name15, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name15,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name15
      call MPI_FILE_OPEN(lbecomm, file_name15, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile15, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile15, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile15, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile15,ierr)
!
! fill from buffer vector
      a15(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 16th population -------------------------------
! open file...
      INQUIRE(FILE=file_name16, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name16,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name16
      call MPI_FILE_OPEN(lbecomm, file_name16, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile16, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile16, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile16, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile16,ierr)
!
! fill from buffer vector
      a16(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 17th population -------------------------------
! open file...
      INQUIRE(FILE=file_name17, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name17,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name17
      call MPI_FILE_OPEN(lbecomm, file_name17, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile17, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile17, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile17, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile17,ierr)
! 
! fill from buffer vector
      a17(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 18th population -------------------------------
! open file...
      INQUIRE(FILE=file_name18, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name18,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name18
      call MPI_FILE_OPEN(lbecomm, file_name18, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile18, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile18, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile18, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile18,ierr)
!
! fill from buffer vector
      a18(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
! ----------------- 19th population -------------------------------
! open file...
      INQUIRE(FILE=file_name19, EXIST=file_exists)
      if (file_exists.neqv..TRUE.) then
          write(*,*) "ERROR: file ",file_name19,"not found!!!!!!"
          stop
      endif
!
      if(myrank == 0) write(6,*) 'Restoring using MPI-IO from ', file_name19
      call MPI_FILE_OPEN(lbecomm, file_name19, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, myfile19, ierr)
!
! define view and read...
      call MPI_FILE_SET_VIEW(myfile19, file_offset, MYMPIREAL, &
                            dump3d, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(myfile19, buffer, buffer_size, MYMPIREAL, &
                         MPI_STATUS_IGNORE, ierr)
!
! close
      call MPI_File_close(myfile19,ierr)
!
! fill from buffer vector
      a19(1:l,1:m,1:n) = buffer(1:l,1:m,1:n)
!
!-----------------------------------------------------------------
! barrier...
      call mpi_barrier(lbecomm,ierr)
! 
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: a01", a01(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a02", a02(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a03", a03(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a04", a04(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a05", a05(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a06", a06(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a07", a07(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a08", a08(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a09", a09(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a10", a10(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a11", a11(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a12", a12(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a13", a13(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a14", a14(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a15", a15(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a16", a16(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a17", a17(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a18", a18(l/2,m/2,n/2)
         write(6,*) "DEBUG1: a19", a19(l/2,m/2,n/2)
      endif
#endif
!
      call mpi_barrier(lbecomm,ierr)
!
! free memory...
      deallocate(buffer)
!
# ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. restore_mpiio"
      endif
# endif
#endif
!
! format
4000    format(i8.8)
!
      return
      end subroutine restore_mpiio
