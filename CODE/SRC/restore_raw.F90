!=======================================================================
!     ****** LBE/restore_raw
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       restore
!     DESCRIPTION
!       restore all microscopic variables (all populations)
!       read from unit 30 (restore.xx.yy.zz.bin, unformatted)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!      
!     NOTES
!       integer itime,i,k
!
!     *****
!=======================================================================
!
      subroutine restore_raw(itime)
!
      use storage
      use timing
      implicit none
!
      character*23 file_name
!
      integer itime,i,k,j
!
      file_name = 'restore.xx.xx.xx.bin'
      write(file_name(9:10),4100) mpicoords(1)
      write(file_name(12:13),4100) mpicoords(2)
      write(file_name(15:16),4100) mpicoords(3)
      open(30,file=file_name,form="unformatted",status='unknown')
!
      write(16,*) 'INFO: I am task', myrank, 'restoring from ', file_name
      write(6,*)  'INFO: I am task', myrank, 'restoring from ', file_name
!
      read(30) itime
!
      read(30) (((a01(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a02(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a03(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a04(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a05(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a06(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a07(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a08(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a09(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a10(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a11(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a12(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a13(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a14(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a15(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a16(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a17(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a18(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      read(30) (((a19(i,j,k),i=0,l1),j=0,m1),k=0,n1)
!
      close(30)
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. restore_raw"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. restore_raw mem =", mem_stop
      endif
# endif
!

4100  format(i2.2)

      return
      end subroutine restore_raw
