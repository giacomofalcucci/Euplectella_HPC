!=======================================================================
!     ****** LBE/vtk_xy_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_xy_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write 2D binary output for VTK with velocity + pressure field
!       write on unit 52 (vtk_yyyy.xxxxxxx.dat) where 
!                                  yyyy is the task id
!                                  xxxxxxx is the timestep
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!       k0      ---> k coordinate of the 2D slice
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       character used:  file_name (19)
!       integer variables used: itime,i,k, x_scale, z_scale
!       max time allowed  99'999'999
!       single precision only, to saving space..
!
!     BUGS
!
!     *****
!=======================================================================
!
        subroutine vtk_xy_bin(itime,k0)
!
        use storage
        implicit none
!
        integer i,j,k0,itime
!	
        character*24 file_name
!
        real(sp) :: u,w,v,den
!
        file_name = 'tec_xy.xxxx.xxxxxxxx.vtk'
!
#ifdef SERIAL
        myrank = 0
#endif
        write(file_name(8:11),3100) myrank
        write(file_name(13:20),4000) itime
!
! first write legal header (ASCII)
        open(52,file=file_name,status='unknown')
!
        write(52,'(A26)')'# vtk DataFile Version 2.0'
        write(52,'(A5)') 'Campo'
        write(52,'(A6)') 'BINARY'
        write(52,'(A25)')'DATASET STRUCTURED_POINTS'
        write(52,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',l,' ',m,' ',1
        write(52,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',offset(1)+1,' ' & 
                                                     ,offset(2)+1,' ' & 
                                                     ,offset(3)+1+k0
        write(52,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(52,'(A10,I10)')'POINT_DATA ',l*m*1
        write(52,'(A23)')'VECTORS velocity float'
        close(52)
!
! then write output (binary)
        open(52,file=file_name,status='old', position='append', & 
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        do j = 1,m
           do i = 1,l
              u = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                 -a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0)-a13(i,j,k0)-a14(i,j,k0)
!
              w = a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a12(i,j,k0) &
                 -a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0)-a17(i,j,k0)-a18(i,j,k0)
!
              v = a04(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a13(i,j,k0)+a18(i,j,k0) &
                 -a02(i,j,k0)-a09(i,j,k0)-a11(i,j,k0)-a15(i,j,k0)-a16(i,j,k0)
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0) &
                   +a05(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0) &
                   +a09(i,j,k0)+a10(i,j,k0)+a11(i,j,k0)+a12(i,j,k0) &
                   +a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0)+a16(i,j,k0) &
                   +a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
              write(52) u/den, w/den, v/den
           end do
        end do
!
#ifdef QQQQ
        stop
!        write(52,'(A21)')'SCALARS u double'
!        write(52,'(A20)')'LOOKUP_TABLE default'
        do j = 1,m
           do i = 1,l

              u = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                 -a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0)-a13(i,j,k0)-a14(i,j,k0)
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                   +a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a10(i,j,k0) &
                   +a11(i,j,k0)+a12(i,j,k0)+a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0) &
                   +a16(i,j,k0)+a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
!              write(52,1004) u/den
           end do
        end do
!
!        write(52,'(A21)')'SCALARS w double'
!        write(52,'(A20)')'LOOKUP_TABLE default'
        do j = 1,m
           do i = 1,l
!
              w = a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a12(i,j,k0) &
                 -a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0)-a17(i,j,k0)-a18(i,j,k0)
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                   +a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a10(i,j,k0) &
                   +a11(i,j,k0)+a12(i,j,k0)+a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0) &
                   +a16(i,j,k0)+a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
!              write(52,1004) w/den
           end do
        end do
!
!        write(52,'(A21)')'SCALARS rho double'
!        write(52,'(A20)')'LOOKUP_TABLE default'
        do j = 1,m
           do i = 1,l
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                   +a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a10(i,j,k0) &
                   +a11(i,j,k0)+a12(i,j,k0)+a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0) &
                   +a16(i,j,k0)+a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
!              write(52,1004) den
           end do
        end do
#endif
!
        close(52)
        write(16,*) "plane xy (vtk,binary) done, plane z=", k0
        write(6,*)  "plane xy (vtk,binary) done, plane z=", k0
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_xy_bin"
        endif
#endif


1004    format(3(e14.6,1x))
3100    format(i4.4)
4000    format(i8.8)

       return 
       end subroutine vtk_xy_bin

