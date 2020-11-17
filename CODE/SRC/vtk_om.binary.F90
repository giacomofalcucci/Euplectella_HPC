!=======================================================================
!     ****** LBE/vtk_om_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_om_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write binary output for VTK with vorticity (stream field to do)
!       write on unit 55 (tec_om.yyyy.xxxxxxx.dat) where 
!                                yyyy is the task id 
!                                     xxxxxxx is the timestep 
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
!       integer variables used: itime,i,k
!       max time allowed  99'999'999
!       max task allowed  9999
!       single precision only, to saving space..
!
!     BUGS
!
!     *****
!=======================================================================
!
        subroutine vtk_om_bin(itime,k0)
!
        use storage
        implicit none
!
        integer i,j,k0,itime
!	
        character*24 file_name
!
        real(sp) ::  u,v,w
        real(sp) ::  den, uv, vv
        real(sp) ::  up1, down1, left1, right1
        real(sp) ::  d_up, d_down, d_left, d_right
        real(sp) ::  vorticity
        real(sp) ::  phi(-1:m1)
!
        file_name = 'tec_om.xxxx.xxxxxxxx.vtk'
!
#ifdef SERIAL
        myrank = 0
#endif
!
        write(file_name( 8:11),3100) myrank
        write(file_name(13:20),4000) itime
        open(55,file=file_name,status='unknown')

!for now only vorticity is computed
#ifdef QQQQQ
!
! pre-compute phi ...
!
        phi(-1) = 0.0
        i = 1
        do j = 0, m
            den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0) & 
                 +a05(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0) & 
                 +a09(i,j,k0)+a10(i,j,k0)+a11(i,j,k0)+a12(i,j,k0) & 
                 +a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0)+a16(i,j,k0) & 
                 +a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
           uv =   a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0) & 
                 +a05(i,j,k0)-a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0) & 
                 -a13(i,j,k0)-a14(i,j,k0)
!
           phi(j) =  phi(j-1) + uv/den
        enddo
#endif
!
        write(55,'(A26)')'# vtk DataFile Version 2.0'
        write(55,'(A5)')'Campo'
        write(55,'(A6)')'BINARY'
        write(55,'(A25)')'DATASET STRUCTURED_POINTS'
        write(55,'(A11,I10,A1,I10,A1,I10)') &
                              'DIMENSIONS ',l-2,' ',m-2,' ',1
        write(55,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',offset(1)+2,' ' &
                                                     ,offset(2)+2,' ' &
                                                     ,offset(3)+1
        write(55,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(55,'(A10,I10)')'POINT_DATA ',(l-2)*(m-2)
        write(55,'(A24)')'SCALARS vorticity float'
        write(55,'(A20)')'LOOKUP_TABLE default'
        close(55)
!
! then write output (binary)
        open(55,file=file_name,status='old', position='append', &
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        do j = 2,m-1
           do i = 2,l-1
!
            d_up= a01(i,j+1,k0)+a02(i,j+1,k0)+a03(i,j+1,k0)+a04(i,j+1,k0) & 
                 +a05(i,j+1,k0)+a06(i,j+1,k0)+a07(i,j+1,k0)+a08(i,j+1,k0) & 
                 +a09(i,j+1,k0)+a10(i,j+1,k0)+a11(i,j+1,k0)+a12(i,j+1,k0) & 
                 +a13(i,j+1,k0)+a14(i,j+1,k0)+a15(i,j+1,k0)+a16(i,j+1,k0) & 
                 +a17(i,j+1,k0)+a18(i,j+1,k0)+a19(i,j+1,k0)
!
           d_down=a01(i,j-1,k0)+a02(i,j-1,k0)+a03(i,j-1,k0)+a04(i,j-1,k0) & 
                 +a05(i,j-1,k0)+a06(i,j-1,k0)+a07(i,j-1,k0)+a08(i,j-1,k0) & 
                 +a09(i,j-1,k0)+a10(i,j-1,k0)+a11(i,j-1,k0)+a12(i,j-1,k0) & 
                 +a13(i,j-1,k0)+a14(i,j-1,k0)+a15(i,j-1,k0)+a16(i,j-1,k0) & 
                 +a17(i,j-1,k0)+a18(i,j-1,k0)+a19(i,j-1,k0)
!
           d_left=a01(i-1,j,k0)+a02(i-1,j,k0)+a03(i-1,j,k0)+a04(i-1,j,k0) & 
                 +a05(i-1,j,k0)+a06(i-1,j,k0)+a07(i-1,j,k0)+a08(i-1,j,k0) & 
                 +a09(i-1,j,k0)+a10(i-1,j,k0)+a11(i-1,j,k0)+a12(i-1,j,k0) & 
                 +a13(i-1,j,k0)+a14(i-1,j,k0)+a15(i-1,j,k0)+a16(i-1,j,k0) & 
                 +a17(i-1,j,k0)+a18(i-1,j,k0)+a19(i-1,j,k0)
!
          d_right=a01(i+1,j,k0)+a02(i+1,j,k0)+a03(i+1,j,k0)+a04(i+1,j,k0) & 
                 +a05(i+1,j,k0)+a06(i+1,j,k0)+a07(i+1,j,k0)+a08(i+1,j,k0) & 
                 +a09(i+1,j,k0)+a10(i+1,j,k0)+a11(i+1,j,k0)+a12(i+1,j,k0) & 
                 +a13(i+1,j,k0)+a14(i+1,j,k0)+a15(i+1,j,k0)+a16(i+1,j,k0) & 
                 +a17(i+1,j,k0)+a18(i+1,j,k0)+a19(i+1,j,k0)
!
           up1 =  a01(i,j+1,k0)+a02(i,j+1,k0)+a03(i,j+1,k0)+a04(i,j+1,k0) & 
                 +a05(i,j+1,k0)-a10(i,j+1,k0)-a11(i,j+1,k0)-a12(i,j+1,k0) & 
                 -a13(i,j+1,k0)-a14(i,j+1,k0)/d_up
!
          down1 = a01(i,j-1,k0)+a02(i,j-1,k0)+a03(i,j-1,k0)+a04(i,j-1,k0) & 
                 +a05(i,j-1,k0)-a10(i,j-1,k0)-a11(i,j-1,k0)-a12(i,j-1,k0) & 
                 -a13(i,j-1,k0)-a14(i,j-1,k0)/d_down
!
          left1 = a03(i-1,j,k0)+a07(i-1,j,k0)+a08(i-1,j,k0)+a09(i-1,j,k0) &
                 +a12(i-1,j,k0)-a01(i-1,j,k0)-a10(i-1,j,k0)-a16(i-1,j,k0) & 
                 -a17(i-1,j,k0)-a18(i-1,j,k0)/d_left
!
          right1= a03(i+1,j,k0)+a07(i+1,j,k0)+a08(i+1,j,k0)+a09(i+1,j,k0) &
                 +a12(i+1,j,k0)-a01(i+1,j,k0)-a10(i+1,j,k0)-a16(i+1,j,k0) & 
                 -a17(i+1,j,k0)-a18(i+1,j,k0)/d_right
!
! delta_x = delta_y = 2.0
              vorticity = 0.5*(up1-down1) - 0.5*(right1-left1)  
              write(55) vorticity
!
           end do
        end do
!
! it works, commented only for comparison with binary one
#ifdef QQQQ
        write(55,'(A25)')'SCALARS stream double'
        write(55,'(A20)')'LOOKUP_TABLE default'
        do j = 2,m-1
           do i = 2,l-1
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0) &
                   +a05(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0) &
                   +a09(i,j,k0)+a10(i,j,k0)+a11(i,j,k0)+a12(i,j,k0) &
                   +a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0)+a16(i,j,k0) &
                   +a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)
!
              vv  = a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0) & 
                   +a12(i,j,k0)-a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0) & 
                   -a17(i,j,k0)-a18(i,j,k0)
!
               phi(j) = phi(j) - vv/den
!
              write(55,1004) phi(j)
           end do
        end do
#endif
!
        close(55)
        write(16,*) "vorticity xy (vtk) done"
        write(6,*) "vorticity xy (vtk) done"
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_om_bin"
        endif
#endif

1004    format((e14.6,1x))
4000    format(i8.8)
3100    format(i4.4)

       return
       end subroutine vtk_om_bin

