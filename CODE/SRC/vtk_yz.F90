!=======================================================================
!     ****** LBE/vtk_yz
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_yz
!     DESCRIPTION
!       Graphic subroutine:
!       write ASCII output for VTK with velocity + pressure field
!       write on unit 52 (vtk_xxxx.xxxxxxx.dat, where xxxxxxx is the timestep)
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
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
        subroutine vtk_yz(itime)
!
        use storage
        implicit none
!
        integer i0,j,k,itime
!	
        character*24 file_name
!
        real(sp) :: u,w,v,den
!
        file_name = 'tec_yz.xxxx.xxxxxxxx.vtk'
!
#ifdef SERIAL
        myrank = 0
#endif
        i0 = l/2

        write(file_name(8:11),3100) myrank
        write(file_name(13:20),4000) itime
        open(52,file=file_name,status='unknown')
!
        write(52,'(A26)')'# vtk DataFile Version 2.0'
        write(52,'(A5)')'Campo'
        write(52,'(A5)')'ASCII'
!        write(52,'(A5)')'BINARY'
        write(52,'(A24)')'DATASET RECTILINEAR_GRID'
        write(52,'(A11,I10,A1,I10,A1,I10)')  'DIMENSIONS ',m,' ',n,' ',1
!
        write(52,'(A14,I10,A7)')'X_COORDINATES ',m,' double'
        do j = 1,m
           write(52, *) j + offset(2)
        enddo
!
        write(52,'(A14,I10,A7)')'Y_COORDINATES ',n,' double'
        do k = 1,n
           write(52, *) k + offset(3)
        enddo
!
        write(52,'(A14,I10,A7)')'Z_COORDINATES ',1,' double'
        write(52, *) 0
!
        write(52,'(A10,I10)')'POINT_DATA ',m*n*1
        write(52,'(A24)')'VECTORS velocity double'
!        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m

              u = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                 -a10(i0,j,k)-a11(i0,j,k)-a12(i0,j,k)-a13(i0,j,k)-a14(i0,j,k)
!
              w = a03(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a12(i0,j,k) &
                 -a01(i0,j,k)-a10(i0,j,k)-a16(i0,j,k)-a17(i0,j,k)-a18(i0,j,k)
!
              v = a04(i0,j,k)+a06(i0,j,k)+a07(i0,j,k)+a13(i0,j,k)+a18(i0,j,k) &
                 -a02(i0,j,k)-a09(i0,j,k)-a11(i0,j,k)-a15(i0,j,k)-a16(i0,j,k)
!
              den = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                   +a06(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a10(i0,j,k) &
                   +a11(i0,j,k)+a12(i0,j,k)+a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k) &
                   +a16(i0,j,k)+a17(i0,j,k)+a18(i0,j,k)+a19(i0,j,k)
!
              write(52,1004) u/den, w/den, v/den
           end do
        end do
!
        write(52,'(A21)')'SCALARS w double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m

              u = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                 -a10(i0,j,k)-a11(i0,j,k)-a12(i0,j,k)-a13(i0,j,k)-a14(i0,j,k)
!
              den = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                   +a06(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a10(i0,j,k) &
                   +a11(i0,j,k)+a12(i0,j,k)+a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k) &
                   +a16(i0,j,k)+a17(i0,j,k)+a18(i0,j,k)+a19(i0,j,k)
!
              write(52,1004) w/den
           end do
        end do
!
        write(52,'(A21)')'SCALARS v double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
!
              v = a04(i0,j,k)+a06(i0,j,k)+a07(i0,j,k)+a13(i0,j,k)+a18(i0,j,k) &
                 -a02(i0,j,k)-a09(i0,j,k)-a11(i0,j,k)-a15(i0,j,k)-a16(i0,j,k)
!
              den = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                   +a06(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a10(i0,j,k) &
                   +a11(i0,j,k)+a12(i0,j,k)+a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k) &
                   +a16(i0,j,k)+a17(i0,j,k)+a18(i0,j,k)+a19(i0,j,k)
!
              write(52,1004) v/den
           end do
        end do
!
        write(52,'(A21)')'SCALARS rho double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
!
              den = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                   +a06(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a10(i0,j,k) &
                   +a11(i0,j,k)+a12(i0,j,k)+a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k) &
                   +a16(i0,j,k)+a17(i0,j,k)+a18(i0,j,k)+a19(i0,j,k)
!
              write(52,1004) den
           end do
        end do
!
        close(52)
        write(16,*) "plane yz (vtk) done"
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_yz"
        endif
#endif


1004    format(3(e14.6,1x))
3100    format(i4.4)
4000    format(i8.8)

       return 
       end subroutine vtk_yz

