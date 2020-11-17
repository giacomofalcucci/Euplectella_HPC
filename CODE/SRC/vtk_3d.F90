!=======================================================================
!     ****** LBE/vtk_3d
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_3d
!     DESCRIPTION
!       Graphic subroutine:
!       write ASCII output for VTK with velocity + pressure field
!       write on unit 52 (vtk_3d.yyyy.xxxxxxx.dat, 
!                         where yyyy is the task id number,
!                         where xxxxxxx is the timestep)
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       character used:  file_name (23)
!       integer variables used: itime
!       max time allowed  99'999'999
!
!     BUGS
!
!     *****
!=======================================================================
!
        subroutine vtk_3d(itime)
!
        use storage
        implicit none
!
        integer i,j, k,itime
!	
        character*24 file_name
!
        real(sp) :: u,v,w,den
!
        file_name = 'tec_3d.xxxx.xxxxxxxx.vtk'
!
#ifdef SERIAL
        myrank = 0
#endif
!
        write(file_name(8:11),3100) myrank
        write(file_name(13:20),4000) itime
        open(52,file=file_name,status='unknown')
!
        write(52,'(A26)')'# vtk DataFile Version 2.0'
        write(52,'(A5)')'Campo'
        write(52,'(A5)')'ASCII'
!        write(52,'(A5)')'BINARY'
        write(52,'(A24)')'DATASET RECTILINEAR_GRID'
        write(52,'(A11,I10,A1,I10,A1,I10)')  'DIMENSIONS ',l,' ',m,' ',n
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
!
        write(52,'(A10,I10)')'POINT_DATA ',l*m*n
        write(52,'(A24)')'VECTORS velocity double'
!        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
              do i = 1,l

              u = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                 -a10(i,j,k)-a11(i,j,k)-a12(i,j,k)-a13(i,j,k)-a14(i,j,k)
!
              w = a03(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a12(i,j,k) &
                 -a01(i,j,k)-a10(i,j,k)-a16(i,j,k)-a17(i,j,k)-a18(i,j,k)
!
              v = a04(i,j,k)+a06(i,j,k)+a07(i,j,k)+a13(i,j,k)+a18(i,j,k) &
                 -a02(i,j,k)-a09(i,j,k)-a11(i,j,k)-a15(i,j,k)-a16(i,j,k)
!
              den =a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                  +a06(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a10(i,j,k) &
                  +a11(i,j,k)+a12(i,j,k)+a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                  +a16(i,j,k)+a17(i,j,k)+a18(i,j,k)+a19(i,j,k)
!
              write(52,1004) u/den, w/den, v/den
              end do
           end do
        end do
!
        write(52,'(A21)')'SCALARS u double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
              do i = 1,l

              u = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                 -a10(i,j,k)-a11(i,j,k)-a12(i,j,k)-a13(i,j,k)-a14(i,j,k)
!
              den = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                   +a06(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a10(i,j,k) &
                   +a11(i,j,k)+a12(i,j,k)+a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                   +a16(i,j,k)+a17(i,j,k)+a18(i,j,k)+a19(i,j,k)
!
              write(52,1004) u/den
              end do
           end do
        end do
!
        write(52,'(A21)')'SCALARS v double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
              do i = 1,l

              w = a03(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a12(i,j,k) &
                 -a01(i,j,k)-a10(i,j,k)-a16(i,j,k)-a17(i,j,k)-a18(i,j,k)
!
              den = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                   +a06(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a10(i,j,k) &
                   +a11(i,j,k)+a12(i,j,k)+a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                   +a16(i,j,k)+a17(i,j,k)+a18(i,j,k)+a19(i,j,k)
!
              write(52,1004) w/den
              end do
           end do
        end do
!
        write(52,'(A21)')'SCALARS w double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
              do i = 1,l
!
              v = a04(i,j,k)+a06(i,j,k)+a07(i,j,k)+a13(i,j,k)+a18(i,j,k) &
                 -a02(i,j,k)-a09(i,j,k)-a11(i,j,k)-a15(i,j,k)-a16(i,j,k)
!
              den = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                   +a06(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a10(i,j,k) &
                   +a11(i,j,k)+a12(i,j,k)+a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                   +a16(i,j,k)+a17(i,j,k)+a18(i,j,k)+a19(i,j,k)
!
              write(52,1004) v/den
              end do
           end do
        end do
!
        write(52,'(A21)')'SCALARS rho double'
        write(52,'(A20)')'LOOKUP_TABLE default'
        do k = 1,n
           do j = 1,m
              do i = 1,l
!
              den = a01(i,j,k)+a02(i,j,k)+a03(i,j,k)+a04(i,j,k)+a05(i,j,k) &
                   +a06(i,j,k)+a07(i,j,k)+a08(i,j,k)+a09(i,j,k)+a10(i,j,k) &
                   +a11(i,j,k)+a12(i,j,k)+a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                   +a16(i,j,k)+a17(i,j,k)+a18(i,j,k)+a19(i,j,k)
!
              write(52,1004) den
              end do
           end do
        end do
!
        close(52)
        write(16,*) "3d (vtk) done"
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_3d"
        endif
#endif


1004    format(3(e14.6,1x))
3100    format(i4.4)
4000    format(i8.8)

       return 
       end subroutine vtk_3d

