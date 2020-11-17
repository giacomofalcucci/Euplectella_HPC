!=======================================================================
!     ****** LBE/cilinder
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       cilinder
!     DESCRIPTION
!       cilinder: obstcle
!     INPUTS
!       icoord --> i center of cilinder (absolute value)
!       kcoord --> k center of cilinder (absolute value)
!       radius --> radius of cilinder (absolute value) 
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       obslete..
!     *****
!=======================================================================
!
        subroutine cilinder(icoord, kcoord, radius)
!
        use storage
        implicit none
!
        integer:: i, k, icoord, kcoord, radius
        integer:: offsetX, offsetZ
        real(mykind):: d2, R2a, R2b, R
!
        R = radius*1.0
        R2a = (R-1)*(R-1)
        R2b = (R+1)*(R+1)
!
        offsetX = mpicoords(1)*l
        offsetZ = mpicoords(2)*n
!
#ifdef QQQQQ
!        write(*,*) "task", myrank, " -->", offsetX, offsetZ, icoord, kcoord,R
!
        do k = 1, n
           do i = 1, l
              d2 = (icoord-(i+offsetX))*(icoord-(i+offsetX))  & 
                  +(kcoord-(k+offsetZ))*(kcoord-(k+offsetZ))
              if((d2.gt.R2a).and.(d2.lt.R2b)) then
                a11(i,k+1) = a04(i-1,k)! v_x = +1
                a13(i,k-1) = a02(i-1,k)
                a14(i,k  ) = a05(i-1,k)
!
                a02(i,k+1) = a13(i+1,k)! v_x = -1
                a04(i,k-1) = a11(i+1,k)
                a05(i,k  ) = a14(i+1,k)
              endif
           end do
        end do
!
!	write(6,*) "end sub. bcond, "
!
#endif
        return
        end subroutine cilinder
