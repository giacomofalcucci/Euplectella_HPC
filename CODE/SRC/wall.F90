!=======================================================================
!     ****** LBE/wall
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       wall
!     DESCRIPTION
!       boundary condition for vertical wall of size = 1 lattice unit
!       the wall is located at distance equal to "icoord" 
!       from the channel inlet and is "ktop-kdown" high
!       A simple bounce-back is used
!     INPUTS
!       icoord --> i location (absolute value)
!       ktop   --> k location of upper part of the wall (absolute value)
!       kdown  --> k location of lower part of the wall (absolute value)
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       obsolete
!       integer variables used: k, icoord, ktop, kdown
!
!     *****
!=======================================================================
!
        subroutine wall(icoord, kdown, ktop)
!
        use storage
        implicit none
!
        integer:: i, k, icoord, ktop, kdown
        integer:: offsetX, offsetZ
!
        offsetX = mpicoords(1)*l
        offsetZ = mpicoords(2)*n
!
#ifdef QQQQQ
! no slip wall 
        do i = 1, l
          if((i+offsetX)==(icoord)) then
            do k = 1, n
              if(((k+offsetZ).lt.(ktop)).and.((k+offsetZ).gt.(kdown))) then
                a11(i,k+1) = a04(i-1,k)  ! v_x = +1
                a13(i,k-1) = a02(i-1,k)
                a14(i,k  ) = a05(i-1,k)
!
                a02(i,k+1) = a13(i+1,k)  ! v_x = -1
                a04(i,k-1) = a11(i+1,k)
                a05(i,k  ) = a14(i+1,k)
              endif
            enddo
          endif
        end do
!
! border fix.... (to fix)
        do i = 1, l
          if((i+offsetX)==(icoord)) then
!
!  fix for ktop  (to fix)
            do k = 1, n
              if((k+offsetZ)==(ktop)) then
                a13(i,k-1) = a02(i-1,k)
                a14(i,k  ) = a05(i-1,k)
                a04(i,k-1) = a11(i+1,k)
                a05(i,k  ) = a14(i+1,k)
              endif
            enddo
!
!  fix for kdown 
            do k = 1, n
              if((k+offsetZ)==(kdown)) then
                a11(i,k+1) = a04(i-1,k)
                a14(i,k  ) = a05(i-1,k)
                a02(i,k+1) = a13(i+1,k)  
                a05(i,k  ) = a14(i+1,k)
              endif
            enddo
!
!  fix for ktop+1 
            do k = 1, n
              if((k+offsetZ)==(ktop)) then
                 a13(i,k)   = a02(i-1,k+1)
                 a06(i,k)   = a15(i  ,k+1)
                 a04(i,k)   = a11(i+1,k+1)
              endif
            enddo
!
!  fix for kdown-1 
            do k = 1, n
              if((k+offsetZ)==(kdown)) then
                a11(i,k)  = a04(i-1,k-1)
                a15(i,k)  = a06(i  ,k-1)
                a02(i,k)  = a13(i+1,k-1)  
              endif
            enddo
          endif
        enddo
!
!	write(6,*) "end sub. bcond, "
!
#endif
        return
        end subroutine wall
