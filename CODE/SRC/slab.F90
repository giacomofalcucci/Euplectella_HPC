!=======================================================================
!     ****** LBE/slab
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       slab
!     DESCRIPTION
!       boundary condition for a  slab of size = 1 lattice unit
!       the slab is located at "kcoord"  from the wall
!       the slab is long iend-istart lattice units
!       A simple bounce-back is used
!     INPUTS
!       kcoord --> k location             (absolute value)
!       istart --> beginnign of the slab  (absolute value)
!       iend   --> end of the slab        (absolute value)
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       obsolete
!       integer variables used: i, kcoord, istart, iend
!
!     *****
!=======================================================================
!
        subroutine slab(kcoord, istart, iend)
!
        use storage
        implicit none
!
        integer:: i, k, istart,iend, kcoord
        integer:: offsetX, offsetZ
!
        offsetX = mpicoords(1)*l
        offsetZ = mpicoords(2)*n
!
#ifdef QQQQQ
! no slip wall 
        do k = 1, n
          if((k+offsetZ)==(kcoord)) then
            do i = 1, l
              if(((i+offsetX).lt.(iend)).and.((i+offsetX).gt.(istart))) then
                a15(i  ,k) = a06(i,k-1)     ! v_z = +1
                a02(i-1,k) = a13(i,k-1)
                a11(i+1,k) = a04(i,k-1)

                a06(i  ,k) = a15(i,k+1)     ! v_z = -1
                a04(i-1,k) = a11(i,k+1)
                a13(i+1,k) = a02(i,k+1)
              endif
            enddo
          endif
        enddo
!
        do k = 1, n
          if((k+offsetZ)==(kcoord)) then
!
!  fix for istart 
            do i = 1, l
              if((i+offsetX)==(istart)) then
                a14(i,k) = a05(i-1,k)
              endif
!
!  fix for iend 
              if((i+offsetX)==(iend)) then
                a05(i  ,k) = a14(i+1,k)
              endif
            enddo
          endif
        enddo
!
!	write(6,*) "end sub. bcond, "
!
#endif
        return
        end subroutine slab
