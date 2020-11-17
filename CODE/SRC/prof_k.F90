!=======================================================================
!     ****** LBE/prof_k
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_k
!     DESCRIPTION
!       Diagnostic subroutine:
!       istantaneous profile along z-direction (i, j fixed)
!       write on unit 60 (prof_k.dat)
!     INPUTS
!       itime   -->  timestep
!       icoord  -->  x coordinate
!       jcoord  -->  y coordinate
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: 
!       real variables used: 
!                            
!
!     *****
!=======================================================================
!
        subroutine prof_k(itime,icoord,jcoord)
!
        use storage
        implicit none
!
        integer:: itime,k
        integer:: icoord,jcoord
!
        real(mykind) :: u(1:n),w(1:n),v(1:n)      ! istantaneous velocity fields
        real(mykind) :: den(1:n)           ! istantaneous density field
!
        do k = 1,n         ! density
           den(k) = +a01(icoord,jcoord,k)+a02(icoord,jcoord,k)+a03(icoord,jcoord,k) &
                    +a04(icoord,jcoord,k)+a05(icoord,jcoord,k)+a06(icoord,jcoord,k) &
                    +a07(icoord,jcoord,k)+a08(icoord,jcoord,k)+a09(icoord,jcoord,k) &
                    +a10(icoord,jcoord,k)+a11(icoord,jcoord,k)+a12(icoord,jcoord,k) &
                    +a13(icoord,jcoord,k)+a14(icoord,jcoord,k)+a15(icoord,jcoord,k) &
                    +a16(icoord,jcoord,k)+a17(icoord,jcoord,k)+a18(icoord,jcoord,k) &
                    +a19(icoord,jcoord,k)
        enddo
 
        do k = 1,n         ! streamwise velocity
           u(k) = +( a01(icoord,jcoord,k)+a02(icoord,jcoord,k)+a03(icoord,jcoord,k) &
                    +a04(icoord,jcoord,k)+a05(icoord,jcoord,k) &
                    -a10(icoord,jcoord,k)-a11(icoord,jcoord,k)-a12(icoord,jcoord,k) &
                    -a13(icoord,jcoord,k)-a14(icoord,jcoord,k) ) / den(k)
        end do

        do k = 1,n         ! spanwise velocity
           w(k) = +( a03(icoord,jcoord,k)+a07(icoord,jcoord,k)+a08(icoord,jcoord,k) &
                    +a09(icoord,jcoord,k)+a12(icoord,jcoord,k) &
                    -a01(icoord,jcoord,k)-a10(icoord,jcoord,k)-a16(icoord,jcoord,k) &
                    -a17(icoord,jcoord,k)-a18(icoord,jcoord,k) ) / den(k)
        end do

        do k = 1,n         ! normal_to_wall velocity
           v(k) = +( a04(icoord,jcoord,k)+a06(icoord,jcoord,k)+a07(icoord,jcoord,k) &
                    +a13(icoord,jcoord,k)+a18(icoord,jcoord,k) &
                    -a02(icoord,jcoord,k)-a09(icoord,jcoord,k)-a11(icoord,jcoord,k) &
                    -a15(icoord,jcoord,k)-a16(icoord,jcoord,k) ) / den(k)
        end do

        write(60,1005) itime

        do k=1,n
           write(60,1002) k+(offset(3)), u(k),w(k),v(k),den(k)
        end do
        write(60,'(a1)') 
        write(60,'(a1)') 

#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_k"
        endif
#endif


!	format
1002    format(i5,4(e14.6,1x))
1005    format("# t=",i7)
        return
        end subroutine prof_k
