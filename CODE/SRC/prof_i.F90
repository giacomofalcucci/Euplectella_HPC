!=======================================================================
!     ****** LBE/prof_i
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_i
!     DESCRIPTION
!       Diagnostic subroutine:
!       istantaneous profile along x-direction (j, k fixed)
!       write on unit 61 (prof_i.dat)
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
        subroutine prof_i(itime,jcoord,kcoord)
!
        use storage
        implicit none
!
        integer:: itime,i
        integer:: jcoord,kcoord
!
        real(mykind) :: u(1:l),w(1:l),v(1:l)      ! istantaneous velocity fields
        real(mykind) :: den(1:l)           ! istantaneous density field
!
        do i = 1,l         ! density
           den(i) = +a01(i,jcoord,kcoord)+a02(i,jcoord,kcoord)+a03(i,jcoord,kcoord) &
                    +a04(i,jcoord,kcoord)+a05(i,jcoord,kcoord)+a06(i,jcoord,kcoord) &
                    +a07(i,jcoord,kcoord)+a08(i,jcoord,kcoord)+a09(i,jcoord,kcoord) &
                    +a10(i,jcoord,kcoord)+a11(i,jcoord,kcoord)+a12(i,jcoord,kcoord) &
                    +a13(i,jcoord,kcoord)+a14(i,jcoord,kcoord)+a15(i,jcoord,kcoord) &
                    +a16(i,jcoord,kcoord)+a17(i,jcoord,kcoord)+a18(i,jcoord,kcoord) &
                    +a19(i,jcoord,kcoord)
        enddo
 
        do i = 1,l         ! streamwise velocity
           u(i) = +( a01(i,jcoord,kcoord)+a02(i,jcoord,kcoord)+a03(i,jcoord,kcoord) &
                    +a04(i,jcoord,kcoord)+a05(i,jcoord,kcoord) &
                    -a10(i,jcoord,kcoord)-a11(i,jcoord,kcoord)-a12(i,jcoord,kcoord) &
                    -a13(i,jcoord,kcoord)-a14(i,jcoord,kcoord) ) / den(i)
        end do

        do i = 1,l         ! spanwise velocity
           w(i) = +( a03(i,jcoord,kcoord)+a07(i,jcoord,kcoord)+a08(i,jcoord,kcoord) &
                    +a09(i,jcoord,kcoord)+a12(i,jcoord,kcoord) &
                    -a01(i,jcoord,kcoord)-a10(i,jcoord,kcoord)-a16(i,jcoord,kcoord) &
                    -a17(i,jcoord,kcoord)-a18(i,jcoord,kcoord) ) / den(i)
        end do

        do i = 1,l         ! normal_to_wall velocity
           v(i) = +( a04(i,jcoord,kcoord)+a06(i,jcoord,kcoord)+a07(i,jcoord,kcoord) &
                    +a13(i,jcoord,kcoord)+a18(i,jcoord,kcoord) &
                    -a02(i,jcoord,kcoord)-a09(i,jcoord,kcoord)-a11(i,jcoord,kcoord) &
                    -a15(i,jcoord,kcoord)-a16(i,jcoord,kcoord) ) / den(i)
        end do

        write(61,1005) itime

        do i=1,l
           write(61,1002) i+(offset(1)), u(i),w(i),v(i),den(i)
        end do
        write(61,'(a1)') 
        write(61,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_i"
        endif
#endif
!
!	format
1002    format(i5,4(e14.6,1x))
1005    format("# t=",i7)
1111    format("#pause")
        return
        end subroutine prof_i
