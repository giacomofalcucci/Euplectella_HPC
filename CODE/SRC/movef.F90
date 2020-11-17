!====================================================
!     ****** LBE/movef
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       movef
!     DESCRIPTION
!       "In place" streaming of populations
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,k
!
!     *****
!====================================================
!
subroutine movef
!
        use storage
#ifdef OPENACCOLD
        use openacc
#endif
!
        implicit none
!
        integer i,j,k
!
!------------------------------------------------------
! Best (?) decompsition for BW
!------------------------------------------------------
!!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(2) & 
!$OMP PARALLEL DEFAULT(SHARED)& 
!$OMP PRIVATE(i,j,k) 
!$OMP DO 
!$acc kernels
!$acc loop independent
         do k = 1, n
!$acc loop independent
            do j = 1, m
!$acc loop independent
               do i = 1, l
                  b01(i,j,k) = a01(i-1,j+1,k  )
                  b02(i,j,k) = a02(i-1,j  ,k+1)
                  b03(i,j,k) = a03(i-1,j-1,k  )
                  b04(i,j,k) = a04(i-1,j  ,k-1)
                  b05(i,j,k) = a05(i-1,j  ,k  )
                  b06(i,j,k) = a06(i  ,j  ,k-1)
                  b07(i,j,k) = a07(i  ,j-1,k-1)
                  b08(i,j,k) = a08(i  ,j-1,k  )
                  b09(i,j,k) = a09(i  ,j-1,k+1)
                  b10(i,j,k) = a10(i+1,j+1,k  )
                  b11(i,j,k) = a11(i+1,j  ,k+1)
                  b12(i,j,k) = a12(i+1,j-1,k  )
                  b13(i,j,k) = a13(i+1,j,  k-1)
                  b14(i,j,k) = a14(i+1,j  ,k  )
                  b15(i,j,k) = a15(i  ,j  ,k+1)
                  b16(i,j,k) = a16(i  ,j+1,k+1)
                  b17(i,j,k) = a17(i  ,j+1,k  )
                  b18(i,j,k) = a18(i  ,j+1,k-1)
               enddo
            enddo
         enddo
!$OMP END PARALLEL
!$acc end kernels


!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. move"
        endif
#endif
!
        return
        end subroutine movef
