!=======================================================================
!     ****** LBE/bcond_obs
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_obs
!     DESCRIPTION
!       obstcle
!     INPUTS
!       
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!     *****
!=======================================================================
!
        subroutine bcond_obs
!
        use timing
        use storage
!
#ifdef SERIAL
! do nothing
#else
        use mpi
#endif
!
        implicit none
!
        integer:: i, j, k 
        integer:: flip,flop
!
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(flip,flop)  &
!$OMP SHARED(imax,jmax,kmax,imin,jmin,kmin)  &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09)  &
!$OMP SHARED(a10,a11,a12,a13,a14,a15,a16,a17,a18)  &
!$OMP SHARED(obs,l,m,n) 
!$OMP DO
!$acc kernels
!$acc loop independent
      do k =kmin, kmax
!      do k = 1, n
!$acc loop independent
         do j = jmin, jmax
!         do j =  1, m
!$acc loop independent
             do i = imin, imax
! s           do i = 1, l
                 if(obs(i,j,k)==1) then
                     a01(i,j,k) = a12(i+1,j-1,k  )
                     a02(i,j,k) = a13(i+1,j  ,k-1)
                     a03(i,j,k) = a10(i+1,j+1,k  )
                     a04(i,j,k) = a11(i+1,j  ,k+1)
                     a05(i,j,k) = a14(i+1,j  ,k  )
                     a06(i,j,k) = a15(i  ,j  ,k+1)
                     a07(i,j,k) = a16(i  ,j+1,k+1)
                     a08(i,j,k) = a17(i  ,j+1,k  )
                     a09(i,j,k) = a18(i  ,j+1,k-1)
                     a10(i,j,k) = a03(i-1,j-1,k  )
                     a11(i,j,k) = a04(i-1,j  ,k-1)
                     a12(i,j,k) = a01(i-1,j+1,k  )
                     a13(i,j,k) = a02(i-1,j  ,k+1)
                     a14(i,j,k) = a05(i-1,j  ,k  )
                     a15(i,j,k) = a06(i  ,j  ,k-1)
                     a16(i,j,k) = a07(i  ,j-1,k-1)
                     a17(i,j,k) = a08(i  ,j-1,k  )
                     a18(i,j,k) = a09(i  ,j-1,k+1)
                 endif
              end do
           end do
        end do
!$acc end kernels
!$OMP END DO
!$OMP END PARALLEL
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_obs"
           write(6,*) " "
        endif
#endif
!
        return
        end subroutine bcond_obs
