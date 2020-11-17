!====================================================
!     ****** LBE/movef_TF
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       movef_TF
!     DESCRIPTION
!       "In place" streaming of populations
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       experimental
!       integer variables used: i,k
!
!     *****
!====================================================
!
subroutine movef_TF
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
!$OMP PARALLEL PRIVATE(i,k)

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b01(i,j,k  ) = a01(i-1,j+1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b02(i,j,k  ) = a02(i-1,j,k+1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b03(i,j,k  ) = a03(i-1,j-1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b04(i,j  ,k) = a04(i-1,j  ,k-1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b05(i,j,k  ) = a05(i-1,j,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b06(i,j  ,k) = a06(i,j  ,k-1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b07(i,j,k) = a07(i,j-1,k-1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b08(i,j,k  ) = a08(i,j-1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b09(i,j,k) = a09(i,j-1,k+1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b10(i,j,k  ) = a10(i+1,j+1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b11(i,j  ,k) = a11(i+1,j  ,k+1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b12(i,j,k  ) = a12(i+1,j-1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b13(i,j  ,k) = a13(i+1,j,  k-1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b14(i,j,k  ) = a14(i+1,j,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b15(i,j  ,k) = a15(i,j  ,k+1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b16(i,j,k) = a16(i,j+1,k+1)
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b17(i,j,k  ) = a17(i,j+1,k  )
               enddo
            enddo
         enddo
!$omp end do nowait

!$omp do
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  b18(i,j,k) = a18(i,j+1,k-1)
               enddo
            enddo
         enddo
!$omp end do
!$omp end parallel

#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. movef_TF"
        endif
#endif
!
        return
        end subroutine movef_TF
