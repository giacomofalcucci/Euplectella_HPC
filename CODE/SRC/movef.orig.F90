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
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(2) & 
!$OMP PRIVATE(i,j,k) 
!$OMP DO 
         do k = 1, n
            do j = 1, m
               do i = l, 1, -1
                  a01(i,j,k  ) = a01(i-1,j+1,k  )
               enddo
            enddo
         enddo
!$OMP END DO

!$OMP DO 
         do k = 1, n
            do j = m, 1, -1
               do i = l, 1, -1
                  a03(i,j,k  ) = a03(i-1,j-1,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n
            do j = m, 1, -1
               do i = l, 1, -1
                  a05(i,j,k  ) = a05(i-1,j,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n
            do j = m, 1, -1
               do i = l, 1, -1
                  a08(i,j,k  ) = a08(i,j-1,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n, 1
            do j = 1, m
               do i = 1, l
                  a10(i,j,k  ) = a10(i+1,j+1,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n, 1
            do j = m, 1, -1
               do i = 1, l
                  a12(i,j,k  ) = a12(i+1,j-1,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n, 1
            do j = 1, m
               do i = 1, l
                  a14(i,j,k  ) = a14(i+1,j,k  )
               enddo
            enddo
         enddo

!$OMP DO 
         do k = 1, n, 1
            do j = 1, m
               do i = 1, l
                  a17(i,j,k  ) = a17(i,j+1,k  )
               enddo
            enddo
         enddo
!$OMP END PARALLEL
!------------------------------------------------------
!$OMP SECTIONS PRIVATE(i,j,k)
!$OMP SECTION 
         do k = 1, n, 1
            do j = m, 1, -1
               do i = l, 1, -1
                  a02(i,j,k  ) = a02(i-1,j,k+1)
               enddo
            enddo
         enddo

!$OMP SECTION 
        do k = n, 1, -1
           do j = m, 1, -1
               do i = l, 1, -1
                  a04(i,j  ,k) = a04(i-1,j  ,k-1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = n, 1, -1
            do j = m, 1, -1
               do i = l, 1, -1
                  a06(i,j  ,k) = a06(i,j  ,k-1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = n, 1, -1
            do j = m, 1, -1
               do i = l, 1, -1
                  a07(i,j,k) = a07(i,j-1,k-1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = 1, n
            do j = m, 1, -1
               do i = l, 1, -1
                  a09(i,j,k) = a09(i,j-1,k+1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = 1, n
            do j = 1, m, 1
               do i = 1, l
                  a11(i,j  ,k) = a11(i+1,j  ,k+1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = n, 1, -1
            do j = 1, m, 1
               do i = 1, l
                  a13(i,j  ,k) = a13(i+1,j,  k-1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = 1, n
            do j = 1, m, 1
               do i = 1, l
                  a15(i,j  ,k) = a15(i,j  ,k+1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = 1, n
            do j = 1, m
               do i = 1, l
                  a16(i,j,k) = a16(i,j+1,k+1)
               enddo
            enddo
         enddo

!$OMP SECTION 
         do k = n, 1, -1
            do j = 1, m
               do i = 1, l
                  a18(i,j,k) = a18(i,j+1,k-1)
               enddo
            enddo
         enddo
!$OMP END SECTIONS

!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. move"
        endif
#endif
!
        return
        end subroutine movef
