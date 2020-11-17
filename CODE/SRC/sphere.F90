!=======================================================================
!     ****** LBE/sphere
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       sphere
!     DESCRIPTION
!       sphere: obstcle
!     INPUTS
!       icoord --> i center of cilinder (absolute value)
!       jcoord --> j center of cilinder (absolute value)
!       kcoord --> k center of cilinder (absolute value)
!       radius --> radius of cilinder (absolute value) 
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       obsolete
!
!     *****
!=======================================================================
!
        subroutine sphere(icoord, jcoord, kcoord)
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
        integer:: i, j, k, icoord, jcoord, kcoord
        integer:: offsetX, offsetY, offsetZ
        real(mykind):: d2, R2a, R2b, R
!
        R = radius*1.0
        R2a = (R-2)*(R-2)
        R2b = (R+2)*(R+2)
!
        offsetX = mpicoords(1)*l
        offsetY = mpicoords(2)*m
        offsetZ = mpicoords(3)*n
!
!        write(*,*) "task", myrank, " -->", offsetX, offsetZ, offsetY, & 
!                                           icoord, jcoord, kcoord, R
!
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(d2)  &
!$OMP FIRSTPRIVATE(offsetZ,offsetY,offsetX,R2a,R2b)  &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09)  &
!$OMP SHARED(a10,a11,a12,a13,a14,a15,a16,a17,a18)  &
!$OMP SHARED(icoord,jcoord,kcoord) &
!$OMP SHARED(l,m,n)
!$OMP DO
        do k = 1, n
           do j = 1, m
              do i = 1, l
                 d2 = (icoord-(i+offsetX))*(icoord-(i+offsetX))  & 
                     +(jcoord-(j+offsetY))*(jcoord-(j+offsetY))  & 
                     +(kcoord-(k+offsetZ))*(kcoord-(k+offsetZ))
                 if((d2.gt.R2a).and.(d2.lt.R2b)) then
                     a12(i,j,k) = a01(i-1,j+1,k  )
                     a13(i,j,k) = a02(i-1,j  ,k+1)
                     a10(i,j,k) = a03(i-1,j-1,k  )
                     a11(i,j,k) = a04(i-1,j  ,k-1)
                     a14(i,j,k) = a05(i-1,j  ,k  )
                     a15(i,j,k) = a06(i  ,j  ,k-1)
                     a16(i,j,k) = a07(i  ,j-1,k-1)
                     a17(i,j,k) = a08(i  ,j-1,k  )
                     a18(i,j,k) = a09(i  ,j-1,k+1)
                     a03(i,j,k) = a10(i+1,j+1,k  )
                     a04(i,j,k) = a11(i+1,j  ,k+1)
                     a01(i,j,k) = a12(i+1,j-1,k  )
                     a02(i,j,k) = a13(i+1,j  ,k-1)
                     a05(i,j,k) = a14(i+1,j  ,k  )
                     a06(i,j,k) = a15(i  ,j  ,k+1)
                     a07(i,j,k) = a16(i  ,j+1,k+1)
                     a08(i,j,k) = a17(i  ,j+1,k  )
                     a09(i,j,k) = a18(i  ,j+1,k-1)
                 endif
              end do
           end do
        end do
!$OMP END DO
!$OMP END PARALLEL
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. sphere"
           write(6,*) " "
        endif
#endif
!
        return
        end subroutine sphere
