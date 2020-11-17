!=====================================================================
!     ****** LBE/bcond_driven
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       ad hoc subroutine for (serial/mpi) driven cavity
!
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
        subroutine bcond_driven
!
        use timing
        use storage
#ifdef SERIAL
! do nothing
#else
        use mpi
#endif
!
        implicit none
!
        integer      :: i,j,k 
        integer      :: tag, ierr
        real(mykind) :: force
        real(mykind) :: den,xj,yj,zj,rho,vsq,z,y,rhoinv
        real(mykind) :: x01,x02,x03,x04,x05,x06
        real(mykind) :: x07,x08,x09,x10,x11,x12
        real(mykind) :: x13,x14,x15,x16,x17,x18
        real(mykind) :: cvsq,crho
        real(mykind) :: cx01,cx02,cx03,cx04,cx05
        real(mykind) :: cx10,cx11,cx12,cx13,cx14
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
!
!!$acc data present(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19) 
!
        force =  u00/(6.0)
!
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(xj,yj,zj)  &
!$OMP PRIVATE(z,y)  &
!$OMP PRIVATE(cx01,cx02,cx03,cx04,cx05)  &
!$OMP PRIVATE(cx10,cx11,cx12,cx13,cx14)  &
!$OMP PRIVATE(rhoinv,crho,cvsq)  &
!$OMP FIRSTPRIVATE(front,rear,right,left,up,down,myrank)  &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09)  &
!$OMP SHARED(a10,a11,a12,a13,a14,a15,a16,a17,a18,a19)  &
!$OMP SHARED(l,m,n,l1,m1,n1,lz,ly) &
!$OMP SHARED(mpicoords) &
!$OMP SHARED(u_inflow) &
!$OMP SHARED(force) 
! ----------------------------------------------
! front wall (x = l)
! noslip
        if(front(1)==1) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do k = 1, n
!$acc loop independent
             do j = 1, m
                a12(l1,j-1,k  ) = a01(l,j,k)
                a13(l1,j  ,k-1) = a02(l,j,k)
                a10(l1,j+1,k  ) = a03(l,j,k)
                a11(l1,j  ,k+1) = a04(l,j,k)
                a14(l1,j  ,k  ) = a05(l,j,k)
             end do
          end do
!$acc end kernels
!$OMP END DO
        endif
!
! ----------------------------------------------
! rear wall (x = 0)
! noslip
        if(rear(1)==1) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do k = 1, n
!$acc loop independent
             do j = 1, m
                a03(0,j-1,k  ) = a10(1,j,k)
                a04(0,j  ,k-1) = a11(1,j,k)
                a01(0,j+1,k  ) = a12(1,j,k)
                a02(0,j  ,k+1) = a13(1,j,k)
                a05(0,j  ,k  ) = a14(1,j,k)
             enddo
          enddo
!$acc end kernels
!$OMP END DO
        endif
!
! ----------------------------------------------
! up (k = n) wall
! lid-wall
        if(up(1)==2) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do j = 1, m
!$acc loop independent
             do i = 1, l
                a15(i  ,j  ,n1) = a06(i,j,n)
                a02(i-1,j  ,n1) = a13(i,j,n) + force
                a09(i  ,j-1,n1) = a18(i,j,n) 
                a11(i+1,j  ,n1) = a04(i,j,n) - force
                a16(i  ,j+1,n1) = a07(i,j,n) 
             enddo
          enddo
!$acc end kernels
!$OMP END DO
        endif
!
! ----------------------------------------------
! down (k = 0) wall: 
! noslip
        if(down(1)==1) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do j = 1, m
!$acc loop independent
             do i = 1, l
                a06(i  ,j  ,0)  = a15(i,j,1)   
                a04(i-1,j  ,0)  = a11(i,j,1)
                a07(i  ,j-1,0)  = a16(i,j,1)
                a13(i+1,j  ,0)  = a02(i,j,1)
                a18(i  ,j+1,0)  = a09(i,j,1)
             enddo
          enddo
!$acc end kernels
!$OMP END DO
        endif
!
! ----------------------------------------------
! left (y = 0) wall: 
! noslip
        if(left(1)==1) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do k = 1, n
!$acc loop independent
             do i = 1, l
                a08(i  ,0,k  )  = a17(i,1,k)
                a12(i+1,0,k  )  = a01(i,1,k)
                a03(i-1,0,k  )  = a10(i,1,k)
                a07(i  ,0,k-1)  = a16(i,1,k)
                a09(i  ,0,k+1)  = a18(i,1,k)
             enddo
          enddo
!$acc end kernels
!$OMP END DO
        endif
!
! ----------------------------------------------
! right (y = m) wall: 
! noslip
        if(right(1)==1) then
!$OMP DO
!$acc kernels
!$acc loop independent
          do k = 1, n
!$acc loop independent
             do i = 1, l
                a10(i+1,m1,k  ) = a03(i,m,k)
                a16(i  ,m1,k+1) = a07(i,m,k)
                a17(i  ,m1,k  ) = a08(i,m,k)
                a18(i  ,m1,k-1) = a09(i,m,k)
                a01(i-1,m1,k  ) = a12(i,m,k)
             enddo
          enddo
!$acc end kernels
!$OMP END DO
        endif
!
!$OMP END PARALLEL
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond.driven"
           write(6,*) " "
        endif
#endif
!
        return
        end subroutine bcond_driven
