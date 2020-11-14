!=====================================================================
!     ****** LBE/bcond_bc
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
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
        subroutine bcond_bc
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
#ifdef DRIVEN
        force =  u00/(6.0)
#endif
!
#ifdef COUETTE 
        force =  u00/(6.0)
#endif
!
#ifdef NOMANAGED
!$acc data present(a01,a02,a03,a04,a05,a06, & 
!$acc              a07,a08,a09,a10,a11,a12, & 
!$acc              a13,a14,a15,a16,a17,a18) 
#endif
!
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
! periodic
        if((front(1)==0).and.(front(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc kernels
!GA!$acc loop independent
          do k = 0, n+1
!GA!$acc loop independent
             do j = 0, m+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do k = 1, n
#ifdef SERIAL
!!# ifdef OBSTACLES
! no idea why.....
             a12(l1,0,k) = a12(1,m,k)
!!# endif
# endif
!GA!$acc loop independent
             do j = 1, m
#endif
                a10(l1,j,k) = a10(1,j,k) 
                a11(l1,j,k) = a11(1,j,k) 
                a12(l1,j,k) = a12(1,j,k) 
                a13(l1,j,k) = a13(1,j,k)
                a14(l1,j,k) = a14(1,j,k)
             end do
          end do
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(front(1)==1) then
!$OMP DO
!GA!$acc data copy(a12,a13,a10,a11,a14, &
!GA!$acc           a01,a02,a03,a04,a05)
!$acc parallel
!$acc loop independent
          do k = 1, n
!!$acc loop independent
             do j = 1, m
                a12(l1,j-1,k  ) = a01(l,j,k)
                a13(l1,j  ,k-1) = a02(l,j,k)
                a10(l1,j+1,k  ) = a03(l,j,k)
                a11(l1,j  ,k+1) = a04(l,j,k)
                a14(l1,j  ,k  ) = a05(l,j,k)
             end do
          end do
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! outflow 
        if(front(1)==5) then
!$OMP DO
!GA!$acc kernels 
!GA!$acc loop independent
          do k = 1, n
!GA!$acc loop independent
             do j = 1, m
!
!                crho = 1.0
!                rhoinv = 1.0/crho
!
                xj = +( a01(l,j,k)+a02(l,j,k)+a03(l,j,k) &
                       +a04(l,j,k)+a05(l,j,k)            &
                       -a10(l,j,k)-a11(l,j,k)-a12(l,j,k) &
                       -a13(l,j,k)-a14(l,j,k) )!*rhoinv

                yj = +( a03(l,j,k)+a07(l,j,k)+a08(l,j,k) &
                       +a09(l,j,k)+a12(l,j,k)            &
                       -a01(l,j,k)-a10(l,j,k)-a16(l,j,k) &
                       -a17(l,j,k)-a18(l,j,k) )!*rhoinv

                zj = +( a04(l,j,k)+a06(l,j,k)+a07(l,j,k) &
                       +a13(l,j,k)+a18(l,j,k)            &
                       -a02(l,j,k)-a09(l,j,k)-a11(l,j,k) &
                       -a15(l,j,k)-a16(l,j,k) )!*rhoinv
!
                cvsq=xj*xj+yj*yj+zj*zj
!
                cx10 = rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                cx11 = rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                cx12 = rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                cx13 = rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                cx14 = rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
!
                a10(l1,j,k) = p2*(1.0+cx10)!*crho
                a11(l1,j,k) = p2*(1.0+cx11)!*crho
                a12(l1,j,k) = p2*(1.0+cx12)!*crho
                a13(l1,j,k) = p2*(1.0+cx13)!*crho
                a14(l1,j,k) = p1*(1.0+cx14)!*crho
          enddo
        enddo
!GA!$acc end kernels
!$OMP END DO
!
! border fix
!GA!$acc kernels 
!GA!$acc loop independent
        do j = 1, m
             a10(l1,j,0) = a10(l1,j,1)
             a11(l1,j,0) = a11(l1,j,1)
             a12(l1,j,0) = a12(l1,j,1)
             a13(l1,j,0) = a13(l1,j,1)
             a14(l1,j,0) = a14(l1,j,1)
!
             a10(l1,j,n1) = a10(l1,j,n)
             a11(l1,j,n1) = a11(l1,j,n)
             a12(l1,j,n1) = a12(l1,j,n)
             a13(l1,j,n1) = a13(l1,j,n)
             a14(l1,j,n1) = a14(l1,j,n)
        enddo
!GA!$acc end kernels
!
!GA!$acc kernels
!GA!$acc loop independent
        do k = 1, n
             a10(l1,0,k) = a10(l1,1,k)
             a11(l1,0,k) = a11(l1,1,k)
             a12(l1,0,k) = a12(l1,1,k)
             a13(l1,0,k) = a13(l1,1,k)
             a14(l1,0,k) = a14(l1,1,k)
!
             a10(l1,m1,k) = a10(l1,m,k)
             a11(l1,m1,k) = a11(l1,m,k)
             a12(l1,m1,k) = a12(l1,m,k)
             a13(l1,m1,k) = a13(l1,m,k)
             a14(l1,m1,k) = a14(l1,m,k)
        enddo
!GA!$acc end kernels
      endif
! ----------------------------------------------
! rear wall (x = 0)
!
! periodic
        if((rear(1)==0).and.(rear(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc loop independent
          do k = 0, n+1
!GA!$acc loop independent
             do j = 0, m+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do k = 1, n
#ifdef SERIAL
!!# ifdef OBSTACLES
             a01(0,m1,k) = a01(l,1,k)
!!# endif
#endif
!GA!$acc loop independent
             do j = 1, m
#endif
                a01(0,j,k)  = a01(l,j,k)   
                a02(0,j,k)  = a02(l,j,k)   
                a03(0,j,k)  = a03(l,j,k)   
                a04(0,j,k)  = a04(l,j,k)
                a05(0,j,k)  = a05(l,j,k)
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(rear(1)==1) then
!$OMP DO
!GA!$acc data copy(a12,a13,a10,a11,a14,a01,a02,a03,a04,a05)
!$acc parallel
!$acc loop independent
          do k = 1, n
!!$acc loop independent
             do j = 1, m
                a03(0,j-1,k  ) = a10(1,j,k)
                a04(0,j  ,k-1) = a11(1,j,k)
                a01(0,j+1,k  ) = a12(1,j,k)
                a02(0,j  ,k+1) = a13(1,j,k)
                a05(0,j  ,k  ) = a14(1,j,k)
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! inflow (fixed velocity)
        if(rear(1)==4) then
!
          xj = u_inflow
          yj = 0.0
          zj = 0.0
!
!$OMP DO
!GA!$acc kernels
!GA!$acc loop independent
          do k = 0, n1
             z = (real(k,mykind) + mpicoords(3)*n-0.5 - 0.5*real(lz,mykind))/(0.5*real(lz,mykind))
!
#ifdef CHANNEL
! poiseuille inflow
             xj = u_inflow  * ( 1 - z*z)              
#endif
!
#ifdef INLET
! poiseuille inflow
             xj = u_inflow  * ( 1 - z*z)              
#endif
!
#ifdef MIXING
             if (z.gt.0.25) then
                xj = 0.2
             else
                xj = 0.02
             endif
#endif
!GA!$acc loop independent
             do j = 0, m1
!
#ifdef DUCT
                y =(real(j)+mpicoords(2)*m-0.5-0.5*float(ly))/(0.5*float(ly))
                xj = u_inflow  * ( 1 - z*z)*( 1 - y*y)
#endif
!
#ifdef CHECK2D
! poiseuille (fixed)
                y =(real(j)+mpicoords(2)*m-0.5-0.5*float(ly))/(0.5*float(ly))
                xj = u_inflow  * ( 1 - y*y)              
#endif
!
                cvsq=xj*xj+yj*yj+zj*zj
!
!               crho = 1.0
!
                cx01 = rf*(xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                cx02 = rf*(xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                cx03 = rf*(xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                cx04 = rf*(xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                cx05 = rf*(xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
!
                a01(0,j,k) = p2*(1.0+cx01)!*crho
                a02(0,j,k) = p2*(1.0+cx02)!*crho
                a03(0,j,k) = p2*(1.0+cx03)!*crho
                a04(0,j,k) = p2*(1.0+cx04)!*crho
                a05(0,j,k) = p1*(1.0+cx05)!*crho
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
!
        endif
!
! inflow (pouiseuille)
!!!        if(rear(1)==4) then
!!! to fix 
!!!          zj = 0.0
!!!          do k = 1, n
!!!             z = (real(k) + mpicoords(2)*n - 0.5 - 0.5*float(lz))/(0.5*float(lz))
!!!             xj = u00 * ( 1 - z*z)                !poiseuille (fixed)
!!!             vsq=xj*xj+zj*zj
!!!          
!!!             x02 = rf*( xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-vsq)
!!!             x04 = rf*( xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-vsq)
!!!             x05 = rf*( xj      )+qf*(3.0*(xj   )*(xj   )-vsq)
!!!
!!!             a02(0,k) = rho*p2*(1.0+x02)
!!!             a04(0,k) = rho*p2*(1.0+x04)
!!!             a05(0,k) = rho*p1*(1.0+x05)
!!!          enddo
!!!          write(6,*) "BC-2: still to fix"
!!!             stop
!!!        endif
!
! ----------------------------------------------
! up (k = n) wall
! periodic
        if((up(1)==0).and.(up(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc kernels
!GA!$acc loop independent
          do j = 0, m+1
!GA!$acc loop independent
             do i = 0, l+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do j = 1, m
!GA!$acc loop independent
             do i = 1, l
#endif
             a02(i,j,n1)  = a02(i,j,1)
             a09(i,j,n1)  = a09(i,j,1)
             a11(i,j,n1)  = a11(i,j,1)
             a15(i,j,n1)  = a15(i,j,1)   
             a16(i,j,n1)  = a16(i,j,1)   
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(up(1)==1) then
!$OMP DO
!GA!$acc data copy(a15,a02,a09,a11,a16,a06,a13,a18,a04,a07)
!$acc parallel
!$acc loop independent
          do j = 1, m
!GA!$acc loop independent
             do i = 1, l
                a15(i  ,j  ,n1) = a06(i,j,n)
                a02(i-1,j  ,n1) = a13(i,j,n) 
                a09(i  ,j-1,n1) = a18(i,j,n) 
                a11(i+1,j  ,n1) = a04(i,j,n) 
                a16(i  ,j+1,n1) = a07(i,j,n) 
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! lid-wall
        if(up(1)==2) then
!$OMP DO
!GA!$acc data copy(a06,a13,a18,a04,a07,a15,a02,a09,a11,a16)
!$acc parallel
!$acc loop independent
          do j = 1, m
!GA!$acc loop independent
             do i = 1, l
                a15(i  ,j  ,n1) = a06(i,j,n)
                a02(i-1,j  ,n1) = a13(i,j,n) + force
                a09(i  ,j-1,n1) = a18(i,j,n) 
                a11(i+1,j  ,n1) = a04(i,j,n) - force
                a16(i  ,j+1,n1) = a07(i,j,n) 
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! free-slip
        if(up(1)==3) then
!
! to understand why...
!$OMP DO
!GA!$acc kernels
!GA!$acc loop independent
          do j = 0, m1
!GA!$acc loop independent
             do i = 0, l1
                a15(i,j,n1) = a06(i,j,n)
                a11(i,j,n1) = a13(i,j,n)
                a16(i,j,n1) = a18(i,j,n)
                a02(i,j,n1) = a04(i,j,n)
                a09(i,j,n1) = a07(i,j,n)
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
! ----------------------------------------------
! down (k = 0) wall: 
! periodic
        if((down(1)==0).and.(down(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc kernels
!GA!$acc loop independent
          do j = 0, m+1
!$acc loop independent
             do i = 0, l+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do j = 1, m
!GA!$acc loop independent
             do i = 1, l
#endif
                a04(i,j,0) = a04(i,j,n)
                a06(i,j,0) = a06(i,j,n)   
                a07(i,j,0) = a07(i,j,n)   
                a13(i,j,0) = a13(i,j,n)
                a18(i,j,0) = a18(i,j,n)
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(down(1)==1) then
!$OMP DO
!GA!$acc data copy(a04,a06,a07,a13,a18,a15,a11,a16,a02,a09)
!$acc parallel
!$acc loop independent
          do j = 1, m
!GA!$acc loop independent
             do i = 1, l
                a06(i  ,j  ,0)  = a15(i,j,1)   
                a04(i-1,j  ,0)  = a11(i,j,1)
                a07(i  ,j-1,0)  = a16(i,j,1)
                a13(i+1,j  ,0)  = a02(i,j,1)
                a18(i  ,j+1,0)  = a09(i,j,1)
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! free-slip
        if(down(1)==3) then
!$OMP DO
          do j = 1, m
             do i = 1, l
                a06(i,j,0)  = a15(i,j,1)
                a13(i,j,0)  = a11(i,j,1)
                a18(i,j,0)  = a16(i,j,1)
                a04(i,j,0)  = a02(i,j,1)
                a07(i,j,0)  = a09(i,j,1)
             enddo
          enddo
!$OMP END DO
        endif
!
! ----------------------------------------------
! left (y = 0) wall: 
! periodic
        if((left(1)==0).and.(left(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc kernels
!GA!$acc loop independent
          do k = 0, n+1
!GA!$acc loop independent
             do i = 0, l+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do k = 1, n
             a03(0,0,k) = a03(l,m,k)
!GA!$acc loop independent
             do i = 1, l
#endif
                a03(i,0,k) = a03(i,m,k)
                a07(i,0,k) = a07(i,m,k)     
                a08(i,0,k) = a08(i,m,k)     
                a09(i,0,k) = a09(i,m,k)
                a12(i,0,k) = a12(i,m,k)
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(left(1)==1) then
!$OMP DO
!GA!$acc data copy(a08,a12,a03,a07,a09,a17,a01,a10,a16,a18) 
!$acc parallel
!$acc loop independent
          do k = 1, n
!GA!$acc loop independent
             do i = 1, l
                a08(i  ,0,k  )  = a17(i,1,k)
                a12(i+1,0,k  )  = a01(i,1,k)
                a03(i-1,0,k  )  = a10(i,1,k)
                a07(i  ,0,k-1)  = a16(i,1,k)
                a09(i  ,0,k+1)  = a18(i,1,k)
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! free-slip
        if(left(1)==3) then
!$OMP DO
          do k = 1, n
             do i = 1, l
                a08(i,0,k)  = a17(i,1,k)
                a03(i,0,k)  = a01(i,1,k)
                a12(i,0,k)  = a10(i,1,k)
                a09(i,0,k)  = a16(i,1,k)
                a07(i,0,k)  = a18(i,1,k)
             enddo
          enddo
!$OMP END DO
        endif

!
! ----------------------------------------------
! right (y = m) wall: 
! periodic
        if((right(1)==0).and.(right(2)==myrank)) then
!$OMP DO
#ifdef PERIODIC
!GA!$acc kernels
!GA!$acc loop independent
          do k = 0, n+1
!GA!$acc loop independent
             do i = 0, l+1
#else
!GA!$acc kernels
!GA!$acc loop independent
          do k = 1, n
             a10(l1,m1,k) = a10(1,1,k)
!GA!$acc loop independent
             do i = 1, l
#endif
                a01(i,m1,k) = a01(i,1,k)
                a10(i,m1,k) = a10(i,1,k)      
                a16(i,m1,k) = a16(i,1,k)      
                a17(i,m1,k) = a17(i,1,k)
                a18(i,m1,k) = a18(i,1,k)
             enddo
          enddo
!GA!$acc end kernels
!$OMP END DO
        endif
!
! noslip
        if(right(1)==1) then
!$OMP DO
!GA!$acc data copy(a10,a16,a17,a01,a18,a03,a07,a08,a09,a12)
!$acc parallel
!$acc loop independent
          do k = 1, n
!GA!$acc loop independent
             do i = 1, l
                a10(i+1,m1,k  ) = a03(i,m,k)
                a16(i  ,m1,k+1) = a07(i,m,k)
                a17(i  ,m1,k  ) = a08(i,m,k)
                a18(i  ,m1,k-1) = a09(i,m,k)
                a01(i-1,m1,k  ) = a12(i,m,k)
             enddo
          enddo
!$acc end parallel
!GA!$acc end data 
!$OMP END DO
        endif
!
! free-slip
        if(right(1)==3) then
!$OMP DO
          do k = 1, n
             do i = 1, l
                a01(i,m1,k) = a03(i,m,k)
                a18(i,m1,k) = a07(i,m,k)     
                a17(i,m1,k) = a08(i,m,k)     
                a16(i,m1,k) = a09(i,m,k)
                a10(i,m1,k) = a12(i,m,k)
             enddo
          enddo
!$OMP END DO
        endif

!$OMP END PARALLEL
!
!
#ifdef NOMANAGED
!$acc end data
#endif
!
! fix for edge for free-slip
! fix 1
        if((left(1)==3).and.(up(1)==3)) then
          do i = 1, l
             a09(i,0,n1) = a18(i,1,n)
          enddo
        endif
!
! fix 2
        if((right(1)==3).and.(up(1)==3)) then
          do i = 1, l
             a16(i,m1,n1) = a07(i,m,n)
          enddo
        endif
!
! fix 3
        if((right(1)==3).and.(down(1)==3)) then
          do i = 1, l
             a18(i,m1,0) = a09(i,m,1)
          enddo
        endif
!
! fix 4
        if((left(1)==3).and.(down(1)==3)) then
          do i = 1, l
             a07(i,0,0) = a16(i,1,1)
          enddo
        endif
!
! fix for free-sleep: to clean as soon as possible
#ifdef WWWWW
#ifdef FREESLIP
        if((front(1)==0)) then
! 
! 4 corner to fix
          a09(l1, 0,n1) = a18(l1,1,n)              ! corner 1
          a16(l1,m1,n1) = a07(l1,m,n)              ! corner 2
          a18(l1,m1, 0) = a09(l1,m,1)              ! corner 3
          a07(l1, 0, 0) = a16(l1,1,1)              ! corner 4
! 
! four (2x2) edge to fix
          do k = 1, n
             a08(l1,0,k) = a17(l1,1,k)
             a03(l1,0,k) = a01(l1,1,k)
             a12(l1,0,k) = a10(l1,1,k)
             a09(l1,0,k) = a16(l1,1,k)
             a07(l1,0,k) = a18(l1,1,k)
!
             a01(l1,m1,k) = a03(l1,m,k)
             a18(l1,m1,k) = a07(l1,m,k)     
             a17(l1,m1,k) = a08(l1,m,k)     
             a16(l1,m1,k) = a09(l1,m,k)
             a10(l1,m1,k) = a12(l1,m,k)
          enddo
!
          do j = 1, m
             a06(l1,j,0)  = a15(l1,j,1)
             a13(l1,j,0)  = a11(l1,j,1)
             a18(l1,j,0)  = a16(l1,j,1)
             a04(l1,j,0)  = a02(l1,j,1)
             a07(l1,j,0)  = a09(l1,j,1)
!
             a15(l1,j,n1) = a06(l1,j,n)
             a11(l1,j,n1) = a13(l1,j,n)
             a16(l1,j,n1) = a18(l1,j,n)
             a02(l1,j,n1) = a04(l1,j,n)
             a09(l1,j,n1) = a07(l1,j,n)
          enddo
        endif
#else
! do nothing
#endif
!
! fix for free-sleep: to clean as soon as possible
#ifdef FREESLIP
        if((rear(1)==0)) then
! 4 corner to fix
          a09(0, 0,n1) = a18(0,1,n)              ! corner 1
          a16(0,m1,n1) = a07(0,m,n)              ! corner 2
          a18(0,m1, 0) = a09(0,m,1)              ! corner 3
          a07(0, 0, 0) = a16(0,1,1)              ! corner 4
! 
! four (2x2) edge to fix
          do k = 1, n
             a08(0,0,k) = a17(0,1,k)
             a03(0,0,k) = a01(0,1,k)
             a12(0,0,k) = a10(0,1,k)
             a09(0,0,k) = a16(0,1,k)
             a07(0,0,k) = a18(0,1,k)
!
             a01(0,m1,k) = a03(0,m,k)
             a18(0,m1,k) = a07(0,m,k)
             a17(0,m1,k) = a08(0,m,k)
             a16(0,m1,k) = a09(0,m,k)
             a10(0,m1,k) = a12(0,m,k)
          enddo
!
          do j = 1, m
             a06(0,j,0)  = a15(0,j,1)
             a13(0,j,0)  = a11(0,j,1)
             a18(0,j,0)  = a16(0,j,1)
             a04(0,j,0)  = a02(0,j,1)
             a07(0,j,0)  = a09(0,j,1)
!
             a15(0,j,n1) = a06(0,j,n)
             a11(0,j,n1) = a13(0,j,n)
             a16(0,j,n1) = a18(0,j,n)
             a02(0,j,n1) = a04(0,j,n)
             a09(0,j,n1) = a07(0,j,n)
          enddo
        endif
#else
! do nothing
#endif
!
! fix for free-sleep: to clean as soon as possible
#ifdef FREESLIP
        if((right(1)==0)) then
! 4 corner to fix
          a11(l1,m1,n1) = a04(l,m1,n)              ! corner 1
          a13(l1,m1, 0) = a02(l,m1,1)              ! corner 2
          a04( 0,m1, 0) = a11(1,m1,1)              ! corner 3
          a02( 0,m1,n1) = a13(1,m1,n)              ! corner 4
! 
! four (2x2) edge to fix
!
          do k = 1, n
             a01(0,m1,k) = a10(1,m1,k)
             a02(0,m1,k) = a11(1,m1,k)
             a03(0,m1,k) = a12(1,m1,k)
             a04(0,m1,k) = a13(1,m1,k)
             a05(0,m1,k) = a14(1,m1,k)
!
             a10(l1,m1,k) = a01(l,m1,k)
             a11(l1,m1,k) = a02(l,m1,k)
             a12(l1,m1,k) = a03(l,m1,k)
             a13(l1,m1,k) = a04(l,m1,k)
             a14(l1,m1,k) = a05(l,m1,k)
          enddo
!
          do i = 1, l
             a04(i,m1,0) = a02(i,m1,1)
             a07(i,m1,0) = a09(i,m1,1)
             a13(i,m1,0) = a11(i,m1,1)
             a06(i,m1,0) = a15(i,m1,1)
             a18(i,m1,0) = a16(i,m1,1)
!
             a02(i,m1,n1) = a04(i,m1,n)
             a09(i,m1,n1) = a07(i,m1,n)
             a11(i,m1,n1) = a13(i,m1,n)
             a15(i,m1,n1) = a06(i,m1,n)
             a16(i,m1,n1) = a18(i,m1,n)
          enddo
        endif
#else
! do nothing
#endif
!
!
! fix for free-sleep: to clean as soon as possible
#ifdef FREESLIP
        if((left(1)==0)) then
! 4 corner to fix
          a11(l1,0,n1) = a04(l,0,n)              ! corner 1
          a13(l1,0, 0) = a02(l,0,1)              ! corner 2
          a04( 0,0, 0) = a11(1,0,1)              ! corner 3
          a02( 0,0,n1) = a13(1,0,n)              ! corner 4
! 
! four (2x2) edge to fix
!
          do k = 1, n
             a01(0,0,k) = a10(1,0,k)
             a02(0,0,k) = a11(1,0,k)
             a03(0,0,k) = a12(1,0,k)
             a04(0,0,k) = a13(1,0,k)
             a05(0,0,k) = a14(1,0,k)
!
             a10(l1,0,k) = a01(l,0,k)
             a11(l1,0,k) = a02(l,0,k)
             a12(l1,0,k) = a03(l,0,k)
             a13(l1,0,k) = a04(l,0,k)
             a14(l1,0,k) = a05(l,0,k)
          enddo
!
          do i = 1, l
             a04(i,0,0) = a02(i,0,1)
             a07(i,0,0) = a09(i,0,1)
             a13(i,0,0) = a11(i,0,1)
             a06(i,0,0) = a15(i,0,1)
             a18(i,0,0) = a16(i,0,1)
!
             a02(i,0,n1) = a04(i,0,n)
             a09(i,0,n1) = a07(i,0,n)
             a11(i,0,n1) = a13(i,0,n)
             a15(i,0,n1) = a06(i,0,n)
             a16(i,0,n1) = a18(i,0,n)
          enddo
        endif
#else
! do nothing
#endif
!
! fix for free-sleep: to clean as soon as possible
#ifdef FREESLIP
        if((up(1)==0)) then
! 4 corner to fix
          a10(l1,m1,n1) = a03(l,m,n1)              ! corner 1
          a12(l1, 0,n1) = a01(l,1,n1)              ! corner 2
          a03( 0, 0,n1) = a10(1,1,n1)              ! corner 3
          a01( 0,m1,n1) = a12(1,m,n1)              ! corner 4
! 
! four (2x2) edge to fix
!
          do i = 1, l
             a03(i,0,n1) = a01(i,1,n1)
             a12(i,0,n1) = a10(i,1,n1)
             a09(i,0,n1) = a16(i,1,n1)
             a08(i,0,n1) = a17(i,1,n1)
             a07(i,0,n1) = a18(i,1,n1)
!
             a01(i,m1,n1) = a03(i,m,n1)
             a10(i,m1,n1) = a12(i,m,n1)
             a16(i,m1,n1) = a09(i,m,n1)
             a17(i,m1,n1) = a08(i,m,n1)
             a18(i,m1,n1) = a07(i,m,n1)
          enddo
!
          do j = 1, m
             a01(0,j,n1) = a10(1,j,n1)
             a02(0,j,n1) = a11(1,j,n1)
             a03(0,j,n1) = a12(1,j,n1)
             a04(0,j,n1) = a13(1,j,n1)
             a05(0,j,n1) = a14(1,j,n1)
!
             a10(l1,j,n1) = a01(l,j,n1)
             a11(l1,j,n1) = a02(l,j,n1)
             a12(l1,j,n1) = a03(l,j,n1)
             a13(l1,j,n1) = a04(l,j,n1)
             a14(l1,j,n1) = a05(l,j,n1)
          enddo
!
        endif
#else
! do nothing
#endif
!
! fix for free-sleep: to clean as soon as possible
#ifdef FREESLIP
        if((down(1)==0)) then
! 4 corner to fix
          a10(l1,m1,0) = a03(l,m,0)              ! corner 1
          a12(l1, 0,0) = a01(l,1,0)              ! corner 2
          a03( 0, 0,0) = a10(1,1,0)              ! corner 3
          a01( 0,m1,0) = a12(1,m,0)              ! corner 4
! 
! four (2x2) edge to fix
!
          do i = 1, l
             a03(i,0,0) = a01(i,1,0)
             a12(i,0,0) = a10(i,1,0)
             a09(i,0,0) = a16(i,1,0)
             a08(i,0,0) = a17(i,1,0)
             a07(i,0,0) = a18(i,1,0)
!
             a01(i,m1,0) = a03(i,m,0)
             a10(i,m1,0) = a12(i,m,0)
             a16(i,m1,0) = a09(i,m,0)
             a17(i,m1,0) = a08(i,m,0)
             a18(i,m1,0) = a07(i,m,0)
          enddo
!
          do j = 1, m
             a01(0,j,0) = a10(1,j,0)
             a02(0,j,0) = a11(1,j,0)
             a03(0,j,0) = a12(1,j,0)
             a04(0,j,0) = a13(1,j,0)
             a05(0,j,0) = a14(1,j,0)
!
             a10(l1,j,0) = a01(l,j,0)
             a11(l1,j,0) = a02(l,j,0)
             a12(l1,j,0) = a03(l,j,0)
             a13(l1,j,0) = a04(l,j,0)
             a14(l1,j,0) = a05(l,j,0)
          enddo
!
        endif
#else
! do nothing
#endif
#endif
!
! fix for serial/openacc with freeslip (hardwired fix)
! (only for upper free-slip wall)
#ifdef SPONGES
# ifdef SERIAL
!GA!$acc kernels 
!GA!$acc loop independent
     do i = 1,l
        a16(i,m1,n1)= a18(i,1,n)
!
        a09(i,0,n1) = a07(i,m,n)
     enddo
!GA!$acc end kernels 
# endif
#endif
!
!#ifdef OBSTACLES
!       call bcond_obs
!       call bcond_comm_noblock_packed_try
!#endif
!  corner fix...
!
! warning: this fix is necessary? 
! it is commented beacuse in this way omp-mpi version produce the same results
! but I'm not sure if it is correct.....
!
#ifdef DRIVEN
# ifdef SERIAL
!!        to port
!!        a11(l1,n1) = a04(l,n) - force
!!        a02(0 ,n1) = a13(1,n) + force
!!        a11(l1,n1) = a04(l,n) 
!!        a02(0 ,n1) = a13(1,n) 
# endif
!!         a13(l1,0 ) = a02(l,1)
!!         a04(0 ,0 ) = a11(1,1)
#endif

#ifdef CHANNEL
# ifdef SERIAL
!!          do j = 1, m
!!             a02(0 ,j,n1) = a02(l,j,1)
!!             a04(0 ,j,0 ) = a04(l,j,n)
!!             a11(l1,j,n1) = a11(1,j,1)
!!             a13(l1,j,0 ) = a13(1,j,n)
!!          enddo
# endif
#endif

#ifdef PERIODIC
# ifdef SERIAL
!!        to port
!!        a02(0 ,n1) = a02(l,1)
!!        a04(0 ,0 ) = a04(l,n)
!!        a11(l1,n1) = a11(1,1)
!!        a13(l1,0 ) = a13(1,n)
# endif
#endif

!!!! to check
#ifdef INLET	
# ifdef SERIAL
!!        to port
!!        a11(l1,n1) = a04(l,n)
!!        a02(0 ,n1) = a13(1,n)
!!        a13(l1,0 ) = a02(l,1)
!!        a04(0 ,0 ) = a11(1,1)
# endif
#endif

!
! set obstacles...
!!        to port
!!#ifdef SLAB
!!        call slab(lz/2,lx/4-lx/8,3*lx/4-lx/8)! in the middle of the box..
!!#endif
!!
!!#ifdef WALL
!!        call wall(lx/4,lz/4,3*lz/4)   
!!#endif
!!
!!#ifdef CILINDER
!!        call cilinder(lx/4,lz/2,lz/8)   ! in the middle of the box..
!!#endif
!!!
!!$acc end data
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_bc = time_bc + real(countA1-countA0)/(count_rate)
        time_bc1 = time_bc1 + (tcountA1-tcountA0)
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond.bc"
           write(6,*) " "
        endif
#endif
        return
        end subroutine bcond_bc
