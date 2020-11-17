!=====================================================================
!     ****** LBE/flux_X
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       flux_X
!     DESCRIPTION
!       diagnostic subroutine:
!       compute sigma (see Artoli article) along 2 planes for drag 
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       fluxX
!     TODO
!
!     NOTES
!       integer variables used: 
!       real variables used: 
!
!     *****
!=====================================================================
!
       subroutine flux_X(itime,istart,istop,jstart,jstop,kstart,kstop,fluxX)
!
       use storage
       implicit none
!
       integer:: itime,i,j,k,ii,ll
       integer:: i0
       integer:: istart, istop       ! i1 < 12
       integer:: jstart, jstop       ! j1 < j2
       integer:: kstart, kstop       ! k1 < k2
       integer:: i1, i2       ! i1 < 12
       integer:: j1, j2       ! j1 < j2
       integer:: k1, k2       ! k1 < k2
       integer, dimension(1:2) :: face
!
       real(mykind) ::  cvsq
       real(mykind) ::  rho,xj,yj,zj,rhoinv
       real(mykind) ::  sigma_xx,sigma_yy,sigma_zz
       real(mykind) ::  sigma_xy,sigma_yz,sigma_xz
       real(mykind) ::  offsetX,offsetY,offsetZ
       real(mykind) ::  fluxX
       real(mykind), dimension(1:19) ::  fe
       real(mykind), dimension(1:19) ::  fn
       real(mykind), dimension(1:2) ::  flux
       real(mykind), dimension(1:2) ::  norm
!
       offsetX = mpicoords(1)*l
       offsetY = mpicoords(2)*m
       offsetZ = mpicoords(3)*n
!
! some normalization
       norm(1)=1.0
       norm(2)=1.0
       flux(1)=0.0
       flux(2)=0.0
! 
! set the right loop values
! ------------------------- X plane --------------------------------
! case A: i1+i2 inside the task
       if((istart.gt.offsetX).and.(istop.le.(offsetX+l))) then
         i1 = istart-offsetX
         i2 = istop -offsetX
!         write(6,*) "A: Here we are", istart, istop, i1, i2, myrank
       endif
!
! case B: only i1 inside the task
       if((istart.gt.offsetX).and.(istop.gt.(offsetX+l))) then
         i1=istart-offsetX
         i2=l
         norm(2)=0.0
!         write(6,*) "B: Here we are", istart, istop, i1, i2, myrank
       endif
!
! case C: only i2 inside the task
       if((istart.le.offsetX).and.(istop.le.(offsetX+l))) then 
         i1=1
         norm(1)=0.0
         i2=istop-offsetX
!         write(6,*) "C: Here we are", istart, istop, i1, i2, myrank
       endif
!
! case D: both outside the task (but to compute)
       if((istart.lt.offsetX).and.(istop.gt.(offsetX+l))) then
!         write(6,*) "case D", myrank
         return
       endif
!
! case E: both outside the task (lower/upper)
       if((istart.gt.(offsetX+l)).or.(istop.lt.offsetX)) then
!         write(6,*) "case E", myrank
         return
       endif
!
! ------------------------- Y plane --------------------------------
! A: j1+j2 inside the task
       if((jstart.gt.offsetY).and.(jstop.le.(offsetY+m))) then
         j1 = jstart-offsetY
         j2 = jstop -offsetY
       endif
!
! case B: only j1 inside the task
       if((jstart.gt.offsetY).and.(jstop.gt.(offsetY+m))) then
         j1=jstart-offsetY
         j2=m
!         write(6,*) "B: Here we are", jstart, jstop, j1, j2, myrank
       endif
!
! case C: only j2 inside the task
       if((jstart.le.offsetY).and.(jstop.le.(offsetY+m))) then 
         j1 =1
         j2 = jstop-offsetY
!         write(6,*) "C: Here we are", jstart, jstop, j1, j2, myrank
       endif
!
! case D: both outside the task (but to compute)
       if((jstart.lt.offsetY).and.(jstop.gt.(offsetY+m))) then
         j1 = 1
         j2 = m
       endif
!
! case E: both outside the task (lower/upper)
       if((jstart.gt.(offsetY+m)).or.(jstop.lt.offsetY)) then
         return
       endif
!
! ------------------------- Z plane --------------------------------
! A: k1+k2 inside the task
       if((kstart.gt.offsetZ).and.(kstop.le.(offsetZ+n))) then
         k1 = kstart-offsetZ
         k2 = kstop -offsetZ
       endif
!
! B: only k1 inside the task
       if((kstart.gt.offsetZ).and.(kstop.gt.(offsetZ+n))) then
         k1 = kstart-offsetZ
         k2 = n
!         write(6,*) "B: Here we are", kstart, kstop, k1, k2
       endif
!
! C: only k2 inside the task
       if((kstart.le.offsetZ).and.(kstop.le.(offsetZ+n))) then 
         k1 = 1
         k2 = kstop-offsetZ
!         write(6,*) "C: Here we are", kstart, kstop, k1, k2
       endif
!
! D: both outside the task (but to compute)
       if((kstart.lt.offsetZ).and.(kstop.gt.offsetZ+n)) then
         k1 = 1
         k2 = n
!         write(6,*) "D: Here we are", kstart, kstop, k1, k2
       endif
!
! E: both outside the task (lower/upper)
       if((kstart.gt.offsetZ+n).or.(kstop.lt.offsetZ)) then
         return
       endif
!
       face(1) = i1
       face(2) = i2
!
       do ii  = 1,2
          i0 =face(ii)
!!$acc parallel loop 
          do k=k1,k2
             do j=j1,j2
!
! 1) compute rho,u,v,w
!
                rho = +a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k) &
                      +a04(i0,j,k)+a05(i0,j,k)+a06(i0,j,k) &
                      +a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k) &
                      +a10(i0,j,k)+a11(i0,j,k)+a12(i0,j,k) &
                      +a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k) &
                      +a16(i0,j,k)+a17(i0,j,k)+a18(i0,j,k) &
                      +a19(i0,j,k)
!
                rhoinv = 1.0/rho
!
                xj = +( a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k) &
                       +a04(i0,j,k)+a05(i0,j,k)             & 
                       -a10(i0,j,k)-a11(i0,j,k)-a12(i0,j,k) &
                       -a13(i0,j,k)-a14(i0,j,k) )*rhoinv
!
                yj = +( a03(i0,j,k)+a07(i0,j,k)+a08(i0,j,k) &
                       +a09(i0,j,k)+a12(i0,j,k)             &
                       -a01(i0,j,k)-a10(i0,j,k)-a16(i0,j,k) &
                       -a17(i0,j,k)-a18(i0,j,k) )*rhoinv
!
                zj = +( a04(i0,j,k)+a06(i0,j,k)+a07(i0,j,k) &
                       +a13(i0,j,k)+a18(i0,j,k)             &
                       -a02(i0,j,k)-a09(i0,j,k)-a11(i0,j,k) &
                       -a15(i0,j,k)-a16(i0,j,k) )*rhoinv
!
! 2) compute f_eq
!
                cvsq=xj*xj+yj*yj+zj*zj
!
                fe(01)=rho*p2*(1.0+ rf*( xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq))
                fe(02)=rho*p2*(1.0+ rf*( xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq))
                fe(03)=rho*p2*(1.0+ rf*( xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq))
                fe(04)=rho*p2*(1.0+ rf*( xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq))
                fe(05)=rho*p1*(1.0+ rf*( xj      )+qf*(3.0*(xj   )*(xj   )-cvsq))
                fe(06)=rho*p1*(1.0+ rf*(       zj)+qf*(3.0*(zj   )*(zj   )-cvsq))
                fe(07)=rho*p2*(1.0+ rf*(    yj+zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq))
                fe(08)=rho*p1*(1.0+ rf*(    yj   )+qf*(3.0*(yj   )*(yj   )-cvsq))
                fe(09)=rho*p2*(1.0+ rf*(    yj-zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq))
                fe(10)=rho*p2*(1.0+ rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq))
                fe(11)=rho*p2*(1.0+ rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq))
                fe(12)=rho*p2*(1.0+ rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq))
                fe(13)=rho*p2*(1.0+ rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq))
                fe(14)=rho*p1*(1.0+ rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq))
                fe(15)=rho*p1*(1.0+ rf*(      -zj)+qf*(3.0*(zj   )*(zj   )-cvsq))
                fe(16)=rho*p2*(1.0+ rf*(   -yj-zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq))
                fe(17)=rho*p1*(1.0+ rf*(   -yj   )+qf*(3.0*(yj   )*(yj   )-cvsq))
                fe(18)=rho*p2*(1.0+ rf*(   -yj+zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq))
                fe(19)=rho*p0*(1.0+ rf*(   0.0   )+qf*(3.0*( 0.0 )*( 0.0 )-cvsq))
!
! 3) compute f_neq
!
                fn(01) = a01(i0,j,k) - fe(01)
                fn(02) = a02(i0,j,k) - fe(02)
                fn(03) = a03(i0,j,k) - fe(03)
                fn(04) = a04(i0,j,k) - fe(04)
                fn(05) = a05(i0,j,k) - fe(05)
                fn(06) = a06(i0,j,k) - fe(06)
                fn(07) = a07(i0,j,k) - fe(07)
                fn(08) = a08(i0,j,k) - fe(08)
                fn(09) = a09(i0,j,k) - fe(09)
                fn(10) = a10(i0,j,k) - fe(10)
                fn(11) = a11(i0,j,k) - fe(11)
                fn(12) = a12(i0,j,k) - fe(12)
                fn(13) = a13(i0,j,k) - fe(13)
                fn(14) = a14(i0,j,k) - fe(14)
                fn(15) = a15(i0,j,k) - fe(15)
                fn(16) = a16(i0,j,k) - fe(16)
                fn(17) = a17(i0,j,k) - fe(17)
                fn(18) = a18(i0,j,k) - fe(18)
!!                fn(19) = a19(i0,j,k) - fe(19)
!
! 4) compute sigma
!
                sigma_xx = -rho/rf
                do ll = 1,18
                   sigma_xx = sigma_xx - (1-0.5*omega)*cx(ll)*cx(ll)*fn(ll)
!!                sigma_xx = -rho/rf -(1-0.5*omega)*cx( 1)*cx( 1)*fn( 1) &
!!                                   -(1-0.5*omega)*cx( 2)*cx( 2)*fn( 2) &
!!                                   -(1-0.5*omega)*cx( 3)*cx( 3)*fn( 3) &
!!                                   -(1-0.5*omega)*cx( 4)*cx( 4)*fn( 4) &
!!                                   -(1-0.5*omega)*cx( 5)*cx( 5)*fn( 5) &
!!                                   -(1-0.5*omega)*cx( 6)*cx( 6)*fn( 6) &
!!                                   -(1-0.5*omega)*cx( 7)*cx( 7)*fn( 7) &
!!                                   -(1-0.5*omega)*cx( 8)*cx( 8)*fn( 8) &
!!                                   -(1-0.5*omega)*cx( 9)*cx( 9)*fn( 9) &
!!                                   -(1-0.5*omega)*cx(10)*cx(10)*fn(10) &
!!                                   -(1-0.5*omega)*cx(11)*cx(11)*fn(11) &
!!                                   -(1-0.5*omega)*cx(12)*cx(12)*fn(12) &
!!                                   -(1-0.5*omega)*cx(13)*cx(13)*fn(13) &
!!                                   -(1-0.5*omega)*cx(14)*cx(14)*fn(14) &
!!                                   -(1-0.5*omega)*cx(15)*cx(15)*fn(15) &
!!                                   -(1-0.5*omega)*cx(16)*cx(16)*fn(16) &
!!                                   -(1-0.5*omega)*cx(17)*cx(17)*fn(17) &
!!                                   -(1-0.5*omega)*cx(18)*cx(18)*fn(18)
                enddo
!
!!                sigma_xy = 0.0       
!!                do ll = 1,18
!!                   sigma_xy = sigma_xy - (1-0.5*omega)*cx(ll)*cy(ll)*fn(ll)
!!                enddo
!!!
!!                sigma_xz = 0.0       
!!                do ll = 1,18
!!                   sigma_xz = sigma_xz - (1-0.5*omega)*cx(ll)*cz(ll)*fn(ll)
!!                enddo
!
!!                sigma_yy = -rho/rf 
!!                do ll = 1,18
!!                   sigma_yy = sigma_yy - (1-0.5*omega)*cy(ll)*cy(ll)*fn(ll)
!!                enddo
!
!!                sigma_yz = 0.0      
!!                do ll = 1,18
!!                   sigma_yz = sigma_yz - (1-0.5*omega)*cy(ll)*cz(ll)*fn(ll)
!!                enddo
!
!!                sigma_zz = -rho/rf 
!!                do ll = 1,18
!!                   sigma_zz = sigma_zz - (1-0.5*omega)*cz(ll)*cz(ll)*fn(ll)
!!                enddo
!
!                flux(ii) = flux(ii)+(sigma_xx + sigma_xy + sigma_xz)
! ask why!!!!
                flux(ii) = flux(ii)+sigma_xx 
             enddo
          enddo
       enddo
!
       fluxX = -flux(1)*norm(1)+flux(2)*norm(2)
!
#ifdef DEBUG_1
        if(myrank == 11) then
           write(6,*) "DEBUG1: Exiting from sub. flux_X", flux(1)*norm(1),flux(2)*norm(2)
           write(6,*) "DEBUG1:", itime, myrank, flux(1)*norm(1),flux(2)*norm(2)
           write(6,*) "DEBUG2:", itime, myrank, i1,i2
           write(6,*) "DEBUG3:", itime, myrank, j1,j2
           write(6,*) "DEBUG4:", itime, myrank, k1,k2
        endif
#endif

      return
      end subroutine flux_X
