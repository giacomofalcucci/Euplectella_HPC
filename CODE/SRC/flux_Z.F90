!=====================================================================
!     ****** LBE/flux_Z
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       diagno
!     DESCRIPTION
!       diagnostic subroutine:
!       compute sigma (see Artoli article)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: 
!       real variables used: 
!
!     *****
!=====================================================================
!
       subroutine flux_Z(itime,i1,i2,j1,j2,k1,k2,fluxZ)
!
       use storage
       implicit none
!
       integer:: itime,i,j,k,ii,ll
       integer:: k0  
       integer:: i1, i2       ! i1 < 12
       integer:: j1, j2       ! j1 < j2
       integer:: k1, k2       ! k1 < k2
       integer, dimension(1:2) :: face
!
       real(mykind) ::  cvsq
       real(mykind) ::  rho,xj,yj,zj,rhoinv
       real(mykind) ::  sigma_xx,sigma_yy,sigma_zz
       real(mykind) ::  sigma_xy,sigma_yz,sigma_xz
       real(mykind) ::  fluxZ
       real(mykind), dimension(1:19) ::  fe
       real(mykind), dimension(1:19) ::  fn
       real(mykind), dimension(1:2) ::  flux
!
       face(1) = k1
       face(2) = k2
!
       do ii  = 1,2
          k0 =face(ii)
          flux(ii) = 0.0
          do j=j1,j2
             do i=i1,i2
!
! 1) compute rho,u,v,w
!
                rho = +a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0) &
                      +a04(i,j,k0)+a05(i,j,k0)+a06(i,j,k0) &
                      +a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0) &
                      +a10(i,j,k0)+a11(i,j,k0)+a12(i,j,k0) &
                      +a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0) &
                      +a16(i,j,k0)+a17(i,j,k0)+a18(i,j,k0) &
                      +a19(i,j,k0)
!
                rhoinv = 1.0/rho
!
                xj = +( a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0) &
                       +a04(i,j,k0)+a05(i,j,k0)             & 
                       -a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0) &
                       -a13(i,j,k0)-a14(i,j,k0) )*rhoinv
!
                yj = +( a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0) &
                       +a09(i,j,k0)+a12(i,j,k0)             &
                       -a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0) &
                       -a17(i,j,k0)-a18(i,j,k0) )*rhoinv
!
                zj = +( a04(i,j,k0)+a06(i,j,k0)+a07(i,j,k0) &
                       +a13(i,j,k0)+a18(i,j,k0)             &
                       -a02(i,j,k0)-a09(i,j,k0)-a11(i,j,k0) &
                       -a15(i,j,k0)-a16(i,j,k0) )*rhoinv
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
                fn(01) = a01(i,j,k0) - fe(01)
                fn(02) = a02(i,j,k0) - fe(02)
                fn(03) = a03(i,j,k0) - fe(03)
                fn(04) = a04(i,j,k0) - fe(04)
                fn(05) = a05(i,j,k0) - fe(05)
                fn(06) = a06(i,j,k0) - fe(06)
                fn(07) = a07(i,j,k0) - fe(07)
                fn(08) = a08(i,j,k0) - fe(08)
                fn(09) = a09(i,j,k0) - fe(09)
                fn(10) = a10(i,j,k0) - fe(10)
                fn(11) = a11(i,j,k0) - fe(11)
                fn(12) = a12(i,j,k0) - fe(12)
                fn(13) = a13(i,j,k0) - fe(13)
                fn(14) = a14(i,j,k0) - fe(14)
                fn(15) = a15(i,j,k0) - fe(15)
                fn(16) = a16(i,j,k0) - fe(16)
                fn(17) = a17(i,j,k0) - fe(17)
                fn(18) = a18(i,j,k0) - fe(18)
                fn(19) = a19(i,j,k0) - fe(19)
!
! 4) compute sigma
!
                sigma_xx = -rho/rf
                do ll = 1,18
                   sigma_xx = sigma_xx - (1-0.5*omega)*cx(ll)*cx(ll)*fn(ll)
                enddo
!
                sigma_xy = 0.0       
                do ll = 1,18
                   sigma_xy = sigma_xy - (1-0.5*omega)*cx(ll)*cy(ll)*fn(ll)
                enddo
!
                sigma_xz = 0.0       
                do ll = 1,18
                   sigma_xz = sigma_xz - (1-0.5*omega)*cx(ll)*cz(ll)*fn(ll)
                enddo
!
                sigma_yy = -rho/rf 
                do ll = 1,18
                   sigma_yy = sigma_yy - (1-0.5*omega)*cy(ll)*cy(ll)*fn(ll)
                enddo
!
                sigma_yz = 0.0      
                do ll = 1,18
                   sigma_yz = sigma_yz - (1-0.5*omega)*cy(ll)*cz(ll)*fn(ll)
                enddo
!
                sigma_zz = -rho/rf 
                do ll = 1,18
                   sigma_zz = sigma_zz - (1-0.5*omega)*cz(ll)*cz(ll)*fn(ll)
                enddo
!
                flux(ii) = flux(ii)+(sigma_zz + sigma_yz + sigma_xz)
             enddo
          enddo
       enddo
!
       fluxZ = -flux(1)+flux(2)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. flux_Z", flux(1),flux(2)
        endif
#endif
!
      return
      end subroutine flux_Z
