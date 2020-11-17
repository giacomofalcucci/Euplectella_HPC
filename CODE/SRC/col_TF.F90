!====================================================
!     ****** LBE/coll_TF
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       coll
!     DESCRIPTION
!       Collision according to bgk style (\omega(f_f^(eq)))
!       forcing is included and is proportional to u_0
!       TWO Field approch, input is different from output...
!     INPUTS
!       itime --> timestep 
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       experimental
!       integer variables used: i,k,itime
!       real variables used: x02,x04,x05,x06,x11,x13,x14,x15,x19
!                            e02,e04,e05,e06,e11,e13,e14,e15,e19
!                            rho,rhoinv,vx,vz,vx2,vz2,vsq
!                            vxpz,vxmz,rp1,rp2
!                            qxpz,qxmz,qx,qz
!                            force,pi
!
!
!     *****
!====================================================
!
        subroutine col_TF(itime)
!
        use storage
        use real_kinds
#ifdef OPENACCOLD
        use openacc
#endif
!
        implicit none
!
        integer i,j,k,itime
!
        real(mykind) dx01,dx02,dx03,dx04,dx05,dx06,dx07
        real(mykind) dx08,dx09,dx10,dx11,dx12,dx13
        real(mykind) dx14,dx15,dx16,dx17,dx18,dx19
        real(mykind) de01,de02,de03,de04,de05,de06,de07
        real(mykind) de08,de09,de10,de11,de12,de13
        real(mykind) de14,de15,de16,de17,de18,de19
        real(mykind) drho,drhoinv,dvx,dvy,dvz
        real(mykind) dvx2,dvy2,dvz2,dvsq
!	real(sp) drho,drhoinv,dvx,dvz,dvx2,dvz2,dvsq	! for check
        real(mykind) drp1,drp2
        real(mykind) dvxpz,dvxmz
        real(mykind) dvxpy,dvxmy
        real(mykind) dvypz,dvymz
        real(mykind) dqx,dqy,dqz
        real(mykind) dqxpz,dqxmz,dqxpy,dqxmy,dqypz,dqymz
        real(mykind) force
        real(mykind) omega1
!
        omega1 = 1 - omega

!	write(6,*) "start sub. coll, ", itime
!
! main loop
!$OMP PARALLEL DEFAULT(NONE)                                      & 
!$OMP PRIVATE(i,j,k)                                              &
!$OMP PRIVATE(dx01,dx02,dx03,dx04,dx05,dx06,dx07,dx08,dx09,dx10)  &
!$OMP PRIVATE(dx11,dx12,dx13,dx14,dx15,dx16,dx17,dx18,dx19)       &
!$OMP PRIVATE(de01,de02,de03,de04,de05,de06,de07,de08,de09,de10)  &
!$OMP PRIVATE(de11,de12,de13,de14,de15,de16,de17,de18,de19)       &
!$OMP PRIVATE(drho,drhoinv,dvx,dvy,dvz,dvx2,dvy2,dvz2,dvsq)       &
!$OMP PRIVATE(dvxpz,dvxmz,dvxpy,dvxmy,dvypz,dvymz,drp1,drp2)      &
!$OMP PRIVATE(dqxpz,dqxmz,dqxpy,dqxmy,dqypz,dqymz,dqx,dqy,dqz)    &
!$OMP PRIVATE(force)                                              &
!$OMP SHARED(omega,omega1,fgrad,l,m,n)                            &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10)             &
!$OMP SHARED(a11,a12,a13,a14,a15,a16,a17,a18,a19)                 & 
!$OMP SHARED(b01,b02,b03,b04,b05,b06,b07,b08,b09,b10)             &
!$OMP SHARED(b11,b12,b13,b14,b15,b16,b17,b18,b19)
!$OMP DO
!$acc data region copy(a02(0:l1,0:n1),a04(0:l1,0:n1),a05(0:l1,0:n1),a06(0:l1,0:n1),a11(0:l1,0:n1),a13(0:l1,0:n1),a14(0:l1,0:n1),a15(0:l1,0:n1),a19(0:l1,0:n1))
!$acc kernels loop private(dx02,dx04,dx05,dx06,dx11,dx13,dx14,dx15,dx19,de02,de04,de05,de06,de11,de13,de14,de15,de19,drho,drhoinv,dvx,dvz,dvx2,dvz2,dvsq, dvxpz,dvxmz,drp1,drp2,dqxpz,dqxmz,dqx,dqz,force)
        do k = 1,n 
           do j = 1,m
!DIR$ SIMD 
!DIR$ IVDEP
           do i = 1,l

        dx01 = b01(i,j,k)
        dx02 = b02(i,j,k)
        dx03 = b03(i,j,k)
        dx04 = b04(i,j,k)
        dx05 = b05(i,j,k)
        dx06 = b06(i,j,k)
        dx07 = b07(i,j,k)
        dx08 = b08(i,j,k)
        dx09 = b09(i,j,k)
        dx10 = b10(i,j,k)
        dx11 = b11(i,j,k)
        dx12 = b12(i,j,k)
        dx13 = b13(i,j,k)
        dx14 = b14(i,j,k)
        dx15 = b15(i,j,k)
        dx16 = b16(i,j,k)
        dx17 = b17(i,j,k)
        dx18 = b18(i,j,k)
        dx19 = a19(i,j,k)

        drho = dx01 + dx02 + dx03 + dx04 + dx05 + dx06 + dx07 + dx08 &
              +dx09 + dx10 + dx11 + dx12 + dx13 + dx14 + dx15 + dx16 &
              +dx17 + dx18 + dx19

        drhoinv = 1.0/drho

        dvx = (dx01+dx02+dx03+dx04+dx05-dx10-dx11-dx12-dx13-dx14)*drhoinv
        dvy = (dx03+dx07+dx08+dx09+dx12-dx01-dx10-dx16-dx17-dx18)*drhoinv
        dvz = (dx04+dx06+dx07+dx13+dx18-dx02-dx09-dx11-dx15-dx16)*drhoinv

        dvx2 = dvx*dvx
        dvy2 = dvy*dvy
        dvz2 = dvz*dvz

        dvsq = dvx2+dvy2+dvz2

        dvxpy = dvx+dvy
        dqxpy = 1.0+qf*(3.0*dvxpy*dvxpy-dvsq)
        dvxmy = dvx-dvy
        dqxmy = 1.0+qf*(3.0*dvxmy*dvxmy-dvsq)
        dvxpz = dvx+dvz
        dqxpz = 1.0+qf*(3.0*dvxpz*dvxpz-dvsq)
        dvxmz = dvx-dvz
        dqxmz = 1.0+qf*(3.0*dvxmz*dvxmz-dvsq)
        dvypz = dvy+dvz
        dqypz = 1.0+qf*(3.0*dvypz*dvypz-dvsq)
        dvymz = dvy-dvz
        dqymz = 1.0+qf*(3.0*dvymz*dvymz-dvsq)
        dqx   = 1.0+qf*(3.0*dvx2       -dvsq)
        dqy   = 1.0+qf*(3.0*dvy2       -dvsq)
        dqz   = 1.0+qf*(3.0*dvz2       -dvsq)
        dvx   = rf*dvx
        dvy   = rf*dvy
        dvz   = rf*dvz
        dvxpy = rf*dvxpy
        dvxmy = rf*dvxmy
        dvxpz = rf*dvxpz
        dvxmz = rf*dvxmz
        dvypz = rf*dvypz
        dvymz = rf*dvymz

        drp1 = drho*p1
        drp2 = drho*p2
!
! forcing term
        force = fgrad*drho
!
! equilibrium distribution
        de01 = drp2*(+dvxmy+dqxmy)
        de02 = drp2*(+dvxmz+dqxmz)
        de03 = drp2*(+dvxpy+dqxpy)
        de04 = drp2*(+dvxpz+dqxpz)
        de05 = drp1*(+dvx  +dqx  )
        de06 = drp1*(+dvz  +dqz  )
        de07 = drp2*(+dvypz+dqypz)
        de08 = drp1*(+dvy  +dqy  )
        de09 = drp2*(+dvymz+dqymz)
        de10 = drp2*(-dvxpy+dqxpy)
        de11 = drp2*(-dvxpz+dqxpz)
        de12 = drp2*(-dvxmy+dqxmy)
        de13 = drp2*(-dvxmz+dqxmz)
        de14 = drp1*(-dvx  +dqx  )
        de15 = drp1*(-dvz  +dqz  )
        de16 = drp2*(-dvypz+dqypz)
        de17 = drp1*(-dvy  +dqy  )
        de18 = drp2*(-dvymz+dqymz)
        de19 = drho*p0*(1.0-qf*dvsq)
!
! loop on populations
        a01(i,j,k) = omega1*dx01+(omega*de01 + force)
        a02(i,j,k) = omega1*dx02+(omega*de02 + force)
        a03(i,j,k) = omega1*dx03+(omega*de03 + force)
        a04(i,j,k) = omega1*dx04+(omega*de04 + force)
        a05(i,j,k) = omega1*dx05+(omega*de05 + force)
        a06(i,j,k) = omega1*dx06+(omega*de06        )
        a07(i,j,k) = omega1*dx07+(omega*de07        )
        a08(i,j,k) = omega1*dx08+(omega*de08        )
        a09(i,j,k) = omega1*dx09+(omega*de09        )
        a10(i,j,k) = omega1*dx10+(omega*de10 - force)
        a11(i,j,k) = omega1*dx11+(omega*de11 - force)
        a12(i,j,k) = omega1*dx12+(omega*de12 - force)
        a13(i,j,k) = omega1*dx13+(omega*de13 - force)
        a14(i,j,k) = omega1*dx14+(omega*de14 - force)
        a15(i,j,k) = omega1*dx15+(omega*de15        )
        a16(i,j,k) = omega1*dx16+(omega*de16        )
        a17(i,j,k) = omega1*dx17+(omega*de17        )
        a18(i,j,k) = omega1*dx18+(omega*de18        )
        a19(i,j,k) = omega1*dx19+(omega*de19        )

              end do
           end do
        end do
!$acc end kernels
!$acc end data region
!$OMP END DO
!$OMP END PARALLEL

#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. coll_TF"
        endif
#endif

        return
        end subroutine col_TF

