!====================================================
!     ****** LBE/coll
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       coll
!     DESCRIPTION
!       Collision according to bgk style (\omega(f_f^(eq)))
!       forcing is included and is proportional to u_0
!     INPUTS
!       itime --> timestep 
!     OUTPUT
!       none
!     TODO
!
!     NOTES
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
        subroutine col(itime)
!
        use storage
        use real_kinds
#ifdef OPENACCOLD
        use openacc
#endif
!
#ifdef _OPENMP
        use omp_lib
#endif
!
        implicit none
!
        integer :: i,j,k,itime
!
#ifdef _OPENMP
        integer:: nthreads, threadid
#endif
!
!
        real(mykind) :: dx01,dx02,dx03,dx04,dx05,dx06,dx07
        real(mykind) :: dx08,dx09,dx10,dx11,dx12,dx13
        real(mykind) :: dx14,dx15,dx16,dx17,dx18,dx19
        real(mykind) :: de01,de02,de03,de04,de05,de06,de07
        real(mykind) :: de08,de09,de10,de11,de12,de13
        real(mykind) :: de14,de15,de16,de17,de18,de19
        real(mykind) :: drho,drhoinv,dvx,dvy,dvz
        real(mykind) :: dvx2,dvy2,dvz2,dvsq
        real(mykind) :: drp1,drp2
        real(mykind) :: dvxpz,dvxmz
        real(mykind) :: dvxpy,dvxmy
        real(mykind) :: dvypz,dvymz
        real(mykind) :: dqx,dqy,dqz
        real(mykind) :: dqxpz,dqxmz,dqxpy,dqxmy,dqypz,dqymz
        real(mykind) :: forcex, forcey, forcez
        real(mykind) :: omega1
!
        omega1 = uno - omega

!	write(6,*) "start sub. coll, ", itime
!
! main loop
!$OMP PARALLEL DO &
!!$OMP DEFAULT(NONE) & 
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(dx01,dx02,dx03,dx04,dx05,dx06,dx07,dx08,dx09,dx10)  &
!$OMP PRIVATE(dx11,dx12,dx13,dx14,dx15,dx16,dx17,dx18,dx19)  &
!$OMP PRIVATE(de01,de02,de03,de04,de05,de06,de07,de08,de09,de10)  &
!$OMP PRIVATE(de11,de12,de13,de14,de15,de16,de17,de18,de19)  &
!$OMP PRIVATE(drho,drhoinv,dvx,dvy,dvz,dvx2,dvy2,dvz2,dvsq)  &
!$OMP PRIVATE(dvxpz,dvxmz,dvxpy,dvxmy,dvypz,dvymz,drp1,drp2)  &
!$OMP PRIVATE(dqxpz,dqxmz,dqxpy,dqxmy,dqypz,dqymz,dqx,dqy,dqz)  &
!$OMP PRIVATE(forcex,forcey,forcez)  &
!$OMP PRIVATE(threadid)  &
!$OMP FIRSTPRIVATE(myrank)  
!!$OMP SHARED(omega,omega1,fgrad,l,m,n,itime)  &
!!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10)  &
!!$OMP SHARED(a11,a12,a13,a14,a15,a16,a17,a18,a19)

!          threadid = OMP_GET_THREAD_NUM()

!!$OMP DO
!$acc kernels 
!$acc loop independent
        do k = 1,n 
!           write(6,*) "I'm", myrank, threadid, k, itime
!$acc loop independent
           do j = 1,m
!!DIR$ SIMD 
!DIR$ IVDEP
!$acc loop independent
           do i = 1,l
        
!           if (obs(i,j,k) == 0) then

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

        drhoinv = uno/drho

        dvx = (dx01+dx02+dx03+dx04+dx05-dx10-dx11-dx12-dx13-dx14)*drhoinv
        dvy = (dx03+dx07+dx08+dx09+dx12-dx01-dx10-dx16-dx17-dx18)*drhoinv
        dvz = (dx04+dx06+dx07+dx13+dx18-dx02-dx09-dx11-dx15-dx16)*drhoinv

        dvx2 = dvx*dvx
        dvy2 = dvy*dvy
        dvz2 = dvz*dvz

        dvsq = dvx2+dvy2+dvz2

        dvxpy = dvx+dvy
        dqxpy = uno+qf*(tre*dvxpy*dvxpy-dvsq)
        dvxmy = dvx-dvy
        dqxmy = uno+qf*(tre*dvxmy*dvxmy-dvsq)
        dvxpz = dvx+dvz
        dqxpz = uno+qf*(tre*dvxpz*dvxpz-dvsq)
        dvxmz = dvx-dvz
        dqxmz = uno+qf*(tre*dvxmz*dvxmz-dvsq)
        dvypz = dvy+dvz
        dqypz = uno+qf*(tre*dvypz*dvypz-dvsq)
        dvymz = dvy-dvz
        dqymz = uno+qf*(tre*dvymz*dvymz-dvsq)
        dqx   = uno+qf*(tre*dvx2       -dvsq)
        dqy   = uno+qf*(tre*dvy2       -dvsq)
        dqz   = uno+qf*(tre*dvz2       -dvsq)
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
!
#ifdef FORCING_Y
        forcex = zero
        forcey = fgrad*drho
        forcez = zero
#else
# ifdef FORCING_Z
        forcex = zero
        forcex = zero
        forcez = fgrad*drho
# else
        forcex = fgrad*drho
        forcey = zero
        forcez = zero
# endif
#endif
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
        de19 = drho*p0*(uno-qf*dvsq)
!
! loop on populations
        a01(i,j,k) = omega1*dx01+(omega*de01 + forcex - forcey         )
        a02(i,j,k) = omega1*dx02+(omega*de02 + forcex          - forcez)
        a03(i,j,k) = omega1*dx03+(omega*de03 + forcex + forcey         )
        a04(i,j,k) = omega1*dx04+(omega*de04 + forcex          + forcez)
        a05(i,j,k) = omega1*dx05+(omega*de05 + forcex                  )
        a06(i,j,k) = omega1*dx06+(omega*de06                   + forcez)
        a07(i,j,k) = omega1*dx07+(omega*de07          + forcey + forcez)
        a08(i,j,k) = omega1*dx08+(omega*de08          + forcey         )
        a09(i,j,k) = omega1*dx09+(omega*de09          + forcey - forcez)
        a10(i,j,k) = omega1*dx10+(omega*de10 - forcex - forcey         )
        a11(i,j,k) = omega1*dx11+(omega*de11 - forcex          - forcez)
        a12(i,j,k) = omega1*dx12+(omega*de12 - forcex + forcey         )
        a13(i,j,k) = omega1*dx13+(omega*de13 - forcex          + forcez)
        a14(i,j,k) = omega1*dx14+(omega*de14 - forcex                  )
        a15(i,j,k) = omega1*dx15+(omega*de15                   - forcez)
        a16(i,j,k) = omega1*dx16+(omega*de16          - forcey - forcez)
        a17(i,j,k) = omega1*dx17+(omega*de17          - forcey         )
        a18(i,j,k) = omega1*dx18+(omega*de18          - forcey + forcez)
        a19(i,j,k) = omega1*dx19+(omega*de19                           )
!
!              endif
!
              end do
           end do
        end do
!$acc end kernels
!!$acc end data region
!!$OMP END DO
!$OMP END PARALLEL DO
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. coll"
        endif
#endif
        return
        end subroutine col
