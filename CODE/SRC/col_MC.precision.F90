!====================================================
!     ****** LBE/coll_MC
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       coll
!     DESCRIPTION
!       Collision according to bgk style (\omega(f_f^(eq)))
!       forcing is included and is proportional to u_0
!       MOVE + COLLIDE version
!     INPUTS
!       itime --> timestep 
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
!                            forcex,pi
!
!
!     *****
!====================================================
!
        subroutine col_MC(itime)
!
        use storage
        use real_kinds
        use timing
!
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
        integer:: i,j,k,itime
!
#ifdef MIXEDPRECISION
# ifdef DOUBLE_P
        real(qp) :: x01,x02,x03,x04,x05,x06,x07,x08,x09,x10
        real(qp) :: x11,x12,x13,x14,x15,x16,x17,x18,x19
        real(qp) :: e01,e02,e03,e04,e05,e06,e07,e08,e09,e10
        real(qp) :: e11,e12,e13,e14,e15,e16,e17,e18,e19
        real(qp) :: rho,rhoinv,vx,vy,vz,vx2,vy2,vz2,vsq
        real(qp) :: vxpy,vxmy,vxpz,vxmz,vypz,vymz,rp1,rp2
        real(qp) :: qxpy,qxmy,qxpz,qxmz,qypz,qymz,qx,qy,qz
        real(qp) :: forcex, forcey, forcez
        real(qp) :: pi
# else
        real(dp) :: x01,x02,x03,x04,x05,x06,x07,x08,x09,x10
        real(dp) :: x11,x12,x13,x14,x15,x16,x17,x18,x19
        real(dp) :: e01,e02,e03,e04,e05,e06,e07,e08,e09,e10
        real(dp) :: e11,e12,e13,e14,e15,e16,e17,e18,e19
        real(dp) :: rho,rhoinv,vx,vy,vz,vx2,vy2,vz2,vsq
        real(dp) :: vxpy,vxmy,vxpz,vxmz,vypz,vymz,rp1,rp2
        real(dp) :: qxpy,qxmy,qxpz,qxmz,qypz,qymz,qx,qy,qz
        real(dp) :: forcex, forcey, forcez
        real(dp) :: pi
# endif
#else
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07,x08,x09,x10
        real(mykind) :: x11,x12,x13,x14,x15,x16,x17,x18,x19
        real(mykind) :: e01,e02,e03,e04,e05,e06,e07,e08,e09,e10
        real(mykind) :: e11,e12,e13,e14,e15,e16,e17,e18,e19
        real(mykind) :: rho,rhoinv,vx,vy,vz,vx2,vy2,vz2,vsq
        real(mykind) :: vxpy,vxmy,vxpz,vxmz,vypz,vymz,rp1,rp2
        real(mykind) :: qxpy,qxmy,qxpz,qxmz,qypz,qymz,qx,qy,qz
        real(mykind) :: forcex, forcey, forcez
        real(mykind) :: pi
#endif
!
        parameter(pi=3.14159265358979)
!
#ifdef ORIGINAL
! skip this subroutine
#else
!
# ifdef DEBUG_3
        real(mykind) :: cte
        character*17 file_nameD
        file_nameD = 'debug.xxx.xxx.log'
        write(file_nameD(7:9),3300) itime
        write(file_nameD(11:13),3300) myrank
        open(41,file=file_nameD, status='unknown')        ! debug file
!
        call probe_global(itime,(3*lx/4),(ly/2),(lz/2-5))
# endif
!
!
!        call time(tcountG0)
!
!        do k = 1,n
! works (fine) with intel compiler!!!
! works (quite well) with IBM compiler !!!
!!IBM* ASSERT (NODEPS)
! works (?) with pgf compiler !!!
!!!!$acc do vector(256)
!!!!pgi$l nodepchk
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(x01,x02,x03,x04,x05,x06,x07,x08,x09,x10)  &
!$OMP PRIVATE(x11,x12,x13,x14,x15,x16,x17,x18,x19)  &
!$OMP PRIVATE(e01,e02,e03,e04,e05,e06,e07,e08,e09,e10)  &
!$OMP PRIVATE(e11,e12,e13,e14,e15,e16,e17,e18,e19)  &
!$OMP PRIVATE(rho,rhoinv,vx,vy,vz,vx2,vy2,vz2,vsq)  &
!$OMP PRIVATE(vxpz,vxmz,vxpy,vxmy,vypz,vymz,rp1,rp2)  &
!$OMP PRIVATE(qxpz,qxmz,qxpy,qxmy,qypz,qymz,qx,qy,qz)  &
!$OMP PRIVATE(forcex,forcey,forcez)  &
!$OMP SHARED(omega,fgrad,l,m,n,itime)  &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10)  &
!$OMP SHARED(a11,a12,a13,a14,a15,a16,a17,a18,a19)      &
!$OMP SHARED(b01,b02,b03,b04,b05,b06,b07,b08,b09,b10)  &
!$OMP SHARED(b11,b12,b13,b14,b15,b16,b17,b18,b19)  
!$OMP DO
!$acc kernels  
!$acc loop independent
        do k = 1,n
!$acc loop independent 
        do j = 1,m
!DIR$ IVDEP
!$acc loop independent 
         do i = 1,l
!
        x01 = a01(i-1,j+1,k  )
        x02 = a02(i-1,j  ,k+1)
        x03 = a03(i-1,j-1,k  )
        x04 = a04(i-1,j  ,k-1)
        x05 = a05(i-1,j  ,k  )
        x06 = a06(i  ,j  ,k-1)
        x07 = a07(i  ,j-1,k-1)
        x08 = a08(i  ,j-1,k  )
        x09 = a09(i  ,j-1,k+1)
        x10 = a10(i+1,j+1,k  )
        x11 = a11(i+1,j  ,k+1)
        x12 = a12(i+1,j-1,k  )
        x13 = a13(i+1,j  ,k-1)
        x14 = a14(i+1,j  ,k  )
        x15 = a15(i  ,j  ,k+1)
        x16 = a16(i  ,j+1,k+1)
        x17 = a17(i  ,j+1,k  )
        x18 = a18(i  ,j+1,k-1)
        x19 = a19(i  ,j  ,k  )
!
        rho = x01+x02+x03+x04+x05+x06+x07+x08+ &
     &        x09+x10+x11+x12+x13+x14+x15+x16+ &
     &        x17+x18+x19
!
        rhoinv = uno/rho
!
        vx = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
        vy = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
        vz = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv
!
# ifdef DEBUG_3
        write(41,3131) i+offset(1),j+offset(2),k+offset(3), obs(i,j,k),   & 
                        vx,vy,vz,rho
!                       x01,x02,x03,x04,x05,x06
!                       x07,x08,x09,x10,x11,x12
!                       x13,x14,x15,x16,x17,x18
# endif
!
        vx2 = vx*vx
        vy2 = vy*vy
        vz2 = vz*vz
!
        vsq = vx2+vy2+vz2
!
        vxpy = vx+vy
        qxpy = +qf*(tre*vxpy*vxpy-vsq)
        vxmy = vx-vy
        qxmy = +qf*(tre*vxmy*vxmy-vsq)
        vxpz = vx+vz
        qxpz = +qf*(tre*vxpz*vxpz-vsq)
        vxmz = vx-vz
        qxmz = +qf*(tre*vxmz*vxmz-vsq)
        vypz = vy+vz
        qypz = +qf*(tre*vypz*vypz-vsq)
        vymz = vy-vz
        qymz = +qf*(tre*vymz*vymz-vsq)
        qx   = +qf*(tre*vx2-vsq)
        qy   = +qf*(tre*vy2-vsq)
        qz   = +qf*(tre*vz2-vsq)
        vx   = rf*vx
        vy   = rf*vy
        vz   = rf*vz
        vxpy = rf*vxpy
        vxmy = rf*vxmy
        vxpz = rf*vxpz
        vxmz = rf*vxmz
        vypz = rf*vypz
        vymz = rf*vymz
!
        rp1 = rho*p1
        rp2 = rho*p2
!
! equilibrium distribution
        e01 = rp2*   (uno+( vxmy+(qxmy+fgrad/(omega*p2))))
        e02 = rp2*   (uno+( vxmz+(qxmz+fgrad/(omega*p2))))
        e03 = rp2*   (uno+( vxpy+(qxpy+fgrad/(omega*p2))))
        e04 = rp2*   (uno+( vxpz+(qxpz+fgrad/(omega*p2))))
        e05 = rp1*   (uno+( vx  +(qx  +fgrad/(omega*p1))))
        e06 = rp1*   (uno+( vz  +qz  ))
        e07 = rp2*   (uno+( vypz+qypz))
        e08 = rp1*   (uno+( vy  +qy  ))
        e09 = rp2*   (uno+( vymz+qymz))
        e10 = rp2*   (uno+(-vxpy+(qxpy-fgrad/(omega*p2))))
        e11 = rp2*   (uno+(-vxpz+(qxpz-fgrad/(omega*p2))))
        e12 = rp2*   (uno+(-vxmy+(qxmy-fgrad/(omega*p2))))
        e13 = rp2*   (uno+(-vxmz+(qxmz-fgrad/(omega*p2))))
        e14 = rp1*   (uno+(-vx  +(qx  -fgrad/(omega*p1))))
        e15 = rp1*   (uno+(-vz  +qz  ))
        e16 = rp2*   (uno+(-vypz+qypz))
        e17 = rp1*   (uno+(-vy  +qy  ))
        e18 = rp2*   (uno+(-vymz+qymz))
        e19 = rho*p0*(uno-qf*vsq)
!
! forcing term
!
# ifdef FORCING_Y
        forcex = zero
        forcey = fgrad*rho
        forcez = zero
# else
#  ifdef FORCING_Z
        forcex = zero
        forcex = zero
        forcez = fgrad*rho
#  else
        forcex = fgrad*rho     ! default value...
        forcey = zero
        forcez = zero
#  endif
# endif
!
! loop on populations
        b01(i,j,k) = x01 - omega*(x01-e01) !+ forcex - forcey     
        b02(i,j,k) = x02 - omega*(x02-e02) !+ forcex          - forcez
        b03(i,j,k) = x03 - omega*(x03-e03) !+ forcex + forcey         
        b04(i,j,k) = x04 - omega*(x04-e04) !+ forcex          + forcez
        b05(i,j,k) = x05 - omega*(x05-e05) !+ forcex
        b06(i,j,k) = x06 - omega*(x06-e06) !                  + forcez
        b07(i,j,k) = x07 - omega*(x07-e07) !         + forcey + forcez
        b08(i,j,k) = x08 - omega*(x08-e08) !         + forcey
        b09(i,j,k) = x09 - omega*(x09-e09) !         + forcey - forcez
        b10(i,j,k) = x10 - omega*(x10-e10) !- forcex - forcey         
        b11(i,j,k) = x11 - omega*(x11-e11) !- forcex          - forcez
        b12(i,j,k) = x12 - omega*(x12-e12) !- forcex + forcey         
        b13(i,j,k) = x13 - omega*(x13-e13) !- forcex          + forcez
        b14(i,j,k) = x14 - omega*(x14-e14) !- forcex 
        b15(i,j,k) = x15 - omega*(x15-e15) !                  - forcez
        b16(i,j,k) = x16 - omega*(x16-e16) !         - forcey - forcez
        b17(i,j,k) = x17 - omega*(x17-e17) !         - forcey         
        b18(i,j,k) = x18 - omega*(x18-e18) !         - forcey + forcez
        a19(i,j,k) = x19 - omega*(x19-e19) !
!
        end do
!!!$OMP END SIMD 
        end do
        end do
!$acc end kernels
!
!
!#ifdef _OPENMP
!        write(6,*) "OpenMP is on...", itime
!#endif
!$OMP END PARALLEL

! fix
!        call time(tcountG0)
!
        c01 => a01
        c02 => a02
        c03 => a03
        c04 => a04
        c05 => a05
        c06 => a06
        c07 => a07
        c08 => a08
        c09 => a09
        c10 => a10
        c11 => a11
        c12 => a12
        c13 => a13
        c14 => a14
        c15 => a15
        c16 => a16
        c17 => a17
        c18 => a18
!!        c19 => a19
!
! new ---> current
!
        a01 => b01
        a02 => b02
        a03 => b03
        a04 => b04
        a05 => b05
        a06 => b06
        a07 => b07
        a08 => b08
        a09 => b09
        a10 => b10
        a11 => b11
        a12 => b12
        a13 => b13
        a14 => b14
        a15 => b15
        a16 => b16
        a17 => b17
        a18 => b18
!!        a19 => b19
!
        b01 => c01
        b02 => c02
        b03 => c03
        b04 => c04
        b05 => c05
        b06 => c06
        b07 => c07
        b08 => c08
        b09 => c09
        b10 => c10
        b11 => c11
        b12 => c12
        b13 => c13
        b14 => c14
        b15 => c15
        b16 => c16
        b17 => c17
        b18 => c18
!!        b19 => c19
!
!
!        call time(tcountG1)
!        write(6,*) tcountG1,tcountG0, tcountG1-tcountG0
!        write(6,*) "end sub. coll, ", itime
#endif
! 
#ifdef DEBUG_3
!
! format
3300  format(i3.3)
!3131  format(3(i,1x),6(e14.6,1x))
3131  format(4(i,1x),4(e14.6,1x))
        close(41)
#endif
!
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. coll_MC"
        endif
#endif
!
        return 
        end subroutine col_MC
