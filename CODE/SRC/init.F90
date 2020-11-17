!=======================================================================
!     ****** LBE/init
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       init
!     DESCRIPTION
!       initial condition for 2d bgk simulation
!       you can choose: poiseuille flow	(OK)
!                       decayinig flow	(OK)
!                       couette flow	(OK)
!                       taylor vortices	(OK)
!                       rest flow	(OK)
!                       double periodic shear flow 
!                       (see Dellar Phys. Rev. E. 64 (2001))
!       write on unit 76 (fort.76) x & z coordinates for check
!     INPUTS
!       opt ---> 0	rest flow
!           ---> 1	poiseuille flow
!           ---> 2	decayng flow
!           ---> 3	Kida vortices
!           ---> 4	couette flow
!           ---> 5	double periodic shear flow
!           ---> other	stop simluation
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,k,opt
!       real variables used: x,z,xj,zj,vsq,rho,pi
!                            x02,x04,x05,x06,x11,x13,x14,x15,x19
!                            kappa,delta
!       u = velocity along x direction (i, streamwise)
!       v = velocity along z direction (k, normal-to-wall)
!
!     *****
!=======================================================================
!
        subroutine init(opt)
!
        use storage
        implicit none
!
        integer i,j,k,opt,ierr
!
! just a check
!        real(sp) ::  x,y,z,xj,yj,zj
        real(mykind) ::  x,y,z,xj,yj,zj
        real(mykind) ::  cvsq,crho,pi,rho
        real(mykind) ::  cx01,cx02,cx03,cx04,cx05,cx06
        real(mykind) ::  cx07,cx08,cx09,cx10,cx11,cx12
        real(mykind) ::  cx13,cx14,cx15,cx16,cx17,cx18
        real(mykind) ::  cx19
        real(mykind) ::  kappa, delta
        real(mykind) ::  zstart, ystart, xstart
!
        integer      :: nn, mm, ll
!
        parameter(pi=3.141592653589793238462643383279)
        parameter(kappa=80)
        parameter(delta=0.05)
!
! check parameter opt
        if((opt.lt.0).or.(opt.gt.5)) then
           write(6,*) "Initial condition out of range[0,5]",opt
           stop
        endif
!!        write(6,*) 'building starting configuration', opt
!
#ifdef SERIAL
        ll = l
        mm = m
        nn = n
        xstart = 0
        ystart = 0
        zstart = 0
#else
        ll = lx
        mm = ly
        nn = lz
        xstart = mpicoords(1)*l
        ystart = mpicoords(2)*m
        zstart = mpicoords(3)*n
#endif
!
! builds the populations
        crho = 1.0_mykind
!
        do k = 1, n
           z = (real(k+zstart) - 0.5 - 0.5*real(nn))/(0.5*real(nn))
!           write(76,*) k, z
        enddo
!	write(76,*) " "
        do j = 1, m
           y = (real(j+ystart) - 0.5 - 0.5*real(mm))/(0.5*real(mm))
!           write(76,*) j, y
        enddo
!	write(76,*) " "
        do i = 1, l
           x = (real(i+xstart) - 0.5 - 0.5*real(ll))/(0.5*real(nn))
!           write(76,*) i, x
        enddo
!
! rest flow 
        xj = 0.0
        yj = 0.0
        zj = 0.0
!
!
! With PGI compiler it creates NaN using OpenMP????? (BUGS)
!
!!$OMP PARALLEL DEFAULT(NONE)                                      & 
!!$OMP PRIVATE(i,j,k)                                              &
!!$OMP PRIVATE(xj,yj,zj)                                           &
!!$OMP PRIVATE(cx01,cx02,cx03,cx04,cx05,cx06,cx07,cx08,cx09,cx10)  &
!!$OMP PRIVATE(cx11,cx12,cx13,cx14,cx15,cx16,cx17,cx18,cx19)       &
!!$OMP PRIVATE(cvsq,crho)                                          & 
!!$OMP PRIVATE(x,y,z)                                              &
!!$OMP SHARED(opt,u0,u00)                                          &
!!$OMP SHARED(l1,m1,n1,ll,mm,nn)                                   &
!!$OMP SHARED(xstart,ystart,zstart)                                &
!!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09)                 &
!!$OMP SHARED(a10,a11,a12,a13,a14,a15,a16,a17,a18,a19)             
!!$OMP DO
        do k = 0, n1
#ifdef PERIODIC
           z = (real(k+zstart,mykind)-0.5_qp)/real(nn,mykind)     ! 0<z<1 (taylor)
#else
           z = (real(k+zstart)-0.5-0.5*real(nn))/(0.5*real(nn))! -1<z<1
#endif
           do j = 0, m1
#ifdef PERIODIC
              y = (real(j+ystart,mykind)-0.5_qp)/real(mm,mykind)  ! 0<x<1 (taylor)
#else
              y = (real(j+ystart)-0.5-0.5*real(mm))/(0.5*real(mm)) ! -1<x<1
#endif
              do i = 0, l1
#ifdef PERIODIC
                 x = (real(i+xstart,mykind)-0.5_qp)/real(ll,mykind)! 0<x<1 (taylor)
#else
                 x = (real(i+xstart)-0.5-0.5*real(ll))/(0.5*real(nn)) ! -1<x<1
#endif

!
                 if(opt.eq.1) then
                   xj = 0.1 * ( 1 - z*z)                    !poiseuille
                   yj = 0.0
                   zj = 0.0
                 endif
!
                 if(opt.eq.2) then
#ifdef HALF_P
#else
                   xj = u0 * sin(2*pi*z)                    !decaying flow
#endif
                   yj = 0.0
                   zj = 0.0
                 endif
!
                 if(opt.eq.3) then
#ifdef HALF_P
#else
                   xj = 0.1*sin(real(2,mykind)*pi*x)*cos(real(2,mykind)*pi*z)          !kida(?) vortices
                   yj = 0.0
                   zj =-0.1*cos(real(2,mykind)*pi*x)*sin(real(2,mykind)*pi*z)
#endif
                 endif
!
! take care of pressure..... (incompressible???)
!
                 if(opt.eq.4) then
                    xj = u00*(1+z)/2                         !couette flow
                    zj = 0.0
                    write(6,*) "INIT-2: still to port!!!!"
                    stop
                 endif
!
                 if(opt.eq.5) then
                    if(z.le.0.5) then                          !double shear
                       xj = 0.2*tanh(kappa*(z-0.25))
                    else
                       xj = 0.2*tanh(kappa*(0.75-z))
                    endif
                    zj = 0.1*delta*sin(2*pi*(x+0.25))
                   write(6,*) "INIT-3: still to port!!!!"
                   stop
                 endif
!
                 cvsq=xj*xj+yj*yj+zj*zj

                 cx01 = rf*( xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx02 = rf*( xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx03 = rf*( xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx04 = rf*( xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx05 = rf*( xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx06 = rf*(       zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx07 = rf*(    yj+zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx08 = rf*(    yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx09 = rf*(    yj-zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx10 = rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx11 = rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx12 = rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx13 = rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx14 = rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx15 = rf*(      -zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx16 = rf*(   -yj-zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx17 = rf*(   -yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx18 = rf*(   -yj+zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx19 = rf*(   0.0   )+qf*(3.0*( 0.0 )*( 0.0 )-cvsq)

                 a01(i,j,k) = crho*p2*(1.0+cx01)
                 a02(i,j,k) = crho*p2*(1.0+cx02)
                 a03(i,j,k) = crho*p2*(1.0+cx03)
                 a04(i,j,k) = crho*p2*(1.0+cx04)
                 a05(i,j,k) = crho*p1*(1.0+cx05)
                 a06(i,j,k) = crho*p1*(1.0+cx06)
                 a07(i,j,k) = crho*p2*(1.0+cx07)
                 a08(i,j,k) = crho*p1*(1.0+cx08)
                 a09(i,j,k) = crho*p2*(1.0+cx09)
                 a10(i,j,k) = crho*p2*(1.0+cx10)
                 a11(i,j,k) = crho*p2*(1.0+cx11)
                 a12(i,j,k) = crho*p2*(1.0+cx12)
                 a13(i,j,k) = crho*p2*(1.0+cx13)
                 a14(i,j,k) = crho*p1*(1.0+cx14)
                 a15(i,j,k) = crho*p1*(1.0+cx15)
                 a16(i,j,k) = crho*p2*(1.0+cx16)
                 a17(i,j,k) = crho*p1*(1.0+cx17)
                 a18(i,j,k) = crho*p2*(1.0+cx18)
                 a19(i,j,k) = crho*p0*(1.0+cx19)

              end do
           end do
        end do
!!$OMP END DO
!!$OMP END PARALLEL
!
        if(myrank == 0) then 
           if (opt == 0) then
              write(16,*) "INFO: initial condition --> rest flow"
           endif
!
           if (opt == 1) then
              write(16,*) "INFO: initial condition --> poiseuille flow"
           endif
!
           if (opt == 2) then
              u0 = 0
              write(16,*) "INFO: initial condition --> decayinig flow"
              write(16,*) "INFO: initial condition --> setting u0 = 0"
           endif
!
           if (opt == 3) then
              u0 = 0
              write(16,*) "INFO: initial condition --> taylor vortices"
              write(16,*) "INFO: initial condition --> setting u0 = 0"
              write(6,*) "WARNING: serial version doesn't work!!!"
              if(lx.NE.ly)  then
                 write(6,*) "ERROR:  the box is not a square!!!!", lx, lz
                 write(16,*) "ERROR:  the box is not a square!!!!", lx, lz
                 stop
              endif
           endif
!
           if (opt == 4) then
              write(16,*) "INFO: initial condition --> couette flow"
           endif
!
           if (opt == 5) then
              write(16,*) "INFO: initial condition --> double periodic shear flow"
           endif
        endif
!
#ifdef SERIAL
! do nothing
#else
        call mpi_barrier(lbecomm,ierr)
#endif
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. init"
        endif
#endif
!
        return
        end subroutine init
