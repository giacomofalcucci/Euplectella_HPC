!=====================================================================
!     ****** LBE/diagno
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       diagno
!     DESCRIPTION
!       diagnostic subroutine:
!       check of conserved quantities and mean value
!       write on unit 63 (diagno.dat)
!       write on unit 16 (bgk.log)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,k
!       real variables used: rtot,xtot,ztot,stot
!                            xj, zj, rho, rhoinv, rdv
!
!     *****
!=====================================================================
!
       subroutine diagno(itime)
!
       use storage
       implicit none
!
       integer:: itime,i,j,k
       integer:: ierr
!
#ifdef MIXEDPRECISION
# ifdef DOUBLE_P
       real(qp):: rtot,xtot,ytot, ztot,stot
       real(qp):: xj, yj, zj, rho, rhoinv, rdv
       real(qp):: loctot(5), glotot(5)
# else
       real(dp):: rtot,xtot,ytot, ztot,stot
       real(dp):: xj, yj, zj, rho, rhoinv, rdv
       real(dp):: loctot(5), glotot(5)
# endif
#else
       real(mykind):: rtot,xtot,ytot, ztot,stot
       real(mykind):: xj, yj, zj, rho, rhoinv, rdv
       real(mykind):: loctot(5), glotot(5)
#endif
!
!
       rdv = 1./(real(l,mykind)*real(m,mykind)*real(n,mykind))
!
       glotot(1) = 0.
       glotot(2) = 0.
       glotot(3) = 0.
       glotot(4) = 0.
       glotot(5) = 0.
!
       rtot = 0.
       xtot = 0.
       ytot = 0.
       ztot = 0.
       stot = 0.
!
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i,j,k)  &
!$OMP PRIVATE(rho,xj,yj,zj,rhoinv)  &
!$OMP SHARED(l,m,n)  &
!$OMP SHARED(obs)  &
!$OMP SHARED(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10)  &
!$OMP SHARED(a11,a12,a13,a14,a15,a16,a17,a18,a19)      &
!$OMP REDUCTION(+:rtot,xtot,ytot,ztot,stot) 
!$OMP DO
       do k=1,n
          do j=1,m
             do i=1,l
                if (obs(i,j,k) == 1) then 
                   rho = 1.0
                   xj = 0.0
                   yj = 0.0
                   zj = 0.0
                else
                   rho = +a01(i,j,k)+a02(i,j,k)+a03(i,j,k) &
                         +a04(i,j,k)+a05(i,j,k)+a06(i,j,k) &
                         +a07(i,j,k)+a08(i,j,k)+a09(i,j,k) &
                         +a10(i,j,k)+a11(i,j,k)+a12(i,j,k) &
                         +a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                         +a16(i,j,k)+a17(i,j,k)+a18(i,j,k) &
                         +a19(i,j,k)

                   rhoinv = 1.0/rho

                   xj = +( a01(i,j,k)+a02(i,j,k)+a03(i,j,k) &
                          +a04(i,j,k)+a05(i,j,k)            & 
                          -a10(i,j,k)-a11(i,j,k)-a12(i,j,k) &
                          -a13(i,j,k)-a14(i,j,k) )*rhoinv

                   yj = +( a03(i,j,k)+a07(i,j,k)+a08(i,j,k) &
                          +a09(i,j,k)+a12(i,j,k)            &
                          -a01(i,j,k)-a10(i,j,k)-a16(i,j,k) &
                          -a17(i,j,k)-a18(i,j,k) )*rhoinv

                   zj = +( a04(i,j,k)+a06(i,j,k)+a07(i,j,k) &
                          +a13(i,j,k)+a18(i,j,k)            &
                          -a02(i,j,k)-a09(i,j,k)-a11(i,j,k) &
                          -a15(i,j,k)-a16(i,j,k) )*rhoinv
                endif
!
                rtot = rtot+rho
                xtot = xtot+xj
                ytot = ytot+yj
                ztot = ztot+zj
                stot = stot+xj*xj+yj*yj+zj*zj
!
             enddo
          enddo
       enddo
!$OMP END PARALLEL
!
#ifdef SERIAL
       rtot = rtot*rdv
       xtot = xtot*rdv
       ytot = ytot*rdv
       ztot = ztot*rdv
       stot = stot*rdv
#else
       loctot(1) = rtot*rdv
       loctot(2) = xtot*rdv
       loctot(3) = ytot*rdv
       loctot(4) = ztot*rdv
       loctot(5) = stot*rdv
!!!       write(6,*) loctot(1), myrank

#ifdef MIXEDPRECISION
# ifdef DOUBLE_P
       write/6,*) "ERROR: still not implemented 
# else
       call mpi_reduce(loctot,glotot,5,MPI_DOUBLE_PRECISION,mpi_sum,0,lbecomm,ierr)
# endif 
#else
       call mpi_reduce(loctot,glotot,5,MYMPIREAL,mpi_sum,0,lbecomm,ierr)
#endif

       rtot = glotot(1)/float(nprocs)
       xtot = glotot(2)/float(nprocs)
       ytot = glotot(3)/float(nprocs)
       ztot = glotot(4)/float(nprocs)
       stot = glotot(5)/float(nprocs)
#endif

       if(myrank.eq.0) then
          write(16,1001) itime
          write(16,1002) rtot
          write(16,1003) xtot,ytot,ztot,stot
          write(63,1004) itime, xtot, ytot, ztot, rtot, stot
          flush(16)
          flush(63)
       endif
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. diagno"
        endif
#endif
!
! formats...
!
1001  format(" Timestep ",i8)
1002  format("       mean rho ",1(e14.6,1x))
1003  format("       mean vel ",4(e14.6,1x))
1004  format(i8,5(e14.6,1x))
!
      return
      end subroutine diagno
