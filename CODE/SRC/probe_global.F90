!=====================================================================
!     ****** LBE/probe_global
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       probe_global (in global variables)
!     DESCRIPTION
!       Diagnostic subroutine:
!       probe popolations and/or velocities in one single point	
!       write on unit 67 (probe_g.dat)
!     INPUTS
!       itime --> timestep
!       i     --> x coordinate
!       j     --> y coordinate
!       k     --> z coordinate
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j,k,itime
!       real variables used: rto,xj,yj,zj
!
!     *****
!=====================================================================
!
        subroutine probe_global(itime,i0,j0,k0)
!
        use storage
        implicit none
!
        integer i0,j0,k0,itime
!
        integer:: offsetX, offsetY, offsetZ
        real(mykind) :: rho, xj, yj, zj
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07,x08,x09,x10
        real(mykind) :: x11,x12,x13,x14,x15,x16,x17,x18,x19
!
        offsetX = mpicoords(1)*l
        offsetY = mpicoords(2)*m
        offsetZ = mpicoords(3)*n

        if((k0.gt.offsetZ).and.(k0.le.(offsetZ+n))) then
          if((j0.gt.offsetY).and.(j0.le.(offsetY+m))) then
            if((i0.gt.offsetX).and.(i0.le.(offsetX+l))) then
!
               x01 = a01(i0-offsetX,j0-offsetY,k0-offsetZ)
               x02 = a02(i0-offsetX,j0-offsetY,k0-offsetZ)
               x03 = a03(i0-offsetX,j0-offsetY,k0-offsetZ)
               x04 = a04(i0-offsetX,j0-offsetY,k0-offsetZ)
               x05 = a05(i0-offsetX,j0-offsetY,k0-offsetZ)
               x06 = a06(i0-offsetX,j0-offsetY,k0-offsetZ)
               x07 = a07(i0-offsetX,j0-offsetY,k0-offsetZ)
               x08 = a08(i0-offsetX,j0-offsetY,k0-offsetZ)
               x09 = a09(i0-offsetX,j0-offsetY,k0-offsetZ)
               x10 = a10(i0-offsetX,j0-offsetY,k0-offsetZ)
               x11 = a11(i0-offsetX,j0-offsetY,k0-offsetZ)
               x12 = a12(i0-offsetX,j0-offsetY,k0-offsetZ)
               x13 = a13(i0-offsetX,j0-offsetY,k0-offsetZ)
               x14 = a14(i0-offsetX,j0-offsetY,k0-offsetZ)
               x15 = a15(i0-offsetX,j0-offsetY,k0-offsetZ)
               x16 = a16(i0-offsetX,j0-offsetY,k0-offsetZ)
               x17 = a17(i0-offsetX,j0-offsetY,k0-offsetZ)
               x18 = a18(i0-offsetX,j0-offsetY,k0-offsetZ)
               x19 = a19(i0-offsetX,j0-offsetY,k0-offsetZ)
!
               rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                    +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                    +x17 + x18 + x19

               xj = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)
               yj = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)
               zj = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)
!
               write(67,1002) itime, xj/rho, yj/rho, zj/rho, rho
!
           endif
         endif
       endif
!
#ifdef DEBUG_2
      if(myrank == 0) then
         write(6,*) "DEBUG2: Exiting from sub. probe_global", i0,j0,k0
!         write(6,*) "DEBUG2: relative index", i0-offsetX 
!         write(6,*) "DEBUG2: relative index", j0-offsetY 
!         write(6,*) "DEBUG2: relative index", k0-offsetZ 
      endif
#endif
!
!format
1002    format(i8,4(e14.6,1x))

       return
       end subroutine probe_global
