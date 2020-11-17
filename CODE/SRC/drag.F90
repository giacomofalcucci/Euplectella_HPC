!=====================================================================
!     ****** LBE/drag
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       drag
!     DESCRIPTION
!       subroutine to check drag, to coustomize as necessary
!	- defined/tested for infinite cylinder
!	- write on file 66 (drag.*****.dat where ***** is taskid)
!     INPUTS
!       itime
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       to normalize for correct drag this is the gnuplot string
!       rep 'drag.0000.dat' $2:($3+$4)/(0.5*0.1*0.1*(2*R*lz))
!       where:
!            u0 --> velocity at inlet
!            R  --> radius of cylinder
!            òz --> height of cylinder
!
!       work correctly only do 1D decomposition or the cykinfer fits 
!       completely one task.
!       
!
!     *****
!=====================================================================
!
        subroutine drag(itime)
!
        use storage
        use timing
        implicit none
!
        integer:: itime
        integer:: xstart,xstop
        integer:: ystart,ystop
        integer:: zstart,zstop
        integer:: delta, i
!
        real(mykind) :: fluxX, fluxY, fluxZ
!
! drag & drop
        radius = 100
!
! cilindro 
        do delta  = 52, 67, 5
!
! absolute values
           xstart = lx/4-radius-delta
           xstop  = lx/4+radius+delta
           ystart = ly/2-radius-delta
           ystop  = ly/2+radius+delta
           zstart = 1
           zstop  = lz
!        
           fluxX = 0.d0
           fluxY = 0.d0
           fluxZ = 0.d0
!
           call flux_X(itime,xstart,xstop,ystart,ystop,zstart,zstop,fluxX)
           call flux_Y(itime,xstart,xstop,ystart,ystop,zstart,zstop,fluxY)
           write(66,*) delta, itime, fluxX, fluxY, myrank
!
        enddo 
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. drag"
           write(6,*) "DEBUG2: ", xstart, xstop, ystart, ystop, zstart, zstop
        endif
#endif
!
        return
        end subroutine drag
