! =====================================================================
!     ****** LBE/hencol
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       hencol
!     DESCRIPTION
!       computes collision parameters and the forcing term
!       according to the following equations:
!       omega = 2.0/(6.0*svisc+1.0)
!       fgrad = 4.0*u0*svisc/(3.0*float(n)*float(n))
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
! =====================================================================
!
        subroutine hencol
!
        use storage
        implicit none
!
        omega = 2.d0/(6.d0*svisc+1.d0)
!        omega = 2.0_16/(6.0_16*svisc+1.0_16)
!        omega = (uno_qp+uno_qp)/((tre_qp*(uno_qp+uno_qp))*(svisc*uno_qp)+uno_qp)
!
! forcing term
!
#ifdef SERIAL
        fgrad = 4.0*u0*svisc/(5.0*real(n,mykind)*real(n,mykind))
#else
        fgrad = 4.0*u0*svisc/(5.0*real(lz,mykind)*real(lz,mykind))
#endif
!
        if(myrank == 0) then
           if(fgrad.le.0.0000001) then
              write( 6,*) "WARNING: volume forcing below 10e-7" 
              write(16,*) "WARNING: volume forcing below 10e-7" 
           endif
        endif
!
#ifdef DEBUG_1
       if(myrank == 0) then
          write(6,*) "DEBUG1: omega, fgrad",omega,fgrad
          write(6,*) "DEBUG1: Exiting from sub. hencol"
       endif
#endif

        return
        end subroutine hencol
