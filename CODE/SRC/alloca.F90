!=====================================================================
!     ****** LBE/alloca
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       alloca
!     DESCRIPTION
!
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine alloca
!
      use storage;
      use timing
!
      implicit none
!
!     define unit vector
!
      cx(  1    ) =  1
      cx(  2    ) =  1
      cx(  3    ) =  1
      cx(  4    ) =  1
      cx(  5    ) =  1
      cx(  6    ) =  0
      cx(  7    ) =  0
      cx(  8    ) =  0
      cx(  9    ) =  0
      cx( 10    ) = -1
      cx( 11    ) = -1
      cx( 12    ) = -1
      cx( 13    ) = -1
      cx( 14    ) = -1
      cx( 15    ) =  0
      cx( 16    ) =  0
      cx( 17    ) =  0
      cx( 18    ) =  0
      cx( 19    ) =  0
!
      cy(  1    ) = -1
      cy(  2    ) =  0
      cy(  3    ) =  1
      cy(  4    ) =  0
      cy(  5    ) =  0
      cy(  6    ) =  0
      cy(  7    ) =  1
      cy(  8    ) =  1
      cy(  9    ) =  1
      cy( 10    ) = -1
      cy( 11    ) =  0
      cy( 12    ) =  1
      cy( 13    ) =  0
      cy( 14    ) =  0
      cy( 15    ) =  0
      cy( 16    ) = -1
      cy( 17    ) = -1
      cy( 18    ) = -1
      cy( 19    ) =  0
!
      cz(  1    ) =  0
      cz(  2    ) = -1
      cz(  3    ) =  0
      cz(  4    ) =  1
      cz(  5    ) =  0
      cz(  6    ) =  1
      cz(  7    ) =  1
      cz(  8    ) =  0
      cz(  9    ) = -1
      cz( 10    ) =  0
      cz( 11    ) = -1
      cz( 12    ) =  0
      cz( 13    ) =  1
      cz( 14    ) =  0
      cz( 15    ) = -1
      cz( 16    ) = -1
      cz( 17    ) =  0
      cz( 18    ) =  1
      cz( 19    ) =  0
!
      allocate(obs(1:l,1:m,1:n))
      obs = 0
!
      allocate(a01(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a02(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a03(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a04(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a05(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a06(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a07(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a08(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a09(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a10(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a11(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a12(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a13(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a14(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a15(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a16(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a17(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a18(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a19(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
!
      allocate(b01(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b02(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b03(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b04(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b05(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b06(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b07(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b08(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b09(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b10(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b11(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b12(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b13(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b14(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b15(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b16(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b17(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b18(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b19(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
!
      a01 = huge(mykind)
      a02 = huge(mykind)
      a03 = huge(mykind)
      a04 = huge(mykind)
      a05 = huge(mykind)
      a06 = huge(mykind)
      a07 = huge(mykind)
      a08 = huge(mykind)
      a09 = huge(mykind)
      a10 = huge(mykind)
      a11 = huge(mykind)
      a12 = huge(mykind)
      a13 = huge(mykind)
      a14 = huge(mykind)
      a15 = huge(mykind)
      a16 = huge(mykind)
      a17 = huge(mykind)
      a18 = huge(mykind)
      a19 = huge(mykind)
!
      b01 = huge(mykind)
      b02 = huge(mykind)
      b03 = huge(mykind)
      b04 = huge(mykind)
      b05 = huge(mykind)
      b06 = huge(mykind)
      b07 = huge(mykind)
      b08 = huge(mykind)
      b09 = huge(mykind)
      b10 = huge(mykind)
      b11 = huge(mykind)
      b12 = huge(mykind)
      b13 = huge(mykind)
      b14 = huge(mykind)
      b15 = huge(mykind)
      b16 = huge(mykind)
      b17 = huge(mykind)
      b18 = huge(mykind)
      b19 = huge(mykind)
!
#ifdef ORIGINAL
! do nothing
#else
      c01 => null()
      c02 => null()
      c03 => null()
      c04 => null()
      c05 => null()
      c06 => null()
      c07 => null()
      c08 => null()
      c09 => null()
      c10 => null()
      c11 => null()
      c12 => null()
      c13 => null()
      c14 => null()
      c15 => null()
      c16 => null()
      c17 => null()
      c18 => null()
      c19 => null()
#endif

!set some variables
      radius = huge(mykind)
!
#ifdef PGI
! do nothing
#else
# ifdef FUGAKU
! do nothing
# else
        mem_start = get_mem()
# endif
#endif

#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. alloca"
        endif
#endif
!
# ifdef MEM_CHECK
        if(myrank == 0) then
           mem_stop = get_mem();
           write(6,*) "MEM_CHECK: after sub. alloca mem =", mem_stop
        endif
# endif
        return
        end subroutine alloca

