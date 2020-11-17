!=====================================================================
!     ****** LBE/collision
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for movef routine...
!       - movef (if -DORIGINAL)
!       - none  else
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
        subroutine collision(itime)
!
        use timing
        use storage
!
#ifdef CUDAFOR
        use cudafor
        use cuda_bgk
#endif
!
        implicit none
!
        integer:: itime, ierr, istat
!
#ifdef CUDAFOR
        type (dim3) :: grid, tBlock

! some stuff....
!!        tBlock = dim3 (192 ,1 ,1)
        tBlock = dim3 (256 ,2 ,1)
        grid = dim3 (ceiling (real(l)/ tBlock % x) , &
                     ceiling (real(n)/ tBlock % y) , 1)
#endif
!
! start timing
        call SYSTEM_CLOCK(countC0, count_rate, count_max)
        call time(tcountC0)
!
#ifdef ORIGINAL
# ifdef CUDAFOR
#  ifdef QQQQQ  
        call copyf_cuda <<< grid, tBlock >>> (l, n,       &
                dev_a02,dev_a04,dev_a05,dev_a06,dev_a11,  &
                dev_a13,dev_a14,dev_a15,                  &
                dev_b02,dev_b04,dev_b05,dev_b06,dev_b11,  &
                dev_b13,dev_b14,dev_b15) 
        istat = cudaThreadSynchronize()

        call col_fused_cuda <<< grid, tBlock >>> (l, n, fgrad, omega,  &
                dev_a02,dev_a04,dev_a05,dev_a06,dev_a11,  &
                dev_a13,dev_a14,dev_a15,dev_a19,          &
                dev_b02,dev_b04,dev_b05,dev_b06,dev_b11,  &
                dev_b13,dev_b14,dev_b15)           
#  endif

        call col_cuda <<< grid, tBlock >>> (l, n, fgrad, omega,  &
                dev_a02,dev_a04,dev_a05,dev_a06,dev_a11,  &
                dev_a13,dev_a14,dev_a15,dev_a19,          &
                dev_b02,dev_b04,dev_b05,dev_b06,dev_b11,  &
                dev_b13,dev_b14,dev_b15)           

!!!!!        istat = cudaThreadSynchronize()
# else
        call col(itime)
!!!        write(*,*) "check ---> in place collision!!!"; stop
# endif
#else
# ifdef CUDAFOR
#  ifdef PERIODIC
        call col_periodic_cuda <<< grid, tBlock >>> (l, n, &
                i_shift, fgrad, omega,  &
                dev_a02,dev_a04,dev_a05,dev_a06,dev_a11,  &
                dev_a13,dev_a14,dev_a15,dev_a19)
        i_shift = mod(i_shift+1,l)
#  else
        write(*,*) "ERROR: M+C version not implemented for CUDAFORT"
        stop
#  endif
# endif
#endif
!
#ifdef TWO_FIELDS
        call col_TF(itime)
#endif
!
#ifdef TWO_FIELDS
! do nothing
#else
# ifdef ORIGINAL
! do nothing
# else
!
#  ifdef TRYBOX
!
! take care that 2*border < l or m or n
! up
       call col_MC_box(itime,1       ,l         ,1       ,m         ,1       ,border    ,0)
! down
       call col_MC_box(itime,1       ,l         ,1       ,m         ,n-border,n         ,0)
! left
       call col_MC_box(itime,1       ,l         ,1       ,border    ,border+1,n-border-1,0)      
! right
       call col_MC_box(itime,1       ,l         ,m-border,m         ,border+1,n-border-1,0)
! front
       call col_MC_box(itime,1       ,border    ,border+1,m-border-1,border+1,n-border-1,0)
! rear
       call col_MC_box(itime,l-border,l         ,border+1,m-border-1,border+1,n-border-1,1)
! bulk
!       call col_MC_box(itime,border+1,l-border-1,border+1,m-border-1,border+1,n-border-1,1)
!
#  else
       call col_MC(itime)
#  endif
!
# endif
#endif
!
! stop timing
        call time(tcountC1)
        call SYSTEM_CLOCK(countC1, count_rate, count_max)
        time_coll = time_coll + real(countC1-countC0)/(count_rate)
        time_coll1 = time_coll1 + (tcountC1-tcountC0)
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. collision"
        endif
#endif
!
        return
        end subroutine collision
