!=====================================================================
!     ****** LBE/boundaries
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for boundaries routine...
!       - call bcond  else
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
        subroutine boundaries
!
        use timing
        use storage
#ifdef CUDAFOR
        use cudafor
        use cuda_bgk
#endif
!
        implicit none
        integer:: ierr, i
        real(mykind):: temp1, temp2
!
#ifdef OBSTACLES
! start timing
        call bcond_comm_noblock_packed_try
        call SYSTEM_CLOCK(countO0, count_rate, count_max)
        call time(tcountO0)
!
        call bcond_obs
!
! stop timing
        call time(tcountO1)
        call SYSTEM_CLOCK(countO1, count_rate, count_max)
        time_obs = time_obs + real(countO1-countO0)/(count_rate)
        time_obs1 = time_obs1 + (tcountO1-tcountO0)
#endif
!
#ifdef CUDAFOR
!
        type (dim3) :: grid, tBlock
        tBlock = dim3 (64 ,4 ,1)
        grid = dim3 (ceiling (real(l)/ tBlock % x) , &
                     ceiling (real(n)/ tBlock % y) , 1)
!
# ifdef PERIODIC
! front wall
!!        a11(1,:) = dev_a11(1,:)
!!        a13(1,:) = dev_a13(1,:)
!!        a14(1,:) = dev_a14(1,:)
! rear wall
!!        a02(l,:)  = dev_a02(l,:)
!!        a04(l,:)  = dev_a04(l,:)
!!        a05(l,:)  = dev_a05(l,:)
! up wall
!!        a02(:,1)  = dev_a02(:,1)
!!        a11(:,1)  = dev_a11(:,1)
!!        a15(:,1)  = dev_a15(:,1)
! down wall: 
!!        a04(:,n) = dev_a04(:,n)
!!        a06(:,n) = dev_a06(:,n)  
!!        a13(:,n) = dev_a13(:,n)
!!!
!!        call bcond
!!!
! front wall
!!!        dev_a11(l+1,:) = a11(l+1,:)
!!!        dev_a13(l+1,:) = a13(l+1,:)
!!!        dev_a14(l+1,:) = a14(l+1,:)
        call bcond_periodic_front <<< grid, tBlock >>> (l, n, dev_a11, dev_a13, dev_a14)
! rear wall
!        dev_a02(0,:)  = a02(0,:)
!        dev_a04(0,:)  = a04(0,:)
!        dev_a05(0,:)  = a05(0,:)
        call bcond_periodic_rear <<< grid, tBlock >>> (l, n, dev_a02, dev_a04, dev_a05)
! up wall
!        dev_a02(:,n+1)  = a02(:,n+1)
!        dev_a11(:,n+1)  = a11(:,n+1)
!        dev_a15(:,n+1)  = a15(:,n+1)
        ierr = cudaMemcpy(dev_a02(:,n+1),dev_a02(:,1),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a11(:,n+1),dev_a11(:,1),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a15(:,n+1),dev_a15(:,1),(l+2), & 
               cudaMemcpyDeviceToDevice)
! down wall: 
!        dev_a04(:,0) = a04(:,0)
!        dev_a06(:,0) = a06(:,0)  
!        dev_a13(:,0) = a13(:,0)
        ierr = cudaMemcpy(dev_a04(:,0),dev_a04(:,n),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a06(:,0),dev_a06(:,n),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a13(:,0),dev_a13(:,n),(l+2), & 
               cudaMemcpyDeviceToDevice)
# else
!
!-----------------------------------------------

        call bcond_noslip_top  <<< grid, tBlock >>> (l, n, u00, &
                          dev_a15, dev_a06,        &
                          dev_a11, dev_a04,        &
                          dev_a02, dev_a13)

#ifdef OLD_STUFF
        call bcond_cuda_top <<< grid, tBlock >>> (l, n, u00,    &
                dev_a13,dev_b02, &
                dev_a04,dev_b11, &
                dev_a06,dev_b15)

        ierr = cudaMemcpy(dev_a02(:,n+1),dev_b02(:,n+1),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a11(:,n+1),dev_b11(:,n+1),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a15(:,n+1),dev_b15(:,n+1),(l+2), & 
               cudaMemcpyDeviceToDevice)
#endif
!
        ierr = cudaDeviceSynchronize () 
!
!-----------------------------------------------
!
        call bcond_noslip_down  <<< grid, tBlock >>> (l, n, &
                          dev_a06, dev_a15,        &
                          dev_a04, dev_a11,        &
                          dev_a13, dev_a02)

#ifdef OLD_STUFF
        call bcond_cuda_down <<< grid, tBlock >>> (l, n,   &
                dev_a15,dev_b06, &
                dev_a11,dev_b04, &
                dev_a02,dev_b13)
!
        ierr = cudaMemcpy(dev_a04(:,0),dev_b04(:,0),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a06(:,0),dev_b06(:,0),(l+2), & 
               cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a13(:,0),dev_b13(:,0),(l+2), & 
              cudaMemcpyDeviceToDevice)
#endif
!
        ierr = cudaDeviceSynchronize () 
!
!-----------------------------------------------
!
!!        call bcond_cuda_rear <<< grid, tBlock >>> (l, n,   &
!!                      dev_a11,dev_bb04, &
!!                      dev_a14,dev_bb05, &
!!                      dev_a13,dev_bb02)
!
!!        ierr = cudaDeviceSynchronize () 
!
        call bcond_noslip_rear  <<< grid, tBlock >>> (l, n, &
                      dev_a04, dev_a11,        &
                      dev_a05, dev_a14,        &
                      dev_a02, dev_a13)
!
           do i = n-1, n
!!           ierr = cudaMemcpy(dev_a04(0,i),dev_bb04(i),1, &
!!!                  cudaMemcpyDeviceToDevice)
           ierr = cudaMemcpy(dev_a04(0,i-1),dev_a11(1,i),1, &
                  cudaMemcpyDeviceToDevice)
           enddo
!
           ierr = cudaMemcpy(dev_a05(0,n),dev_a14(1,n),1, &
                  cudaMemcpyDeviceToDevice)
!
!!!           ierr = cudaMemcpy(dev_a05(0,n),dev_bb05(n),1, &
!!!                 cudaMemcpyDeviceToDevice)
           ierr = cudaMemcpy(dev_a05(0,n+1),dev_a14(1,n+1),1, &
                  cudaMemcpyDeviceToDevice)
!
!!           ierr = cudaMemcpy(dev_a05(0,n+1),dev_bb05(n+1),1, &
!!                  cudaMemcpyDeviceToDevice)
!
!
!!!        call copy_rear <<< grid, tBlock >>> (l, n, &
!!!                      dev_bb02,dev_a02,          &
!!!                      dev_bb04,dev_a04,          &
!!!                      dev_bb05,dev_a05)
!
!-----------------------------------------------
!
!!        call bcond_cuda_front <<< grid, tBlock >>> (l, n,   &
!!                      dev_a02,dev_bb13, &
!!                      dev_a05,dev_bb14, &
!!                      dev_a04,dev_bb11)
!!        ierr = cudaDeviceSynchronize () 
!
        call bcond_noslip_front  <<< grid, tBlock >>> (l, n, &
                          dev_a11, dev_a04,        &
                          dev_a14, dev_a05,        &
                          dev_a13, dev_a02)
!        ierr = cudaDeviceSynchronize () 
!
! fix (to understand...)
!
        do i = n-1, n
           ierr = cudaMemcpy(dev_a13(l+1,i),dev_a02(l,i+1),1, &
                  cudaMemcpyDeviceToDevice)

!           ierr = cudaMemcpy(dev_a13(l+1,i),dev_bb13(i),1, &
!                  cudaMemcpyDeviceToDevice)
        enddo
!
!!        ierr = cudaMemcpy(dev_a14(l+1,n),dev_bb14(n),1, &
!!                  cudaMemcpyDeviceToDevice)
        ierr = cudaMemcpy(dev_a14(l+1,n),dev_a05(l,n),1, &
                  cudaMemcpyDeviceToDevice)
!
!!!        call copy_front <<< grid, tBlock >>> (l, n, &
!!!                      dev_bb11,dev_a11,          &
!!!                      dev_bb13,dev_a13,          &
!!!                      dev_bb14,dev_a14)
!
        ierr = cudaDeviceSynchronize () 
!
!-----------------------------------------------

# endif
#else
! 
# ifdef SENDRECV
#  ifdef SENDRECV_TRY
! standard sendrecv
!        call bcond_comm
!
! sendrecv + section (1 for every sendrecv)
        call bcond_comm_openmp
!        call bcond_comm_packed2
!        write(6,*) "ERROR: Wrong comm"
!        stop
#  else
        call bcond_comm_packed
#  endif
# endif
!
# ifdef NOBLOCK
#  ifdef NOBLOCK_TRY
        call bcond_comm_noblock
!        write(6,*) "ERROR: Wrong comm"
!        stop
#  else
        call bcond_comm_noblock_packed
!        call bcond_comm_noblock_packed_overlap
#  endif
# endif
!
# ifdef DUPLEX
        write(6,*) "ERROR: Wrong comm"
        call bcond_comm_duplex 
        stop
# endif
!        call bcond_bc
#endif
!
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. boundaries"
        endif
#endif
!
        return
      end subroutine boundaries
