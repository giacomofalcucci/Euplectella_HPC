! =====================================================================
!     ****** LBE/storage.f90
!
!     COPYRIGHT
!       (c) 2000-2008 by CASPUR/G.Amati
!     NAME
!       storage
!     DESCRIPTION
!       module for storage
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables defined:  l,n,l1,n1
!       real variables defined: lx, ly, dt, p0, p1, p2, rf, qf
!                               svisc, u0, omega,fgrad
!                               a02, a04, a05, a06, a11, a13, a14, a15, a19
!
!     *****
! =====================================================================
!
        module storage
!
        use real_kinds
!
#ifdef SERIAL
!do nothing
#else
        use mpi
#endif
!
        integer:: lx, ly, lz              ! global size        (MPI stuff)
        integer:: proc_x, proc_y, proc_z  ! task decomposition (MPI stuff)
!
        integer:: l, m, n                    ! local (task) size
        integer:: l1, m1, n1
        integer:: up(2),down(2),left(2)
        integer:: front(2),rear(2),right(2)
!
        integer, parameter::  mpid=3      ! mpi dimension
        integer, parameter::  zeroi=0      ! hint to help compiler...
        integer, parameter::  unoi=1       ! hint to help compiler...
!
        integer, parameter::  border=32    ! hint to help compiler...
!
#if defined PGI || defined ARM
        real(dp), parameter::  zero_qp=0.d0      ! hint to help compiler...
        real(dp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(dp), parameter::  tre_qp=3.d0      ! hint to help compiler...
!
        real(dp), parameter :: rf_qp = 3.d0
        real(dp), parameter :: qf_qp = 1.5d0
!
! 1/3
        real(dp), parameter :: p0_qp = uno_qp/rf_qp
! 1/18
        real(dp), parameter :: p1_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+uno_qp))
! 1/36
        real(dp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp)) 
#else
        real(qp), parameter::  zero_qp=0.d0      ! hint to help compiler...
        real(qp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(qp), parameter::  tre_qp=3.d0      ! hint to help compiler...
!
        real(qp), parameter :: rf_qp = 3.d0
        real(qp), parameter :: qf_qp = 1.5d0
!
! 1/3
        real(qp), parameter :: p0_qp = uno_qp/rf_qp
! 1/18
        real(qp), parameter :: p1_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+uno_qp))
! 1/36
        real(qp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp)) 
#endif


        integer:: nprocs, myrank, lbecomm, localcomm
        integer:: rear_task, front_task
        integer:: left_task, right_task
        integer:: down_task, up_task
        integer:: imax, imin                    ! obstacle stuff
        integer:: jmax, jmin                    ! obstacle stuff
        integer:: kmax, kmin                    ! obstacle stuff
        integer:: nobs                          ! #of obstacles per task
        integer:: xyplane, xzplane, yzplane, myxrank, yzcomm
        integer:: prgrid(mpid)
        integer:: mpicoords(mpid)
        integer:: gsizes(3),lsizes(3),start_idx(3)
        integer:: offset(3)
        integer:: dump3d
        integer:: buffer_size
        integer:: i_shift
        integer:: ipad,jpad,kpad
        integer:: flag1, flag2, flag3
!
        real(mykind), dimension(1:19) :: cx,cy,cz
!
#ifdef SERIAL
!do nothing
#else
        integer(kind=MPI_OFFSET_KIND):: file_offset
#endif
!
        logical remdims(mpid)
        logical periodic(mpid)
        logical rreorder
!
#ifdef ORIGINAL
        real(mykind), dimension(:,:,:), allocatable :: a01,a02,a03,a04,a05
        real(mykind), dimension(:,:,:), allocatable :: a06,a07,a08,a09,a10
        real(mykind), dimension(:,:,:), allocatable :: a11,a12,a13,a14,a15
        real(mykind), dimension(:,:,:), allocatable :: a16,a17,a18,a19
!
#else
        real(mykind), dimension(:,:,:), pointer :: a01,a02,a03,a04,a05
        real(mykind), dimension(:,:,:), pointer :: a06,a07,a08,a09,a10
        real(mykind), dimension(:,:,:), pointer :: a11,a12,a13,a14,a15
        real(mykind), dimension(:,:,:), pointer :: a16,a17,a18,a19
#endif

#ifdef ORIGINAL
        real(mykind), dimension(:,:,:), allocatable :: b01,b02,b03,b04,b05
        real(mykind), dimension(:,:,:), allocatable :: b06,b07,b08,b09,b10
        real(mykind), dimension(:,:,:), allocatable :: b11,b12,b13,b14,b15
        real(mykind), dimension(:,:,:), allocatable :: b16,b17,b18,b19
#else
        real(mykind), dimension(:,:,:), pointer :: b01,b02,b03,b04,b05,b06, &
                                                   b07,b08,b09,b10,b11,b12, &
                                                   b13,b14,b15,b16,b17,b18, &
                                                   b19
        real(mykind), dimension(:,:,:), pointer :: c01,c02,c03,c04,c05,c06, &
                                                   c07,c08,c09,c10,c11,c12, &
                                                   c13,c14,c15,c16,c17,c18, &
                                                   c19
#endif
!
        integer, dimension(:,:,:), pointer :: obs
!
        real(mykind):: svisc, u0, u00, fgrad
        real(mykind):: u0x, u0y, u0z
        real(mykind):: u_inflow
        integer:: mydev, ndev              ! openacc variables

!
! correct casting
#ifdef MIXEDPRECISION
# ifdef DOUBLE_P
        real(qp) :: omega
        real(qp) :: radius
        real(qp), parameter :: zero = zero_qp
        real(qp), parameter :: uno  = uno_qp
        real(qp), parameter :: tre  = tre_qp
! 
        real(qp), parameter :: rf = rf_qp
        real(qp), parameter :: qf = qf_qp
        real(qp), parameter :: p0 = p0_qp
        real(qp), parameter :: p1 = p1_qp
        real(qp), parameter :: p2 = p2_qp
# else
#  ifdef HALF_P
        real(sp) :: omega
        real(sp) :: radius
        real(sp), parameter :: zero = zero_qp
        real(sp), parameter :: uno  = uno_qp
        real(sp), parameter :: tre  = tre_qp
!
        real(sp), parameter :: rf = rf_qp
        real(sp), parameter :: qf = qf_qp
        real(sp), parameter :: p0 = p0_qp
        real(sp), parameter :: p1 = p1_qp
        real(sp), parameter :: p2 = p2_qp
#  else
        real(dp) :: omega
        real(dp) :: radius
        real(dp), parameter :: zero = zero_qp
        real(dp), parameter :: uno  = uno_qp
        real(dp), parameter :: tre  = tre_qp
!
        real(dp), parameter :: rf = rf_qp
        real(dp), parameter :: qf = qf_qp
        real(dp), parameter :: p0 = p0_qp
        real(dp), parameter :: p1 = p1_qp
        real(dp), parameter :: p2 = p2_qp
#  endif
# endif
#else
        real(mykind) :: omega
        real(mykind) :: radius
        real(mykind), parameter :: zero = zero_qp
        real(mykind), parameter :: uno  = uno_qp
        real(mykind), parameter :: tre  = tre_qp
! 
        real(mykind), parameter :: rf = rf_qp
        real(mykind), parameter :: qf = qf_qp
        real(mykind), parameter :: p0 = p0_qp
        real(mykind), parameter :: p1 = p1_qp
        real(mykind), parameter :: p2 = p2_qp
#endif

        end module  storage
