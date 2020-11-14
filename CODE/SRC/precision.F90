module real_kinds

#ifdef SERIAL
!do nothing
#else
  use mpi
#endif
!
#ifdef HALF_P
  real*2 , parameter :: qq = 1.0
  integer, parameter :: hp = kind(qq)
#endif
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
  integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))
!
#ifdef DOUBLE_P
  integer, parameter :: mykind = dp
#else
# ifdef QUAD_P
  integer, parameter :: mykind = qp
# else
#  ifdef HALF_P
  integer, parameter :: mykind = hp 
#  else
  integer, parameter :: mykind = sp ! default
#  endif
# endif
#endif


#ifdef SERIAL
! do nothing
#else
# ifdef DOUBLE_P
  integer, parameter :: MYMPIREAL = MPI_DOUBLE_PRECISION
# else
  integer, parameter :: MYMPIREAL = MPI_REAL
# endif
!
# ifdef HALF_P
  Write(6,*) "MPI version + half precision still to implement"
  stop
# endif
!
# ifdef QUAD_P
  Write(6,*) "MPI version + quad precision still to implement"
  stop
# endif
!
#endif

end module real_kinds
