!*******************************************************************************
!     MD 6.2.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Author:            Don Berry
!     Last modified by:  Don Berry, 2012-Apr-06
!  ---------------------------------------------------------------------
!
!*******************************************************************************


module  md_comm
   use md_globals

   character*3       ::  mp_type = 'F90'

   integer, private  ::  ierr

   interface MPI_send
      module procedure MPI_send_i
      module procedure MPI_send_d
   end interface

   interface MPI_recv
      module procedure MPI_recv_i
      module procedure MPI_recv_d
   end interface

   interface MPI_bcast
      module procedure MPI_bcast_i1xn
      module procedure MPI_bcast_d3xn
      module procedure MPI_bcast_species1xn
      module procedure MPI_bcast_d1xn
   end interface

   interface MPI_reduce
      module procedure MPI_reduce_da
   end interface

   interface MPI_allreduce
      module procedure MPI_allreduce_d
      module procedure MPI_allreduce_da
      module procedure MPI_allreduce_da2
   end interface


!*******************************************************************************
!*******************************************************************************

CONTAINS


!*******************************************************************************
!  In this subroutine, process 0 broadcasts all the run parameters to the other
!  processes. Since there is no need for such a broadcast in the serial code, we
!  just return.

   subroutine bcast_parms
   implicit none
   return
   end subroutine bcast_parms



!*******************************************************************************
!  Stub for MPI_reduce routine. Works only for MPI_DOUBLE_PRECISION datatype.

   subroutine MPI_reduce_da(x,xx,count,datatype,op,root,comm,ierror)
   implicit none
   real(dble)  x(:,:)
   real(dble)  xx(:,:)
   integer     count, datatype, op, root, comm, ierror
   xx = x
   ierror=0
   return
   end subroutine MPI_reduce_da



!*******************************************************************************
!  Stub for MPI_allreduce routine. Works only for scalar MPI_DOUBLE_PRECISION.

   subroutine MPI_allreduce_d(x,xx,count,data_type,op,comm,ierror)
   implicit none
   real(dble)  x
   real(dble)  xx
   integer     count, data_type, op, comm, ierror
   xx = x
   ierror=0
   return
   end subroutine MPI_allreduce_d



!*******************************************************************************
!  Stub for MPI_allreduce routine. Works only for array of MPI_DOUBLE_PRECISION.

   subroutine MPI_allreduce_da(x,xx,count,data_type,op,comm,ierror)
   implicit none
   real(dble)  x(:)
   real(dble)  xx(:)
   integer     count, data_type, op, comm, ierror
   xx = x
   ierror=0
   return
   end subroutine MPI_allreduce_da



!*******************************************************************************
!  Stub for MPI_allreduce routine. Works only for array of MPI_DOUBLE_PRECISION.

   subroutine MPI_allreduce_da2(x,xx,count,data_type,op,comm,ierror)
   implicit none
   real(dble)  x(:,:)
   real(dble)  xx(:,:)
   integer     count, data_type, op, comm, ierror
   xx = x
   ierror=0
   return
   end subroutine MPI_allreduce_da2



!*******************************************************************************
!  Stubs for MPI_init, MPI_finalize, MPI_comm_size and MPI_comm_rank routines.

   subroutine MPI_init(ierror)
   implicit none
   integer    ierror
   ierror=0
   return
   end subroutine MPI_init

   subroutine MPI_finalize(ierror)
   implicit none
   integer    ierror
   ierror=0
   return
   end subroutine MPI_finalize

   subroutine MPI_comm_size(comm,nprocs,ierror)
   implicit none
   integer    comm
   integer    nprocs
   integer    ierror
   nprocs=1
   ierror=0
   return
   end subroutine MPI_comm_size

   subroutine MPI_comm_rank(comm,myrank,ierror)
   implicit none
   integer    comm
   integer    myrank
   integer    ierror
   myrank=0
   ierror=0
   return
   end subroutine MPI_comm_rank



!*******************************************************************************
!  Stubs for MPI_send routines.

   subroutine MPI_send_i(buff,count,datatype,dest,tag,comm,ierror)
   implicit none
   integer    buff(*)
   integer    count, datatype, dest, tag, comm, ierror
   ierror = 0
   return
   end subroutine MPI_send_i

   subroutine MPI_send_d(buff,count,datatype,dest,tag,comm,ierror)
   implicit none
   real(dble)   buff(*)
   integer    count, datatype, dest, tag, comm, ierror
   ierror = 0
   return
   end subroutine MPI_send_d



!*******************************************************************************
!  Stubs for MPI_recv routines.

   subroutine MPI_recv_i(buff,count,datatype,source,tag,comm,status,ierror)
   implicit none
   include 'mpif.h'
   integer    buff(*)
   integer    count, datatype, source, tag, comm, status(MPI_STATUS_SIZE), ierror
   ierror = 0
   return
   end subroutine MPI_recv_i

   subroutine MPI_recv_d(buff,count,datatype,source,tag,comm,status,ierror)
   implicit none
   include 'mpif.h'
   real(dble)   buff(*)
   integer    count, datatype, source, tag, comm, status(MPI_STATUS_SIZE), ierror
   ierror = 0
   return
   end subroutine MPI_recv_d



!*******************************************************************************
!  Stubs for MPI_bcast routine.

   subroutine MPI_bcast_i1xn(buf, count, datatype, root, comm, ierror)
   implicit  none
   integer     buf(*)
   integer     count, datatype, root, comm, ierror
   ierror=0
   return
   end subroutine MPI_bcast_i1xn

   subroutine MPI_bcast_d3xn(buf, count, datatype, root, comm, ierror)
   implicit  none
   real(dble)  buf(3,*)
   integer     count, datatype, root, comm, ierror
   ierror=0
   return
   end subroutine MPI_bcast_d3xn

   subroutine MPI_bcast_species1xn(buf, count, datatype, root, comm, ierror)
   use md_types
   implicit none
   type(species)  buf(*)
   integer        count, datatype, root, comm, ierror
   ierror=0
   return
   end subroutine MPI_bcast_species1xn

   subroutine MPI_bcast_d1xn(buf, count, datatype, root, comm, ierror)
   implicit none
   real(dble)  buf(*)
   integer     count, datatype, root, comm, ierror
   ierror=0
   return
   end subroutine MPI_bcast_d1xn




!*******************************************************************************
!  Stub for MPI_barrier.

   subroutine MPI_barrier(comm,ierror)
   implicit none
   integer    comm
   integer    ierror
   ierror = 0
   return
   end subroutine MPI_barrier



end module md_comm
