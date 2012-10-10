!*******************************************************************************
!    MD 6.3.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-09
! ---------------------------------------------------------------------
!
!*******************************************************************************


module  md_comm
      use md_globals

      character*7   ::  mp_type = 'F90+MPI'

      integer, private      ::  ierr

!*******************************************************************************
!*******************************************************************************

CONTAINS

!*******************************************************************************
!  This subroutine broadcasts all the run parameters from process 0 to the
!  other processes.

      subroutine bcast_parms
      use  md_types
      use  md_globals
      implicit none
      include  'mpif.h'

      integer, parameter :: NRPARM=80
      integer, parameter :: NIPARM=80
      integer, parameter :: NLPARM=20

      integer  i

! Arrays for broadcasting parameters to MPI processes.
      real(dble)  xparm(NRPARM)
      integer     iparm(NIPARM)
      logical     lparm(NLPARM)

   10 if(myrank.eq.0) then
         xparm(1)     = tstart
         xparm(2)     = dt
         xparm(3)     = tend
         xparm(10)    = zi
         xparm(11)    = ai
         xparm(20)    = rho
         xparm(21)    = rmax
         xparm(22:24) = xl
         xparm(25:30) = strnfac
         xparm(31)    = kT
         xparm(32)    = xmass
         xparm(40)    = qmin
         xparm(41)    = dq
         xparm(50)    = xmuc
         xparm(51)    = frp
         xparm(52)    = rccut
         xparm(60)    = aa
         xparm(61)    = bb
         xparm(62)    = cc
         xparm(63)    = xpacket
         xparm(64)    = rncut
         xparm(70)    = bfield
         xparm(71)    = efield
         xparm(72)    = q0
         xparm(73)    = w0
         xparm(74)    = tref
         call MPI_Bcast(xparm, NRPARM, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      else
         call MPI_Bcast(xparm, NRPARM, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         tstart     = xparm(1)
         dt         = xparm(2)
         tend       = xparm(3)
         zi         = xparm(10)
         ai         = xparm(11)
         rho        = xparm(20)
         rmax       = xparm(21)
         xl         = xparm(22:24)
         strnfac    = xparm(25:30)
         kT         = xparm(31)
         xmass      = xparm(32)
         qmin       = xparm(40)
         dq         = xparm(41)
         xmuc       = xparm(50)
         frp        = xparm(51)
         rccut      = xparm(52)
         aa         = xparm(60)
         bb         = xparm(61)
         cc         = xparm(62)
         xpacket    = xparm(63)
         rncut      = xparm(64)
         bfield     = xparm(70)
         efield     = xparm(71)
         q0         = xparm(72)
         w0         = xparm(73)
         tref       = xparm(74)
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)

   20  if(myrank.eq.0) then
         iparm(1)   = nwgroup
         iparm(2)   = nwsteps
         iparm(3)   = ngroup
         iparm(4)   = ntot
         iparm(5)   = nind
         iparm(10)  = ntnorm
         iparm(11)  = ncom
         iparm(12)  = nptensor
         iparm(13)  = nckpt
         iparm(14)  = nout
         iparm(20)  = n
         iparm(21)  = ni
         iparm(30)  = ngofr
         iparm(31)  = nbin
         iparm(32)  = nsbin
        !iparm(33)  = gtype
        !iparm(34)  = gspec
         iparm(40)  = irnd
         iparm(41)  = iseed
         call MPI_Bcast(iparm, NIPARM, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       else
         call MPI_Bcast(iparm, NIPARM, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         nwgroup    = iparm(1)
         nwsteps    = iparm(2)
         ngroup     = iparm(3)
         ntot       = iparm(4)
         nind       = iparm(5)
         ntnorm     = iparm(10) 
         ncom       = iparm(11)
         nptensor   = iparm(12)
         nckpt      = iparm(13)
         nout       = iparm(14)
         n          = iparm(20)
         ni         = iparm(21)
         ngofr      = iparm(30)
         nbin       = iparm(31)
         nsbin      = iparm(32)
        !gtype      = iparm(33)
        !gspec      = iparm(34)
         irnd       = iparm(40)
         iseed      = iparm(41)
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)

   30 if(myrank.eq.0) then
         lparm(1) = append
         call MPI_Bcast(lparm, NLPARM, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      else
         call MPI_Bcast(lparm, NLPARM, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         append      = lparm(1)
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)

   40 continue
         call MPI_Bcast(sim_type, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(coulomb, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(nuclear, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(ftype, 6, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         call MPI_barrier(MPI_COMM_WORLD,ierr)

      return
      end subroutine bcast_parms



      end module md_comm
