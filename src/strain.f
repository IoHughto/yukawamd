!*******************************************************************************
!
!     MD 6.2.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Authors:           Don Berry
!     Last modified by:  Don Berry, 2012-Apr-04
!  ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Calculate strain factors. These are the factors by which the box edge lengths
! will be multiplied each time step.  Only diagonal components of the strain
! rate tensor are allowed, not shearing components. If the user specifies one or
! two components, the remaining two or one are calculated so that volume of the
! box remains constant. If he specifies all three components, we assume he knows
! what he's doing.
!
      subroutine strainfactors
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      integer      i0,i1,i2

! All strain rates start out defaulted to a no-strain value. This is just a par-
! ameter that would never be selected as a realistic, physical strain rate.  The
! user then defines one, two or all three strain rates in the input parameter
! file. Thus there are eight possibilities.
      i0=0; i1=0; i2=0
      if(deps(1).ne.XNOSTRAIN) i0=1
      if(deps(2).ne.XNOSTRAIN) i1=2
      if(deps(3).ne.XNOSTRAIN) i2=4

      select case (i2+i1+i0)
        case(0)         !no strain
          strnfac(1)=1.d0
          strnfac(2)=1.d0
          strnfac(3)=1.d0
          deps(1)=0.d0
          deps(2)=0.d0
          deps(3)=0.d0
        case(1)         !xx strain
          strnfac(1)=1.d0+deps(1)*dt
          strnfac(2)=1.d0/sqrt(strnfac(1))
          strnfac(3)=1.d0/sqrt(strnfac(1))
          deps(2)=(strnfac(2)-1.d0)/dt
          deps(3)=(strnfac(3)-1.d0)/dt
        case(2)         !yy strain
          strnfac(2)=1.d0+deps(2)*dt
          strnfac(1)=1.d0/sqrt(strnfac(2))
          strnfac(3)=1.d0/sqrt(strnfac(2))
          deps(1)=(strnfac(1)-1.d0)/dt
          deps(3)=(strnfac(3)-1.d0)/dt
        case(3)         !xx, yy strains
          strnfac(1)=1.d0+deps(1)*dt
          strnfac(2)=1.d0+deps(2)*dt
          strnfac(3)=1.d0/(strnfac(1)*strnfac(2))
          deps(3)=(strnfac(3)-1.d0)/dt
        case(4)         !zz strain
          strnfac(3)=1.d0+deps(3)*dt
          strnfac(1)=1.d0/sqrt(strnfac(3))
          strnfac(2)=1.d0/sqrt(strnfac(3))
          deps(1)=(strnfac(1)-1.d0)/dt
          deps(2)=(strnfac(2)-1.d0)/dt
        case(5)         !xx, zz strains
          strnfac(1)=1.d0+deps(1)*dt
          strnfac(3)=1.d0+deps(3)*dt
          strnfac(2)=1.d0/(strnfac(1)*strnfac(3))
          deps(2)=(strnfac(2)-1.d0)/dt
        case(6)         !yy, zz strains
          strnfac(2)=1.d0+deps(2)*dt
          strnfac(3)=1.d0+deps(3)*dt
          strnfac(1)=1.d0/(strnfac(2)*strnfac(3))
          deps(1)=(strnfac(1)-1.d0)/dt
        case(7)         !xx, yy, zz strains
          strnfac(1)=1.d0+deps(1)*dt
          strnfac(2)=1.d0+deps(2)*dt
          strnfac(3)=1.d0+deps(3)*dt
      end select

      return
      end subroutine strainfactors



!*******************************************************************************
!  Apply strain to the simulation box. Only diagonal components of strain
!  tensor are allowed, not shearing components.
!
      subroutine strain
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      call starttimer()   !DKB-perf (strain)

!  ---------------------------------------------------------
!  Apply stretch/compress factor to each dimension.
      xl(1) = xl(1)*strnfac(1)
      xl(2) = xl(2)*strnfac(2)
      xl(3) = xl(3)*strnfac(3)

!  ---------------------------------------------------------
!  Find particles that have been exlcluded from the box due to compression,
!  and wrap them around to periodic image inside.
      !$omp parallel do schedule(runtime)
      do i=0,n-1
        do k=1,3
          if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k) 
        enddo
      enddo
      !$omp end parallel do
 
      call stoptimer(1,t_strain,ts_strain,n_strain)  !DKB-perf (strain)
      return
      end subroutine strain
