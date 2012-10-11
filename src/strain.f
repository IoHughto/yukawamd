!*******************************************************************************
!
!     MD 6.2.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Authors:           Don Berry
!     Last modified by:  Don Berry, 2012-Jul-09
!  ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Apply strain to the simulation box.  Only normal strains (tension and compres-
! sion) are allowed.  If the user specifies one or two components, the remaining
! two or one are calculated so that volume of the box remains constant. If he
! specifies all three components, we assume he knows what he's doing.
!
      subroutine strain
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      integer      i0,i1,i2

! All strain rates default to a no-strain value. This is just a parameter that
! would never be selected as a realistic, physical strain rate.  The user then
! defines one, two or all three strain rates in the input parameter file. Thus
! there are eight possibilities. Note that xl0(1:3) are the gage dimensions,
! while xl(1:3) are the current dimensions of the strained volume.
      i0=0; i1=0; i2=0
      if(deps(1).ne.XNOSTRAIN) i0=1
      if(deps(2).ne.XNOSTRAIN) i1=2
      if(deps(3).ne.XNOSTRAIN) i2=4

      select case (i2+i1+i0)
        case(0)         !no strain
          continue
        case(1)         !x strain
          xl(1) = xl(1)+deps(1)*dt*xl0(1)
          xl(2) = xl0(2)*sqrt(xl0(1)/xl(1))
          xl(3) = xl0(3)*sqrt(xl0(1)/xl(1))
        case(2)         !y strain
          xl(2) = xl(2)+deps(2)*dt*xl0(2)
          xl(1) = xl0(1)*sqrt(xl0(2)/xl(2))
          xl(3) = xl0(3)*sqrt(xl0(2)/xl(2))
        case(3)         !x, y strains
          xl(1) = xl(1)+deps(1)*dt*xl0(1)
          xl(2) = xl(2)+deps(2)*dt*xl0(2)
          xl(3) = xl0(3)*(xl0(1)/xl(1))*(xl0(2)/xl(2))
        case(4)         !z strain
          xl(3) = xl(3)+deps(3)*dt*xl0(3)
          xl(1) = xl0(1)*sqrt(xl0(3)/xl(3))
          xl(2) = xl0(2)*sqrt(xl0(3)/xl(3))
        case(5)         !x, z strains
          xl(1) = xl(1)+deps(1)*dt*xl0(1)
          xl(3) = xl(3)+deps(3)*dt*xl0(3)
          xl(2) = xl0(2)*(xl0(1)/xl(1))*(xl0(3)/xl(3))
        case(6)         !y, z strains
          xl(2) = xl(2)+deps(2)*dt*xl0(2)
          xl(3) = xl(3)+deps(3)*dt*xl0(3)
          xl(1) = xl0(1)*(xl0(2)/xl(2))*(xl0(3)/xl(3))
        case(7)         !x, y, z strains
          xl(1) = xl(1)+deps(1)*dt*xl0(1)
          xl(2) = xl(2)+deps(2)*dt*xl0(2)
          xl(3) = xl(3)+deps(3)*dt*xl0(3)
      end select

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

      return
      end subroutine strain
