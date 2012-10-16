!*******************************************************************************
!
!     MD 6.3.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Authors:           Don Berry, Joe Hughto
!     Last modified by:  Joe Hughto, 2012-Oct-15
!  ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Apply strain to the simulation box. If the user specifies one or two 
! components of the normal strain, the remaining two or one are calculated so 
! that volume of the box remains constant. If he specifies all three components, 
! we assume he knows what he's doing.  Only one component of the shear strain 
! may be specified.
!
      subroutine strain
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      integer      i0,i1,i2,i3,i4,i5,ix,iy,iz

! All strain rates default to a no-strain value. This is just a parameter that
! would never be selected as a realistic, physical strain rate.  The user then
! defines one, two or all three strain rates in the input parameter file. Thus
! there are eight possibilities. Note that xl0(1:3) are the gage dimensions,
! while xl(1:3) are the current dimensions of the strained volume.  Note that
! shear strains will be applied virtually and xl(1:3) will not be changed.
      i0=0; i1=0; i2=0; i3=0; i4=0; i5=0
      if(deps(1).ne.XNOSTRAIN) i0=1
      if(deps(2).ne.XNOSTRAIN) i1=2
      if(deps(3).ne.XNOSTRAIN) i2=4
      if(deps(4).ne.XNOSTRAIN) i3=1
      if(deps(5).ne.XNOSTRAIN) i4=2
      if(deps(6).ne.XNOSTRAIN) i5=4

      select case (i2+i1+i0)
      case(0)                   !no strain
         continue
      case(1)                   !x strain
         xl(1) = xl(1)+deps(1)*dt*xl0(1)
         xl(2) = xl0(2)*sqrt(xl0(1)/xl(1))
         xl(3) = xl0(3)*sqrt(xl0(1)/xl(1))
      case(2)                   !y strain
         xl(2) = xl(2)+deps(2)*dt*xl0(2)
         xl(1) = xl0(1)*sqrt(xl0(2)/xl(2))
         xl(3) = xl0(3)*sqrt(xl0(2)/xl(2))
      case(3)                   !x, y strains
         xl(1) = xl(1)+deps(1)*dt*xl0(1)
         xl(2) = xl(2)+deps(2)*dt*xl0(2)
         xl(3) = xl0(3)*(xl0(1)/xl(1))*(xl0(2)/xl(2))
      case(4)                   !z strain
         xl(3) = xl(3)+deps(3)*dt*xl0(3)
         xl(1) = xl0(1)*sqrt(xl0(3)/xl(3))
         xl(2) = xl0(2)*sqrt(xl0(3)/xl(3))
      case(5)                   !x, z strains
         xl(1) = xl(1)+deps(1)*dt*xl0(1)
         xl(3) = xl(3)+deps(3)*dt*xl0(3)
         xl(2) = xl0(2)*(xl0(1)/xl(1))*(xl0(3)/xl(3))
      case(6)                   !y, z strains
         xl(2) = xl(2)+deps(2)*dt*xl0(2)
         xl(3) = xl(3)+deps(3)*dt*xl0(3)
         xl(1) = xl0(1)*(xl0(2)/xl(2))*(xl0(3)/xl(3))
      case(7)                   !x, y, z strains
         xl(1) = xl(1)+deps(1)*dt*xl0(1)
         xl(2) = xl(2)+deps(2)*dt*xl0(2)
         xl(3) = xl(3)+deps(3)*dt*xl0(3)
      end select
      
!  ---------------------------------------------------------
!  Find particles that have been exlcluded from the box due to compression,
!  and wrap them around to periodic image inside.
      if(i3+i4+i5.eq.0) then
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            do k=1,3
               if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k) 
            enddo
         enddo
         !$omp end parallel do
      else
         sfac=0.d0
         ix=1
         iy=2
         iz=3
         select case (i3+i4+i5)
         case(1)
            sfac=deps(4)*(time-tstart)/2.d0
         case(2)
            sfac=deps(5)*(time-tstart)/2.d0
            ix=2
            iy=3
            iz=1
         case(4)
            sfac=deps(6)*(time-tstart)/2.d0
            ix=3
            iy=1
            iz=2
         end select
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            if(x(ix,i).gt.xl(ix)*(1+sfac*sfac)) then
               x(ix,i)=x(ix,i)-xl(ix)*(1+sfac*sfac)
            end if
            if(x(iy,i).lt.sfac*x(iz,i)) then
               x(iy,i)=x(iy,i)+xl(iy)
               x(iz,i)=x(iz,i)+xl(iy)*sfac
            end if
            if(x(iy,i).gt.sfac*x(iz,i)+xl(iy)*(1-sfac*sfac)) then
               x(iy,i)=x(iy,i)-xl(iy)
               x(iz,i)=x(iz,i)-xl(iy)*sfac
            end if
            if(x(iz,i).lt.sfac*x(iy,i)) then
               x(iz,i)=x(iz,i)+xl(iz)
               x(iy,i)=x(iy,i)+xl(iz)*sfac
            end if
            if(x(iz,i).gt.sfac*x(iy,i)+xl(iz)*(1-sfac*sfac)) then
               x(iz,i)=x(iz,i)-xl(iz)
               x(iy,i)=x(iy,i)-xl(iz)*sfac
            end if
         enddo
         !$omp end parallel do
      end if


      return
      end subroutine strain
