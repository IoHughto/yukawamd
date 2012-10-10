!*******************************************************************************
!
!     MD 6.3.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Authors:           Don Berry, Joe Hughto
!     Last modified by:  Don Berry, 2012-Oct-09
!  ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Calculate strain factors. These are the factors by which the box edge lengths
! will be multiplied each time step.  If the user specifies one or two of the 
! diagonal components of the strain rate tensor, the remaining two or one are 
! calculated so that volume of the box remains constant. If he specifies all 
! three components, we assume he knows what he's doing.  Only one of the off-
! diagonal components may be selected.  If more than one is chosen, the code 
! will simply not shear.  No error message is given.
!
      subroutine strainfactors
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      integer  i0,i1,i2,i3,i4,i5

! All strain rates start out defaulted to a no-strain value. This is just a par-
! ameter that would never be selected as a realistic, physical strain rate.  The
! user then defines one, two or all three strain rates in the input parameter
! file. Thus there are eight possibilities.
      i0=0; i1=0; i2=0; i3=0; i4=0; i5=0
      if(deps(1).ne.XNOS) i0=1
      if(deps(2).ne.XNOS) i1=2
      if(deps(3).ne.XNOS) i2=4
      if(deps(4).ne.XNOS) i3=1
      if(deps(5).ne.XNOS) i4=2
      if(deps(6).ne.XNOS) i5=4

      if(i0+i1+i2.eq.0 .or. i3+i4+i5.eq.0) then
         select case (i2+i1+i0)
         case(0)                !no strain
            strnfac(1)=1.d0
            strnfac(2)=1.d0
            strnfac(3)=1.d0
            deps(1)=0.d0
            deps(2)=0.d0
            deps(3)=0.d0
         case(1)                !xx strain
            strnfac(1)=1.d0+deps(1)*dt
            strnfac(2)=1.d0/sqrt(strnfac(1))
            strnfac(3)=1.d0/sqrt(strnfac(1))
            deps(2)=(strnfac(2)-1.d0)/dt
            deps(3)=(strnfac(3)-1.d0)/dt
         case(2)                !yy strain
            strnfac(2)=1.d0+deps(2)*dt
            strnfac(1)=1.d0/sqrt(strnfac(2))
            strnfac(3)=1.d0/sqrt(strnfac(2))
            deps(1)=(strnfac(1)-1.d0)/dt
            deps(3)=(strnfac(3)-1.d0)/dt
         case(3)                !xx, yy strains
            strnfac(1)=1.d0+deps(1)*dt
            strnfac(2)=1.d0+deps(2)*dt
            strnfac(3)=1.d0/(strnfac(1)*strnfac(2))
            deps(3)=(strnfac(3)-1.d0)/dt
         case(4)                !zz strain
            strnfac(3)=1.d0+deps(3)*dt
            strnfac(1)=1.d0/sqrt(strnfac(3))
            strnfac(2)=1.d0/sqrt(strnfac(3))
            deps(1)=(strnfac(1)-1.d0)/dt
            deps(2)=(strnfac(2)-1.d0)/dt
         case(5)                !xx, zz strains
            strnfac(1)=1.d0+deps(1)*dt
            strnfac(3)=1.d0+deps(3)*dt
            strnfac(2)=1.d0/(strnfac(1)*strnfac(3))
            deps(2)=(strnfac(2)-1.d0)/dt
         case(6)                !yy, zz strains
            strnfac(2)=1.d0+deps(2)*dt
            strnfac(3)=1.d0+deps(3)*dt
            strnfac(1)=1.d0/(strnfac(2)*strnfac(3))
            deps(1)=(strnfac(1)-1.d0)/dt
         case(7)                !xx, yy, zz strains
            strnfac(1)=1.d0+deps(1)*dt
            strnfac(2)=1.d0+deps(2)*dt
            strnfac(3)=1.d0+deps(3)*dt
         end select
         
         select case (i3+i4+i5)
         case(0)
            strnfac(4)=0.d0
            strnfac(5)=0.d0
            strnfac(6)=0.d0
         case(1)
            strnfac(4)=deps(4)*dt
            strnfac(5)=0.d0
            strnfac(6)=0.d0
         case(2)
            strnfac(4)=0.d0
            strnfac(5)=deps(5)*dt
            strnfac(6)=0.d0
         case(4)
            strnfac(4)=0.d0
            strnfac(5)=0.d0
            strnfac(6)=deps(6)*dt
         end select
      end if
      
      return
      end subroutine strainfactors



!*******************************************************************************
!  Apply strain to the simulation box. 
!
      subroutine strain
      use  md_types
      use  md_globals
      implicit real(dble)(a-h,o-z)
      include  'perf.h'
      integer  i0,i1,i2,i3,i4,i5

      call starttimer()   !DKB-perf (strain)
      i0=0; i1=0; i2=0; i3=0; i4=0; i5=0
      if(deps(1).ne.XNOS) i0=1
      if(deps(2).ne.XNOS) i1=2
      if(deps(3).ne.XNOS) i2=4
      if(deps(4).ne.XNOS) i3=1
      if(deps(5).ne.XNOS) i4=2
      if(deps(6).ne.XNOS) i5=4

!  ---------------------------------------------------------
!  Apply stretch/compress factor to each dimension.
      if(i3+i4+i5.eq.0) then
         xl(1) = xl(1)*strnfac(1)
         xl(2) = xl(2)*strnfac(2)
         xl(3) = xl(3)*strnfac(3)
      else 
         select case (i3+i4+i5)
         case(1)
            xl(1)=xl(1)*(1.d0+strnfac(4)*strnfac(4)/4.d0)
         case(2)
            xl(2)=xl(2)*(1.d0+strnfac(5)*strnfac(5)/4.d0)
         case(4)
            xl(3)=xl(3)*(1.d0+strnfac(6)*strnfac(6)/4.d0)
         end select
      end if

!  ---------------------------------------------------------
!  Find particles that have been exlcluded from the box due to compression,
!  and wrap them around to periodic image inside.
      select case (i3+i4+i5)
      case(0)
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            do k=1,3
               if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k) 
            enddo
         enddo
         !$omp end parallel do
      case(1)
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            if(x(1,i).gt.xl(1)) then
               x(1,i)=x(1,i)-xl(1)
            end if
            if(x(2,i).lt.strnfac(4)/2.d0*x(3,i)) then
               x(2,i)=x(2,i)+xl(2)
               x(3,i)=x(3,i)+xl(2)*strnfac(4)/2.d0
            end if
            if(x(2,i).gt.strnfac(4)/2.d0*x(3,i)+xl(2)*(1-0.25*strnfac(4)**2.d0) then
               x(2,i)=x(2,i)-xl(2)
               x(3,i)=x(3,i)-xl(2)*strnfac(4)/2.d0
            end if
            if(x(3,i).lt.strnfac(4)/2.d0*x(2,i)) then
               x(3,i)=x(3,i)+xl(3)
               x(2,i)=x(2,i)+xl(3)*strnfac(4)/2.d0
            end if
            if(x(3,i).gt.strnfac(4)/2.d0*x(2,i)+xl(3)*(1-0.25*strnfac(4)**2.d0) then
               x(3,i)=x(3,i)-xl(3)
               x(2,i)=x(2,i)-xl(3)*strnfac(4)/2.d0
            end if
         enddo
         !$omp end parallel do
      case(2)
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            if(x(2,i).gt.xl(2)) then
               x(2,i)=x(2,i)-xl(2)
            end if
            if(x(3,i).lt.strnfac(5)/2.d0*x(1,i)) then
               x(3,i)=x(3,i)+xl(3)
               x(1,i)=x(1,i)+xl(3)*strnfac(5)/2.d0
            end if
            if(x(3,i).gt.strnfac(5)/2.d0*x(1,i)+xl(3)*(1-0.25*strnfac(5)**2.d0) then
               x(3,i)=x(3,i)-xl(3)
               x(1,i)=x(1,i)-xl(3)*strnfac(5)/2.d0
            end if
            if(x(1,i).lt.strnfac(5)/2.d0*x(3,i)) then
               x(1,i)=x(1,i)+xl(1)
               x(3,i)=x(3,i)+xl(1)*strnfac(5)/2.d0
            end if
            if(x(1,i).gt.strnfac(5)/2.d0*x(3,i)+xl(1)*(1-0.25*strnfac(5)**2.d0) then
               x(1,i)=x(1,i)-xl(1)
               x(3,i)=x(3,i)-xl(1)*strnfac(5)/2.d0
            end if
         enddo
         !$omp end parallel do
      case(4)
         !$omp parallel do schedule(runtime)
         do i=0,n-1
            if(x(3,i).gt.xl(3)) then
               x(3,i)=x(3,i)-xl(3)
            end if
            if(x(1,i).lt.strnfac(6)/2.d0*x(2,i)) then
               x(1,i)=x(1,i)+xl(1)
               x(2,i)=x(2,i)+xl(1)*strnfac(6)/2.d0
            end if
            if(x(1,i).gt.strnfac(6)/2.d0*x(2,i)+xl(1)*(1-0.25*strnfac(6)**2.d0) then
               x(1,i)=x(1,i)-xl(1)
               x(2,i)=x(2,i)-xl(1)*strnfac(6)/2.d0
            end if
            if(x(2,i).lt.strnfac(6)/2.d0*x(1,i)) then
               x(2,i)=x(2,i)+xl(2)
               x(1,i)=x(1,i)+xl(2)*strnfac(6)/2.d0
            end if
            if(x(2,i).gt.strnfac(6)/2.d0*x(1,i)+xl(2)*(1-0.25*strnfac(6)**2.d0) then
               x(2,i)=x(2,i)-xl(2)
               x(1,i)=x(1,i)-xl(2)*strnfac(6)/2.d0
            end if
         enddo
         !$omp end parallel do
      end select
      
      call stoptimer(1,t_strain,ts_strain,n_strain)  !DKB-perf (strain)
      return
      end subroutine strain
