!*******************************************************************************
!
!    MD 6.3.0
! ---------------------------------------------------------------------
!    Copyright 2011, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-16
! ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Determine a file type based on a file name extension of up to six characters.
! If there is no extension return 'none'. If it is longer than six characters
! return 'long'.
!
      character(6) function filetype(filename)
      character*(*) filename
      integer   i,j

      j=len_trim(filename)
      i=scan(filename,'.',.true.)
      if(i.eq.0) then
         filetype='none'
      else if(j-i.gt.6) then
         filetype='long'
      else
         filetype=filename(i+1:j)
      endif
      return
         
      end function filetype
   


!*******************************************************************************
! Calculate average kinetic energy per particle.
!
      subroutine ttot(eka)
      use  md_globals
      implicit real(dble) (a-h,o-z)

      real(dble)   ekx

      ekx=0.
     !$omp parallel do reduction(+:ekx)
      do i=0,n-1
         ekx = ekx + aii(i)*(v(1,i)**2+v(2,i)**2+v(3,i)**2)
      enddo
     !$omp end parallel do
      eka=(0.5*xmass*ekx)/float(n)

      return
      end subroutine ttot



!*******************************************************************************
! By the equipartition theorem, the average kinetic energy per particle is
! related to the temperature by  eka=(3/2)kT. This subroutine scales particle
! velocities so that eka matches the specified temperature kT.
!
      subroutine tnorm
      use  md_globals
      implicit real(dble)(a-h,o-z)
        
      call ttot(eka)
      fac=sqrt(eka/(1.5*kT))
     !$omp parallel do
      do i=0,n-1
        v(:,i)=v(:,i)/fac
      enddo
     !$omp end parallel do
      return
      end subroutine tnorm



!*******************************************************************************
! Cancel center-of-mass motion.
!
!MD_6.2.0:
! Energy lost due center-of-mass motion cancellation is now added back as
! thermal energy, by calling tnorm after the cancellation.

      subroutine cancel_vcom
      use  md_globals
      implicit real(dble) (a-h,o-z)
      real(dble)  mtot      ! total mass
      real(dble)  ptot(3)   ! total momentum
      real(dble)  vcom(3)   ! center-of-mass velocity

      mtot = 0.0d0
      ptot = 0.0d0
     !$omp parallel do reduction(+:mtot,ptot)
      do i=0,n-1
        mtot = mtot+aii(i)
        ptot(:) = ptot(:)+aii(i)*v(:,i)
      enddo
      !$omp end parallel do
      vcom(:) = ptot(:)/mtot
     !$omp parallel do
      do i=0,n-1
         v(:,i) = v(:,i) - vcom(:)
      enddo
      !$omp end parallel do
      call tnorm

      end subroutine cancel_vcom



!*******************************************************************************
! Calculate distance between particles i and j. Return result in r.
!
      subroutine pdist(i,j,r)
      use  md_globals
      implicit real(dble) (a-h,o-z)
      integer i3, i4, i5, ix, iy, iz

      i3=0;i4=0;i5=0
      if(deps(4).ne.XNOSTRAIN) i3=1
      if(deps(5).ne.XNOSTRAIN) i4=2
      if(deps(6).ne.XNOSTRAIN) i5=4

      r=0.
      if(i3+i4+i5.eq.0) then
         do k=1,3
            yy=abs(x(k,i)-x(k,j))
            yy=min(yy,xl(k)-yy)   
            r=r+yy*yy
         enddo
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
         yy=abs(x(ix,i)-x(ix,j))
         yy=min(yy,xl(ix)*(1.d0+sfac*sfac-yy))
         r=r+yy*yy
         yyy=abs(x(iy,i)-x(iy,j))
         zzz=abs(x(iz,i)-x(iz,j))
         yy=yyy*yyy+zzz*zzz
         yy=min(yy,(yyy+xl(iy))**2.d0+(zzz+sfac*xl(iy))**2.d0)
         yy=min(yy,(yyy+sfac*xl(iz))**2.d0+(zzz+xl(iz))**2.d0)
         yy=min(yy,(yyy-xl(iy))**2.d0+(zzz-sfac*xl(iy))**2.d0)
         yy=min(yy,(yyy-sfac*xl(iz))**2.d0+(zzz-xl(iz))**2.d0)
         yy=min(yy,(yyy+xl(iy)+xl(iz)*sfac)**2.d0+(zzz+sfac*xl(iy)+xl(iz))**2.d0)
         yy=min(yy,(yyy+xl(iy)-xl(iz)*sfac)**2.d0+(zzz+sfac*xl(iy)-xl(iz))**2.d0)
         yy=min(yy,(yyy-xl(iy)+xl(iz)*sfac)**2.d0+(zzz-sfac*xl(iy)+xl(iz))**2.d0)
         yy=min(yy,(yyy-xl(iy)-xl(iz)*sfac)**2.d0+(zzz-sfac*xl(iy)-xl(iz))**2.d0)
         r=r+yy
      end if
      r=sqrt(r)!+1.e-15

      return
      end subroutine pdist

!*******************************************************************************
! Calculate distance between particles i and j. Return result in r.
!
      subroutine pvec(i,j,r,xx)
      use  md_globals
      implicit real(dble) (a-h,o-z)
      real(dble)   xx(3)     !relative position vector from i-th to j-th ions
      real(dble)   halfl(3)    !simulation box half edge lengths
      integer i3, i4, i5, ix, iy, iz

      i3=0;i4=0;i5=0
      if(deps(4).ne.XNOSTRAIN) i3=1
      if(deps(5).ne.XNOSTRAIN) i4=2
      if(deps(6).ne.XNOSTRAIN) i5=4

      halfl=0.5*xl

      r=0.d0
      if(i3+i4+i5.eq.0) then
         do k=1,3
            xx(k)=x(k,i)-x(k,j)
            if(xx(k).gt.+halfl(k)) xx(k)=xx(k)-xl(k)
            if(xx(k).lt.-halfl(k)) xx(k)=xx(k)+xl(k)
            r=r+xx(k)*xx(k)
         enddo
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

         halfl(ix)=halfl(ix)*(1+sfac*sfac)
         xx(ix)=x(ix,i)-x(ix,j)
         if(xx(ix).gt.+halfl(ix)) xx(ix)=xx(ix)-xl(ix)
         if(xx(ix).lt.-halfl(ix)) xx(ix)=xx(ix)+xl(ix)
         r=r+xx(ix)*xx(ix)

         xx(iy)=x(iy,i)-x(iy,j)
         xx(iz)=x(iz,i)-x(iz,j)
         yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         yyy=xx(iy)
         zzz=xx(iz)
         if(yy.gt.(yyy+xl(iy))**2.d0+(zzz+sfac*xl(iy))**2.d0) then
            xx(iy)=yyy+xl(iy)
            xx(iz)=zzz+sfac*xl(iy)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy+sfac*xl(iz))**2.d0+(zzz+xl(iz))**2.d0) then
            xx(iy)=yyy+sfac*xl(iz)
            xx(iz)=zzz+xl(iz)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy-xl(iy))**2.d0+(zzz-sfac*xl(iy))**2.d0) then
            xx(iy)=yyy-xl(iy)
            xx(iz)=zzz-sfac*xl(iy)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy-sfac*xl(iz))**2.d0+(zzz-xl(iz))**2.d0) then
            xx(iy)=yyy-sfac*xl(iy)
            xx(iz)=zzz-xl(iy)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy+xl(iy)+xl(iz)*sfac)**2.d0+(zzz+sfac*xl(iy)+xl(iz))**2.d0) then
            xx(iy)=yyy+xl(iy)+xl(iz)*sfac
            xx(iz)=zzz+sfac*xl(iy)+xl(iz)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy+xl(iy)-xl(iz)*sfac)**2.d0+(zzz+sfac*xl(iy)-xl(iz))**2.d0) then
            xx(iy)=yyy+xl(iy)-xl(iz)*sfac
            xx(iz)=zzz+sfac*xl(iy)-xl(iz)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy-xl(iy)+xl(iz)*sfac)**2.d0+(zzz-sfac*xl(iy)+xl(iz))**2.d0) then
            xx(iy)=yyy-xl(iy)+xl(iz)*sfac
            xx(iz)=zzz-sfac*xl(iy)+xl(iz)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         if(yy.gt.(yyy-xl(iy)-xl(iz)*sfac)**2.d0+(zzz-sfac*xl(iy)-xl(iz))**2.d0) then
            xx(iy)=yyy-xl(iy)-xl(iz)*sfac
            xx(iz)=zzz-sfac*xl(iy)-xl(iz)
            yy=xx(iy)*xx(iy)+xx(iz)*xx(iz)
         end if
         r=r+yy
      end if
      r=sqrt(r)!+1.e-15

      return
      end subroutine pvec
