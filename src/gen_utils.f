!*******************************************************************************
!
!    MD 6.3.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-09
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
! Calculate distance between particles i and j using sheared BCs if relevant. 
! Return result in rr.
!
      subroutine pdist(i,j,rr)
      use  md_globals
      implicit real(dble) (a-h,o-z)
      integer i3, i4, i5

      i3=0; i4=0; i5=0
      if(deps(4).ne.XNOS) i3=1
      if(deps(5).ne.XNOS) i4=2
      if(deps(6).ne.XNOS) i5=4

      rr=0.
      select case (i3+i4+i5)
      case(0)
         do k=1,3
            yy=abs(x(k,i)-x(k,j))
            yy=min(yy,xl(k)-yy)   
            rr=rr+yy*yy
         enddo
      case(1)
         yy=abs(x(1,i)-x(1,j))
         yy=min(yy,xl(1)-yy)
         rr=rr+yy*yy
         yy=abs(x(2,i)-x(2,j))
         zz=abs(x(3,i)-x(3,j))
         rrr=yy*yy+zz*zz
         rrr=min(rrr,(yy+xl(2))**2.d0+(zz+strnfac(4)/2.d0*xl(3))**2.d0
      end select
         rr=sqrt(rr)
      return
      end subroutine pdist
