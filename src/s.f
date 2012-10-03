!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Mar-30
! ---------------------------------------------------------------------
!
!  Calculate the static structure factor S(q) for the particle species specified
!  by gspec. This is done by numerically evaluating the integral
!
!         S(q) = 1 + rhox \int d^3r (g(r)-1)exp(iq.r),
!
!  where q and r here are 3-vectors, and rhox is the density of the particle
!  species gspec. We do this for each group of g(r) data, and for each value of
!  q given by  q = qmin, qmin+dq, qmin+2*dq, ... qmin+nsbin*dq.
!
!  Only MPI process 0 performs this subroutine. Other MPI processes immediately
!  return.
!*******************************************************************************

!DKB-todo : Not sure rectangular boxes are handled correctly. Check statements
!DKB-todo+:    referring to xl and xlmin.

      subroutine s
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)

      real(dble), allocatable::  s2(:), smean(:), r(:)

      if(myrank.ne.0) return
      if(nsbin.eq.0) return

!  Allocate arrays for static structure function.
      allocate(ss(nsbin,ngroup))  ! for storing static structure function
      allocate(cs(ngroup))

      nbmax=nbin/sqrt(3.)
      xlmin=min(xl(1),xl(2),xl(3))
      dd=xlmin/float(nbin)*sqrt(3.)/2
!  Calculate for the given particle species at a density of rhox.
      allocate(r(nbin))
      ngspec=0
      do i=0,n-1
         ngspec=ngspec+1
      enddo
      rhox=float(ngspec)/float(n)*rho
      do k=1,nbmax
         r(k)=(float(k)-.5)*dd
      enddo
       
!  Evaluate the integral for S(q) for each group of measurements (i=1...ngroup),
!  for the values q=qmin+j*dq, j=1...nsbin.
      do i=1,ngroup
         q=qmin-dq
         do j=1,nsbin
            q=q+dq
            ss(j,i)=0.
            do k=1,nbmax
               rq=q*r(k)
               ss(j,i)=ss(j,i)+(r(k)**2)*(sin(rq)/rq)*(gg(k,i,0)-1.)
            enddo
         ss(j,i)=4.d0*PI*dd*rhox*ss(j,i)+1.
         enddo
      enddo
      deallocate(r)

!  For each q, compute mean(S(q)) and std.dev.(S(q)) over the groups of data.
      allocate(s2(nsbin))
      allocate(smean(nsbin))
      open(UNIT=9,FILE='s.out'//suffix,STATUS='UNKNOWN')
      q=qmin-dq
      do j=1,nsbin
         q=q+dq
         smean(j)=0.
         s2(j)=0.
         do i=1,ngroup
            smean(j)=smean(j)+ss(j,i)
            s2(j)=s2(j)+ss(j,i)**2
         enddo
         smean(j)=smean(j)/float(ngroup)
         s2(j)=s2(j)/float(ngroup)
         s2(j)=sqrt(abs(s2(j)-smean(j)**2)/float(ngroup))
         write(9,10000) q,smean(j),s2(j)
      enddo
10000 format(f15.8,2x,1pe15.8,2x,1pe15.8)
      close(9)
      deallocate(s2)
      deallocate(smean)

      deallocate(ss)
      deallocate(cs)

      return
      end subroutine s
