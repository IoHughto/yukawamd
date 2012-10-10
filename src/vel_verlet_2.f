!*******************************************************************************
!
!    MD 6.3.0
! ---------------------------------------------------------------------
!    Copyright 2011, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-09
! ---------------------------------------------------------------------
!
!*******************************************************************************
! Integrate Newton's law using a velocity Verlet scheme altered for a uniform
! magnetic field along the z-direction. This scheme first appeared in MD_5.2.2
! for nuclear pasta runs (neutron/proton systems).
!
! In the MPI version, each process updates all particles' positions and
! coordinates, so that no communication is needed.
!
      subroutine newton(do_measurements)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble)(a-h,o-z)
      include  'perf.h'
      integer  i0,i1,i2,i3,i4,i5

      logical      do_measurements
      real(dble)   vp(3)

      call starttimer()   !DKB-perf (newton)
      i0=0; i1=0; i2=0; i3=0; i4=0; i5=0
      if(deps(1).ne.XNOS) i0=1
      if(deps(2).ne.XNOS) i1=2
      if(deps(3).ne.XNOS) i2=4
      if(deps(4).ne.XNOS) i3=1
      if(deps(5).ne.XNOS) i4=2
      if(deps(6).ne.XNOS) i5=4

! ---------------------------------------------------------
! Integrate velocities and accelerations a full time step to get
! coordinates at t+dt.
      afac=0.5*(dt**2) 
      !$omp parallel
      !$omp do schedule(runtime)
      do i=0,n-1
         select case (i3+i4+i5)
         case(0)
            do k=1,3
               x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
               if(x(k,i).lt.0.0)   x(k,i)=x(k,i)+xl(k) !periodic boundary conditions
               if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k)
            enddo
         case(1)
            do k=1,3
               x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
            enddo
            if(x(1,i).lt.0.0)   x(1,i)=x(1,i)+xl(1) 
            if(x(1,i).gt.xl(1)) x(1,i)=x(1,i)-xl(1)
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
         case(2)
            do k=1,3
               x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
            enddo
            if(x(2,i).lt.0.0)   x(2,i)=x(2,i)+xl(2) 
            if(x(2,i).gt.xl(2)) x(2,i)=x(2,i)-xl(2)
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
         case(4)
            do k=1,3
               x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
            enddo
            if(x(3,i).lt.0.0)   x(3,i)=x(3,i)+xl(3) 
            if(x(3,i).gt.xl(3)) x(3,i)=x(3,i)-xl(3)
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
         end select
      enddo
      !$omp end do

! ---------------------------------------------------------
! Save old velocities. Integrate old accelerations a full time step to
! get an estimate of velocities at time t+dt.
      !$omp do schedule(runtime)
      do i=0,n-1
        vold(:,i)=v(:,i)
        v(:,i)=v(:,i)+dt*a(:,i)
      enddo
      !$omp end do
      !$omp end parallel

! ---------------------------------------------------------
! Calculate acceleration at t+dt due to two-particle interactions.
      call accel(do_measurements)

! Add acceleration due to a uniform B-field along z, if B-field is non-zero.
      if(bfield.ne.0.0d0) then
        !$omp parallel do private(omega) schedule(runtime)
        do i=0,n-1
          omega=(3.218404d-20)*zii(i)*bfield/aii(i)
          a(1,i) = a(1,i)+v(2,i)*omega
          a(2,i) = a(2,i)-v(1,i)*omega
        enddo   
        !$omp end parallel do
      endif

! ---------------------------------------------------------
! Using old velocities, integrate new accelerations a full time step to
! get another estimate of velocities at t+dt. Set velocities at t+dt to
! the average of the two velocity estimates.
      !$omp parallel do private(vp) schedule(runtime)
      do i=0,n-1
        vp(:)=vold(:,i)+dt*a(:,i)
        v(:,i)=0.5*(v(:,i)+vp(:))
      enddo

! ---------------------------------------------------------
! If measurements were called for, complete calculation of the pressure tensor,
! and calculate average kinetic energy.
      if(do_measurements) then
        vxx=0.0d0; vxy=0.0d0; vxz=0.0d0; vyy=0.0d0; vyz=0.0d0; vzz=0.0d0
        !$omp parallel do reduction(+:vxx,vxy,vxz,vyy,vyz,vzz) schedule(runtime)
        do i=0,n-1
          vxx=vxx+aii(i)*v(1,i)*v(1,i)
          vxy=vxy+aii(i)*v(1,i)*v(2,i)
          vxz=vxz+aii(i)*v(1,i)*v(3,i)
          vyy=vyy+aii(i)*v(2,i)*v(2,i)
          vyz=vyz+aii(i)*v(2,i)*v(3,i)
          vzz=vzz+aii(i)*v(3,i)*v(3,i)
        enddo
        !$omp end parallel do
        ek=(0.5d0*xmass)*(vxx+vyy+vzz)/float(n)
        pp(1,1)=pp(1,1)+xmass*vxx
        pp(1,2)=pp(1,2)+xmass*vxy
        pp(1,3)=pp(1,3)+xmass*vxz
        pp(2,2)=pp(2,2)+xmass*vyy
        pp(2,3)=pp(2,3)+xmass*vyz
        pp(3,3)=pp(3,3)+xmass*vzz
        pp(2,1)=pp(1,2)
        pp(3,1)=pp(1,3)
        pp(3,2)=pp(2,3)
        pp(:,:)=pp(:,:)/(xl(1)*xl(2)*xl(3))
      else
        ek=0.0d0
        pp(:,:)=0.0d0
      endif

! ---------------------------------------------------------
! Update the total integration time
      time=time+dt

     !call dump_data    !DKB-debug

      call stoptimer(1,t_newton,ts_newton,n_newton)  !DKB-perf (newton)
      return
      end subroutine newton
