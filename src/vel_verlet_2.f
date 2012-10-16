!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2011, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2011-Dec-29
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
      integer i3,i4,i5,ix,iy,iz,inbox
      include  'perf.h'

      logical      do_measurements
      real(dble)   vp(3)

      call starttimer()   !DKB-perf (newton)

      i3=0;i4=0;i5=0
      if(deps(4).ne.XNOSTRAIN) i3=1
      if(deps(5).ne.XNOSTRAIN) i4=2
      if(deps(6).ne.XNOSTRAIN) i5=4

! ---------------------------------------------------------
! Integrate velocities and accelerations a full time step to get
! coordinates at t+dt.
      afac=0.5*(dt**2) 
      !$omp parallel
      !$omp do schedule(runtime)
      do i=0,n-1
         do k=1,3
            x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
         enddo
      enddo
      !$omp end do
      if(i3+i4+i5.eq.0) then
         !$omp do schedule(runtime)
         do i=0,n-1
            do k=1,3
               if(x(k,i).lt.0.0)   x(k,i)=x(k,i)+xl(k) !periodic boundary conditions
               if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k)
            enddo
         enddo
         !$omp end do
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
         !$omp do schedule(runtime)
         do i=0,n-1
            inbox=0
            do while(inbox.eq.0)
               inbox=1
               if(x(ix,i).gt.xl(ix)*(1+sfac*sfac)) then
                  x(ix,i)=x(ix,i)-xl(ix)*(1+sfac*sfac)
                  inbox=0
               end if
               if(x(ix,i).lt.0.0) then
                  x(ix,i)=x(ix,i)+xl(ix)*(1+sfac*sfac)
                  inbox=0
               end if
               if(x(iy,i).lt.sfac*x(iz,i)) then
                  x(iy,i)=x(iy,i)+xl(iy)
                  x(iz,i)=x(iz,i)+xl(iy)*sfac
                  inbox=0
               end if
               if(x(iy,i).gt.sfac*x(iz,i)+xl(iy)*(1-sfac*sfac)) then
                  x(iy,i)=x(iy,i)-xl(iy)
                  x(iz,i)=x(iz,i)-xl(iy)*sfac
                  inbox=0
               end if
               if(x(iz,i).lt.sfac*x(iy,i)) then
                  x(iz,i)=x(iz,i)+xl(iz)
                  x(iy,i)=x(iy,i)+xl(iz)*sfac
                  inbox=0
               end if
               if(x(iz,i).gt.sfac*x(iy,i)+xl(iz)*(1-sfac*sfac)) then
                  x(iz,i)=x(iz,i)-xl(iz)
                  x(iy,i)=x(iy,i)-xl(iz)*sfac
                  inbox=0
               end if
            enddo
         enddo
         !$omp end do
      end if

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
