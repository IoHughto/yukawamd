!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2011, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2011-Nov-29
! ---------------------------------------------------------------------
!
!*******************************************************************************
! Integrate Newton's law using the velocity Verlet scheme as amended by
! Spreiter and Walter for a uniform magnetic field along the z-direction.
! See Q.Spreiter, M.Walter, J.Comp.Phys. 152:102 (1999). This subroutine
! implements the simple method they call "inversion of ez X v", which is
! specified by eqs. (13)-(15) of their paper.
!
! In the MPI version, each process updates all particles' positions and
! coordinates, so that no communication is needed.
!
      subroutine newton
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble)(a-h,o-z)
      include  'perf.h'

      call starttimer()   !DKB-perf (newton)

! ---------------------------------------------------------
! Integrate a full time step to get coordinates at t+dt.
      afac=0.5*(dt**2) 
      !$omp parallel do schedule(runtime)
      do i=0,n-1
        do k=1,3
          x(k,i)=x(k,i)+v(k,i)*dt+afac*a(k,i)
          if(x(k,i).lt.0)     x(k,i)=x(k,i)+xl(k)  !periodic boundary conditions
          if(x(k,i).gt.xl(k)) x(k,i)=x(k,i)-xl(k) 
        enddo
      enddo
      !$omp end parallel do

! ---------------------------------------------------------
! First stage of Spreiter-Walter velocity update: use velocities and
! accelerations at time t.
      halfdt=0.5*dt
      !$omp parallel do private(omega,omegap) schedule(runtime)
      do i=0,n-1
        omega=3.218404e-20*zii(i)*bfield/aii(i)
        omegap=halfdt*omega
        v(1,i)=(1.0-omegap**2)*v(1,i) + 2.0*omegap*v(2,i) + halfdt*(a(1,i)+omegap*a(2,i))
        v(2,i)=(1.0-omegap**2)*v(2,i) - 2.0*omegap*v(1,i) + halfdt*(a(2,i)-omegap*a(1,i))
        v(3,i)=v(3,i) + halfdt*a(3,i)
      enddo
      !$omp end parallel do

! ---------------------------------------------------------
! Calculate new acceleration at t+dt
      call accel

! ---------------------------------------------------------
! Second stage of Spreiter-Walter velocity update: use accelerations
! just calculated for time t+dt.  Also, complete calculation of the
! pressure tensor.
      vxx=0.0; vxy=0.0d0; vxz=0.0d0; vyy=0.0d0; vyz=0.0d0; vzz=0.0d0
      !$omp parallel do private(omega,omegap) reduction(+:vxx,vxy,vxz,vyy,vyz,vzz) schedule(runtime)
      do i=0,n-1
        omega=3.218404e-20*zii(i)*bfield/aii(i)
        omegap=halfdt*omega
        v(1,i)=(v(1,i)+halfdt*(a(1,i)+omegap*a(2,i)))/(1.0+omegap**2)
        v(2,i)=(v(2,i)+halfdt*(a(2,i)-omegap*a(1,i)))/(1.0+omegap**2)
        v(3,i)=v(3,i) + halfdt*a(3,i)
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

! ---------------------------------------------------------
! Update the total integration time
      time=time+dt

     !call dump_data    !DKB-debug

      call stoptimer(1,t_newton,ts_newton,n_newton)  !DKB-perf (newton)
      return
      end subroutine newton
