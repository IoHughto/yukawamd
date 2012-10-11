!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Author:            Don Berry
!    Last modified by:  Don Berry, 2012-May-17
! ---------------------------------------------------------------------
!
! This file contains subroutines for calculating potential energy, accelerations,
! pressure tensor, and pressure for pure-ion simulations (single species of ions)
! of systems interacting via a screened coulomb potential. The total number of
! ions in the system is n, the ions' charge number is zi, and their mass number
! is ai. The screening length is 1/xmuc. All subroutines use the direct particle-
! particle method on a general purpose computer.
!
! This file can be compiled as a serial, OpenMP, MPI, or MPI+OpenMP code. Serial
! and OpenMP variants use stub MPI routines in md_comm_ser.f, and definitions of
! MPI constants in mpif_stubs.h.
!
! Target particles are distributed over MPI processes. Source particles are
! distributed over OpenMP threads.
!
 
 
!*******************************************************************************
! Calculate average potential energy per ion. Calculate instantaneous pressure.
!
      subroutine vtot_ion_pure(evavg,pressure)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'perf.h'
      include 'mpif.h'

      real(dble)   evavg     !avg potential energy per ion
      real(dble)   pressure  !pressure
      real(dble)   evx       !potential energy of a particle pair
      real(dble)   xtmp(2)   !for accumulating potential energy and virial
      real(dble)   ytmp(2)   !(ditto)
      real(dble)   xx        !temporary variable, for distance calculation
      real(dble)   r2,r      !distance squared and distance between two ions
      real(dble)   rccut2    !cutoff radius squared
      real(dble)   evcut     !potential at cutoff radius

      call starttimer()   !DKB-perf (vtot)
      call starttimer()   !DKB-perf (calc_v)

      rccut2=rccut**2
      if( rccut .le. sqrt(xl(1)**2+xl(2)**2+xl(3)**2)/2.0 ) then
         evcut=exp(-xmuc*rccut)/rccut
      else
         evcut=0.0d0
      endif

      ytmp=0.0d0
      !$omp parallel private(xx,r2,r,evx)
      do 100 i=myrank,n-2,nprocs
         !$omp do reduction(+:ytmp), schedule(runtime)
         do 90 j=i+1,n-1
            r2=0.
            do k=1,3
               xx=abs(x(k,i)-x(k,j))
               xx=min(xx,xl(k)-xx)
               r2=r2+xx*xx
            enddo
            if(r2.le.rccut2) then
               r=sqrt(r2)
               evx = exp(-xmuc*r)/r                 !potential energy of pair ij
               ytmp(1) = ytmp(1) + (evx-evcut)      !accumulate into total for i
               ytmp(2) = ytmp(2) + (1.+xmuc*r)*evx  !contribution to virial
            endif
   90    continue
         !$omp end do
  100 continue
      !$omp end parallel
      ytmp(1) = frp*zi*zi*vc*ytmp(1)
      ytmp(2) = frp*zi*zi*vc*ytmp(2)
      call stoptimer(1,t_calc_v,ts_calc_v,n_calc_v)  !DKB-perf (calc_v)

      call MPI_allreduce(ytmp,xtmp,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      evavg = xtmp(1)/float(n)
      pressure = rho*(kT+xtmp(2)/(3.0*n))
      call stoptimer(1,t_vtot,ts_vtot,n_vtot)        !DKB-perf (vtot)
      return

      end subroutine vtot_ion_pure



!*******************************************************************************
! Calculate accelerations due to two-body screened Coulomb force. Use Newton's
! 3rd law, so only n*(n-1)/2 interactions must be calculated. If do_measurements
! =.true., also calculate the system potential energy, pressure, and pressure
! tensor. If only the potential energy and pressure are needed, subroutine
! vtot_ion_pure would be more efficient.
!
      subroutine accel_ion_pure(do_measurements)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'perf.h'
      include 'mpif.h'

      logical      do_measurements

      real(dble)   xx(3)     !relative position vector from i-th to j-th ions
      real(dble)   halfl(3)  !half edge lengths
      real(dble)   r2,r      !distance-squared and distance between two ions
      real(dble)   rccut2    !cutoff radius squared
      real(dble)   evcut     !potential at cutoff radius
      real(dble)   fc        !r-dependence of screened Coulomb interaction
      real(dble)   fcx(3)    !force of particle j on i
      real(dble)   fi(3)     !for accumulating force of j particles on i
      real(dble), allocatable :: fj(:,:)  !for reaction of i on j

      call starttimer()   !DKB-perf (accel)
      call starttimer()   !DKB-perf (calc_a)

      allocate(fj(3,0:n+2))   !locations n, n+1, n+2 are used in MPI_allreduce.
      fj=0.0
      halfl=0.5*xl

      rccut2=rccut**2
      if( rccut .le. sqrt(xl(1)**2+xl(2)**2+xl(3)**2)/2.0 ) then
        evcut=exp(-xmuc*rccut)/rccut
      else
        evcut=0.0d0
      endif

      evx=0.0d0
      virx=0.0d0
      fi(:)=0.0d0
      pxx=0.0d0; pxy=0.0d0; pxz=0.0d0; pyy=0.0d0; pyz=0.0d0; pzz=0.0d0

!-------------------------------------------------------------------------------
! If measurements were called for, use this longer section, which  computes
! potential energy, virial and pressure tensor in addition to all the forces.
      if(do_measurements) then

      !$omp parallel private(xx,r2,r,fc,fcx)
      do 100 i=myrank,n-2,nprocs
        !$omp do reduction(+:evx,virx,fi,pxx,pxy,pxz,pyy,pyz,pzz)
        do 90 j=i+1,n-1
          r2=0.0d0
          do k=1,3
            xx(k)=x(k,i)-x(k,j)
            if(xx(k).gt.+halfl(k)) xx(k)=xx(k)-xl(k)
            if(xx(k).lt.-halfl(k)) xx(k)=xx(k)+xl(k)
            r2=r2+xx(k)*xx(k)
          enddo
          if(r2.le.rccut2) then
            r=sqrt(r2)
            fc = exp(-xmuc*r)/(r*r2)
            evx = evx + (fc*r2-evcut)
            fc = (1.+xmuc*r)*fc
            virx = virx + fc*r2
            fcx(:)  = fc*xx(:)
            fi(:)   = fi(:)   + fcx(:)      !action of j on i
            fj(:,j) = fj(:,j) - fcx(:)      !reaction of i on j
            pxx = pxx + xx(1)*fcx(1)
            pxy = pxy + xx(1)*fcx(2)
            pxz = pxz + xx(1)*fcx(3)
            pyy = pyy + xx(2)*fcx(2)
            pyz = pyz + xx(2)*fcx(3)
            pzz = pzz + xx(3)*fcx(3)
          endif
   90   continue
        !$omp end do
        !$omp single
        fj(:,i) = fj(:,i)+fi(:)
        fi(:)=0.0d0
        !$omp end single
  100 continue
      !$omp end parallel

!-------------------------------------------------------------------------------
! If measurements were not called for, use this section, which is a quicker
! calculation of just forces.
      else

      !$omp parallel private(xx,r2,r,fc,fcx)
      do 200 i=myrank,n-2,nprocs
        !$omp do reduction(+:fi)
        do 190 j=i+1,n-1
          r2=0.0d0
          do k=1,3
            xx(k)=x(k,i)-x(k,j)
            if(xx(k).gt.+halfl(k)) xx(k)=xx(k)-xl(k)
            if(xx(k).lt.-halfl(k)) xx(k)=xx(k)+xl(k)
            r2=r2+xx(k)*xx(k)
          enddo
          if(r2.le.rccut2) then
            r=sqrt(r2)
            fc = (1.+xmuc*r)*exp(-xmuc*r)/(r*r2)
            fcx(:)  = fc*xx(:)
            fi(:)   = fi(:)   + fcx(:)      !action of j on i
            fj(:,j) = fj(:,j) - fcx(:)      !reaction of i on j
          endif
  190   continue
        !$omp end do
        !$omp single
        fj(:,i) = fj(:,i)+fi(:)
        fi(:)=0.0d0
        !$omp end single
  200 continue
      !$omp end parallel

!-------------------------------------------------------------------------------
      endif
      call stoptimer(1,t_calc_a,ts_calc_a,n_calc_a)  !DKB-perf (calc_a)

! Pack evx, virx and the pressure tensor into locations n, n+1 and n+2 of fj
! so we need to do only one MPI_allreduce. Note that evx, virx and the pressure
! tensor will all be zero if only the force calculation was done.
      fj(1,n) = evx
      fj(2,n) = virx
      fj(3,n) = 0.0d0      !not used
      fj(1,n+1) = pxx
      fj(2,n+1) = pxy
      fj(3,n+1) = pxz
      fj(1,n+2) = pyy
      fj(2,n+2) = pyz
      fj(3,n+2) = pzz
      call MPI_allreduce(fj,a,3*(n+3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$omp parallel do schedule(runtime)
      do i=0,n-1
         a(:,i) = frp*zi*zi*vc*a(:,i)/(ai*xmass)
      enddo
      !$omp end parallel do
      if(do_measurements) then
        ev      = frp*zi*zi*vc*a(1,n)
        vir     = frp*zi*zi*vc*a(2,n)
        pp(1,1) = frp*zi*zi*vc*a(1,n+1)
        pp(1,2) = frp*zi*zi*vc*a(2,n+1)
        pp(1,3) = frp*zi*zi*vc*a(3,n+1)
        pp(2,2) = frp*zi*zi*vc*a(1,n+2)
        pp(2,3) = frp*zi*zi*vc*a(2,n+2)
        pp(3,3) = frp*zi*zi*vc*a(3,n+2)
        ev = ev/float(n)           !average potential energy per ion
        px = rho*(kT+vir/(3.0*n))  !(instantaneous) pressure
      else
        ev  = 0.0d0
        vir = 0.0d0
        pp  = 0.0d0
        px  = 0.0d0
      endif

      deallocate(fj)
      call stoptimer(1,t_accel,ts_accel,n_accel)     !DKB-perf (accel)

      return
      end subroutine accel_ion_pure
