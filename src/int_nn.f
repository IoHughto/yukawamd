!*******************************************************************************
!    MD 6.3.0
! -----------------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Charles J. Horowitz, Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-16
! -----------------------------------------------------------------------------
!
! This file contains subroutines for calculating potential energy, accelera-
! tions, pressure and pressure tensor of a classical N-body system of nucleons
! interacting via the two-body central potential described in C.J.Horowitz,
! et al., Phys.Rev.C 69, 045804 (2004). The potential consists of an effective
! part modeling the strong nuclear force, and a screened coulomb part.
!
! The nuclear interaction is turned on by setting global character variable
! nuclear='HPP' in module md_globals. Otherwise it is off.
!
! The screened Coulomb interaction is turned on by setting global character
! variable coulomb='screened-coulomb' in module md_globals. Otherwise it is off.
!
! These subroutines use the direct particle-particle method on a general purpose
! computer, but can be compiled for serial, OpenMP, MPI, or MPI+OpenMP codes.
! Serial and OpenMP variants use stub MPI routines in md_comm_ser.f, and defini-
! tions of MPI constants in mpif_stubs.h.
!
! Target particles are distributed over MPI processes. Source particles are
! distributed over OpenMP threads.
!  -----------------------------------------------------------------------------
!
!


!*******************************************************************************
! Calculate average potential energy per nucleon, and instantaneous pressure.
! Result are returned in evavg and pressure.
!
      subroutine vtot_nn(evavg,pressure)
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'perf.h'
      include 'mpif.h'

      real(dble)   evavg       !avg potential energy per nucleon
      real(dble)   pressure    !pressure
      real(dble)   expfac      !(scratch variable)
      real(dble)   xtmp(2)     !for accumulating potential energy and virial
      real(dble)   ytmp(2)     !scratch variable for potential energy and pressure
      real(dble)   xx          !temporary variable, for distance calculation
      real(dble)   r           !distance squared and distance between two nucleons
      real(dble)   xpacketi    !1/Lambda
      real(dble)   xpacket2i   !1/(2*Lambda)
      logical      coulomb_on  !turn coulomb interaction on/off
      logical      nuclear_on  !turn nuclear interaction on/off
      real(dble)   xbbcc       ! = bb+cc for n-n, p-p interactions,
                               ! = bb-cc for n-p interaction
      real(dble)   evcut       !screened Coulomb potential at rccut

      call starttimer()   !DKB-perf (vtot)
      call starttimer()   !DKB-perf (calc_v)

      coulomb_on = coulomb.eq.'screened-coulomb'
      nuclear_on = nuclear.eq.'HPP'

      xpacketi = 1./xpacket
      xpacket2i = 1./(2.0*xpacket)

      if( rccut .le. sqrt(xl(1)**2+xl(2)**2+xl(3)**2)/2.0 ) then
        evcut=vc*exp(-xmuc*rccut)/rccut
      else
        evcut=0.0d0
      endif

      ytmp=0.0d0
      xtmp=0.0d0
      !$omp parallel private(xx,r,expfac)
      do 100 i=myrank,n-2,nprocs
        !$omp do reduction(+:xtmp)
         do 90 j=i+1,n-1
            r=0.d0
            call pdist(i,j,r)
            if(nuclear_on .and. r.le.rncut) then
               expfac=exp(-r*r*xpacket2i)
               if( zii(i).eq.zii(j) ) then !n-n, p-p
                  xbbcc = bb+cc
               else             !n-p
                  xbbcc = bb-cc
               endif
               xtmp(1) = xtmp(1) + (aa*expfac+xbbcc)*expfac !potential
               xtmp(2) = xtmp(2) + (r*r*xpacketi)*(2.*aa*expfac+xbbcc)*expfac !virial
            endif
            if(coulomb_on .and. r.le.rccut) then
               if((zii(i).eq.1.d0).and.(zii(j).eq.1.d0)) then !p-p
                  expfac = vc*exp(-xmuc*r)/r
                  xtmp(1) = xtmp(1) + (expfac-evcut) !potential
                  xtmp(2) = xtmp(2) + (1.+xmuc*r)*expfac !virial
               endif
            endif
 90      continue
        !$omp end do
 100  continue
      !$omp end parallel
      call stoptimer(1,t_calc_v,ts_calc_v,n_calc_v) !DKB-perf (calc_v)

      call MPI_allreduce(xtmp,ytmp,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      evavg=ytmp(1)/float(n)	 
      pressure = rho*(kT+ytmp(2)/(3.0*n))
      call stoptimer(1,t_vtot,ts_vtot,n_vtot)        !DKB-perf (vtot)
      return
      end subroutine vtot_nn



!*******************************************************************************
! Calculate acceleration of each nucleon due to two-body forces from all others.
! Use Newton's 3rd law, so only n*(n-1)/2 interactions must be calculated. If
! do_measurements=.true., also calculate the system potential energy, pressure,
! and pressure tensor. If only the potential energy and pressure are needed, sub-
! routine vtot_nn would be more efficient.
!
      subroutine accel_nn(do_measurements)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'perf.h'
      include 'mpif.h'

      logical      do_measurements !.true.=calculate potential energy, pressure,
                                   !   pressure tensor

      real(dble)   xx(3)       !position vector of nucleon i relative to j
      real(dble)   r           !distance squared and distance between two particles
      real(dble)   evcut       !screened Coulomb potential at rccut
      real(dble)   evx         !scratch for calcluating evavg
      real(dble)   vir,virx    !for accumulating virial
      real(dble)   pxx,pxy,pxz,pyy,pyz,pzz   !scratch for pressure tensor
      real(dble)   xpacketi    !1/Lambda
      real(dble)   xpacket2i   !1/(2*Lambda)
      real(dble)   fn          !r-dependence of nuclear interaction
      real(dble)   fc          !r-dependence of screened Coulomb interaction
      real(dble)   fncx(3)     !combined nuclear and coulomb force
      real(dble)   fi(3)       !for accumulating force of j particles on i
      real(dble), allocatable :: fj(:,:)  !for reaction of i on j particles
      real(dble)   expfac      !exponential factor
      logical      coulomb_on  !turns on/off coulomb interaction
      logical      nuclear_on  !turns on/off nuclear interaction
      logical      fzero       !set .true. unless i and j interact
      real(dble)   xbbcc       ! = bb+cc for n-n, p-p interactions,
                               ! = bb-cc for n-p interaction

      call starttimer()   !DKB-perf (accel)
      call starttimer()   !DKB-perf (calc_a)

      coulomb_on = coulomb.eq.'screened-coulomb'
      nuclear_on = nuclear.eq.'HPP'

      allocate(fj(3,0:n+2))   !locations n, n+1, n+2 are used in MPI_allreduce.
      fj=0.0
      xpacketi = 1./xpacket
      xpacket2i = 1./(2.*xpacket)

      if( rccut .le. sqrt(xl(1)**2+xl(2)**2+xl(3)**2)/2.0 ) then
        evcut=vc*exp(-xmuc*rccut)/rccut
      else
        evcut=0.0d0
      endif

      evx=0.0d0
      virx=0.0d0
      fi(:)=0.0d0
      pxx=0.0; pxy=0.0d0; pxz=0.0d0; pyy=0.0d0; pyz=0.0d0; pzz=0.0d0

!-------------------------------------------------------------------------------
! If measurements were called for, use this longer calculation, which  computes
! potential energy, virial and pressure tensor in addition to all the forces.
      if(do_measurements) then

         !$omp parallel private(xx,r,expfac,xbbcc,fn,fc,fncx,fzero)
         do 100 i=myrank,n-2,nprocs
            !$omp do reduction(+:evx,virx,fi,pxx,pxy,pxz,pyy,pyz,pzz)
            do 90 j=i+1,n-1
               fn=0.0
               fc=0.0
               fzero=.true.
               r=0.d0
               call pvec(i,j,r,xx)
               if(nuclear_on .and. r.le.rncut) then
                  expfac=exp(-r*r*xpacket2i)
                  if( zii(i).eq.zii(j) ) then !n-n, p-p
                     xbbcc = bb+cc
                  else          !n-p
                     xbbcc = bb-cc
                  endif
                  evx = evx + (aa*expfac+xbbcc)*expfac
                  fn = xpacketi*(2.*aa*expfac+xbbcc)*expfac
                  virx = virx + fn*r*r
                  fzero=.false.
               endif
               if(coulomb_on .and. r.le.rccut) then
                  if((zii(i).eq.1.d0).and.(zii(j).eq.1.d0)) then !p-p
                     fc=vc*exp(-xmuc*r)/(r*r*r)
                     evx = evx + (fc*r*r-evcut)
                     fc = (1.+xmuc*r)*fc
                     virx = virx + fc*r*r
                     fzero=.false.
                  endif
               endif
               if(.not.fzero) then
                  fncx(:)=(fn+fc)*xx(:)
                  fi(:)   = fi(:)   + fncx(:) !action of j on i
                  fj(:,j) = fj(:,j) - fncx(:) !reaction of i on j
                  pxx = pxx + xx(1)*fncx(1)
                  pxy = pxy + xx(1)*fncx(2)
                  pxz = pxz + xx(1)*fncx(3)
                  pyy = pyy + xx(2)*fncx(2)
                  pyz = pyz + xx(2)*fncx(3)
                  pzz = pzz + xx(3)*fncx(3)
               endif
 90         continue
            !$omp end do
            !$omp single
            fj(:,i)=fj(:,i)+fi(:)
            fi(:)=0.0d0
            !$omp end single
 100     continue
         !$omp end parallel

!-------------------------------------------------------------------------------
!     If measurements were not called for, use this quicker calculation of just
!     forces.
      else

         !$omp parallel private(xx,r,expfac,xbbcc,fn,fc,fncx,fzero)
         do 200 i=myrank,n-2,nprocs
            !$omp do reduction(+:fi)
            do 190 j=i+1,n-1
               fn=0.0
               fc=0.0
               fzero=.true.
               r=0.d0
               call pvec(i,j,r,xx)
               if(nuclear_on .and. r.le.rncut) then
                  expfac=exp(-r*r*xpacket2i)
                  if( zii(i).eq.zii(j) ) then !n-n, p-p
                     xbbcc = bb+cc
                  else          !n-p
                     xbbcc = bb-cc
                  endif
                  fn = xpacketi*(2.*aa*expfac+xbbcc)*expfac
                  fzero=.false.
               endif
               if(coulomb_on .and. r.le.rccut) then
                  if((zii(i).eq.1.d0).and.(zii(j).eq.1.d0)) then !p-p
                     fc = vc*(1.+xmuc*r)*exp(-xmuc*r)/(r*r*r)
                     fzero=.false.
                  endif
               endif
               if(.not.fzero) then
                  fncx(:)=(fn+fc)*xx(:)
                  fi(:)   = fi(:)   + fncx(:) !action of j on i
                  fj(:,j) = fj(:,j) - fncx(:) !reaction of i on j
               endif
 190        continue
            !$omp end do
            !$omp single
            fj(:,i) = fj(:,i)+fi(:)
            fi(:)=0.0d0
            !$omp end single
 200     continue
         !$omp end parallel

!-------------------------------------------------------------------------------
      endif
      call stoptimer(1,t_calc_a,ts_calc_a,n_calc_a) !DKB-perf (calc_a)


! Pack evx, virx and the pressure tensor into locations n, n+1 and n+2 of fj,
! so we need to do only one MPI_allreduce. Note that evx, virx and the pressure
! tensor will all be zero if only the force calculation was done.
      fj(1,n) = evx
      fj(2,n) = virx
      fj(3,n) = 0.0d0       !not used
      fj(1,n+1) = pxx
      fj(2,n+1) = pxy
      fj(3,n+1) = pxz
      fj(1,n+2) = pyy
      fj(2,n+2) = pyz
      fj(3,n+2) = pzz
      call MPI_allreduce(fj,a,3*(n+3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$omp parallel do
      do i=0,n-1 
         a(:,i) = a(:,i)/xmass
      enddo   
      !$omp end parallel do
      if(do_measurements) then
        ev      = a(1,n)
        vir     = a(2,n)
        pp(1,1) = a(1,n+1)
        pp(1,2) = a(2,n+1)
        pp(1,3) = a(3,n+1)
        pp(2,2) = a(1,n+2)
        pp(2,3) = a(2,n+2)
        pp(3,3) = a(3,n+2)
        ev = ev/float(n)           !average potential energy per nucleon
        px = rho*(kT+vir/(3.0*n))  !(instantaneous) pressure
      else
        ev  = 0.0d0
        vir = 0.0d0
        px  = 0.0d0
        pp  = 0.0d0
      endif

      deallocate(fj)
      call stoptimer(1,t_accel,ts_accel,n_accel)     !DKB-perf (accel)

      return
      end subroutine accel_nn
