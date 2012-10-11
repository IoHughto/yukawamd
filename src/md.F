!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Jul-12
! ---------------------------------------------------------------------
!
!*******************************************************************************

      program md
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'perf.h'   !DKB-perf
      include 'mpif.h'

! Mean, 2nd moment, and std.dev. of potential energy and pressure.
      real(dble)   evave,ev2ave,devave
      real(dble)   pave,p2ave,dp2ave

      logical      do_measurements

      character*8     date
      character*10    daytime
      character*5     timezone
      character*256   line

! Initialize timers to monitor program performance.
      call perf_init      !DKB-perf
      call starttimer()   !DKB-perf (md)

! Initialize MPI message passing system. (Stubbed out in the non-MPI code.)
      call MPI_init(ierror)
      call MPI_comm_size(MPI_COMM_WORLD,nprocs,ierror)
      call MPI_comm_rank(MPI_COMM_WORLD,myrank,ierror)

! Read run parameters. Read charges and masses. Initalize positions and velocities.
      call input1

! Calculate simulation end time. Use it in naming the log and measurement files.
! Since simulation times may go into the 1e10 range, we need to be careful about
! integer overflow (for 32-bit integers). Hence the complicated write statement.
! Copy temporary log file to permanent one, and delete temporary.
      tend = tstart + dt*(nwgroup*nwsteps + ngroup*ntot*nind)
      if(myrank.eq.0) then
        write(suffix,22) nint(tend/1000000.d0), nint(mod(tend,1000000.d0))
   22   format(i5.5,i6.6)
        open(8,FILE='md.'//trim(suffix)//'.log', FORM='FORMATTED',STATUS='UNKNOWN')
        rewind(7)
        do
          read(7,'(a)',end=28) line
          write(8,'(a)') trim(line)
        enddo
   28   continue
        close(7,status='DELETE')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierror)

! Initalize array for correlation function g(r).
!DKB-todo: g(r) disabled in this version.
     !if(ngofr.gt.0) call g(0)

! Calculate potential energy, pressure, pressure tensor and accelerations for
! initial configuration. This is all done in the accel routine if its argument
! is set to .true. Also call ttot to calculate kinetic energy.
      call accel(.true.)
     !call dump_data    !DKB-debug
      call ttot(ek)
      e=ev+ek
      if(myrank.eq.0) then
        write(6,10050) '  Initial Vtot, Ttot, E, P = ',ev,ek,e,px
        write(8,10050) '  Initial Vtot, Ttot, E, P = ',ev,ek,e,px
10050   format(/,a,4(2x,es16.9))
      endif


! ----------------------------------------------------------------------
! Do warm up steps.
      if(nwgroup.ge.1 .and. myrank.eq.0) then
        write(6,10080)
        write(8,10080)
10080   format(/'  Warmup phase')
      endif
      istep=0
      do i=1,nwgroup
        do j=1,nwsteps
          call newton(.false.)
          istep=istep+1
          if(mdmod(istep,nckpt).eq.0) call checkpoint
        enddo
        call cancel_vcom
        call tnorm
        if(myrank.eq.0) then
          write(6,10082) istep,time
          write(8,10082) istep,time
10082     format(2x,i7,1x,f14.2,'   warm up')
        endif
      enddo
     !call accel(.true.)  !DKB-debug
     !call dump_data      !DKB-debug
      if(istep.gt.0 .and. mdmod(istep,nckpt).ne.0) call checkpoint

! ----------------------------------------------------------------------
! Do measurement steps. This is the main MD loop.
      if(ngroup.ge.1 .and. myrank.eq.0) then
        write(6,10090)
        write(8,10090)
10090   format(                     &
          '  Measurement phase',/,  &
          'e   step     time            Vtot              Ttot               E              pressure',/, &
          'p   step     time            P_xx              P_xy              P_xz              P_yy              P_yz              P_zz')
        write(suffix,22) nint(tend/1000000.d0), nint(mod(tend,1000000.d0))
        open(9,FILE='md.'//trim(suffix)//'.meas',FORM='UNFORMATTED',STATUS='UNKNOWN')
      endif
      evave=0.;  ev2ave=0.
      pave=0.;   p2ave=0.
      pp=0.
      istep=0
      do 100 ig=1,ngroup
        eva(ig)=0.;  ev2a(ig)=0.
        pa(ig)=0.;   p2a(ig)=0.
        do 40 i=1,ntot*nind
          if(mdmod(istep,ncom).eq.0) call cancel_vcom
          if(mdmod(istep,ntnorm).eq.0) call tnorm
          call strain
          do_measurements = mdmod(istep+1,nind).eq.0 .or. mdmod(istep+1,nptensor).eq.0
          call newton(do_measurements)
          istep=istep+1
          if(mdmod(istep,nind).eq.0) then
            eva(ig)  = eva(ig)+ev
            ev2a(ig) = ev2a(ig)+ev**2
            pa(ig)  = pa(ig)+px
            p2a(ig) = p2a(ig)+px**2
            e=ek+ev
            if(myrank.eq.0) then
              write(6,10101) 'e',istep,time,ev,ek,e,px
              write(8,10101) 'e',istep,time,ev,ek,e,px
            endif
          endif
          if(mdmod(istep,nptensor).eq.0) then
            if(myrank.eq.0) then
              write(6,10101) 'p',istep,time,pp(1,1),pp(1,2),pp(1,3),pp(2,2),pp(2,3),pp(3,3)
              write(8,10101) 'p',istep,time,pp(1,1),pp(1,2),pp(1,3),pp(2,2),pp(2,3),pp(3,3)
            endif
          endif
          if(do_measurements) then
            if(myrank.eq.0) then
              write(9) istep,time, ev,ek,px, pp(1,1),pp(1,2),pp(1,3),pp(2,2),pp(2,3),pp(3,3)
            endif
          endif
         !if(mdmod(istep,ngofr).eq.0) call g(ig)   !DKB-todo: g(r) disabled in this version
          if(mdmod(istep,nout).eq.0)  call save_config(ftype,append)
          if(mdmod(istep,nckpt).eq.0) call checkpoint
   40   continue
        call statt(eva(ig),ev2a(ig),ntot,dev(ig),evave,ev2ave)
        call statt(pa(ig),p2a(ig),ntot,dp(ig),pave,p2ave)
        if(myrank.eq.0) then
          write(6,10110) ig,time, eva(ig),dev(ig), pa(ig),dp(ig)
          write(8,10110) ig,time, eva(ig),dev(ig), pa(ig),dp(ig)
10110     format('s',i4,1x,f14.2,  2x,1pe16.9,2x,1pe16.9, 2x,16x,2x,1pe16.9,2x,1pe16.9)
        endif
  100 continue
10101 format(a1,1x,i7,1x,f14.2,6(2x,es16.9))
! -----------------------------------------------------------------------

      x1=0.
      x2=0. 
      call statt(evave,ev2ave,ngroup,devave,x1,x2)
      call statt(pave,p2ave,ngroup,dpave,x1,x2)

      if(myrank.eq.0) then
        write(6,10200) '  average potential = ', evave, ' +/-', devave
        write(8,10200) '  average potential = ', evave, ' +/-', devave
        write(6,10200) '  average pressure  = ', pave, ' +/-', dpave
        write(8,10200) '  average pressure  = ', pave, ' +/-', dpave
10200   format(2(a,1pe16.9))
      endif
      if(istep.gt.0 .and. mdmod(istep,nckpt).ne.0) call checkpoint
     !DKB-todo: g(r) disabled in this version
     !if(ngofr.gt.0) then
     !  call g(-1)
     !  call s
     !  call g(-2)
     !endif
      call accel(.true.)
      call ttot(ek)
      e=ev+ek
      if(myrank.eq.0) then
        stretch=xl/xl0
        aspect=stretch*aspect
        rho=float(n)/(xl(1)*xl(2)*xl(3))
        write(6,10206) xl, aspect, stretch, rho
        write(8,10206) xl, aspect, stretch, rho
10206   format('  Final dimensions  =',3(2x,f17.9), ' fm', /   &
               '  Final aspect      =',3(2x,f17.9),        /   & 
               '  Final stretch     =',3(2x,f17.9),        /   &
               '  Final density     =',1x,es15.8, ' fm-3' )
        write(6,10210) '  Final Vtot, Ttot, E, P   = ',ev,ek,e,px
        write(8,10210) '  Final Vtot, Ttot, E, P   = ',ev,ek,e,px
10210   format(a,4(2x,es16.9))
      endif

      call stoptimer(1,t_md,ts_md,n_md)  !DKB-perf (md)

      call date_and_time(date,daytime,timezone)
      if(myrank.eq.0) then
        write(6,10020) date,daytime,timezone
        write(8,10020) date,daytime,timezone
10020   format('  End date/time   = ',a8,2x,a10,2x,a5)
      endif
      call perf_report    !DKB-perf

      if(myrank.eq.0) then
        close(8)
        close(9)
      endif

      call md_exit(0,'normal exit')
      end program md



!*******************************************************************************
      subroutine statt(xa,xa2,n,dx,xave,x2ave)
      use  md_types
      implicit real(dble)(a-h,o-z)
      xa=xa/float(n)
      xa2=xa2/float(n)
      if(n.gt.1) then
        dx=sqrt(abs(xa2-xa*xa)/float(n-1))
      else
        dx=0.0d0
      endif
      xave=xave+xa
      x2ave=x2ave+xa*xa
      return
      end subroutine statt


      
!*******************************************************************************
!  This is similar to the fortran mod function, except that it returns -1 if
!  q=0. This allows the calling program to use q=0 as a special value.

      integer function mdmod(p,q)
      integer  p,q
      if(q.eq.0) then
        mdmod = -1
      else
        mdmod = mod(p,q)
      endif
      return
      end function mdmod
