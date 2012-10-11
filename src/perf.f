!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-May-10
! ---------------------------------------------------------------------
!
! Performance measurement routines.
!DKB-todo: g(r) disabled in this version
!DKB-todo: Separate calculations for virial and pressure removed in MD_6.2.0
!
!*******************************************************************************

!*******************************************************************************
!  Initialize performance monitor counters and timers.
!

      subroutine perf_init
      implicit none
      include  'perf.h'

      double precision   t0

      call inittime()

      t_md       = 0.0;   ts_md       = 0.0
      t_newton   = 0.0;   ts_newton   = 0.0
      t_accel    = 0.0;   ts_accel    = 0.0
      t_calc_a   = 0.0;   ts_calc_a   = 0.0
     !t_vtot     = 0.0;   ts_vtot     = 0.0
     !t_calc_v   = 0.0;   ts_calc_v   = 0.0
     !t_g        = 0.0;   ts_g        = 0.0
     !t_pressure = 0.0;   ts_pressure = 0.0
     !t_calc_vir = 0.0;   ts_calc_vir = 0.0
     !t_strain   = 0.0;   ts_strain   = 0.0

      n_md       = 0
      n_newton   = 0
      n_accel    = 0
      n_calc_a   = 0
     !n_vtot     = 0
     !n_calc_v   = 0
     !n_g        = 0
     !n_pressure = 0
     !n_calc_vir = 0
     !n_strain   = 0

      return
      end subroutine perf_init




!*******************************************************************************
!  Output a performance report.
!

      subroutine perf_report
      use  md_comm
      implicit none
      include  'perf.h'
      include  'mpif.h'

      integer             m
      double precision    tavg    !average time per call
      double precision    tself   !self time
      double precision    tsavg   !average self time per call
      integer             ierror
      integer             iproc
      integer             status(MPI_STATUS_SIZE)
      integer             ibuff(10)
      double precision    xbuff(20)

      call MPI_comm_size(MPI_COMM_WORLD,nprocs,ierror)
      call MPI_comm_rank(MPI_COMM_WORLD,myrank,ierror)

      ! Pack up all performance data for shipping to MPI process 0.
      xbuff(1)  = t_md;         xbuff(2)  = ts_md;         ibuff(1)  = n_md
      xbuff(3)  = t_newton;     xbuff(4)  = ts_newton;     ibuff(2)  = n_newton
      xbuff(5)  = t_accel;      xbuff(6)  = ts_accel;      ibuff(3)  = n_accel
      xbuff(7)  = t_calc_a;     xbuff(8)  = ts_calc_a;     ibuff(4)  = n_calc_a
     !xbuff(9)  = t_vtot;       xbuff(10) = ts_vtot;       ibuff(5)  = n_vtot
     !xbuff(11) = t_calc_v;     xbuff(12) = ts_calc_v;     ibuff(6)  = n_calc_v
     !xbuff(13) = t_g;          xbuff(14) = ts_g;          ibuff(7)  = n_g
     !xbuff(15) = t_pressure;   xbuff(16) = ts_pressure;   ibuff(8)  = n_pressure
     !xbuff(17) = t_calc_vir;   xbuff(18) = ts_calc_vir;   ibuff(9)  = n_calc_vir
     !xbuff(19) = t_strain;     xbuff(20) = ts_strain;     ibuff(10) = n_strain


!===============================================================================
      select case(myrank)

      !-------------------------------------------------------------------------
      case(0)
         do 100 iproc = 0,nprocs-1
            t_md       = xbuff(1);   ts_md       = xbuff(2);   n_md       = ibuff(1)
            t_newton   = xbuff(3);   ts_newton   = xbuff(4);   n_newton   = ibuff(2)
            t_accel    = xbuff(5);   ts_accel    = xbuff(6);   n_accel    = ibuff(3)
            t_calc_a   = xbuff(7);   ts_calc_a   = xbuff(8);   n_calc_a   = ibuff(4)
           !t_vtot     = xbuff(9);   ts_vtot     = xbuff(10);  n_vtot     = ibuff(5)
           !t_calc_v   = xbuff(11);  ts_calc_v   = xbuff(12);  n_calc_v   = ibuff(6)
           !t_g        = xbuff(13);  ts_g        = xbuff(14);  n_g        = ibuff(7)
           !t_pressure = xbuff(15);  ts_pressure = xbuff(16);  n_pressure = ibuff(8)
           !t_calc_vir = xbuff(17);  ts_calc_vir = xbuff(18);  n_calc_vir = ibuff(9)
           !t_strain   = xbuff(19);  ts_strain   = xbuff(20);  n_strain   = ibuff(10)
            do m=6,8,2

               write(m,1000) iproc

               tavg  = t_md/max(n_md,1)
               tsavg = ts_md/max(n_md,1)
               write(m,1010) 'md:       ', t_md, n_md, tavg, ts_md, tsavg

               tavg  = t_newton/max(n_newton,1)
               tsavg = ts_newton/max(n_newton,1)
               write(m,1010) 'newton:   ', t_newton, n_newton, tavg, ts_newton, tsavg

               tavg  = t_accel/max(n_accel,1)
               tsavg = ts_accel/max(n_accel,1)
               write(m,1010) 'accel:    ', t_accel, n_accel, tavg, ts_accel, tsavg

               tavg  = t_calc_a/max(n_calc_a,1)
               tsavg = ts_calc_a/max(n_calc_a,1)
               write(m,1010) 'calc_a:   ', t_calc_a, n_calc_a, tavg, ts_calc_a, tsavg

              !tavg  = t_vtot/max(n_vtot,1)
              !tsavg = ts_vtot/max(n_vtot,1)
              !write(m,1010) 'vtot:     ', t_vtot, n_vtot, tavg, ts_vtot, tsavg

              !tavg  = t_calc_v/max(n_calc_v,1)
              !tsavg = ts_calc_v/max(n_calc_v,1)
              !write(m,1010) 'calc_v:   ', t_calc_v, n_calc_v, tavg, ts_calc_v, tsavg

              !tavg  = t_g/max(n_g,1)
              !tsavg = ts_g/max(n_g,1)
              !write(m,1010) 'g:        ', t_g, n_g, tavg, ts_g, tsavg

              !tavg  = t_pressure/max(n_pressure,1)
              !tsavg = ts_pressure/max(n_pressure,1)
              !write(m,1010) 'pressure: ', t_pressure, n_pressure, tavg, ts_pressure, tsavg

              !tavg  = t_calc_vir/max(n_calc_vir,1)
              !tsavg = ts_calc_vir/max(n_calc_vir,1)
              !write(m,1010) 'calc_vir: ', t_calc_vir, n_calc_vir, tavg, ts_calc_vir, tsavg

              !tavg  = t_strain/max(n_strain,1)
              !tsavg = ts_strain/max(n_strain,1)
              !write(m,1010) 'strain:   ', t_strain, n_strain, tavg, ts_strain, tsavg

            enddo

            if(iproc.lt.nprocs-1) then
               call MPI_recv(xbuff,20,MPI_DOUBLE_PRECISION,iproc+1,1,MPI_COMM_WORLD,status,ierror)
               call MPI_recv(ibuff,10,MPI_INTEGER,iproc+1,1,MPI_COMM_WORLD,status,ierror)
            endif

  100    continue

      !-------------------------------------------------------------------------
      case default
        !call MPI_recv(token,1,MPI_INTEGER,myrank-1,1,MPI_COMM_WORLD,status,ierror)
         call MPI_send(xbuff,20,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierror)
         call MPI_send(ibuff,10,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierror)

      end select
!===============================================================================



 1000 format(/,'MPI process ',i4,'   Performance statistics:',    &
      /,'               total                             self       self',   &
      /,'routine         time   ncalls       avg          time        avg')
 1010 format(a,f11.2,i8,f14.5,f11.2,f14.5)

      return
      end subroutine perf_report
