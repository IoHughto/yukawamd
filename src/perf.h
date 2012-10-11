!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2011, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2011-Nov-28
! ---------------------------------------------------------------------
!
! Timers and counters for doing performance analysis
!
!*******************************************************************************


      real*8   t_md,        ts_md,       &
               t_newton,    ts_newton,   &
               t_accel,     ts_accel,    &
               t_calc_a,    ts_calc_a,   &
               t_vtot,      ts_vtot,     &
               t_calc_v,    ts_calc_v,   &
               t_g,         ts_g,        &
               t_pressure,  ts_pressure, &
               t_calc_vir,  ts_calc_vir, &
               t_strain,    ts_strain

      common /timers/    &
         t_md,        ts_md,           &
         t_newton,    ts_newton,       &
         t_accel,     ts_accel,        &
         t_calc_a,    ts_calc_a,       &
         t_vtot,      ts_vtot,         &
         t_calc_v,    ts_calc_v,       &
         t_g,         ts_g,            &
         t_pressure,  ts_pressure,     &
         t_calc_vir,  ts_calc_vir,     &
         t_strain,    ts_strain

      integer  n_md,        &
               n_newton,    &
               n_accel,     &
               n_calc_a,    &
               n_vtot,      &
               n_calc_v,    &
               n_g,         &
               n_pressure,  &
               n_calc_vir,  &
               n_strain

      common /counters/  &
         n_md,                         &
         n_newton,                     &
         n_accel,                      &
         n_calc_a,                     &
         n_vtot,                       &
         n_calc_v,                     &
         n_g,                          &
         n_pressure,                   &
         n_calc_vir,                   &
         n_strain
