!*******************************************************************************
!     MD 6.2.0
!  -----------------------------------------------------------------------------
!     Copyright 2011, The Trustees of Indiana University
!     Original author:   Don Berry
!     Last modified by:  Don Berry, 2011-Dec-29
!  -----------------------------------------------------------------------------
!
!  Subroutines in this file control which potential energy and force subroutines
!  to call. The choice is made by global variable sim_type, which defines the
!  type of simulation being run.
!
!*******************************************************************************


      subroutine vtot(evavg,pressure)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)

      real(dble)   evavg
      real(dble)   pressure

      select case(trim(sim_type))

      case('nucleon')
        call vtot_nn(evavg)

      case('pure-ion')
        call vtot_ion_pure(evavg)

      case('ion-mixture')
        call vtot_ion_mix(evavg,pressure)

      case('B-field')  !DKB-todo: This was for test only.
        evavg=0.0d0
        pressure=0.0d0

      case default
        write(6,"('sim_type=***',a,'***')") trim(sim_type)
        write(8,"('sim_type=***',a,'***')") trim(sim_type)
        call md_exit(1,'*** Program not compiled to do this simulation type ***')

      end select

      end subroutine vtot



!*******************************************************************************
      subroutine accel(do_measurements)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)

      logical   do_measurements

      select case(trim(sim_type))

      case('nucleon')
        call accel_nn(do_measurements)

      case('pure-ion')
        call accel_ion_pure(do_measurements)

      case('ion-mixture')
        call accel_ion_mix(do_measurements)

      case('B-field')  !DKB-todo: This was for test only.
        a=0.0d0

      case default
        write(6,"('sim_type=***',a,'***')") trim(sim_type)
        write(8,"('sim_type=***',a,'***')") trim(sim_type)
        call md_exit(1,'*** Program not compiled to do this simulation type ***')

      end select

      end subroutine accel
