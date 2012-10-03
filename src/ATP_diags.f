!*******************************************************************************
!     MD 6.2.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Authors:           Don Berry
!     Last modified by:  Don Berry, 2012-May-23
!  ---------------------------------------------------------------------
!
!  ATP diagnostics to check correctness and consistency of run parameters
!  after call to subroutine input1.
!
!*******************************************************************************


      subroutine ATP_input1
      use  md_types
      use  md_globals
      use  md_comm
      implicit real*8 (a-h,o-z)
      include  'mpif.h'


      integer        atpn   !logical unit number for diagnostic output
      character*999  atpf   !file name for diagnotsic output

10010 format(a,' =  ***',a,'***')
10020 format(a,' =  ',i11)
10030 format(a,' =  ',f15.3)
10032 format(a,' =  ',es22.15)
10040 format(a,' =  ',l1)
10050 format(a,' =  ',i6,f8.3,f8.3)
10060 format('   spec(',i4,'): ',i6,f8.3,f8.3)


      call MPI_comm_rank(MPI_COMM_WORLD,myrank,ierror)
      atpn = 100+myrank
      write(atpf,"('ATP_diags.',i5.5)") atpn
      open(unit=atpn,file=atpf,form='formatted',status='unknown')


      write(atpn,10020) 'myrank      ', myrank

! Simulation time parameters:
      write(atpn,10030) 'tstart      ', tstart   !simulation start time (fm/c)
      write(atpn,10030) 'dt          ', dt       !simulation time step (fm/c)
      write(atpn,10030) 'time        ', time     !simulation clock (fm/c)
      write(atpn,10030) 'tend        ', tend     !simulation end time (fm/c)

! Timestep control parameters:
      write(atpn,10020) 'nwgroup     ', nwgroup   !warmup groups
      write(atpn,10020) 'nwsteps     ', nwsteps   !warmup steps per group
      write(atpn,10020) 'ngroup      ', ngroup    !measurement groups
      write(atpn,10020) 'ntot        ', ntot      !measurements per group
      write(atpn,10020) 'nind        ', nind      !steps between measurements
      write(atpn,10020) 'ncom        ', ncom      !steps between center-of-mass motion cancellation
      write(atpn,10020) 'ntnorm      ', ntnorm    !steps between temperature normalizations

!  Output control parameters:
      write(atpn,10020) 'nckpt       ', nckpt     !micro steps between checkpoints
      write(atpn,10020) 'nout        ', nout      !micro steps between config ouputs
      write(atpn,10010) 'ftype       ', ftype
      write(atpn,10040) 'append      ', append

! For pure-ion simulations:
      write(atpn,10020) 'ni          ', ni        !number of ions
      write(atpn,10030) 'zi          ', zi        !charge and mass numbers for pure-ion simulations
      write(atpn,10030) 'ai          ', ai        !charge and mass numbers for pure-ion simulations

! For ion-mixture simulations:
      write(atpn,10010) 'za_file     ', trim(za_file) !ion data file name

! Total number of particles:
      write(atpn,10020) 'n           ', n         !total number of particles

! Special values for type 'species':
      write(atpn,10050) 'proton      ', proton
      write(atpn,10050) 'neutron     ', neutron
      write(atpn,10050) 'nullspecies ', nullspecies
      write(atpn,10050) 'allspecies  ', allspecies

      write(atpn,10020) 'iseed       ', iseed     !seed for random number generator
      write(atpn,10020) 'irnd        ', irnd      !selects which random number generator to use

      write(atpn,10032) 'rho         ', rho       !particle density (particles/fm^3)
      write(atpn,10032) 'xl          ', xl        !edge length of simulation box (fm)
      write(atpn,10032) 't           ', t         !temperature (MeV)
      write(atpn,10032) 'xmass       ', xmass     !nucleon mass (MeV)
      write(atpn,10032) 'rmax        ', rmax      !radius of nucleus, for doing large-nucleus sims.

! List particle species:
      write(atpn,'(//,a)') 'List of spec array:'
      k=1
      do while(za_list(k)%n .gt. 0)
         write(atpn,10060) za_list(k)%n,za_list(k)%z,za_list(k)%a
         k=k+1
      enddo


      call flush(atpn)
      close(atpn)

      return


      end subroutine ATP_input1
