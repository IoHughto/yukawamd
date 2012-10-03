!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Author:            Don Berry
!    Last modified by:  Don Berry, 2012-May-23
! ---------------------------------------------------------------------
!
!*******************************************************************************


      module md_globals

      use md_types

      save

      character*10, parameter :: code_name='MD'
      character*8, parameter  :: code_version='6.2.0'

      integer, parameter    :: MAXSPEC=1000       !max number of species allowed
      real(dble), parameter :: XNOSTRAIN=1000.d0  !indicates no strain

      real(dble), parameter :: HBARC=197.327d0           !hbar*c
      real(dble), parameter :: ALPHA=1.d0/137.036d0      !fine structure constant
      real(dble), parameter :: ALPHAI=137.036d0
      real(dble), parameter :: ME=0.510999d0             !electron mass (MeV)
      real(dble), parameter :: PI=3.14159265358979324d0

!===============================================================================
!                     PARAMETERS DEFINING THE SIMULATION
!
! Many of these are user-definable in the runmd.in file. Others are calculated
! from the runmd.in parameters. Some runmd.in parameters have default values.

      character*256 :: runmdin='runmd.in' !run parameter file
      character*256 :: start='md.in.xv82' !initial configuration file, for restarts
      character*40  :: suffix=''  !for making file names.
      real(dble)    :: tstart=-1. !simulation start time (fm/c)
      real(dble)    :: time       !simulation clock (fm/c)
      real(dble)    :: dt=0.0     !simulation time step (fm/c)
      real(dble)    :: tend       !simulation end time (fm/c)
      integer       :: iseed=0    !seed for random number generator
      integer       :: irnd       !selects which random number generator to use

! For ion simulations.
      integer       :: nspec=0      !number of particle species
      integer       :: ni=0         !number of ions for pure-ion simulations
      real(dble)    :: zi=0., ai=0. !charge and mass numbers for pure-ion simulations
      character*256 :: spec_file='' !ion data file ("species file")
      type(species),save :: spec_list(MAXSPEC) !ion data list ("species list")

! For all simulations.
      integer       :: n=0          !total number of particles
      real(dble)    :: ztot=0.0     !total charge of system
      real(dble)    :: rho=0.0      !particle density (particles/fm^3)
      real(dble)    :: xl0=0.0      !fundamental edge length
      real(dble)    :: aspect(3)=(/1.0,1.0,1.0/)     !simulation box edge aspect ratio
      real(dble)    :: xl(3)=(/0.0,0.0,0.0/)         !simulation box edge lengths (fm)
      real(dble)    :: deps(3)=(/XNOSTRAIN,XNOSTRAIN,XNOSTRAIN/) !xx, yy and zz strain rates
      real(dble)    :: strnfac(3)=(/1.0,1.0,1.0/)    !strain factors in each dimension
      real(dble)    :: kT=0.0       !temperature (MeV)
      real(dble)    :: xmass=931.00 !nucleon mass (MeV)
      real(dble)    :: rmax=0.0     !radius of nucleus, for doing large-nucleus sims.

      type(species), parameter  :: proton  = species(0,1.,1.)
      type(species), parameter  :: neutron = species(0,0.,1.)
      type(species), parameter  :: nullspecies = species(0,0.,0.)
      type(species), parameter  :: allspecies = species(0,-1.,-1.)

! Time step numbers controlling behavior of the program. Unless otherwise stated,
! a value of 0 means never do the indicated action.
      integer       :: istep=0     !current step number
      integer       :: nwgroup=0   !number of groups of warmup steps
      integer       :: nwsteps=0   !warmup steps per group
      integer       :: ngroup=0    !number of measurement groups
      integer       :: ntot=0      !measurements per group
      integer       :: nind=0      !steps between measurements
      integer       :: nptensor=0  !steps between pressure tensor measurements
      integer       :: ncom=0      !steps between center-of-mass velocity cancellations
      integer       :: ntnorm=0    !steps between temperature normalizations
      integer       :: ngofr=0     !steps between g(r) calculations
      integer       :: nckpt=0     !steps between checkpoints. A checkpoint file
                                   !   is always output at the end of a run.
      integer       :: nout=0      !steps between configuration ouputs

! File type for output configuration files:
      character*6   :: ftype = 'x4b'
                         !Revision B formats are recommended:
                         !x4b  = positions only, real*4
                         !xv4b = positions and velocities, real*4
                         !x8b  = positions only, real*8
                         !xv8b = positions and velocities, real*8
      logical       :: append = .true.
                         !.true.  = append all configs to a single md.traj file
                         !.false. = configs go to individial md.out files

!-------------------------------------------------------------------------------
! Parameters for two-particle correlation function, and static structure factor.
      type(species),save :: gspec    !particle species to calculate g(r) for
      integer            :: nbin     !number of bins for g(r)
      integer            :: nsbin    !number of bins for S(q)
      real(dble)         :: qmin,dq  !minimum q and bin width for S(q)


!-------------------------------------------------------------------------------
! Random number generator stuff.
      character*6  ranname(0:4)
      common  /rndtype/irnd,ranname


!-------------------------------------------------------------------------------
! Type of simulation to be performed.
      character*20 ::  sim_type = ''

! Parameters defining the Coulomb interaction.
      character*20  ::  coulomb = ''     !type of coulomb interaction
      real(dble)    ::  xlambda = -1.0d0 !screening length
      real(dble)    ::  rccut=1.0d50     !cutoff radius
      real(dble)    ::  xmuc = 0.0d0     !inverse screening length
      real(dble)    ::  vc=ALPHA*HBARC   !Coulomb coupling constant
      real(dble)    ::  frp = 1.0d0      !ion form factor, f(Rp/xlambda)

! Parameters defining the nuclear interaction.
      character*20 ::  nuclear = ''      !type of nuclear interaction
      real(dble)   ::  xpacket = 1.25
      real(dble)   ::  aa = 110.00
      real(dble)   ::  bb = -26.00
      real(dble)   ::  cc =  24.00
      real(dble)   ::  rncut=0.0d0       !cutoff radius

! Parameters defining external electric and magnetic fields.
! Note that this version does not do external E fields, but we leave the
! associated parameters (efield, q0, w0, tref) defined here for future use.
      real(dble)   ::  bfield = 0.0d0  !uniform B field in z direction (Gauss)
      real(dble)   ::  efield = 0.0d0  !amplitude of oscillating E field (Mev/fm)
      real(dble)   ::  q0              !wave number of E field (1/fm)
      real(dble)   ::  w0              !frequency of E field (1/fm)
      real(dble)   ::  tref            !time reference for E field (fm/c)

!-------------------------------------------------------------------------------
! MPI variables.
      integer      ::  myrank = 0      !MPI rank, or process, number
      integer      ::  nprocs = 1      !number of MPI processes



!===============================================================================
!                            PARTICLE DATA
!
! Position, velocitie, acceleration, particle type, charge, mass.
      real(dble), allocatable, target ::   x(:,:)
      real(dble), allocatable         ::   v(:,:),vold(:,:)
      real(dble), allocatable         ::   a(:,:)
!DKB-todo: Do we still need the ctype array?
      character*6,allocatable         ::   ctype(:)  !character code for particle type
      real(dble), allocatable         ::   zii(:)    !Z for each ion
      real(dble), allocatable         ::   aii(:)    !A for each ion

! Energies: ek=kinetic, ev=potential, e=total
      real(dble)   ek, ev, e            !mean, per particle, for current configuration
      real(dble),allocatable :: eva(:)  !mean ev per measurement group
      real(dble),allocatable :: ev2a(:) !2nd moment of ev, per measurement group
      real(dble),allocatable :: dev(:)  !std.dev. of ev, per measurement group

! Pressure:
      real(dble)   px                   !pressure
      real(dble),allocatable :: pa(:)   !mean per measurement group
      real(dble),allocatable :: p2a(:)  !2nd moment, per measurement group
      real(dble),allocatable :: dp(:)   !std.dev., per measurement group
      real(dble)   ::  pp(3,3)          !pressure tensor

!-------------------------------------------------------------------------------
! Arrays needed for computing two-particle correlation function.
      real(dble), allocatable    :: gg(:,:,:)
      real(dble), allocatable    :: cgg(:)

!-------------------------------------------------------------------------------
! Parameters and arrays needed for computing the static structure factor.
      real(dble), allocatable    :: ss(:,:)
      real(dble), allocatable    :: cs(:)


      end module md_globals
