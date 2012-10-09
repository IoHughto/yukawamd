!*******************************************************************************
!    MD 6.3.0
! ----------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-09
! ----------------------------------------------------------------------
!
! Read the parameter file for the simulation. Use input parameters to calculate
! further simulation parameters. Write all parameters to standard out, and to a
! run log file md.xxxxxxxxxx.log, where "xxxxxxxxxxx" is an 11-digit integer
! specifying the end time of the simulation. It is calculated from the start
! time, time step, nwgroup, nwsteps, nind, ntot, ngroup.
!
!*******************************************************************************


      subroutine input1
      use  md_types
      use  md_globals
      use  md_comm
      implicit real*8 (a-h,o-z)
      include  'mpif.h'

! All the input parameters are specified in the following namelist, which will
! be used to read them from the parameter file. All of these parameters are
! defined in module md_globals.
      namelist  /run_parameters/                         &
                 sim_type,                               &
                 start,                                  &
                 tstart, dt,                             &
                 nwgroup, nwsteps,                       &
                 ngroup, ntot, nind,                     &
                 nptensor,                               &
                 ncom,                                   &
                 ntnorm,                                 &
                 ngofr,                                  &
                !gspec,                                  &
                 nbin, nsbin, qmin, dq,                  &
                 nckpt,                                  &
                 nout, ftype, append,                    &
                 spec_file, spec_list,                   &
                 rho, aspect,                            &
                 kT,                                     &
                 deps,                                   &
                 rccut,                                  &
                 xmass,                                  &
                 coulomb, xlambda, frp,                  &
                 bfield,                                 &
                 efield,q0,w0,tref,                      &
                 nuclear, xpacket, aa, bb, cc,           &
                 rncut,                                  &
                 irnd, iseed, rmax

      integer  inyu(4)
      integer  omp_get_num_threads
      integer  omp_get_max_threads
      character*8     date
      character*10    daytime
      character*5     timezone
      character*256   line
      character*256   string


! MPI process (or rank) 0 starts at label 10. It reads the run parameter file
! and writes info to stdout and the log file. All other MPI processes start at
! label 1010 and wait for process 0 to broadcast the parameters.

   10 if(myrank.eq.0) then
!===============================================================================
! Read the run parameter file.
!DKB-todo: Had to set irnd here instead of in md_globals to avoid multiple def.
      irnd=0
!DKB-todo: handle file-file-not-found error
      open(10,FILE=runmdin,FORM='FORMATTED',STATUS='OLD')
      read(10,NML=run_parameters)
      close(10)

      open(7,FILE='md.aaa.log', FORM='FORMATTED',STATUS='UNKNOWN')

!-------------------------------------------------------------------------------
! Write program version, and what facilities we are using.
      write(6,30,advance='no') trim(code_name)//'_'//trim(code_version),mp_type
      write(7,30,advance='no') trim(code_name)//'_'//trim(code_version),mp_type
   30 format(2x,a,1x,a)
#ifdef _OPENMP
      write(6,34,advance='no')
      write(7,34,advance='no')
   34 format('+OpenMP')
#endif
      write(6,36)
      write(7,36)
   36 format()
      call date_and_time(date,daytime,timezone)
      write(6,37) date,daytime,timezone
      write(7,37) date,daytime,timezone
   37 format('  Start date/time = ',a8,2x,a10,2x,a5)

      write(6,40) nprocs
      write(7,40) nprocs
   40 format('  Number of processes:  ',i6)
#ifdef _OPENMP
     !nthreads = omp_get_num_threads()
      nthreads = omp_get_max_threads()
      write(6,44) nthreads
      write(7,44) nthreads
   44 format('  OpenMP threads/process: ',i4)
#endif

      write(6,48)
      write(7,48)
   48 format('  ----------------------------------------------------------')

!-------------------------------------------------------------------------------
! Write what kind of starting configuration.
!     write(6,50) trim(start)
!     write(7,50) trim(start)
!  50 format('  start = ',a)

!-------------------------------------------------------------------------------
! Type of simulation to run.
      write(6,90) trim(sim_type)
      write(7,90) trim(sim_type)
   90 format(/,'  SIMULATION TYPE:  ',a /  &
               '  ------------------------')

!-------------------------------------------------------------------------------
! There are two ways to specify the particle species list:
! (1) Number of particles and (Z,A) of each species may be listed in a "spec_list"
!     directly in the runmd.in file. These are read into the spec_list array by the
!     namelist facility.
! (2) The same data may be read from a "spec_file".
! Here we just partially handle the spec_list case. We write out the list of species
! to the log file, and count the number of species in the list and the total number
! of particles. Setting up the the zii and aii arrays will be done in initconfig.
      if(maxval(spec_list(:)%n).gt.0) then
         spec_file='internal'
         write(6,104)
         write(7,104)
  104    format('  Species  #particles    Z        A')
         k=1
         n=0
         do while (spec_list(k).ne.nullspecies)
            i=spec_list(k)%n
            if(i.gt.0) then
               write(6,105) k,i,spec_list(k)%z,spec_list(k)%a
               write(7,105) k,i,spec_list(k)%z,spec_list(k)%a
  105          format(3x,i4,4x,i9,2x,f7.2,2x,f7.2)
               n=n+i
            endif
            k=k+1
         enddo
         nspec=k-1
         write(6,106) n
         write(7,106) n
  106    format('   Total:  ',i9)
      endif

!-------------------------------------------------------------------------------
! Variable irnd identifies which random number generator to use. A random number
! generator is used only for generating an initial configuration, and to initial-
! ize velocities for a Maxwell distribution. Random numbers are not needed when
! starting from configuration files containing both positions and velocities.
! Default value is irnd=0, indicating no RNG is used.
      if(irnd.gt.0) then
        write(6,210) irnd,trim(ranname(irnd)),iseed
        write(7,210) irnd,trim(ranname(irnd)),iseed
  210   format(/'  irnd  = ',i9,2x,'(',a,')'  /  &
                '  iseed = ',i9)
      endif

!Set random number seed. Note iseed should be < 4096.
      do i=1,4
        inyu(i)=iseed/i
      end do
      call setrn(inyu)
      call sgrnd(iseed)

!-------------------------------------------------------------------------------
! Read charges and masses, and read or generate initial configuration.
      call initconfig

      rhop=(ztot*rho)/n

! If we are doing an ion simulation, and there is only one ion species, reset
! sim_type to 'pure-ion' so that the pure-ion force and potential routines will
! be used. These might be faster than the ion-mixture routines.
      if(sim_type.eq.'ion') then
         if(nspec.eq.1) then
           sim_type='pure-ion'
           ni=n
           zi=zii(1)
           ai=aii(1)
         else
           sim_type='ion-mixture'
           ni=0
           zi=0.d0
           ai=0.d0
         endif
      endif

!-------------------------------------------------------------------------------
!  Particle density, box edge aspect ratio, edge lengths, strain rates.
      write(6,124) rho
      write(7,124) rho
  124 format(/,'  rho          = ', es15.8, ' fm-3')
      call strainfactors
      write(6,125) xl0, aspect, xl, deps
      write(7,125) xl0, aspect, xl, deps
  125 format('  L0                       = ', 2x,f17.8, ' fm'      /    &
             '  edge aspect              = ', 6(6x,f13.8)          /    &
             '  Lx,Ly,Lz,Lyz,Lzx,Lxy     = ', 6(2x,f17.8), ' fm'   /    &
             '  strainrates              = ', 6(2x,es17.9), ' (c/fm)')

!  System temperature.
      write(6,130) kT
      write(7,130) kT
  130 format('  kT           = ', f10.6,   ' MeV')

!  Nucleon mass.
      write(6,110) xmass
      write(7,110) xmass
  110 format('  nucleon mass = ', f8.3, ' MeV')

!-------------------------------------------------------------------------------
! Two-particle correlation function, and static structure function.
!DKB-todo : Calculation of g(r) and S(q) disabled in this version until we figure
!DKB-todo+:    out how to do species-specific calculation.
      if(ngofr.gt.0) then
         ngofr = 0
         write(6,142)
         write(7,142)
  142    format('  *** g(r) and S(q) disabled in this version ***')
!        write(6,140) gspec
!        write(7,140) gspec
! 140    format('  Will compute g(r) and S(q) for species (N,Z)=(',i3,',',i3,')')
!        write(6,141) nbin,nsbin,qmin,dq
!        write(7,141) nbin,nsbin,qmin,dq
! 141    format('      Nbin      Nsbin      Qmin     dQ (fm-1)',/,            &
!               4x,i6,5x,i6,2x,1pe11.4,2x,1pe11.4)
      endif

!-------------------------------------------------------------------------------
      write(6,150)
      write(7,150)
  150 format(/ '  INTERACTIONS' / '  ------------')
!-------------------------------------------------------------------------------

! Coulomb interactions
! screened-coulomb = coulomb interaction with exponential screening
      if(coulomb.eq.'screened-coulomb') then
         write(6,151) trim(coulomb)
         write(7,151) trim(coulomb)
  151    format('  ',a)
         !If a screening length has been specified, use that. Otherwise
         !calculate Thomas-Fermi screening length.
         if(xlambda.gt.0) then
            xmuc=1.d0/xlambda
            string= '(defined)'
         else
            rhoe=rhop
            xkf=(3.d0*(PI**2)*rhoe)**(1.d0/3.d0)*HBARC
            ef=sqrt(xkf**2+ME**2)
           !ef=xkf                   !could be used in ultra-relativistic limit
            qtf=2.d0*sqrt((ALPHA/PI)*xkf*ef)
            xmuc=qtf/HBARC
            xlambda=1.d0/xmuc
            string= '(Thomas-Fermi)'
         endif
         write(6,152) xlambda,trim(string),frp,rccut
         write(7,152) xlambda,trim(string),frp,rccut
  152    format('     screening length = ', f13.8, ' fm  ', a,  / &
                '     ion form factor  = ', f13.8,              / &
                '     rccut            = ', f13.8, ' fm')
      endif

!-------------------------------------------------------------------------------
! Nuclear interactions.
! HPP = Horowitz-Perez-Piekarewicz interaction, Phys.Rev.C69, 045804 (2004)
      if(nuclear.eq.'HPP') then
        !If the user has left rncut at its default value of zero, this indicates
        !we should calculate it. We choose it so that the force at rncut will
        !be 2**-53 of Fmin=0.0001. The value Fmin should only rarely be achieved
        !in any 2-particle interaction. Since in double precision arithmetic
        !1+(2**-53)=1, then Fmin times this value will only rarely register when
        !summing up two-body forces.
         if(rncut.eq.0.0d0) then
          !xepsilon = (sqrt(2.d0*xpacket)/(4.d0*aa))*0.010d0*(2.d0**-53)
           xepsilon = (sqrt(2.d0*xpacket)/(2.d0*(-bb+cc)))*0.0001d0*(2.d0**(-53))
           rncut = sqrt(-log((-bb-cc)/(2.d0*aa)))
           xdeltar=1.d0
           do
             do while(rncut*exp(-rncut**2).gt.xepsilon)
               rncut0=rncut
               rncut=rncut+xdeltar
             enddo
             if(xdeltar.lt.0.11d0) exit
             xdeltar=xdeltar/10.d0
             rncut=rncut0
           enddo
           rncut = rncut*sqrt(2.d0*xpacket)
         endif
         write(6,160) trim(nuclear)
         write(7,160) trim(nuclear)
  160    format('  ',a)
         write(6,162) aa,bb,cc,xpacket,rncut
         write(7,162) aa,bb,cc,xpacket,rncut
  162    format('     aa       = ', f8.3, ' MeV' /  &
                '     bb       = ', f8.3, ' MeV' /  &
                '     cc       = ', f8.3, ' MeV' /  &
                '     Lambda   = ', f8.3, ' fm^2' / &
                '     rncut    = ', f8.3, ' fm')
      endif

!-------------------------------------------------------------------------------
! External uniform B field.
      if(bfield.gt.0.0d0) then
         write(6,170) bfield
         write(7,170) bfield
  170    format('  B-field:' /                         &
                '     B0   = ', 1pe10.3, ' Gauss')
      endif

!-------------------------------------------------------------------------------
! External oscillating E field.
      if(efield.gt.0.0d0) then
         efield=0.0d0
         write(6,180) efield
         write(7,180) efield
  180    format('  *** external E field disabled in this version ***')
! 180    format('  E-field:'                             /  &
!               '     eE0(MeV/fm) = ', es13.4, ' MeV/fm' /  &
!               '     q0(fm-1)    = ', es13.4, ' fm-1 '  /  &
!               '     w0(fm-1)    = ', es13.4, ' fm-1'   /  &
!               '     tref (fm)   = ', f13.3,  ' fm' )
      endif

!-------------------------------------------------------------------------------
! Write out timestep numbers:
!   - warmup groups, steps per warmup group
!   - measurement groups, measurements per group, steps between measurements
!   - steps between temperature normalizations during measurement phase
!   - steps between center-of-mass motion cancellations during meas. phase
!   - steps between pressure tensor measurements
      write(6,60) nwgroup,nwsteps,ngroup,ntot,nind,ntnorm,ncom,nptensor
      write(7,60) nwgroup,nwsteps,ngroup,ntot,nind,ntnorm,ncom,nptensor
   60 format(/,'  TIME STEP INTERVALS' /                     &
               '  -------------------' /                     &
               '  nwgroup, nwsteps:      ', i8,1x,i8 /       &
               '  ngroup, ntot, nind:    ', i8,1x,i8,1x,i8 / &
               '  ntnorm:                ', i8 /             &
               '  ncom:                  ', i8 /             &
               '  nptensor:              ', i8)
              !'  ngofr:                 ', i8)

! Allocate arrays for measurment statistics.
      allocate(eva(ngroup),ev2a(ngroup),dev(ngroup))
      allocate(pa(ngroup),p2a(ngroup),dp(ngroup))

! Info about configuration and checkpoint files.
      write(6,70) nout, ftype, append, nckpt
      write(7,70) nout, ftype, append, nckpt
   70 format('  nout, ftype, append:   ', i8,7x,a5,5x,l1 / &
             '  nckpt:                 ', i8)

!-------------------------------------------------------------------------------

! MPI process 0 broadcasts all parameters to other processes. Then it broadcasts
! charge, mass and type arrays, and coordinate and velocity arrays.
      call bcast_parms
      call MPI_bcast(zii, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(aii, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(x, 3*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(v, 3*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)



!===============================================================================
! All MPI processes other than 0 come here to participate in the broadcasts.
! Subroutine bcast_parms handles the broadcast of run parameters.

 1010 else

      call bcast_parms
      allocate(zii(0:n-1))       !charges
      allocate(aii(0:n-1))       !masses
      allocate(x(3,0:n-1))       !positions
      allocate(v(3,0:n-1))       !velocities
      allocate(vold(3,0:n-1))    !old velocities
      allocate(a(3,0:n+2))       !accelerations. Indices n,n+1,n+2 are used
                                 !   by MPI_allreduce in force routine
      allocate(eva(ngroup),ev2a(ngroup),dev(ngroup))
      allocate(pa(ngroup),p2a(ngroup),dp(ngroup))

      call MPI_bcast(zii, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(aii, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(x, 3*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_bcast(v, 3*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

!===============================================================================
      endif

#ifdef ATP
      call ATP_input1
#endif

      return
      end subroutine input1
