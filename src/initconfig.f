!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Jul-13
! ---------------------------------------------------------------------
!
! Initconfig sets the charge, mass, and initial position and velocity for all
! particles. Charges and masses are read either from a list in the run parameter
! file ("spec_list"), or from a file ("spec_file"). The initial configuration
! is determined by the "start" parameter.  This parameter can be either a type
! of initial configuration to construct, or a file containing the configuration.
! The types of configurations initconfig can construct are,
!
!   'random'  -- particles coordinates are initialized for random start
!   'nuclear' -- particles coordinates are initialized in small spherical volume
!
! Initconfig can read a number of file types, identified by their extension.
!
! Old formatted file types used prior in versions prior to 6.2.
!    onp    -- positions. This is a format used in MD_0.
!    md1    -- positions and velocities. An old format used in MD_1.
!    xfa    -- positions, Rev.A. (Revision A)
!    xvfa   -- positions and velocities, Rev.A
!    ZAfa   -- charge, mass and number of particles of each species, Rev.A
!    xvZAfa -- positions, veolcities, charges masses, Rev.A
!
! Old Fortran unformatted file types, used in versions prior to 6.2:
!    md0    -- real*4 positions and velocities. An old format used in MD_0.
!    x4a    -- real*4 positions, Revision A ("Rev.A")
!    xv4a   -- real*4 positions and velocities, Rev.A
!    x8a    -- real*8 positions, Rev.A
!    xv8a   -- real*8 positions and velocities, Rev.A
!    xva8a  -- real*8 positions, velocities, accelerations, Rev.A
!
! Revision B Fortran unformatted file types:
!    x4b    -- real*4 positions. Rev.B
!    x8b    -- real*8 positions. Rev.B
!    xv4b   -- real*4 positions and velocities. Rev.B
!    xv8b   -- real*8 positions and velocities. Rev.B
!    xva8b  -- real*8 positions, velocities, accelerations. Rev.B
!
! Revision B Fortran formatted file types:
!    xfb    -- positions, Rev. B
!    xvfb   -- positions and velocities, Rev. B
!    xZAfb  -- positions, charges, masses. Rev.B
!    xvZAfb -- positions, velocities, charges, masses. Rev.B
!    ZAfb   -- charge, mass and number of particles of each species
!    xvZAfb -- same as ZAfb, but also includes positions and velocities
!
! For starts that do not define velocities, we generate them according to a Max-
! well distribution corresponding to the temperature kT specified in the run
! parameter file. Center of mass velocity is then calculated and subtracted from
! all particles, and velocities are normalized so that the average kinetic energy
! exactly corresponds to kT.
!
! NOTE: For MPI programs, only process 0 should call this subroutine.
!
!*******************************************************************************


      subroutine initconfig
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'mpif.h'

      character*256 mdoutfile     !output file name
      character*6   xftype        !input file name extension
      character*20  xtype         !type of initial configuration
      integer       irtn          !return code from file reader routines
                                  ! bit field
                                  !  0- 3 : File status
                                  !       :  =   0  : opened successfully
                                  !       :  =   1  : unknown file type
                                  !  4- 7 : positions 
                                  !       :  =   1  : got_x
                                  !  8-11 : velocities
                                  !       :  =   1  : got_v
                                  ! 12-15 : accelerations
                                  !       :  =   1  : got_a
                                  ! 16-19 : charges
                                  !       :  =   1  : got_Z
                                  ! 20-23 : masses
                                  !       :  =   1  : got_A
      integer       nni           !number of particles of a given species
      real(dble)    xzi,xai       !charge and mass
      logical       override_gage !use gage dimensions from parameter file
      real(dble)    xxl0(3)       !temporary, user-specified gage dimensions

      character(6)  filetype      !function that returns file name extension



!===============================================================================
! Read ion charges and masses.
!
! If the charges and masses are in the input configuration file, they will be
! read when that file is read instead of here.
!
! ------------------------------------------------------------------------------
! Set xtype to type of species file. Type 'internal' is for a list of species
! given directly in the runmd.in file.
      select case(trim(spec_file))
        case('internal')
          xtype=trim(spec_file)
        case default
          xftype=filetype(spec_file)
          xtype=xftype
      end select
      irtn=0

! ------------------------------------------------------------------------------
      select case(trim(xtype))

      case ('internal')
      !------------------
      !Get charges and masses from a spec_list read from the input parameter
      !file. Number of species and total number of particles should have been
      !determined when the spec_list was read.
        allocate(zii(0:n-1))
        allocate(aii(0:n-1))
        k=0
        j=0
        ztot=0.0   !counts total charge
        do
          k=k+1
          if(spec_list(k).eq.nullspecies) goto 103
          nni=spec_list(k)%n
          xzi=spec_list(k)%z
          xai=spec_list(k)%a
          ztot=ztot+nni*xzi
          do i=1,nni
            zii(j)=xzi
            aii(j)=xai
            j=j+1
          enddo
        enddo
        irtn=65536+1048576
  103   continue

      case ('ZAfa')
      !------------------
        write(6,10110) trim(spec_file)
        write(7,10110) trim(spec_file)
        call read_xvfa(spec_file,xftype,.true.,irtn)
        write(6,10120) n
        write(7,10120) n

      case ('ZAfb')
      !------------------
        write(6,10110) trim(spec_file)
        write(7,10110) trim(spec_file)
        call read_ZAfb(spec_file,xftype,.true.,irtn)
        write(6,10120) n
        write(7,10120) n

      case default
      !------------------
        continue

      end select

10110 format('  spec file:  ',a)
10120 format('     n = ',i9)



!===============================================================================
! Read or construct initial configuration.
!
! Some configuration files also contain charges and masses. If so, charge and
! mass arrays will be allocated and filled.
!
! ------------------------------------------------------------------------------
! Set xtype to either a type of initial configuration to generate, or type of
! file to read.
      select case(trim(start))
        case('random','nuclear')
          xtype=trim(start)
        case default
          xftype=filetype(start)
          xtype=xftype
      end select
      irtn=0

! If gage dimensions xl0(1:3) have been specified in the parameter file, they
! should override any read from a configuration file.  Rev.A and earlier con-
! figuration files do not contain either gage or current dimensions. For these
! files gage dimensions will be set to the calculated dimensions of the box
! unless the user overrides them. The same applies to 'random' and 'nuclear'
! starts.
      override_gage = xl0(1).gt.0.0d0

      write(6,10200)
      write(7,10200)
10200 format(/,'  INITIAL CONFIGURATION',/  &
               '  ---------------------')

! ------------------------------------------------------------------------------
      select case(trim(xtype))

        case('random')
        !------------------
        !Distribute particles randomly in a rectangular box.
        !DKB-todo: The following does not spread particles uniformly in the
        !DKB-todo: volume if the edge lengths are different.
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) aspect=1.d0
          xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
          xl0=aspect*xlf
          xl=xl0
          write(6,10300)
          write(7,10300)
10300     format('  Generating random coordinates in box')
          if(tstart.lt.0.) call md_exit(1,'*** initconfig:  tstart not specified ***')
          allocate(x(3,0:n-1))       !positions
          do i=0,n-1
            x(1,i)=xl(1)*ran1(iseed)
            x(2,i)=xl(2)*ran1(iseed)
            x(3,i)=xl(3)*ran1(iseed)
          enddo
          irtn=irtn+16

        case('nuclear')
        !------------------
        !Initialize particle coordinates in spherical volume of radius rmax cen-
        !tered in the rectangular simulation volume. Limit rmax to half the min-
        !imum of the box edge lengths. This type of simulation is meant for sim-
        !ulating the nucleons in a nucleus, such as Pb-208. A cubical simulation
        !volume is probably the only meaningful kind; nevertheless, we allow dif-
        !ferent edge lengths.
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) aspect=1.d0
          xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
          xl0=aspect*xlf
          xl=xl0
          xlmin=min(xl(1),xl(2),xl(3))
          rmax=min( rmax, xlmin*0.5 )
          write(6,10320) rmax
          write(7,10320) rmax
10320     format('  Generating random coordinates in sphere' /   &
                 '  radius     = ',f8.4,' fm')
          if(tstart.lt.0.) call md_exit(1,'*** initconfig:  tstart not specified ***')
          allocate(x(3,0:n-1))       !positions
          allocate(v(3,0:n-1))       !velocities
          do i=0,n-1
            xrr=ran1(iseed)**.3333333*rmax
            costheta = 2.*ran1(iseed)-1.
            sintheta = sqrt(1.d0-costheta*costheta)
            xphi = 2.d0*PI*ran1(iseed)
            x(1,i) = xrr*sintheta*cos(xphi)+0.5*xl(1)
            x(2,i) = xrr*sintheta*sin(xphi)+0.5*xl(2)
            x(3,i) = xrr*costheta+0.5*xl(3)
          enddo
          irtn=irtn+16

        case('onp','md1','xfa','xvfa')
        !--------------------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xvfa(start,xftype,.true.,irtn)
          if(tstart.lt.0.) tstart=time
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) aspect=1.d0
          xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
          xl0=aspect*xlf
          xl=stretch*xl0

        case('md0','x4a','xv4a','x8a','xv8a','xva8a')
        !----------------------------------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xva(start,xftype,.true.,irtn)
          if(tstart.lt.0.) tstart=time
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) aspect=1.d0
          xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
          xl0=aspect*xlf
          xl=stretch*xl0

        case('xvZAfa')
        !------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xvfa(start,xftype,.true.,irtn)
          write(6,10500) n
          write(7,10500) n
          if(tstart.lt.0.) tstart=time
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) aspect=1.d0
          xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
          xl0=aspect*xlf
          xl=stretch*xl0

        case('x4b','x8b','xv4b','xv8b')
        !--------------------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xvb(start,xftype,.true.,irtn)
          if(tstart.lt.0.) tstart=time
          rho=float(n)/(xl(1)*xl(2)*xl(3))
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) then
            xlf=minval(xl)
            aspect=xl/xlf
            xl0=xl
          else
            xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
            xl0=aspect*xlf
          endif
          stretch=xl/xl0

        case('xfb','xvfb')
        !-------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xvfb(start,xftype,.true.,irtn)
          if(tstart.lt.0.) tstart=time
          rho=float(n)/(xl(1)*xl(2)*xl(3))
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) then
            xlf=minval(xl)
            aspect=xl/xlf
            xl0=xl
          else
            xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
            xl0=aspect*xlf
          endif
          stretch=xl/xl0

        case('xZAfb','xvZAfb')
        !------------------
          write(6,10400) trim(start)
          write(7,10400) trim(start)
          call read_xvfb(start,xftype,.true.,irtn)
          write(6,10500) n
          write(7,10500) n
          if(tstart.lt.0.) tstart=time
          rho=float(n)/(xl(1)*xl(2)*xl(3))
          if(aspect(1)*aspect(2)*aspect(3).eq.0.d0) then
            xlf=minval(xl)
            aspect=xl/xlf
            xl0=xl
          else
            xlf=(float(n)/(aspect(1)*aspect(2)*aspect(3)*rho))**(1.d0/3.d0)
            xl0=aspect*xlf
          endif
          stretch=xl/xl0

        case default
        !------------------
          call md_exit(1,'*** initconfig:  input configuration type not supported ***')

      end select

10400     format('  file          = ',a)
10500     format('  This file also includes charges and masses.', / &
                 '     n          = ',i10)

!===============================================================================
! If velocities need to be generated, initalize them from a Maxwell distribution
! at temperature kT.  Ideally, the system will have zero center-of-mass motion,
! and satisfy eka=(3/2)kT, but the random numbers may not work out that way. So
! we call a subroutine that enforces these two conditions.
! -------------
      if(iand(irtn,256).eq.0) then
        write(6,10600)
        write(7,10600)
10600   format('  Generating velocities from Maxwell distribution')
        allocate(v(3,0:n-1))
        do i=0,n-1
          velfac=sqrt(kT/(aii(i)*xmass))
          v(1,i)=velfac*gasdev(iseed)
          v(2,i)=velfac*gasdev(iseed)
          v(3,i)=velfac*gasdev(iseed)
        enddo
        call cancel_vcom   !cancel center of mass motion; normalize to kT
      endif

! ------------------------------------------------------------------------------
! Start time and time step for simulation (fm/c).
      time=tstart
      write(6,10700) tstart,dt
      write(7,10700) tstart,dt
10700 format('  start time    =  ', f14.2,' fm/c' /       &
             '  time step     =        ', f8.2,' fm/c')

! ------------------------------------------------------------------------------
! If the input configuration was not from an xv8b file, write out the initial
! configuration to such a file.
      if(xtype.ne.'xv8b') then
        write(mdoutfile,110) int(time/1000000.d0), int(mod(time,1000000.d0))
  110   format('md.',i5.5,i6.6,'.xv8b')
        call write_xvb(mdoutfile,'xv8b  ',.false.,.true.)
      endif

! ------------------------------------------------------------------------------
! Allocate other particle arrays that will be needed.
      allocate(vold(3,0:n-1))      !old velocities
      if(.not.allocated(a)) then
        allocate(a(3,0:n+2))       !accelerations. Indices n,n+1,n+2 are used
                                   !   by MPI_allreduce in force routine
      endif

      return
      end subroutine initconfig
