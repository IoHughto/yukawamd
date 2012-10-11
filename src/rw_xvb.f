!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Jul-12
! ---------------------------------------------------------------------
!
! This file contains I/O routines for revision b unformatted MD configuration
! files. These files are known in general as 'xv' files, and share a similar
! layout. Depending on input arguments, the routines read/write either an MD
! configuration (positions only) or a phase space point (positions and veloci-
! ties). There is also a file type that saves positions, velocities and accel-
! erations. Argument xftype determines what to write out, and in what precision.
! Argument xappend determines whether to write each configuration to a separate
! file, or append configurations to a single trajectory file.
!
! Arguments:
!
!   character*6  xftype  -- specifies what kind of data the file contains:
!       'x4b'  = real positions
!       'xv4b' = real positions, velocities.
!       'x8b'  = real(dble) positions
!       'xv8b' = real(dble) positions, velocities.
!       'xva8b'= real(dble) positions, velocities, accelerations.
!
!   logical      xappend  -- 
!       .true . = append all configs to a single trajectory file
!       .false. = write each config to a separate file
!
!*******************************************************************************



!*******************************************************************************
      subroutine read_xvb(file,xftype,xclose,irtn)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !file to read
      character*6     xftype         !file type. If unknown, should be 'unknwn'
      logical         xclose         !.true.  = close file after reading
                                     !.false. = leave file open
      integer         irtn           !return code

      character*6     xfiletype
      character*10    xcode_name
      character*8     xcode_version
      character*8     date
      character*10    daytime
      character*5     timezone
      character*20    xsim_type
      integer         i
      real            xtime

      real, allocatable  :: x4(:,:)  !particle coordinates, single precision
      real, allocatable  :: v4(:,:)  !particle velocities, single precision

      irtn=0
      if(xftype.eq.'unknwn') then
        i = scan(file,'.',.true.)
        if(i.eq.0 .or. (len_trim(file)-i).gt.6) then
          irtn=1
          return
        endif
        xftype=file(i+1:)
      endif

      open(13,FILE=trim(file),STATUS='OLD',FORM='UNFORMATTED',POSITION='ASIS')
      read(13) xfiletype
      read(13) xcode_name, xcode_version
      read(13) date,daytime,timezone
      read(13) xsim_type
      read(13) time,xl(1),xl(2),xl(3), ev,ek, px, pp, n  !for current configuration

!DKB-debug:
!      write(6,10010) xfiletype
!10010 format('rw_xvb: filetype   =',a)
!      write(6,10012) xcode_name,xcode_version
!10012 format('rw_xvb: code_name  =',a,'  code_version = ',a)
!      write(6,10014) date,daytime,timezone
!10014 format('rw_xvb: date   =',a,2x,a,2x,a)
!      write(6,10016) xsim_type
!10016 format('rw_xvb: sim_type    =',a)
!      write(6,10018) xlf, aspect, time, xl, ev,ek,px, pp, n
!10018 format('rw_xvb: xlf     =',es17.8  / &
!             'rw_xvb: aspect  =',3es17.8 / &
!             'rw_xvb: time    =',f12.3   / &
!             'rw_xvb: xl      =',3es17.8 / &
!             'rw_xvb: ev,ek,px=',3es17.8 / &
!             'rw_xvb: pp      =',3es17.8 / &
!             'rw_xvb:          ',3es17.8 / &
!             'rw_xvb:          ',3es17.8 / &
!             'rw_xvb: n       =',i10)

      select case (trim(xftype))

        case('x4b')
          allocate(x4(3,0:n-1))
          read(13) (x4(1,i),x4(2,i),x4(3,i), i=0,n-1)
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          x = x4
          deallocate(x4)
          irtn=irtn+16

        case('xv4b')
          allocate(x4(3,0:n-1))
          allocate(v4(3,0:n-1))
          read(13) (x4(1,i),x4(2,i),x4(3,i),v4(1,i),v4(2,i),v4(3,i), i=0,n-1)
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          x = x4
          v = v4
          deallocate(x4)
          deallocate(v4)
          irtn=irtn+16+256

        case('x8b')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          read(13) (x(1,i),x(2,i),x(3,i), i=0,n-1)
          irtn=irtn+16

        case('xv8b')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
          irtn=irtn+16+256

        case('xva8b')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
         !Indices n,n+1,n+2 of a are used by MPI_allreduce in force routine
          if(.not.allocated(a)) allocate(a(3,0:n+2))
          read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),a(1,i),a(2,i),a(3,i), i=0,n-1)
          irtn=irtn+16+256+4096

      end select

      if(xclose) close(13)

      return
      end subroutine read_xvb



!*******************************************************************************
      subroutine write_xvb(file,xftype,xappend,xclose)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !name of file to write
      character*6     xftype         !file type
      logical         xappend        !.true.  = append to existing file
                                     !.false. = write new file
      logical         xclose         !.true.  = close file after writing
                                     !.false. = leave file open
      character*8     date
      character*10    daytime
      character*5     timezone
      integer         i
      logical         xopnd          !.true. = unit 14 is already opened

      inquire(14,opened=xopnd)
      if(.not.xopnd) then
        if(xappend) then
          open(14,FILE=file,STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
        else
          open(14,FILE=file,STATUS='UNKNOWN',FORM='UNFORMATTED')
        endif
      endif
      call date_and_time(date,daytime,timezone)
      write(14) xftype
      write(14) code_name,code_version
      write(14) date,daytime,timezone
      write(14) sim_type
      write(14) time, xl(1),xl(2),xl(3), ev,ek, px, pp, n

      select case (trim(xftype))

        case('x4b')
          write(14) (real(x(1,i)),real(x(2,i)),real(x(3,i)), i=0,n-1)

        case('xv4b')
          write(14) (real(x(1,i)),real(x(2,i)),real(x(3,i)),            &
                       real(v(1,i)),real(v(2,i)),real(v(3,i)), i=0,n-1)

        case('x8b')
          write(14) (x(1,i),x(2,i),x(3,i), i=0,n-1)

        case('xv8b')
          write(14) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)

        case('xva8b')
          write(14) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),a(1,i),a(2,i),a(3,i), i=0,n-1)

      end select

      if(xclose) close(14)

      return
      end subroutine write_xvb
