!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-May-16
! ---------------------------------------------------------------------
!
! Reader and writer routines for revision A unformatted MD configuration files.
! These files were used in versions of the MD code prior to MD_6.2.0, and are
! known in general as 'xv' files. Readers and writers for these files are grouped
! into one reader routine and one writer routine because they share a similar
! layout. Depending on input arguments, the routines read/write either an MD con-
! figuration (positions only) or a phase space point (positions and velocities).
! There is also a file type for positions, velocities and accelerations. Argument
! filetype determines what to write out, and in what precision. Argument xappend
! determines whether to write each configuration to a separate file, or append
! configurations to a single trajectory file.
!
! Arguments:
!
!   character*6  xfiletype  -- specifies what kind of data the file contains
!       'md0'  = real positions, velocities. no time stamp.
!       'x4a'  = timestamp. real positions
!       'xv4a' = timestamp. real positions, velocities.
!       'x8a'  = timestamp. real(dble) positions
!       'xv8a' = timestamp. real(dble) positions, velocities.
!       'xva8a'= timestamp. real(dble) positions, velocities, accelerations.
!
!   logical      xappend  -- 
!       .true . = append all configs to a single trajectory file
!       .false. = write each config to a separate file
!
!*******************************************************************************



!*******************************************************************************
      subroutine read_xva(file,xftype,xclose,irtn)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !file to read
      character*6     xftype         !file type. If unknown, should be 'unknwn'
      logical         xclose         !.true.  = close file after reading
                                     !.false. = leave file open
      integer         irtn           !return code

      character*6     xfiletype
      real            time4
      integer         size
      integer         i,j

      real, allocatable  :: x4(:,:)  !particle coordinates, single precision
      real, allocatable  :: v4(:,:)  !particle velocities, single precision

      irtn=0
      if(xftype.eq.'unknwn') then
        j = len_trim(file)
        i = scan(file,'.',.true.)
        if(i.eq.0 .or. (j-i).gt.6) then
          irtn=irtn+1
          return
        endif
        xftype=file(i+1:j)
      endif
      open(13,FILE=file,STATUS='OLD',FORM='UNFORMATTED',POSITION='ASIS')

      select case (trim(xftype))

        case('md0')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+4)    !record 1: (4+4 bytes for record markers)
            n=n/(6*8)       !record 1: 3 real*8 coordinates and velocities
          endif
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
          irtn=irtn+16+256

        case('x4a')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+4+4)  !record 1: (4+4 byte record markers) + (real*4 timestamp)
            n=n-(4+4)       !record 2: (4+4 byte record markers)
            n=n/(3*4)       !record 2: each particle has 3 real*4 positions
          endif
          allocate(x4(3,0:n-1))
          read(13) time4
          time=time4
          read(13) (x4(1,i),x4(2,i),x4(3,i), i=0,n-1)
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          x = x4
          deallocate(x4)
          irtn=irtn+16

        case('xv4a')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+4+4)  !record 1: (4+4 byte record markers) + (real*4 timestamp)
            n=n-(4+4)       !record 2: (4+4 byte record markers)
            n=n/(6*4)       !record 2: each particle has 3 real*4 positions and velocities
          endif
          allocate(x4(3,0:n-1))
          allocate(v4(3,0:n-1))
          read(13) time4
          time=time4
          read(13) (x4(1,i),x4(2,i),x4(3,i),v4(1,i),v4(2,i),v4(3,i), i=0,n-1)
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          x = x4
          v = v4
          deallocate(x4)
          deallocate(v4)
          irtn=irtn+16+256

        case('x8a')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+8+4)  !record 1: (4+4 byte record markers) + (real*8 timestamp)
            n=n-(4+4)       !record 2: (4_4 byte record markers)
            n=n/(3*8)       !record 2: 3 real*8 coordinates
          endif
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          read(13) time
          read(13) (x(1,i),x(2,i),x(3,i), i=0, n-1)
          irtn=irtn+16

        case('xv8a')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+8+4)  !record 1: (4+4 byte record markers) + (real*8 timestamp)
            n=n-(4+4)       !record 2: (4_4 byte record markers)
            n=n/(6*8)       !record 2: 3 real*8 coordinates and velocities
          endif
          write(6,*) 'xv8a: ',n
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          read(13) time
          read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
          irtn=irtn+16+256

        case('xva8a')
          if(n.eq.0) then
            call file_stat(trim(file),size)
            n=size-(4+8+4)  !record 1: (4+4 byte record markers) + (real*8 timestamp)
            n=n-(4+4)       !record 2: (4_4 byte record markers)
            n=n/(9*8)       !record 2: 3 real*8 coordinates, velocities, accelerations
          endif
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
         !Indices n,n+1,n+2 of a are used by MPI_allreduce in force routine
          if(.not.allocated(a)) allocate(a(3,0:n+2))
          read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),a(1,i),a(2,i),a(3,i), i=0,n-1)
          irtn=irtn+16+256+4096

      end select

      if(xclose) close(13)

      return
      end subroutine read_xva



!*******************************************************************************
      subroutine write_xva(file,xftype,xappend,xclose)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !name of file to write
      character*6     xftype         !file type
      logical         xappend        !.true.  = append to existing file
                                     !.false. = write new file
      logical         xclose         !.true.  = close file after writing
                                     !.false. = leave file open
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

      select case (trim(xftype))

        case('md0')
          write(14) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)

        case('x4a')
          write(14) real(time)
          write(14) (real(x(1,i)),real(x(2,i)),real(x(3,i)), i=0,n-1)

        case('xv4a')
          write(14) real(time)
          write(14) (real(x(1,i)),real(x(2,i)),real(x(3,i)),            &
                       real(v(1,i)),real(v(2,i)),real(v(3,i)), i=0,n-1)

        case('x8a')
          write(14) time
          write(14) (x(1,i),x(2,i),x(3,i), i=0,n-1)

        case('xv8a')
          write(14) time
          write(14) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)

        case('xva8a')
          write(14) time
          write(14) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),a(1,i),a(2,i),a(3,i), i=0,n-1)

      end select

      if(xclose) close(14)

      return
      end subroutine write_xva
