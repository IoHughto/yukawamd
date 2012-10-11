!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Jul-09
! ---------------------------------------------------------------------
!
!*******************************************************************************

      subroutine read_ZAfb(file,xftype,xclose,irtn)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !file to read
      character*6     xftype         !file type. If unknown, should be 'unknwn'
      logical         xclose         !.true.  = close file after reading
                                     !.false. = leave file open
      integer         irtn           !return code

      real(dble)      dummy
      real(dble)      xzi,xai
      integer         i,j,k
      character*256   line

      irtn=0
      if(xftype.eq.'unknwn') then
        j = len_trim(file)
        i = scan(file,'.',.true.)
        if(i.eq.0 .or. (j-i).gt.6) then
          irtn=1
          return
        endif
        xftype=file(i+1:j)
      endif
!Open the file, read and discard any comments (lines beginning with #), and
!read the total number of particles. This must be by itself on the first non-
!comment line.
      open(13,FILE=trim(file),FORM='FORMATTED',STATUS='OLD')
      line(1:1)='#'
      do while(line(1:1).eq.'#')
        read(13,'(a)') line
      enddo
      read(line,*) n

!Allocate arrays for n particles.
      if(.not.allocated(zii)) allocate(zii(0:n-1))  !charge
      if(.not.allocated(aii)) allocate(aii(0:n-1))  !mass
      i=0
      ztot=0.0d0
      do
        read(13,*,end=500) k, xzi, xai
        do j=1,k
          zii(i)=xzi
          aii(i)=xai
          i=i+1
        enddo
        ztot=ztot+k*xzi
      enddo
      irtn=irtn+65536+1048576
  500 continue

      if(xclose) close(13)
      return
      end subroutine read_ZAfb



!*******************************************************************************
      subroutine write_ZAfb(file,xftype,xappend,xclose)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !name of file to write
      character*6     xftype         !file type
      logical         xappend        !.true.  = append to existing file
                                     !.false. = write new file
      logical         xclose         !.true.  = close file after writing
                                     !.false. = leave file open
      integer         i,j,k
      logical         xopnd          !.true. = unit 14 is already opened

      inquire(14,opened=xopnd)
      if(.not.xopnd) then
        if(xappend) then
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        else
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED')
        endif
      endif

      write(14,10020) n
10020 format(i10)
      i=0
      do while(i.lt.n)
        zi=zii(i)
        ai=aii(i)
        k=0
        do while(zii(i).eq.zi .and. aii(i).eq.ai)
          k=k+1
          i=i+1
          if(i.eq.n) exit
        enddo
        write(14,10024) k,zi,ai
10024   format(i8,1x,f8.3,1x,f8.3)
      enddo

      if(xclose) close(14)

      return
      end subroutine write_ZAfb



!*******************************************************************************
      subroutine read_xvfb(file,xftype,xclose,irtn)
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
      real(dble)      dummy
      integer         i,j,k
      character*256   line

      irtn=0
      if(xftype.eq.'unknwn') then
        j = len_trim(file)
        i = scan(file,'.',.true.)
        if(i.eq.0 .or. (j-i).gt.6) then
          irtn=1
          return
        endif
        xftype=file(i+1:j)
      endif
      open(13,FILE=file,STATUS='OLD',FORM='FORMATTED',POSITION='ASIS')
      read(13,*) xfiletype
      read(13,*) xcode_name, xcode_version
      read(13,*) date,daytime,timezone
      read(13,*) xsim_type
      read(13,*) xl0(1),xl0(2),xl0(3)
      read(13,*) time,xl(1),xl(2),xl(3), ev,ek, px, pp, n
      rho=n/(xl(1)*xl(2)*xl(3))

      select case (trim(xftype))

        case('xfb')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          do i=0,n-1
            read(13,*,END=300) x(1,i),x(2,i),x(3,i)
          enddo
          irtn=irtn+16
  300     continue

        case('xvfb')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          do i=0,n-1
            read(13,*,END=400) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
          enddo
          irtn=irtn+16+256
  400     continue

        case('xZAfb')
          if(.not.allocated(x)) allocate(x(3,0:n-1))       !position
          if(.not.allocated(zii)) allocate(zii(0:n-1))       !charge
          if(.not.allocated(aii)) allocate(aii(0:n-1))       !mass
          ztot=0.0d0
          do i=0,n-1
            read(13,*,END=500) x(1,i),x(2,i),x(3,i),zii(i),aii(i)
            ztot=ztot+zii(i)
          enddo
          irtn=irtn+16+65536+1048576
  500     continue

        case('xvZAfb')
          if(.not.allocated(x)) allocate(x(3,0:n-1))       !position
          if(.not.allocated(v)) allocate(v(3,0:n-1))       !velocity
          if(.not.allocated(zii)) allocate(zii(0:n-1))       !charge
          if(.not.allocated(aii)) allocate(aii(0:n-1))       !mass
          ztot=0.0d0
          do i=0,n-1
            read(13,*,END=600) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),zii(i),aii(i)
            ztot=ztot+zii(i)
          enddo
          irtn=irtn+16+256+65536+1048576
  600     continue

      end select

      if(xclose) close(13)

      end subroutine read_xvfb



!*******************************************************************************
      subroutine write_xvfb(file,xftype,xappend,xclose)
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
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        else
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED')
        endif
      endif
      call date_and_time(date,daytime,timezone)
      if(trim(xftype).ne.'XYZ') then
        write(14,10010) trim(xftype)
        write(14,10010) trim(code_name), trim(code_version)
        write(14,10010) trim(date),trim(daytime),trim(timezone)
        write(14,10010) trim(sim_type)
        write(14,10050) xl0(1),xl0(2),xl0(3)
        write(14,10050) time
        write(14,10050) xl(1),xl(2),xl(3)
        write(14,10054) ev,ek,px
        write(14,10058) pp
        write(14,10060) n
      endif
10010 format('# ',a,2x,a,2x,a)
10050 format('#',/,'# ',3(1x,f12.6))
10054 format('#',/,'# ',3(1x,es15.8))
10058 format('#',/,2('# ',3(1x,es15.8),/),'# ',3(1x,es15.8))
10060 format('# ',i10)

      select case (trim(xftype))

        case('xfb')
          write(14,10114) (x(1,i),x(2,i),x(3,i), i=0,n-1)
10114     format(es24.16,es24.16,es24.16)

        case('xvfb')
          write(14,10124) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
10124     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16)

        case('xZAfb')
          write(14,10134) (x(1,i),x(2,i),x(3,i),zii(i),aii(i), i=0,n-1)
10134     format(es24.16,es24.16,es24.16,1x,f8.3,1x,f8.3)

        case('xvZAfb')
          write(14,10144) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),zii(i),aii(i), i=0,n-1)
10144     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16,1x,f8.3,1x,f8.3)

        case('XYZ')
          write(14,10160) n
10160     format(i10)
          write(14,10062) xl(1),xl(2),xl(3)
         !write(14,"('[XYZ file]')")
          write(14,10164) (ctype(i),x(1,i),x(2,i),x(3,i), i=0,n-1)
10062     format('# ',3(1x,f12.6))
10164     format(1x,a6,1x,f12.6,1x,f12.6,1x,f12.6)

      end select

      if(xclose) close(14)

      return
      end subroutine write_xvfb
