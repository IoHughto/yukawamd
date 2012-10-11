!*******************************************************************************
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-May-16
! ---------------------------------------------------------------------
!
!*******************************************************************************

      subroutine read_xvfa(file,xftype,xclose,irtn)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !file to read
      character*6     xftype         !file type. If unknown, should be 'unknwn'
      logical         xclose         !.true.  = close file after reading
                                     !.false. = leave file open
      integer         irtn           !return code

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

      select case (trim(xftype))

        case('onp')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          do i=0,n-1
            read(13,*,END=100) x(1,i),x(2,i),x(3,i),dummy
          enddo
          irtn=irtn+16
  100     continue

        case('md1')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          do i=0,n-1
            read(13,*,END=200) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
          enddo
          irtn=irtn+16+256
  200     continue

        case('xfa')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          read(13,*) time
          do i=0,n-1
            read(13,*,END=300) x(1,i),x(2,i),x(3,i)
          enddo
          irtn=irtn+16
  300     continue

        case('xvfa')
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          read(13,*) time
          do i=0,n-1
            read(13,*,END=400) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
          enddo
          irtn=irtn+16+256
  400     continue

        case('ZAfa')
          i=0
          line(1:1)='#'
          do while(line(1:1).eq.'#')
             read(13,'(a)') line
          enddo
          read(line,20010) k
          if(n==0) n=k
20010     format(i10)
          if(.not.allocated(zii)) allocate(zii(0:n-1))
          if(.not.allocated(aii)) allocate(aii(0:n-1))
          do i=0,n-1
            read(13,*) k, zii(k-1), aii(k-1)
          enddo
          irtn=irtn+65536+1048576

        case('xvZAfa')
          i=0
          line(1:1)='#'
          do while(line(1:1).eq.'#')
             read(13,'(a)') line
          enddo
          read(line,20010) k
          if(n==0) n=k
          if(.not.allocated(zii)) allocate(zii(0:n-1))
          if(.not.allocated(aii)) allocate(aii(0:n-1))
          if(.not.allocated(x)) allocate(x(3,0:n-1))
          if(.not.allocated(v)) allocate(v(3,0:n-1))
          read(13,*) time
          do i=0,n-1
            read(13,*,END=500) x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),zii(i),aii(i)
          enddo
          irtn=irtn+16+256+65536+1048576
  500     continue

      end select

      if(xclose) close(13)

      end subroutine read_xvfa



!*******************************************************************************
      subroutine write_xvfa(file,xftype,xappend,xclose)
      use  md_types
      use  md_globals
      implicit none

      character*(*)   file           !name of file to write
      character*6     xftype         !file type
      logical         xappend        !.true.  = append to existing file
                                     !.false. = write new file
      logical         xclose         !.true.  = close file after writing
                                     !.false. = leave file open
      integer         i,j
      logical         xopnd          !.true. = unit 14 is already opened

      inquire(14,opened=xopnd)
      if(.not.xopnd) then
        if(xappend) then
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        else
          open(14,FILE=file,STATUS='UNKNOWN',FORM='FORMATTED')
        endif
      endif

      select case (trim(xftype))

        case('md1')
          write(14,10014) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
10014     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16)

        case('xfa')
          write(14,10020) time
10020     format(f15.3)
          write(14,10024)  (x(1,i),x(2,i),x(3,i), i=0,n-1)
10024     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16)

        case('xvfa')
          write(14,10030) time
10030     format(f15.3)
          write(14,10034) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=0,n-1)
10034     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16)

        case('ZAfa')
          write(14,10040) n
10040     format(i10)
          write(14,10044) (i,zii(i),aii(i), i=0,n-1)
10044     format(i8,1x,f8.3,1x,f8.3)

        case('xvZAfa')
          write(14,10050) time
10050     format(f15.3)
          write(14,10054) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i),zii(i),aii(i), i=0,n-1)
10054     format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16,1x,f8.3,1x,f8.3)

      end select

      if(xclose) close(14)

      return
      end subroutine write_xvfa
