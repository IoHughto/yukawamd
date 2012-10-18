
!*******************************************************************************
!      subroutine read_xv8(file,n,time,x,v)
!      implicit none
!      character*(*)   file           !file to read
!      integer         n
!      character*6     xfiletype
!      character*10    xcode_name
!      character*8     xcode_version
!      character*8     date
!      character*10    daytime
!      character*5     timezone
!      character*20    xsim_type
!      integer         i
!      real            xtime

!      real*8          time,xl(3),ev,ek,px,pp,rho
!      real*8          x(3,n)  !positions
!      real*8          v(3,n)  !velocities


!      open(UNIT=13,FILE=trim(file),STATUS='OLD',FORM='UNFORMATTED',POSITION='ASIS')
!      read(13) xfiletype
!      read(13) xcode_name, xcode_version
!      read(13) date,daytime,timezone
!      read(13) xsim_type
!      read(13) time,xl(1),xl(2),xl(3), ev,ek, px, pp, n
!      read(13) (x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i), i=1,n)

!      return
!      end subroutine read_xv8


      subroutine read_xv8(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      character*6     xfiletype
      character*10    xcode_name
      character*8     xcode_version
      character*8     date
      character*10    daytime
      character*5     timezone
      character*20    xsim_type
      real            xtime

      real*8          time,xl(3),ev,ek,px,pp,rho
      integer        i
      integer        size
      open(UNIT=13, FILE=trim(file), STATUS='OLD', FORM='UNFORMATTED')
      read(13) xfiletype
      read(13) xcode_name, xcode_version
      read(13) date,daytime,timezone
      read(13) xsim_type
      read(13) time,xl(1),xl(2),xl(3), ev,ek, px, pp, n
      read(13) (x8(1,i),x8(2,i),x8(3,i),v8(1,i),v8(2,i),v8(3,i), i=1,n)
      close(13)
      return
      end subroutine read_xv8
