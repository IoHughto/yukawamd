!*******************************************************************************
!  This file contains readers and writers for various files used in the MD
!  series of programs.
!*******************************************************************************


!*******************************************************************************
!  File type: x4
!  ------------------
!  Read/write time-stamped, unformatted real*4 configuration file containing two
!  records. The first contains a simulation timestamp; the second contains car-
!  tesian coordinates of all particles.
!
!     Record 1:
!     ----------------
!     time4                  timestamp, real*4
!
!     Record 2:
!     ----------------
!     x(1) y(1) z(1)         particle 1 coordinates, real*4
!     x(2) y(2) z(2)         particle 2 coordinates, real*4
!       .    .    .
!       .    .    .
!       .    .    .
!     x(n) y(n) z(n)         particle n position, real*4
!
!  The calling program sends read_x4 the number of particles n that it wants
!  from the file. The first n will be read and returned. If n is set to 0, all
!  particles will be returned, and n will be set to the number of particles.
!  When writing, the program must send n to write_x4, so it will know how many
!  particles to write.

      subroutine read_x4(file,n,time4,x4)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !particle coordinates
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(4+4+4)  !record 1: (8 bytes for record size) + (real*4 timestamp)
        n=n-(4+4)       !record 2: (8 bytes for record size)
        n=n/(3*4)       !record 2: each particle has 3 real*4 coordinates
      endif
      open(UNIT=13, FILE=file, STATUS='OLD', FORM='UNFORMATTED')
      read(13) time4
      read(13) (x4(1,i),x4(2,i),x4(3,i), i=1,n)
      close(13)
      return
      end subroutine read_x4

      subroutine write_x4(file,n,time4,x4)
      implicit none
      character*(*)  file     !file to write
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      integer        i
      open(UNIT=14, FILE=file, STATUS='UNKNOWN', FORM='UNFORMATTED')
      write(14) time4
      write(14) (x4(1,i),x4(2,i),x4(3,i), i=1,n)
      close(14)
      return
      end subroutine write_x4


!*******************************************************************************
!  File type: xv4
!  ------------------
!  Read/write time-stamped, unformatted real*4 configuration file containing two
!  records. The first contains a simulation timestamp; the second contains posi-
!  tions and velocities of all particles.
!
!     Record 1:
!     ----------------
!     time4                 timestamp, real*4
!
!     Record 2:
!     ----------------
!     x(1) y(1) z(1) vx(1) vy(1) vz(1)   particle 1 position, velocity, real*4
!     x(2) y(2) z(2) vx(2) vy(2) vz(2)   particle 2 position, velocity, real*4
!       .    .    .     .     .     .
!       .    .    .     .     .     .
!       .    .    .     .     .     .
!     x(n) y(n) z(n) vx(n) vy(n) vz(n)   particle n position, velocity, real*4
!
!  The calling program sends read_xv4 the number of particles n that it wants
!  from the file. The first n will be read and returned. If n is set to 0, all
!  particles will be returned, and n will be set to the number of particles.
!  When writing, the program must send n to write_xv4, so it will know how many
!  particles to write.

      subroutine read_xv4(file,n,time4,x4,v4)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      real*4         v4(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+4)    !record 1: (8 byte record size) + (real*4 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*4)       !record 2: each particle has 3 real*4 positions and velocities
      endif
      open(UNIT=13, FILE=file, STATUS='OLD', FORM='UNFORMATTED')
      read(13) time4
      read(13) (x4(1,i),x4(2,i),x4(3,i),v4(1,i),v4(2,i),v4(3,i), i=1,n)
      close(13)
      return
      end subroutine read_xv4

      subroutine open_xv4a(file,n,time4,x4,v4)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      real*4         v4(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+4)    !record 1: (8 byte record size) + (real*4 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*4)       !record 2: each particle has 3 real*4 positions and velocities
      endif
      open(UNIT=13, FILE=file, STATUS='OLD', FORM='UNFORMATTED')
      return
      end subroutine open_xv4a

      subroutine read_xv4a(file,n,time4,x4,v4)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      real*4         v4(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+4)    !record 1: (8 byte record size) + (real*4 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*4)       !record 2: each particle has 3 real*4 positions and velocities
      endif
      read(13) time4
      read(13) (x4(1,i),x4(2,i),x4(3,i),v4(1,i),v4(2,i),v4(3,i), i=1,n)
      return
      end subroutine read_xv4a

      subroutine close_xv4a(file,n,time4,x4,v4)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      real*4         v4(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+4)    !record 1: (8 byte record size) + (real*4 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*4)       !record 2: each particle has 3 real*4 positions and velocities
      endif
      close(13)
      return
      end subroutine close_xv4a

      subroutine write_xv4(file,n,time4,x4,v4)
      implicit none
      character*(*)  file     !file to write
      integer        n        !number of particles
      real*4         time4    !time stamp
      real*4         x4(3,n)  !positions
      real*4         v4(3,n)  !velocities
      integer        i
      open(UNIT=14, FILE=file, STATUS='UNKNOWN', FORM='UNFORMATTED')
      write(14) time4
      write(14) (x4(1,i),x4(2,i),x4(3,i),v4(1,i),v4(2,i),v4(3,i), i=1,n)
      close(14)
      return
      end subroutine write_xv4


!*******************************************************************************
!  File type: x8
!  ------------------
!  Read/write time-stamped, unformatted real*8 configuration file containing two
!  records. The first contains a simulation timestamp; the second contains posi-
!  tions of all particles.
!
!     Record 1:
!     ----------------
!     time8            timestamp, real*8
!
!     Record 2:
!     ----------------
!     x(1) y(1) z(1)   particle 1 position, real*8
!     x(2) y(2) z(2)   particle 2 position, real*8
!       .    .    .  
!       .    .    .  
!       .    .    .  
!     x(n) y(n) z(n)   particle n position, real*8
!
!  The calling program sends read_x8 the number of particles n that it wants
!  from the file. The first n will be read and returned. If n is set to 0, all
!  particles will be returned, and n will be set to the number of particles.
!  When writing, the program must send n to write_x8, so it will know how many
!  particles to write.

      subroutine read_x8(file,n,time8,x8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+8)    !record 1: (8 byte record size) + (real*8 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(3*8)       !record 2: each particle has 3 real*8 coordinates
      endif
      open(UNIT=13, FILE=trim(file), STATUS='OLD', FORM='UNFORMATTED')
      read(13) time8
      read(13) (x8(1,i),x8(2,i),x8(3,i), i=1,n)
      close(13)
      return
      end subroutine read_x8

      subroutine write_x8(file,n,time8,x8)
      implicit none
      character*(*)  file     !file to write
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      integer        i
      open(UNIT=14, FILE=file, STATUS='UNKNOWN', FORM='UNFORMATTED')
      write(14) time8
      write(14) (x8(1,i),x8(2,i),x8(3,i), i=1,n)
      close(14)
      return
      end subroutine write_x8


!*******************************************************************************
!  File type: xv8
!  ------------------
!  Read/write time-stamped, unformatted real*8 configuration file containing two
!  records. The first contains a simulation timestamp; the second contains posi-
!  tions and velocities of all particles.
!
!     Record 1:
!     ----------------
!     time8                 timestamp, real*8
!
!     Record 2:
!     ----------------
!     x(1) y(1) z(1) vx(1) vy(1) vz(1)   particle 1 position, velocity, real*8
!     x(2) y(2) z(2) vx(2) vy(2) vz(2)   particle 2 position, velocity, real*8
!       .    .    .     .     .     .
!       .    .    .     .     .     .
!       .    .    .     .     .     .
!     x(n) y(n) z(n) vx(n) vy(n) vz(n)   particle n position, velocity, real*8
!
!  The calling program sends read_xv8 the number of particles n that it wants
!  from the file. The first n will be read and returned. If n is set to 0, all
!  particles will be returned, and n will be set to the number of particles.
!  When writing, the program must send n to write_xv8, so it will know how many
!  particles to write.

      subroutine read_xv8(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+8)    !record 1: (8 byte record size) + (real*8 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*8)       !record 2: 3 real*8 coordinates and velocities
      endif
      open(UNIT=13, FILE=trim(file), STATUS='OLD', FORM='UNFORMATTED')
      read(13) time8
      read(13) (x8(1,i),x8(2,i),x8(3,i),v8(1,i),v8(2,i),v8(3,i), i=1,n)
      close(13)
      return
      end subroutine read_xv8

      subroutine open_xv8a(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+8)    !record 1: (8 byte record size) + (real*8 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*8)       !record 2: 3 real*8 coordinates and velocities
      endif
      open(UNIT=13, FILE=trim(file), STATUS='OLD', FORM='UNFORMATTED')
      return
      end subroutine open_xv8a

      subroutine read_xv8a(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+8)    !record 1: (8 byte record size) + (real*8 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*8)       !record 2: 3 real*8 coordinates and velocities
      endif
      read(13) time8
      read(13) (x8(1,i),x8(2,i),x8(3,i),v8(1,i),v8(2,i),v8(3,i), i=1,n)
      return
      end subroutine read_xv8a

      subroutine close_xv8a(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      integer        size
      if(n.eq.0) then
        call file_stat(trim(file),size)
        n=size-(8+8)    !record 1: (8 byte record size) + (real*8 timestamp)
        n=n-8           !record 2: (8 byte record size)
        n=n/(6*8)       !record 2: 3 real*8 coordinates and velocities
      endif
      close(13)
      return
      end subroutine close_xv8a

      subroutine write_xv8(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to write
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      open(UNIT=14, FILE=file, STATUS='UNKNOWN', FORM='UNFORMATTED')
      write(14) time8
      write(14) (x8(1,i),x8(2,i),x8(3,i),v8(1,i),v8(2,i),v8(3,i), i=1,n)
      close(14)
      return
      end subroutine write_xv8


!*******************************************************************************
!  File type: xvf
!  ------------------
!  Read/write time-stamped, formatted configuration file containing:
!     x v
!  This would be an md.ckpt file, or a double precision md.out file. The calling
!  program should send the number of particles n in the file. But if this is
!  unknown, set n=0, and the reader should be able to figure it out. It will
!  simply read until reaching an end-of-file, and set n to the number of particles
!  read. This number is returned to the calling program.

      subroutine read_xvf(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to read
      integer        n        !number of particles
      real*8         time8    !time stamp
      real*8         x8(3,n)  !positions
      real*8         v8(3,n)  !velocities
      integer        i
      open(UNIT=13, FILE=file, STATUS='OLD', FORM='FORMATTED')
      read(13,*) time8
      i=1
      do
        read(13,*,END=200) x8(1,i), x8(2,i), x8(3,i), v8(1,i), v8(2,i), v8(3,i)
        i=i+1
        if(n>0 .and. i==n+1) goto 200
      enddo
  200 continue
      i=i-1
      close(13)
      if(n==0) n=i
      end subroutine read_xvf

      subroutine write_xvf(file,n,time8,x8,v8)
      implicit none
      character*(*)  file     !file to write  (in)
      integer        n        !number of particles  (in)
      real*8         time8    !time stamp  (in)
      real*8         x8(3,n)  !positions  (in)
      real*8         v8(3,n)  !velocities  (in)
      integer        i
      open(UNIT=14, FILE=file, STATUS='UNKNOWN', FORM='FORMATTED')
      write(14,10010) time8
10010 format(f15.3)
      write(14,10020)  (x8(1,i),x8(2,i),x8(3,i),v8(1,i),v8(2,i),v8(3,i), i=1,n)
10020 format(es24.16,es24.16,es24.16,es24.16,es24.16,es24.16)
      close(14)
      return
      end subroutine write_xvf
      


!*******************************************************************************
!  File type: xvZAf
!  ------------------
!  Read/write a time-stamped, formatted configuration file containing:
!    x,  v,  Z,  A
!  where Z is the charge (atomic number), and A is the mass (atomic weight). The
!  calling program should send the number of particles n in the file. But if this
!  is unknown, set n=0, and the reader should be able to figure it out. It will
!  simply read until reaching an end-of-file, and set n to the number of particles
!  read. This number is returned to the calling program.

      subroutine read_xvzaf(file,n,time8,x8,v8,z,a)
      implicit none
      character*(*)  file     !file to read  (in)
      integer        n        !number of particles  (in/out)
      real*8         time8    !time stamp  (out)
      real*8         x8(3,*)  !positions  (out)
      real*8         v8(3,*)  !velocities  (out)
      real*8         z(*)     !charges  (out)
      real*8         a(*)     !masses  (out)
      character*256  line     !for reading lines from file
      integer        i

      open(UNIT=13, FILE=file, STATUS='OLD', FORM='FORMATTED')
     !i=0
     !line(1:1)='#'
     !do while(line(1:1).eq.'#')
     !   read(13,'(a)') line
     !enddo
     !read(line,20010) n
      read(13,*) time8
      i=1
      do
        read(13,*,end=200) x8(1,i), x8(2,i), x8(3,i), v8(1,i), v8(2,i), v8(3,i), z(i), a(i)
       !write(12,10006) x8(1,i), x8(2,i), x8(3,i), v8(1,i), v8(2,i), v8(3,i), z(i), a(i)
       !write(12,10006) z(i), a(i)
        i=i+1
        if(n>0 .and. i==n+1) goto 200
      enddo
  200 continue
      i=i-1
      close(13)
      if(n==0) n=i
     !open(UNIT=12, File="fort.12", STATUS="UNKNOWN", FORM="FORMATTED")
     !do i=1,n
     !  write(12,10006) x8(1,i), x8(2,i), x8(3,i), v8(1,i), v8(2,i), v8(3,i), z(i), a(i)
     !enddo
     !close(12)
10006 format(3f13.6, 3x, 3f11.6, 3x, f8.2, f8.2)
10010 format(3es26.17, 3es26.17, f11.5, f11.5)
      return
      end


      subroutine write_xvzaf(file,n,time8,x8,v8,z,a)
      implicit none
      character*(*)  file     !file to write  (in)
      integer        n        !number of particles  (in)
      real*8         time8    !time stamp  (in)
      real*8         x8(3,n)  !positions  (in)
      real*8         v8(3,n)  !velocities  (in)
      real*8         z(n)     !charges  (in)
      real*8         a(n)     !masses  (in)
      integer        i
      open(UNIT=13, FILE=file, STATUS='UNKNOWN', FORM='FORMATTED')
      write(13) time8
      do i=1,n
        write(13,10010) x8(1,i), x8(2,i), x8(3,i), v8(1,i), v8(2,i), v8(3,i), z(i), a(i) 
      enddo
      close(13)
10010 format(3es24.17, 3es24.17, 1x,f7.2, 1x,f7.2)
      return
      end



!*******************************************************************************
!  File type: ZAf
!  ------------------
!  Read/write a formatted ion data file which contains the charge and mass of
!  each ion in a corresponding configuration file:
!    Z,  A
!  where Z is the charge (atomic number), and A is the mass (atomic weight).
!  This file type contains the number of particles n as the first line, so the
!  subroutine can easily determine how many data lines to read. If the calling
!  program sets n=0, the whole file will be read in and the number n returned.
!  If the calling program sets n>0, only the first n particles will be read in.
!  A ZAf file is meant to go along with a configuration file containing the
!  position and (possibly also) the velocity of each ion.


      subroutine read_zaf(file,n,z,a)
      implicit none
      character*(*)  file     !file to read (in)
      integer        n        !number of particles (out)
     !integer        z(*)     !charges (out)
     !integer        a(*)     !masses (out)
      real*8         z(*)     !charges (out)
      real*8         a(*)     !masses (out)
      character*256  line     !for reading lines from file
      integer        i,k

      open(UNIT=13, FILE=file, STATUS='OLD', FORM='FORMATTED')
      i=0
      line(1:1)='#'
      do while(line(1:1).eq.'#')
         read(13,'(a)') line
      enddo
      read(line,20010) k
      if(n==0) n=k
20010 format(i10)
      do i=1,n
       !read(13,20020) k, z(k), a(k)
        read(13,*) k, z(k), a(k)
      enddo
!20020 format(i10, i6, i6)
      close(13)
      return
      end

      subroutine write_zaf(file,n,z,a)
      implicit none
      character*(*)  file     !file to write  (in)
      integer        n        !number of particles  (in)
     !integer        z(n)     !charges  (in)
     !integer        a(n)     !masses  (in)
      real*8         z(n)     !charges  (in)
      real*8         a(n)     !masses  (in)
      integer        i
      open(UNIT=13, FILE=file, STATUS='UNKNOWN', FORM='FORMATTED')
      write(13,10010) n
10010 format(i10)
      do i=1,n
        write(13,10020) i, z(i), a(i)
      enddo
!10020 format(i10, i6, i6)
10020 format(i10,1x,f7.2,1x,f7.2)
      close(13)
      return
      end
