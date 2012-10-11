!*******************************************************************************
!
!    MD 6.2.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry
!    Last modified by:  Don Berry, 2012-Jul-12
! ---------------------------------------------------------------------
!
!*******************************************************************************


!*******************************************************************************
! Checkpoint the run. This is not a true checkpoint, as it does not save the
! run parameters, nor the statistics arrays. It saves the positions and veloci-
! ties in double precision binary. Actual output is done in subroutine write_xvb.
!DKB-todo : Fix checkpoint subroutine so it really does a complete checkpointing of the job.
!
      subroutine checkpoint
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o,z)
      include 'mpif.h'

      character*256   mdoutfile
      character*256   mdckptfile

!  Only MPI process 0 writes to the checkpoint file.
      if(myrank.eq.0) then
        write(mdckptfile,110) nint(time/1000000.d0), nint(mod(time,1000000.d0)), '.ckpt'
        write(mdoutfile,110) nint(time/1000000.d0), nint(mod(time,1000000.d0)), '.xv8b'
  110   format('md.',i5.5,i6.6,a)
        call write_xvb(mdoutfile,'xv8b  ',.false.,.true.)
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierror)
      return
      end subroutine checkpoint


!*******************************************************************************
! Save configuration to a file. Depending on input arguments, this routine saves
! a configuration (positions only) or a phase space point (positions and veloci-
! ties). Argument filetype determines what to write out, and in what precision.
! Argument xappend determines whether to write each configuration to a separate
! file, or append configurations to a single trajectory file. Actual output is
! done in subroutine write_xvb.
!
! Arguments:
!
!   character*6  xfiletype  -- specifies what to save in configuration file
!       'x4b'  = real positions
!       'xv4b' = real positions, real velocities.
!       'x8b'  = real(dble) positions
!       'xv8b' = real(dble) positions, real(dble) velocities.
!
!   logical      xappend  -- 
!       .true . = append all configs to a single md.traj file
!       .false. = write each config to its own md.out file

      subroutine save_config(xfiletype,xappend)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include 'mpif.h'

      character*6     xfiletype
      logical         xappend

      character*256   mdoutfile

! Only MPI process 0 writes to the file.
      if(myrank.eq.0) then
        if(xappend) then
          write(mdoutfile,100) nint(tend/1000000.d0), nint(mod(tend,1000000.d0)), trim(xfiletype)
  100     format('md.traj.',i5.5,i6.6,'.',a)
        else
          write(mdoutfile,110) nint(time/1000000.d0), nint(mod(time,1000000.d0)), trim(xfiletype)
  110     format('md.',i5.5,i6.6,'.',a)
        endif
        !-----------------------------------------------------------------------
        select case(trim(xfiletype))

          case('onp','md1','xfa','xvfa')
          !------------------
            call write_xvfa(start,xfiletype,.true.,.true.)

          case('md0','x4a','xv4a','x8a','xv8a','xva8a')
          !----------------------------------------------
            call write_xva(mdoutfile,xfiletype,xappend,.true.)

          case('xvZAfa')
          !------------------
            call write_xvfa(start,xfiletype,xappend,.true.)

          case('x4b','x8b','xv4b','xv8b')
          !--------------------------------
            call write_xvb(mdoutfile,xfiletype,xappend,.true.)

          case('xfb','xvfb')
          !-------------------
            call write_xvfb(start,xfiletype,xappend,.true.)

          case('xZAfb','xvZAfb')
          !------------------
            call write_xvfb(start,xfiletype,xappend,.true.)

        end select
        !-----------------------------------------------------------------------
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierror)
      return
      end subroutine save_config


!*******************************************************************************
! Dump x, v, and a arrays to a file. This subroutine is meant for debugging.
! Normally one does not need to save the accelerations. Dumped data files are
! suffixed with the current simulation time, and an additional sequence number,
! so that up to 100 dump files may be saved during each time step.
!
      subroutine dump_data
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o,z)
      include 'mpif.h'

      character*256   md_dump_file
      real(dble), save :: xtime=0.0d0
      integer, save    :: nseq=0

      if(xtime.ne.time) then
         xtime=time
         nseq=0
      endif

! Only MPI process 0 writes to the md_dump file.
      if(myrank.eq.0) then
        write(md_dump_file,110) int(time/1000000.d0), int(mod(time,1000000.d0)), nseq, 'xva8b'
  110   format('md.dump.',i5.5,i6.6,'.',i2.2,'.',a)
        call write_xvb(md_dump_file,'xva8b ',.false.,.true.)
      endif

      nseq=nseq+1
      call MPI_barrier(MPI_COMM_WORLD,ierror)

      return
      end subroutine dump_data



!*******************************************************************************
! This subroutine should be called for both normal completion of the program,
! as well as when a fatal error occurs.
!
      subroutine md_exit(icode,string)
      use md_globals
      use md_comm
      implicit real(dble) (a-h,o-z)

      integer        icode
      character*(*)  string

      if(myrank.eq.0) then
        write(6,10010) icode,string
10010   format('exit code = ',i4,2x,a)
      endif
      
! Deallocate particle arrays.
      deallocate(a, v, vold, x)
      if(allocated(aii)) deallocate(aii)      !mass numbers A
      if(allocated(zii)) deallocate(zii)      !charge numbers Z

! Deallocate statistics arrays.
      deallocate(eva,ev2a,dev)
      deallocate(pa,p2a,dp)

      call MPI_finalize(ierror)

      stop
      end subroutine md_exit
