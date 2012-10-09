!*******************************************************************************
!    MD 6.3.0
! ---------------------------------------------------------------------
!    Copyright 2012, The Trustees of Indiana University
!    Authors:           Don Berry, Joe Hughto
!    Last modified by:  Joe Hughto, 2012-Oct-09
! ---------------------------------------------------------------------
!
!  Two-particle correlation function, g(r).
!
!DKB-todo : Subroutine g is not ready for use. It needs to be fixed so that
!DKB-todo+:    it can compute g(r) for a given pair of particle species.
!DKB-todo : Not sure rectangular boxes are handled correctly. Check statements
!DKB-todo+:    referring to xl and xlmin.
!  Compute the two-particle correlation function g(r) for a given configuration
!  and particle species.
!
!  Arguments:
!    integer  ig  --  (in) measurement group to which config belongs
!*******************************************************************************

      subroutine g(ig)
      use  md_types
      use  md_globals
      use  md_comm
      implicit real(dble) (a-h,o-z)
      include  'perf.h'
      include  'mpif.h'

      real(dble), allocatable ::   g2(:), gmean(:)
      real(dble), allocatable ::   ggx(:,:)
      integer    nthd     ! number of OpenMP threads
      save nthd

      integer  omp_get_max_threads, omp_get_thread_num

      if (nbin.eq.0) return

!-------------------------------------------------------------------------------
!  If ig=0, then just initalize arrays and return.
      nthd=1
   !$ nthd = omp_get_max_threads()
      if(ig.eq.0) then
        allocate(gg(nbin,ngroup,0:nthd-1))   ! for storing two-particle correlation
        allocate(cgg(ngroup))                !   data during the simulation
        gg  = 0.
        cgg = 0.
      endif

!-------------------------------------------------------------------------------
!  If ig>0, then take data. Variable gspec tells which particle species we are
!  calculating g(r) for. Array element cgg(ig) counts how many times g(r) has
!  been calculated for measurement group ig.
      if(ig.gt.0) then
         call starttimer()   !DKB-perf (g)
         !$omp parallel private(k,fac2,rr,nbinx,ii)
         k = 0
      !$ k = omp_get_thread_num()
         nbinx = nbin
         xlmin=min(xl(1),xl(2),xl(3))
         fac2=2.*float(nbin)/sqrt(3.)/xlmin
         !$omp do schedule(runtime)
         do i=myrank,n-2,nprocs
            do j=i+1,n-1
               call pdist(i,j,rr)
               ii = min( nbinx, int(fac2*rr+1) )
               gg(ii,ig,k)=gg(ii,ig,k)+1.
            enddo
         enddo
         !$omp end do
         !$omp end parallel
         cgg(ig)=cgg(ig)+1
         call stoptimer(1,t_g,ts_g,n_g)  !DKB-perf (g)
      endif

!-------------------------------------------------------------------------------
!  If ig= -1 then finish the calculation of g(r).
      if(ig.eq.-1) then

!  Accumulate the results of all OpenMP threads into the array section used by
!  thread 0.
         do k=1,nthd-1
            do j=1,nbin
               do i=1,ngroup
                  gg(j,i,0) = gg(j,i,0)+gg(j,i,k)
               enddo
            enddo
         enddo

!  If this is an MPI program, gather the results from each process.
         allocate(ggx(nbin,ngroup))
         call MPI_reduce(gg(:,:,0),ggx,ngroup*nbin,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                         0,MPI_COMM_WORLD,ierror)
         gg(:,:,0)=ggx
         deallocate(ggx)

!  Normalize. 
         if(myrank.eq.0) then
            ngspec=0
            do i=0,n-1
               ngspec=ngspec+1
            enddo
            xlmin=min(xl(1),xl(2),xl(3))
            dd=xlmin/float(nbin)*sqrt(3.)/2.
            fac=float(n)/(float(ngspec)**2*rho*dd*2.*PI)
            do j=1,nbin
               rr=(j-.5)*dd
               do i=1,ngroup
                  gg(j,i,0)=fac*gg(j,i,0)/(cgg(i)+1.e-30)/rr**2
               enddo
            enddo

!  For each bin j, compute mean and std.dev. of g(j,i) over all measurement groups i.
            open(UNIT=7,FILE='g.out'//suffix,STATUS='UNKNOWN')
            allocate(gmean(nbin))
            allocate(g2(nbin))
            do j=1,nbin
               rr=(j-.5)*dd
               gmean(j)=0.
               g2(j)=0.
               do i=1,ngroup
                  gmean(j)=gmean(j)+gg(j,i,0)
                  g2(j)=g2(j)+gg(j,i,0)**2
               enddo
               gmean(j)=gmean(j)/float(ngroup)
               g2(j)=g2(j)/float(ngroup)
               g2(j)=sqrt(abs(g2(j)-gmean(j)**2)/ngroup)
               write(7,10000) rr,gmean(j),g2(j)
            enddo
            deallocate(gmean)
            deallocate(g2)
            close(7)
         endif
      endif
10000 format(f15.8,2x,1pe15.8,2x,1pe15.8)

!-------------------------------------------------------------------------------
      if(ig.eq.-2) then
        deallocate(cgg)
        deallocate(gg)
      endif

      return
      end subroutine g
