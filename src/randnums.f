


!*******************************************************************************
!  Returns a normally distributed random number with zero mean and unit
!  variance. Uses ran1(idum) as the source of uniform random numbers.
!
      DOUBLE PRECISION FUNCTION gasdev(idum)
      INTEGER idum

      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if(iset.eq.0) then
   10   v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.) goto 10
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software #Q2Z3*(121~4:.


!*******************************************************************************
!  Dummy subroutine to use either ran11 or rannyu random number generator.
!
      function ran1(idum)
      implicit double precision (a-h,o-z)
      character*6  ranname(0:4)
      data  ranname/'N/A','ran11','rannyu','ran3','grnd'/
      common /rndtype/irnd,ranname
      if (irnd.eq.1) then
         ran1=ran11(idum)
      else if (irnd.eq.2) then
         ran1=rannyu()
      else if (irnd.eq.3) then
         ran1=ran3(idum)
      else
         ran1=grnd()
      endif
      return
      end
 

!*******************************************************************************
!
!  Changes:
!  2005-Mar-16: (D.Berry) - Convert to Fortran 95 syntax.

      FUNCTION RAN11(IDUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, SAVE      :: R(97)
      INTEGER, SAVE               :: IX1, IX2, IX3, J, IFF
      INTEGER, PARAMETER          :: M1=259200, IA1=7141, IC1=54773
      DOUBLE PRECISION, PARAMETER :: RM1=3.8580247E-6
      INTEGER, PARAMETER          :: M2=134456, IA2=8121, IC2=28411
      DOUBLE PRECISION, PARAMETER :: RM2=7.4373773E-6
      INTEGER, PARAMETER          :: M3=243000, IA3=4561, IC3=51349
      DATA IFF /0/
      IF (IDUM<0 .OR. IFF==0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
        ENDDO
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J>97 .OR. J<1) THEN
        WRITE(6,*) 'Error in ran11: j=', j, '     Max value is 97.'
        STOP
      ENDIF
      RAN11=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END



!*******************************************************************************
!  This a linear congruential 48 bit Random number generator (RNG) from Robert
!  Edwards (NYU) with increment of 1. The seed is kept in a  integer array of
!  length 4. The values in each element must be less than 4096.  By default this
!  RNG is in double precision. You could change the double's to real's and make
!  it single precision.
!
!  Related subroutines:
!        setrn(lseed)  --  sets the seed
!        savern(lseed) --  extracts the seed
!
!  More notes about rannyu:
!     linear congruential with modulus m = 2**48, increment c = 1,
!     multiplier a = (2**36)*m1 + (2**24)*m2 + (2**12)*m3 + m4. 
!     The multiplier is stored in common (see subroutine setrn)
!     and is set to a = 31167285 (recommended by Knuth, vol. 2,
!     2nd ed., p. 102).

      function rannyu()
      double precision rannyu, twom12
      parameter (twom12 = 1/4096.0e0)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4

!     i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1
      i1 = l1*m4 + l2*m3 + l3
!     i2 = l2*m4 + l3*m3 + l4*m2
      i2 = l2*m4 + l3*m3 + l4
      i3 = l3*m4 + l4*m3
      i4 = l4*m4  +  1
      l4 = mod(i4, 4096)
      i3 = i3 + i4/4096
      l3 = mod(i3, 4096)
      i2 = i2 + i3/4096
      l2 = mod(i2, 4096)
      l1 = mod(i1 + i2/4096, 4096)
      rannyu = twom12*(float(l1) + twom12*(float(l2) + twom12*(float(l3)+    &
               twom12*(float(l4)))))
      return
      end


!*******************************************************************************
!  This function sets the seed for Robert Edwards' (NYU) linear congruential
!  48-bit random number generator.
!
!  Multiplier is 31167285 = (2**24) + 3513*(2**12) + 821.
!     (Recommended by Knuth, vol. 2, 2nd ed., p. 102.)
!  Generator is linear congruential with odd increment and maximal period, so
!  seed is unrestricted: it can be either even or odd.

      subroutine setrn(iseed)
      common /rnyucm/ m(4),l(4)
      INTEGER iseed(4)
 
!cc   data m /   0,   1,3513, 821/
!cc   data l /   0,   0,   0,   1/

      m(1) = 0
      m(2) = 1
      m(3) = 3513
      m(4) = 821

      do 10 i = 1, 4
         l(i) = iseed(i)
10    continue
      return
      end


!*******************************************************************************
!  This function extracts the seed for Robert Edwards' (NYU) linear congruential
!  48-bit random number generator.
!
      subroutine savern(iseed)
      common /rnyucm/ m(4),l(4)
      INTEGER iseed(4)
 
      do 10 i = 1, 4
         iseed(i) = l(i)
10    continue

      return
      end


!*******************************************************************************
!  Random number generator ran3 from Numerical Recipes. It returns a uniform
!  random deviate between 0.0 and 1.0. Set idum to any negative value to
!  initialize or reinitialize the sequence.
!
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER  i,iff,ii,inext,inextp,k
      INTEGER  mj,mk,ma(55)
!     REAL     mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/

      if(idum.lt.0 .or. iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ) mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ) ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ) mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software #Q2Z3*(121~4:.



!*******************************************************************************
      subroutine sgrnd(seed)
      implicit integer(a-z)
!
! Period parameters
      parameter(N     =  624)
      dimension mt(0:N-1)      !the array for the state vector
      common /block/mti,mt
      save   /block/

!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
      mt(0)= iand(seed,-1)
      do mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
      enddo

      return
      end


!*******************************************************************************
      double precision function grnd()
      implicit integer(a-z)

! Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)    !constant vector a
      parameter(UMASK = -2147483648)    !most significant w-r bits
      parameter(LMASK =  2147483647)    !least significant r bits
! Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)

      dimension mt(0:N-1)       !the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/            !mti==N+1 means mt[N] is not initialized

      dimension mag01(0:1)
      save mag01
      data mag01/0, MATA/       !mag01(x) = x * MATA for x=0,1

      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)

      if(mti.ge.N) then         !generate N words at one time.
        if(mti.eq.N+1) then     !if sgrnd() has not been called,
          call sgrnd(4357)      !a default initial seed is used
        endif

        do kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        do kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif

      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))

      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif

      return
      end
