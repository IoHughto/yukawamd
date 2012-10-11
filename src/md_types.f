!*******************************************************************************
!     MD 6.2.0
!  ---------------------------------------------------------------------
!     Copyright 2012, The Trustees of Indiana University
!     Author:            Don Berry
!     Last modified by:  Don Berry, 2012-May-23
!  ---------------------------------------------------------------------
!
!  08-Mar-2007 (Don Berry) -- Changed components of type species so that charge
!    and mass are now real numbers instead of integers.
!
!*******************************************************************************


module md_types

   integer, parameter  :: dble=kind(1.0d0)

!  Create a data type to define ion species.
   type species
      integer    ::  n = 0     ! number of ions
      real(dble) ::  z = 0.    ! charge number
      real(dble) ::  a = 0.    ! mass number
   end type species


   interface operator (.eq.)
      module procedure species_eq
   end interface

   interface operator (.ne.)
      module procedure species_neq
   end interface


!===============================================================================
   CONTAINS

   logical function species_eq(s1,s2)
      implicit none
      type(species), intent(in)  ::  s1
      type(species), intent(in)  ::  s2
      species_eq = (s1%z .eq. s2%z) .and. (s1%a .eq. s2%a)
      return
   end function species_eq

   logical function species_neq(s1,s2)
      implicit none
      type(species), intent(in)  ::  s1
      type(species), intent(in)  ::  s2
      species_neq = (s1%z .ne. s2%z) .or. (s1%a .ne. s2%a)
      return
   end function species_neq


end module md_types
