!
!  This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
!
!  (C) Copyright 2022 CERFACS
!
!  AQ is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  any later version.
!
!  AQ is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  A copy of the GNU Lesser General Public License is distributed
!  along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
!

module aq_geom_iter_mod

use iso_c_binding
use kinds
use aq_geom_mod

implicit none

private
public :: aq_geom_iter
public :: aq_geom_iter_registry
public :: aq_geom_iter_setup,aq_geom_iter_clone,aq_geom_iter_equals,aq_geom_iter_current,aq_geom_iter_next
! ------------------------------------------------------------------------------
type :: aq_geom_iter
  type(aq_geom),pointer :: geom => null() !< Geometry
  integer :: ilon = 1                     !< Longitude index
  integer :: ilat = 1                     !< Latitude index
end type aq_geom_iter

#define LISTED_TYPE aq_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_geom_iter_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup for the AQ model's geometry iterator
subroutine aq_geom_iter_setup(self,geom,ind)

! Passed variables
type(aq_geom_iter),intent(inout) :: self !< Geometry iterator
type(aq_geom),pointer,intent(in) :: geom !< Geometry
integer,intent(in) :: ind                !< Index

! Associate geometry
self%geom => geom

! Define ilon/ilat
self%ilat = (ind-1)/geom%nx+1
self%ilon = ind-(self%ilat-1)*geom%nx

end subroutine aq_geom_iter_setup
! ------------------------------------------------------------------------------
!> Clone for the AQ model's geometry iterator
subroutine aq_geom_iter_clone(self,other)

! Passed variables
type(aq_geom_iter),intent(inout) :: self !< Geometry iterator
type(aq_geom_iter),intent(in) :: other   !< Other geometry iterator

! Associate geometry
self%geom => other%geom

! Copy ilon/ilat
self%ilon = other%ilon
self%ilat = other%ilat

end subroutine aq_geom_iter_clone
! ------------------------------------------------------------------------------
!> Check for the AQ model's geometry iterator equality
subroutine aq_geom_iter_equals(self,other,equals)

! Passed variables
type(aq_geom_iter),intent(in) :: self  !< Geometry iterator
type(aq_geom_iter),intent(in) :: other !< Other geometry iterator
integer,intent(out) :: equals          !< Equality flag

! Initialization
equals = 0

! Check equality
if (associated(self%geom,other%geom).and.(self%ilon==other%ilon).and.(self%ilat==other%ilat)) equals = 1

end subroutine aq_geom_iter_equals
! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
subroutine aq_geom_iter_current(self,lat,lon)

! Passed variables
type(aq_geom_iter),intent(in) :: self !< Geometry iterator
real(kind_real),intent(out) :: lat    !< Latitude
real(kind_real),intent(out) :: lon    !< Longitude

! Check ilon/ilat
if (self%ilon*self%ilat>self%geom%nx*self%geom%ny) call abor1_ftn('aq_geom_iter_current: iterator out of bounds')

! Get lat/lon
!APlat = self%geom%lat(self%ilon,self%ilat)
!APlon = self%geom%lon(self%ilon,self%ilat)
lat = 0.0_kind_real
lon = 0.0_kind_real

end subroutine aq_geom_iter_current
! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
subroutine aq_geom_iter_next(self)

! Passed variables
type(aq_geom_iter),intent(inout) :: self !< Geometry iterator

! Update ilon/ilat
if (self%ilon==self%geom%nx) then
  self%ilon = 1
  self%ilat = self%ilat+1
else
  self%ilon = self%ilon+1
endif

end subroutine aq_geom_iter_next
! ------------------------------------------------------------------------------
end module aq_geom_iter_mod
