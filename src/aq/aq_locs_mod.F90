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

module aq_locs_mod

use atlas_module, only: atlas_field
use iso_c_binding
use kinds
use datetime_mod

implicit none
private
public :: aq_locs
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 1 !< Random seed (for reproducibility)

type :: aq_locs
  type(c_ptr), private :: ptr
contains
  procedure, public :: nlocs => locs_nlocs
  procedure, public :: lonlat => locs_lonlat
  procedure, public :: altitude => locs_altitude
  procedure, public :: times => locs_times
end type aq_locs

interface aq_locs
   module procedure ctor_from_ptr
end interface

#include "aq_locs_interface.f"

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
function ctor_from_ptr(ptr) result(this)
  type(aq_locs)    :: this
  type(c_ptr), intent(in) :: ptr

  this%ptr = ptr
end function ctor_from_ptr

! ------------------------------------------------------------------------------
function locs_nlocs(self)
  class(aq_locs), intent(in) :: self  ! locations object
  integer :: locs_nlocs
  locs_nlocs = aq_locs_nlocs_c(self%ptr)
end function

! ------------------------------------------------------------------------------
function locs_lonlat(self) result(field)
  class(aq_locs), intent(in) :: self  ! locations object
  type(atlas_Field) :: field
  field = atlas_Field(aq_locs_lonlat_c(self%ptr))
  call field%return()
end function

! ------------------------------------------------------------------------------
function locs_altitude(self) result(field)
  class(aq_locs), intent(in) :: self  ! locations object
  type(atlas_Field) :: field
  field = atlas_Field(aq_locs_altitude_c(self%ptr))
  call field%return()
end function

! ------------------------------------------------------------------------------
function locs_times(self, jj) result(dt)
  class(aq_locs), intent(in) :: self  ! locations object
  integer, intent(in) :: jj
  type(datetime) :: dt
  integer(c_size_t) :: idx
  idx = jj - 1
  call c_f_datetime(aq_locs_times_c(self%ptr, idx), dt)
end function

! ------------------------------------------------------------------------------
end module aq_locs_mod
