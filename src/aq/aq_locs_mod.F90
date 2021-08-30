! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

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
