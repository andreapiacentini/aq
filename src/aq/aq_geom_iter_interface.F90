! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2021-2022 CERFACS.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_geom_iter_interface

use iso_c_binding
use kinds
use aq_geom_iter_mod
use aq_geom_interface
use aq_geom_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup geometry iterator
subroutine aq_geom_iter_setup_c(c_key_self,c_key_geom,c_index) bind(c,name='aq_geom_iter_setup_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry iterator
integer(c_int),intent(in) :: c_key_geom    !< Geometry
integer(c_int),intent(in) :: c_index       !< Index

! Local variables
type(aq_geom_iter),pointer :: self
type(aq_geom),pointer :: geom

! Interface
call aq_geom_iter_registry%init()
call aq_geom_iter_registry%add(c_key_self)
call aq_geom_iter_registry%get(c_key_self,self)
call aq_geom_registry%get(c_key_geom,geom)

! Call Fortran
call aq_geom_iter_setup(self,geom,c_index)

end subroutine aq_geom_iter_setup_c
! ------------------------------------------------------------------------------
!> Clone geometry iterator
subroutine aq_geom_iter_clone_c(c_key_self,c_key_other) bind(c,name='aq_geom_iter_clone_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry iterator
integer(c_int),intent(in) :: c_key_other   !< Other geometry iterator

! Local variables
type(aq_geom_iter),pointer :: self,other

! Interface
call aq_geom_iter_registry%get(c_key_other,other)
call aq_geom_iter_registry%init()
call aq_geom_iter_registry%add(c_key_self)
call aq_geom_iter_registry%get(c_key_self,self)

! Call Fortran
call aq_geom_iter_clone(self,other)

end subroutine aq_geom_iter_clone_c
! ------------------------------------------------------------------------------
!> Delete geometry iterator
subroutine aq_geom_iter_delete_c(c_key_self) bind(c,name='aq_geom_iter_delete_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry iterator

! Clear interface
call aq_geom_iter_registry%remove(c_key_self)

end subroutine aq_geom_iter_delete_c
! ------------------------------------------------------------------------------
!> Check geometry iterator equality
subroutine aq_geom_iter_equals_c(c_key_self,c_key_other,c_equals) bind(c,name='aq_geom_iter_equals_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry iterator
integer(c_int),intent(in) :: c_key_other   !< Other geometry iterator
integer(c_int),intent(inout) :: c_equals   !< Equality flag

! Local variables
type(aq_geom_iter),pointer :: self,other

! Interface
call aq_geom_iter_registry%get(c_key_self,self)
call aq_geom_iter_registry%get(c_key_other,other)

! Call Fortran
call aq_geom_iter_equals(self,other,c_equals)

end subroutine aq_geom_iter_equals_c
! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
subroutine aq_geom_iter_current_c(c_key_self,c_lat,c_lon) bind(c,name='aq_geom_iter_current_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Geometry iterator
real(c_double),intent(inout) :: c_lat   !< Latitude
real(c_double),intent(inout) :: c_lon   !< Longitude

! Local variables
type(aq_geom_iter),pointer :: self

! Interface
call aq_geom_iter_registry%get(c_key_self,self)

! Call Fortran
call aq_geom_iter_current(self,c_lon,c_lat)

end subroutine aq_geom_iter_current_c
! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
subroutine aq_geom_iter_next_c(c_key_self) bind(c,name='aq_geom_iter_next_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Geometry iterator

! Local variables
type(aq_geom_iter),pointer :: self

! Interface
call aq_geom_iter_registry%get(c_key_self,self)

! Call Fortran
call aq_geom_iter_next(self)

end subroutine aq_geom_iter_next_c
! ------------------------------------------------------------------------------
end module aq_geom_iter_interface
