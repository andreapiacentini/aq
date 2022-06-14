! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_geom_interface

use atlas_module, only: atlas_fieldset,atlas_structuredgrid,atlas_functionspace_structuredcolumns
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use fckit_mpi_module,only: fckit_mpi_comm
use kinds
use iso_c_binding
!AP use aq_projection_mod
use aq_geom_mod
use aq_constants_mod, only : aq_strlen

implicit none

private
public :: aq_geom_registry

#define LISTED_TYPE aq_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!> Setup geometry
subroutine aq_geom_setup_c(c_key_self,c_conf,c_comm,c_agrid,c_afs,c_afs_surf) bind(c,name='aq_geom_setup_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
type(c_ptr),intent(in),value :: c_conf     !< Configuration
type(c_ptr),intent(in),value :: c_comm
type(c_ptr),intent(in),value :: c_agrid
type(c_ptr),intent(in),value :: c_afs !< ATLAS function space pointer
type(c_ptr),intent(in),value :: c_afs_surf !< ATLAS function space pointer

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_mpi_comm)      :: f_comm
type(aq_geom),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)
call aq_geom_registry%init()
call aq_geom_registry%add(c_key_self)
call aq_geom_registry%get(c_key_self,self)
self%grid = atlas_structuredgrid(c_agrid)
self%fs = atlas_functionspace_structuredcolumns(c_afs)
self%fs_surf = atlas_functionspace_structuredcolumns(c_afs_surf)

! Call Fortran
call self%create(f_conf,f_comm)

end subroutine aq_geom_setup_c
! ------------------------------------------------------------------------------
!> Fill geometry extra fields
subroutine aq_geom_fill_extra_fields_c(c_key_self,c_afieldset) &
 & bind(c,name='aq_geom_fill_extra_fields_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self     !< Geometry
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(aq_geom),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call aq_geom_registry%get(c_key_self,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call self%fill_extra_fields(afieldset)

end subroutine aq_geom_fill_extra_fields_c
! ------------------------------------------------------------------------------
!> Clone geometry
subroutine aq_geom_clone_c(c_key_self,c_key_other) bind(c,name='aq_geom_clone_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
integer(c_int),intent(in) :: c_key_other   !< Other geometry

! Local variables
type(aq_geom),pointer :: self,other

! Interface
call aq_geom_registry%get(c_key_other,other)
call aq_geom_registry%init()
call aq_geom_registry%add(c_key_self)
call aq_geom_registry%get(c_key_self ,self )

! Call Fortran
call self%clone(other)

end subroutine aq_geom_clone_c
! ------------------------------------------------------------------------------
!> Delete geometry
subroutine aq_geom_delete_c(c_key_self) bind(c,name='aq_geom_delete_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry

! Local variables
type(aq_geom),pointer :: self

! Interface
call aq_geom_registry%get(c_key_self,self)

! Call Fortran
call self%delete()

! Clear interface
call aq_geom_registry%remove(c_key_self)

end subroutine aq_geom_delete_c
! ------------------------------------------------------------------------------
!> Get geometry info
subroutine aq_geom_info_c(c_key_self,c_nx,c_ny,c_nz,c_deltax,c_deltay,c_mod_levels, &
   & c_orientation_ptr, c_domname_ptr, c_model_ptr) bind(c,name='aq_geom_info_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Geometry
integer(c_int),intent(inout) :: c_nx     !< Number of points in the zonal direction
integer(c_int),intent(inout) :: c_ny     !< Number of points in the meridional direction
integer(c_int),intent(inout) :: c_nz     !< Number of vertical levels
real(c_double),intent(inout) :: c_deltax !< Zonal cell size
real(c_double),intent(inout) :: c_deltay !< Meridional cell size
integer(c_int),intent(inout) :: c_mod_levels !< Number of model levels
character(kind=c_char), target, intent(inout) :: c_orientation_ptr(*)
character(kind=c_char), target, intent(inout) :: c_domname_ptr(*)
character(kind=c_char), target, intent(inout) :: c_model_ptr(*)

! Local variables
integer :: i, length
character(kind=c_char), pointer :: c_orientation(:)
character(kind=c_char), pointer :: c_domname(:)
character(kind=c_char), pointer :: c_model(:)
type(aq_geom),pointer :: self
real(c_double) :: xmin, ymin
integer(c_int) :: halo
character(len=:), allocatable :: domname, orientation, model

! Interface
call aq_geom_registry%get(c_key_self,self)

! Call Fortran
call self%info(c_nx,c_ny,c_nz,c_deltax,c_deltay,xmin, ymin, c_mod_levels, halo, domname, orientation, model)
length = len_trim(orientation)
call c_f_pointer(c_loc(c_orientation_ptr), c_orientation, [aq_strlen])
do i = 1, length
   c_orientation(i)=orientation(i:i)
end do
c_orientation(length+1) = C_NULL_CHAR
length = len_trim(domname)
call c_f_pointer(c_loc(c_domname_ptr), c_domname, [aq_strlen])
do i = 1, length
   c_domname(i)=domname(i:i)
end do
c_domname(length+1) = C_NULL_CHAR
length = len_trim(model)
call c_f_pointer(c_loc(c_model_ptr), c_model, [aq_strlen])
do i = 1, length
   c_model(i)=model(i:i)
end do
c_model(length+1) = C_NULL_CHAR

end subroutine aq_geom_info_c

! ------------------------------------------------------------------------------
end module aq_geom_interface
