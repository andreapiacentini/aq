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

module aq_interpolator_interface

   use iso_c_binding
   use kinds
   use oops_variables_mod
   use aq_interpolator_mod
   use aq_fields_mod
   use aq_fields_interface, only: aq_fields_registry
   use aq_geom_mod
   use aq_geom_interface, only: aq_geom_registry

   implicit none
   public :: aq_interpolator_registry

#define LISTED_TYPE aq_interpolator

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

   !> Global registry
   type(registry_t) :: aq_interpolator_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

   !>  Create an interpolator from a geometry and obs locations lons and lats
   subroutine aq_interpolator_create_c(c_key_self, c_key_geom, c_nlocs, c_lats, c_lons) bind(c,name='aq_interpolator_create_f90')

      implicit none

      ! Passed variables
      integer(c_int),intent(inout) :: c_key_self       !< Interpolator
      integer(c_int),intent(in)    :: c_key_geom       !< Geometry
      integer(c_int), intent(in)   :: c_nlocs
      real(c_double), intent(in)   :: c_lats(c_nlocs)
      real(c_double), intent(in)   :: c_lons(c_nlocs)

      ! Local variables
      type(aq_interpolator),pointer :: self
      type(aq_geom),pointer :: geom

      ! Interface
      call aq_interpolator_registry%init()
      call aq_interpolator_registry%add(c_key_self)
      call aq_interpolator_registry%get(c_key_self,self)
      call aq_geom_registry%get(c_key_geom,geom)

      ! Call Fortran
      call self%create(geom,c_nlocs,c_lats,c_lons)

   end subroutine aq_interpolator_create_c
   ! ------------------------------------------------------------------------------
   !> Delete interpolator
   subroutine aq_interpolator_delete_c(c_key_self) bind(c,name='aq_interpolator_delete_f90')

      implicit none

      ! Passed variables
      integer(c_int),intent(inout) :: c_key_self !< Interpolator

      ! Local variables
      type(aq_interpolator),pointer :: self

      ! Interface
      call aq_interpolator_registry%get(c_key_self,self)

      ! Call Fortran
      call self%delete()

      ! Clear interface
      call aq_interpolator_registry%remove(c_key_self)

   end subroutine aq_interpolator_delete_c

   ! ------------------------------------------------------------------------------
   subroutine aq_interpolator_apply_c(c_key_self, c_key_fld, c_vars, c_nvals, c_mask, c_vals) &
 & bind (c,name='aq_interpolator_apply_f90')

      implicit none
      ! Passed variables
      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(in)    :: c_key_fld
      type(c_ptr),value,intent(in)  :: c_vars
      integer(c_int), intent(in)    :: c_nvals
      integer(c_int), intent(in)    :: c_mask(c_nvals)
      real(c_double), intent(inout) :: c_vals(c_nvals)

      ! Local variables
      type(aq_interpolator),pointer :: self
      type(aq_fields),pointer :: fld
      type(oops_variables) :: vars

      call aq_interpolator_registry%get(c_key_self, self)
      call aq_fields_registry%get(c_key_fld, fld)
      vars = oops_variables(c_vars)

      call self%apply(fld, vars, c_mask, c_vals)

   end subroutine aq_interpolator_apply_c
   ! ------------------------------------------------------------------------------
   subroutine aq_interpolator_applyAD_c(c_key_self, c_key_fld, c_vars, c_nvals, c_mask, c_vals) &
 & bind (c,name='aq_interpolator_applyAD_f90')

      implicit none
      ! Passed variables
      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(in)    :: c_key_fld
      type(c_ptr),value,intent(in)  :: c_vars
      integer(c_int), intent(in)    :: c_nvals
      integer(c_int), intent(in)    :: c_mask(c_nvals)
      real(c_double), intent(in)    :: c_vals(c_nvals)

      ! Local variables
      type(aq_interpolator),pointer :: self
      type(aq_fields),pointer :: fld
      type(oops_variables) :: vars

      call aq_interpolator_registry%get(c_key_self, self)
      call aq_fields_registry%get(c_key_fld, fld)
      vars = oops_variables(c_vars)

      call self%applyAD(fld, vars, c_mask, c_vals)

   end subroutine aq_interpolator_applyAD_c

end module aq_interpolator_interface
