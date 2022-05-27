! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_fields_interface

use atlas_module, only: atlas_fieldset, atlas_field
use datetime_mod
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use oops_variables_mod
use aq_fields_mod
use aq_geom_mod
use aq_geom_interface, only: aq_geom_registry
use aq_geom_iter_mod
!AP use aq_locs_mod

implicit none

private
public :: aq_fields_registry

#define LISTED_TYPE aq_fields

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_fields_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!> Create fields from geometry and variables
subroutine aq_fields_create_c(c_key_self,c_key_geom,c_vars) bind(c,name='aq_fields_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self       !< Fields
integer(c_int),intent(in) :: c_key_geom          !< Geometry
type(c_ptr),value,intent(in) :: c_vars           !< List of variables

! Local variables
type(aq_fields),pointer :: self
type(aq_geom),pointer :: geom
type(oops_variables) :: vars

! Interface
call aq_fields_registry%init()
call aq_fields_registry%add(c_key_self)
call aq_fields_registry%get(c_key_self,self)
call aq_geom_registry%get(c_key_geom,geom)
vars = oops_variables(c_vars)

! Call Fortran
call self%create(geom,vars)

end subroutine aq_fields_create_c
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine aq_fields_create_from_other_c(c_key_self,c_key_other) bind(c,name='aq_fields_create_from_other_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self  !< Fields
integer(c_int),intent(in)    :: c_key_other !< Other fields

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: other

! Interface
call aq_fields_registry%get(c_key_other,other)
call aq_fields_registry%init()
call aq_fields_registry%add(c_key_self)
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%create_from(other)

end subroutine aq_fields_create_from_other_c
! ------------------------------------------------------------------------------
!> Delete fields
subroutine aq_fields_delete_c(c_key_self) bind(c,name='aq_fields_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Fields

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%delete()

! Clear interface
call aq_fields_registry%remove(c_key_self)

end subroutine aq_fields_delete_c
! ------------------------------------------------------------------------------
!> Set fields to zero
subroutine aq_fields_zero_c(c_key_self) bind(c,name='aq_fields_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%zero()

end subroutine aq_fields_zero_c
! ------------------------------------------------------------------------------
!> Set fields to ones
subroutine aq_fields_ones_c(c_key_self) bind(c,name='aq_fields_ones_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%ones()

end subroutine aq_fields_ones_c
! ------------------------------------------------------------------------------
!> Set fields to Diracs
subroutine aq_fields_dirac_c(c_key_self,c_conf) bind(c,name='aq_fields_dirac_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_fields),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%dirac(f_conf)

end subroutine aq_fields_dirac_c
! ------------------------------------------------------------------------------
!> Generate random fields
subroutine aq_fields_random_c(c_key_self) bind(c,name='aq_fields_random_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%random()

end subroutine aq_fields_random_c
! ------------------------------------------------------------------------------
!> Copy fields
subroutine aq_fields_copy_c(c_key_self,c_key_other) bind(c,name='aq_fields_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Fields
integer(c_int),intent(in) :: c_key_other !< Other fields

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: other

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_other,other)

! Call Fortran
call self%copy(other)

end subroutine aq_fields_copy_c
! ------------------------------------------------------------------------------
!> Add fields
subroutine aq_fields_self_add_c(c_key_self,c_key_rhs) bind(c,name='aq_fields_self_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: rhs

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%self_add(rhs)

end subroutine aq_fields_self_add_c
! ------------------------------------------------------------------------------
!> Subtract fields
subroutine aq_fields_self_sub_c(c_key_self,c_key_rhs) bind(c,name='aq_fields_self_sub_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: rhs

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%self_sub(rhs)

end subroutine aq_fields_self_sub_c
! ------------------------------------------------------------------------------
!> Multiply fields by a scalar
subroutine aq_fields_self_mul_c(c_key_self,c_zz) bind(c,name='aq_fields_self_mul_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%self_mul(c_zz)

end subroutine aq_fields_self_mul_c
! ------------------------------------------------------------------------------
!> Apply axpy operator to fields
subroutine aq_fields_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='aq_fields_axpy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: rhs

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%axpy(c_zz,rhs)

end subroutine aq_fields_axpy_c
! ------------------------------------------------------------------------------
!> Schur product of fields
subroutine aq_fields_self_schur_c(c_key_self,c_key_rhs) bind(c,name='aq_fields_self_schur_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: rhs

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%schur(rhs)

end subroutine aq_fields_self_schur_c
! ------------------------------------------------------------------------------
!> Compute dot product for fields
subroutine aq_fields_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='aq_fields_dot_prod_f90')

implicit none

! Passed variables
integer(c_int),intent(in)    :: c_key_fld1 !< First fields
integer(c_int),intent(in)    :: c_key_fld2 !< Second fields
real(c_double),intent(inout) :: c_prod     !< Dot product

! Local variables
type(aq_fields),pointer :: fld1,fld2

! Interface
call aq_fields_registry%get(c_key_fld1,fld1)
call aq_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call fld1%dot_prod_with(fld2,c_prod)

end subroutine aq_fields_dot_prod_c
! ------------------------------------------------------------------------------
!> Add increment to fields
subroutine aq_fields_add_incr_c(c_key_self,c_key_rhs) bind(c,name='aq_fields_add_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(aq_fields),pointer :: self
type(aq_fields),pointer :: rhs

! Interface
call aq_fields_registry%get(c_key_self,self)
call aq_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%add_incr(rhs)

end subroutine aq_fields_add_incr_c
! ------------------------------------------------------------------------------
!> Compute increment from the difference of two fields
subroutine aq_fields_diff_incr_c(c_key_lhs,c_key_fld1,c_key_fld2) bind(c,name='aq_fields_diff_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_lhs  !< Left-hand side
integer(c_int),intent(in) :: c_key_fld1 !< First fields
integer(c_int),intent(in) :: c_key_fld2 !< Second fields

! Local variables
type(aq_fields),pointer :: lhs
type(aq_fields),pointer :: fld1
type(aq_fields),pointer :: fld2

! Interface
call aq_fields_registry%get(c_key_lhs,lhs)
call aq_fields_registry%get(c_key_fld1,fld1)
call aq_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call lhs%diff_incr(fld1,fld2)

end subroutine aq_fields_diff_incr_c
! ------------------------------------------------------------------------------
! !> Change fields resolution
! subroutine aq_fields_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='aq_fields_change_resol_f90')

! implicit none

! ! Passed variables
! integer(c_int),intent(in) :: c_key_fld !< Fields
! integer(c_int),intent(in) :: c_key_rhs !< Right-hand side

! ! Local variables
! type(aq_fields),pointer :: fld,rhs

! ! Interface
! call aq_fields_registry%get(c_key_fld,fld)
! call aq_fields_registry%get(c_key_rhs,rhs)

! ! Call Fortran
! call aq_fields_change_resol(fld,rhs)

! end subroutine aq_fields_change_resol_c
! ------------------------------------------------------------------------------
!> Get fields info
subroutine aq_fields_info_c(c_key_fld,c_conf) bind(c,name='aq_fields_info_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_fields),pointer :: fld

! Interface
f_conf = fckit_configuration(c_conf)
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%info(f_conf)

end subroutine aq_fields_info_c
! ------------------------------------------------------------------------------
!> Read fields from file
subroutine aq_fields_read_file_c(c_key_fld,c_conf,c_dt) bind(c,name='aq_fields_read_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration
type(c_ptr),value,intent(in) :: c_dt    !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(aq_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call aq_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call fld%read(f_conf,fdate)

end subroutine aq_fields_read_file_c
! ------------------------------------------------------------------------------
!> Write fields to file
subroutine aq_fields_write_file_c(c_key_fld,c_conf,c_dt) bind(c,name='aq_fields_write_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
type(c_ptr),value,intent(in) :: c_conf !< Configuration
type(c_ptr),value,intent(in) :: c_dt   !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(aq_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call aq_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call fld%write(f_conf,fdate)

end subroutine aq_fields_write_file_c
! ------------------------------------------------------------------------------
!> Analytic initialization of fields
subroutine aq_fields_analytic_init_c(c_key_fld,c_conf) bind(c,name='aq_fields_analytic_init_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_fields),pointer :: fld

! Interface
f_conf = fckit_configuration(c_conf)
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
 call fld%analytic_IC(f_conf)

end subroutine aq_fields_analytic_init_c
! ------------------------------------------------------------------------------
!> Fields statistics
subroutine aq_fields_gpnorm_c(c_key_fld,vpresent,vmin,vmax,vrms) bind(c,name='aq_fields_gpnorm_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld      !< Fields
integer(c_int),intent(inout) :: vpresent(6) !< Variables presence flag
real(c_double),intent(inout) :: vmin(6)     !< Variables minimum
real(c_double),intent(inout) :: vmax(6)     !< Variables maximum
real(c_double),intent(inout) :: vrms(6)     !< Variables RMS

! Local variables
type(aq_fields),pointer :: fld

! Interface
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
! call fld%stats(vpresent,vmin,vmax,vrms)

end subroutine aq_fields_gpnorm_c
! ------------------------------------------------------------------------------
!> Fields RMS
subroutine aq_fields_rms_c(c_key_fld,prms) bind(c,name='aq_fields_rms_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
real(c_double),intent(inout) :: prms

! Local variables
type(aq_fields),pointer :: fld

! Interface
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%norm(prms)

end subroutine aq_fields_rms_c
! ------------------------------------------------------------------------------
!> Create ATLAS fields
subroutine aq_fields_set_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='aq_fields_set_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(aq_fields),pointer :: fld
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

! Interface
call aq_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call fld%set_atlas(vars,afieldset)

end subroutine aq_fields_set_atlas_c
! ------------------------------------------------------------------------------
!> Convert fields to ATLAS
subroutine aq_fields_to_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='aq_fields_to_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(aq_fields),pointer :: fld
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

! Interface
call aq_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call fld%to_atlas(vars,afieldset)

end subroutine aq_fields_to_atlas_c
! ------------------------------------------------------------------------------
!AQ !> Get fields from ATLAS
!AQ subroutine aq_fields_from_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='aq_fields_from_atlas_f90')
!AQ
!AQ implicit none
!AQ
!AQ ! Passed variables
!AQ integer(c_int),intent(in) :: c_key_fld           !< Fields
!AQ type(c_ptr),value,intent(in) :: c_vars           !< List of variables
!AQ type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer
!AQ
!AQ ! Local variables
!AQ type(aq_fields),pointer :: fld
!AQ type(oops_variables) :: vars
!AQ type(atlas_fieldset) :: afieldset
!AQ
!AQ ! Interface
!AQ call aq_fields_registry%get(c_key_fld,fld)
!AQ vars = oops_variables(c_vars)
!AQ afieldset = atlas_fieldset(c_afieldset)
!AQ
!AQ ! Call Fortran
!AQ ! call aq_fields_from_atlas(fld,vars,afieldset)
!AQ
!AQ end subroutine aq_fields_from_atlas_c
! ------------------------------------------------------------------------------
!> Get points from fields
subroutine aq_fields_getpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='aq_fields_getpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(inout) :: c_vals(c_nval) !< Values

! Local variables
type(aq_fields),pointer :: fld
type(aq_geom_iter),pointer :: iter

! Interface
call aq_fields_registry%get(c_key_fld,fld)
call aq_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
! call aq_fields_getpoint(fld,iter,c_nval,c_vals)

end subroutine aq_fields_getpoint_c
! ------------------------------------------------------------------------------
!> Set points for the fields
subroutine aq_fields_setpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='aq_fields_setpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(in) :: c_vals(c_nval)    !< Values

! Local variables
type(aq_fields),pointer :: fld
type(aq_geom_iter),pointer :: iter

! Interface
call aq_fields_registry%get(c_key_fld,fld)
call aq_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
! call aq_fields_setpoint(fld,iter,c_nval,c_vals)

end subroutine aq_fields_setpoint_c
! ------------------------------------------------------------------------------
!> Serial Size
subroutine aq_fields_serialsize_c(c_key_fld,c_vsize) bind(c,name='aq_fields_serialsize_f90')

implicit none

! Passed variables
integer(c_int), intent(in)  :: c_key_fld !< Fields
integer(c_int), intent(out) :: c_vsize   !< Size

! Local variables
type(aq_fields),pointer :: fld

! Interface
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
c_vsize = fld%serialsize()

end subroutine aq_fields_serialsize_c
! ------------------------------------------------------------------------------
!> Serialize fields
subroutine aq_fields_serialize_c(c_key_fld,c_vsize,c_vect_fld) bind(c,name='aq_fields_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld            !< Fields
integer(c_int),intent(in) :: c_vsize              !< Size
real(c_double),intent(out) :: c_vect_fld(c_vsize) !< Vector

! Local variables
type(aq_fields),pointer :: fld

! Interface
call aq_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%serialize(c_vect_fld)

end subroutine aq_fields_serialize_c
! ------------------------------------------------------------------------------
!> Deserialize fields
subroutine aq_fields_deserialize_c(c_key_self,c_vsize,c_vect_fld,c_index) bind(c,name='aq_fields_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Fields
integer(c_int),intent(in) :: c_vsize             !< Size
real(c_double),intent(in) :: c_vect_fld(c_vsize) !< Vector
integer(c_int), intent(inout):: c_index          !< Index

! Local variables
type(aq_fields),pointer :: self

! Interface
call aq_fields_registry%get(c_key_self,self)

! Call Fortran
call self%deserialize(c_vect_fld,c_index)

end subroutine aq_fields_deserialize_c
! ------------------------------------------------------------------------------
end module aq_fields_interface
