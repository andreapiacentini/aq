! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2021 UCAR.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_obsvec_interface

use iso_c_binding
use aq_obsvec_mod
use fckit_configuration_module, only: fckit_configuration

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup observation vector
subroutine aq_obsvec_setup_c(c_key_self,nlev,nobs) bind(c,name='aq_obsvec_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation vector
integer(c_int),intent(in) :: nlev          !< Number of levels
integer(c_int),intent(in) :: nobs          !< Number of observations

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%init()
call aq_obsvec_registry%add(c_key_self)
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_setup(self,nlev,nobs)

end subroutine aq_obsvec_setup_c
! ------------------------------------------------------------------------------
!> Clone observation vector
subroutine aq_obsvec_clone_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_clone_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation vector
integer(c_int),intent(in) :: c_key_other   !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_other,other)
call aq_obsvec_registry%init()
call aq_obsvec_registry%add(c_key_self)
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_clone(self,other)

end subroutine aq_obsvec_clone_c
! ------------------------------------------------------------------------------
!> Delete observation vector
subroutine aq_obsvec_delete_c(c_key_self) bind(c,name='aq_obsvec_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation vector

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_delete(self)

! Clear interface
call aq_obsvec_registry%remove(c_key_self)

end subroutine aq_obsvec_delete_c
! ------------------------------------------------------------------------------
!> Copy observation vector
subroutine aq_obsvec_copy_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_copy(self,other)

end subroutine aq_obsvec_copy_c
! ------------------------------------------------------------------------------
!> Set observation vector to zero
subroutine aq_obsvec_zero_c(c_key_self) bind(c,name='aq_obsvec_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_zero(self)

end subroutine aq_obsvec_zero_c
! ------------------------------------------------------------------------------
!> Set i-th value of the observation vector to missing value
subroutine aq_obsvec_settomissing_ith_c(c_key_self, i) bind(c,name='aq_obsvec_settomissing_ith_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
integer(c_int),intent(in) :: i          !< index of value to be set to missing value
! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
! increase index by 1 (C indices start with 0; Fortran indices start with 1)
call aq_obsvec_settomissing_ith(self, i+1)

end subroutine aq_obsvec_settomissing_ith_c
! ------------------------------------------------------------------------------
!> Set observation vector to ones
subroutine aq_obsvec_ones_c(c_key_self) bind(c,name='aq_obsvec_ones_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_ones(self)

end subroutine aq_obsvec_ones_c
! ------------------------------------------------------------------------------
!> Mask self observation vector (set values to missing where mask is set)
subroutine aq_obsvec_mask_c(c_key_self,c_key_mask) bind(c,name='aq_obsvec_mask_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_mask  !< Mask

! Local variables
type(aq_obsvec),pointer :: self,mask

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_mask,mask)

! Call Fortran
call aq_obsvec_mask(self,mask)

end subroutine aq_obsvec_mask_c
! ------------------------------------------------------------------------------
!> Mask self observation vector (set values to missing where mask is a missing value)
subroutine aq_obsvec_mask_with_missing_c(c_key_self,c_key_mask) bind(c,name='aq_obsvec_mask_with_missing_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_mask  !< Mask

! Local variables
type(aq_obsvec),pointer :: self,mask

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_mask,mask)

! Call Fortran
call aq_obsvec_mask_with_missing(self,mask)

end subroutine aq_obsvec_mask_with_missing_c
! ------------------------------------------------------------------------------
!> Mask self observation vector (set values to missing where mask is set)
subroutine aq_obsvec_threshold_check_c(c_key_self,c_key_other,c_key_mask,c_key_config) bind(c,name='aq_obsvec_threshold_check_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other  !< Observation vector
integer(c_int),intent(in) :: c_key_mask  !< Mask
type(c_ptr),value,intent(in) :: c_key_config !< Filter configuration 

! Local variables
type(aq_obsvec),pointer :: self,other,mask
type(fckit_configuration) :: f_conf

! Interface
f_conf = fckit_configuration(c_key_config)
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)
call aq_obsvec_registry%get(c_key_mask,mask)

! Call Fortran
call aq_obsvec_threshold_check(self,other,mask,f_conf)

end subroutine aq_obsvec_threshold_check_c
! ------------------------------------------------------------------------------
!> Multiply observation vector with a scalar
subroutine aq_obsvec_mul_scal_c(c_key_self,zz) bind(c,name='aq_obsvec_mul_scal_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
real(c_double),intent(in) :: zz         !< Multiplier

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_mul_scal(self,zz)

end subroutine aq_obsvec_mul_scal_c
! ------------------------------------------------------------------------------
!> Add observation vector
subroutine aq_obsvec_add_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_add(self,other)

end subroutine aq_obsvec_add_c
! ------------------------------------------------------------------------------
!> Subtract observation vector
subroutine aq_obsvec_sub_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_sub_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_sub(self,other)

end subroutine aq_obsvec_sub_c
! ------------------------------------------------------------------------------
!> Multiply observation vector
subroutine aq_obsvec_mul_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_mul_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_mul(self,other)

end subroutine aq_obsvec_mul_c
! ------------------------------------------------------------------------------
!> Divide observation vector
subroutine aq_obsvec_div_c(c_key_self,c_key_other) bind(c,name='aq_obsvec_div_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_div(self,other)

end subroutine aq_obsvec_div_c
! ------------------------------------------------------------------------------
!> Apply axpy on observation vector
subroutine aq_obsvec_axpy_c(c_key_self,zz,c_key_other) bind(c,name='aq_obsvec_axpy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Observation vector
real(c_double),intent(in) :: zz          !< Multiplier
integer(c_int),intent(in) :: c_key_other !< Other observation vector

! Local variables
type(aq_obsvec),pointer :: self,other

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_other,other)

! Call Fortran
call aq_obsvec_axpy(self,zz,other)

end subroutine aq_obsvec_axpy_c
! ------------------------------------------------------------------------------
!> Invert observation vector
subroutine aq_obsvec_invert_c(c_key_self) bind(c,name='aq_obsvec_invert_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_invert(self)

end subroutine aq_obsvec_invert_c
! ------------------------------------------------------------------------------
!> Generate random observation vector
subroutine aq_obsvec_random_c(c_odb,c_self) bind(c,name='aq_obsvec_random_f90')

implicit none

! Passed variables
type(c_ptr),intent(in) :: c_odb     !< Observation data base
integer(c_int),intent(in) :: c_self !< Observation vector

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_self,self)

! Call Fortran
call aq_obsvec_random(c_odb,self)

end subroutine aq_obsvec_random_c
! ------------------------------------------------------------------------------
!> Compute dot product between observation vectors
subroutine aq_obsvec_dotprod_c(c_key_obsvec1,c_key_obsvec2,zz) bind(c,name='aq_obsvec_dotprod_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_obsvec1 !< Observation vector 1
integer(c_int),intent(in) :: c_key_obsvec2 !< Observation vector 2
real(c_double),intent(inout) :: zz         !< Dot product

! Local variables
type(aq_obsvec),pointer :: obsvec1,obsvec2

! Interface
call aq_obsvec_registry%get(c_key_obsvec1,obsvec1)
call aq_obsvec_registry%get(c_key_obsvec2,obsvec2)

! Call Fortran
call aq_obsvec_dotprod(obsvec1,obsvec2,zz)

end subroutine aq_obsvec_dotprod_c
! ------------------------------------------------------------------------------
!> Compute observation vector statistics
subroutine aq_obsvec_stats_c(c_key_self,zmin,zmax,zavg) bind(c,name='aq_obsvec_stats_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
real(c_double),intent(inout) :: zmin    !< Minimum
real(c_double),intent(inout) :: zmax    !< Maximum
real(c_double),intent(inout) :: zavg    !< Average

! Local variables
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_stats(self,zmin,zmax,zavg)

end subroutine aq_obsvec_stats_c
! ------------------------------------------------------------------------------
!> Get number of observations (not missing)
subroutine aq_obsvec_nobs_c(c_key_self,kobs) bind(c,name='aq_obsvec_nobs_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
integer(c_int),intent(inout) :: kobs    !< Observation vector size

! Local vector
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_nobs(self,kobs)

end subroutine aq_obsvec_nobs_c
! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine aq_obsvec_size_c(c_key_self,kobs) bind(c,name='aq_obsvec_size_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
integer(c_int),intent(inout) :: kobs    !< Observation vector size

! Local vector
type(aq_obsvec),pointer :: self

! Interface
call aq_obsvec_registry%get(c_key_self,self)

! Call Fortran
call aq_obsvec_size(self,kobs)

end subroutine aq_obsvec_size_c
! ------------------------------------------------------------------------------
!> Get observation vector size (only non-masked observations)
subroutine aq_obsvec_nobs_withmask_c(c_key_self,c_key_mask,kobs) bind(c,name='aq_obsvec_nobs_withmask_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
integer(c_int),intent(in) :: c_key_mask !< Mask
integer(c_int),intent(inout) :: kobs    !< Observation vector size

! Local vector
type(aq_obsvec),pointer :: self, mask

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_mask,mask)

! Call Fortran
call aq_obsvec_nobs_withmask(self,mask,kobs)

end subroutine aq_obsvec_nobs_withmask_c

! ------------------------------------------------------------------------------
!> Get all non-masked out observation values
subroutine aq_obsvec_get_withmask_c(c_key_self,c_key_mask,vals,nvals) bind(c,name='aq_obsvec_get_withmask_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Observation vector
integer(c_int),intent(in) :: c_key_mask !< Mask
integer(c_int),intent(in) :: nvals      !< number of obs
real(c_double),intent(out),dimension(nvals) :: vals  !< ob. values

! Local vector
type(aq_obsvec),pointer :: self, mask

! Interface
call aq_obsvec_registry%get(c_key_self,self)
call aq_obsvec_registry%get(c_key_mask,mask)

! Call Fortran
call aq_obsvec_get_withmask(self,mask,vals,nvals)

end subroutine aq_obsvec_get_withmask_c

! ------------------------------------------------------------------------------
end module aq_obsvec_interface
