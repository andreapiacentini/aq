! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module aq_getvalues_interface

use datetime_mod
use iso_c_binding
use kinds
use aq_fields_mod
use aq_fields_interface
use aq_getvalues_mod
use aq_geom_mod
use aq_gom_mod
use aq_locs_mod
use interp_matrix_structure_mod, only : csr_format
use matrix_manipulations, only: deallocate_operator

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine aq_getvalues_interp_c(c_locs,c_key_fld,c_t1,c_t2,c_key_gom) bind(c,name='aq_getvalues_interp_f90')

implicit none

! Passed variables
type(c_ptr),value,intent(in) :: c_locs          !< locations
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(aq_locs) :: locs
type(aq_fields),pointer :: fld
type(aq_gom), pointer :: gom
type(datetime) :: t1, t2

! Interface
locs = aq_locs(c_locs)
call aq_fields_registry%get(c_key_fld,fld)
call aq_gom_registry%get(c_key_gom,gom)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
! Call Fortran
call aq_getvalues_interp(locs,fld,t1,t2,gom)

end subroutine aq_getvalues_interp_c
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine aq_getvalues_build_c(c_locs,c_key_fld,c_t1,c_t2,c_key_hmat) bind(c,name='aq_getvalues_build_f90')

implicit none

! Passed variables
type(c_ptr),value,intent(in) :: c_locs          !< locations
integer(c_int),intent(in) :: c_key_fld          !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2      !< times
integer(c_int),intent(in) :: c_key_hmat         !< Interpolation matrix

! Local variables
type(aq_locs) :: locs
type(aq_fields),pointer :: fld
type(csr_format), pointer :: hmat
type(datetime) :: t1, t2

! Interface
locs = aq_locs(c_locs)
call aq_fields_registry%get(c_key_fld,fld)
call aq_hmat_registry%get(c_key_hmat,hmat)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

! Call Fortran
call aq_getvalues_build(locs,fld,t1,t2,hmat)

end subroutine aq_getvalues_build_c
! ------------------------------------------------------------------------------
!> Interpolation from fields, tangent linear
subroutine aq_getvalues_interp_tl_c(c_locs, c_key_fld,c_t1,c_t2,c_key_hmat,c_key_gom) bind(c,name='aq_getvalues_interp_tl_f90')

implicit none

! Passed variables
type(c_ptr),value,intent(in) :: c_locs           !< Locations
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_hmat          !< Interpolation matrix
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(aq_locs) :: locs
type(aq_fields),pointer :: fld
type(datetime) :: t1, t2
type(csr_format), pointer :: hmat
type(aq_gom), pointer :: gom

! Interface
locs = aq_locs(c_locs)
call aq_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
call aq_hmat_registry%get(c_key_hmat,hmat)
call aq_gom_registry%get(c_key_gom,gom)

! Call Fortran
call aq_getvalues_interp_tl(locs,fld,t1,t2,hmat,gom)

end subroutine aq_getvalues_interp_tl_c
! ------------------------------------------------------------------------------
!> Interpolation from fields, adjoint
subroutine aq_getvalues_interp_ad_c(c_locs,c_key_fld,c_t1,c_t2,c_key_hmat,c_key_gom) bind(c,name='aq_getvalues_interp_ad_f90')

implicit none

! Passed variables
type(c_ptr),value,intent(in) :: c_locs          !< locations
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_hmat          !< Interpolation matrix
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(aq_locs) :: locs
type(aq_fields),pointer :: fld
type(datetime) :: t1, t2
type(csr_format), pointer :: hmat
type(aq_gom), pointer :: gom

! Interface
locs = aq_locs(c_locs)
call aq_fields_registry%get(c_key_fld,fld)
call aq_hmat_registry%get(c_key_hmat,hmat)
call aq_gom_registry%get(c_key_gom,gom)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

! Call Fortran
call aq_getvalues_interp_ad(locs,fld,t1,t2,hmat,gom)

end subroutine aq_getvalues_interp_ad_c
! ------------------------------------------------------------------------------
subroutine aq_getvalues_setup_c(c_key_hmat) bind(c,name='aq_getvalues_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_hmat          !< Interpolation matrix

! Local variables
type(csr_format), pointer :: hmat

! Interface
call aq_hmat_registry%init()
call aq_hmat_registry%add(c_key_hmat)
call aq_hmat_registry%get(c_key_hmat,hmat)

end subroutine aq_getvalues_setup_c
! ------------------------------------------------------------------------------
subroutine aq_getvalues_delete_c(c_key_hmat) bind(c,name='aq_getvalues_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_hmat          !< Interpolation matrix

! Local variables
type(csr_format), pointer :: hmat

! Interface
call aq_hmat_registry%get(c_key_hmat,hmat)

! Call Fortran
call deallocate_operator(hmat)

! Clear interface
call aq_hmat_registry%remove(c_key_hmat)

end subroutine aq_getvalues_delete_c
! ------------------------------------------------------------------------------
end module aq_getvalues_interface
