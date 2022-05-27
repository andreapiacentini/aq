! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_geovals_interface

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,only: fckit_mpi_comm
use iso_c_binding
use aq_geom_mod
use aq_geovals_mod
use aq_locs_mod
use oops_variables_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup GeoVals
subroutine aq_geovals_setup_c(c_key_self,c_locs,c_vars,c_comm) bind(c,name='aq_geovals_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GeoVals
type(c_ptr),value,intent(in) :: c_locs     !< Locations
type(c_ptr),value,intent(in) :: c_vars     !< Variables
type(c_ptr),intent(in),value :: c_comm

! Local variables
type(aq_geovals),pointer :: self
type(aq_locs) :: locs

! Interface
call aq_geovals_registry%init()
call aq_geovals_registry%add(c_key_self)
call aq_geovals_registry%get(c_key_self,self)
locs = aq_locs(c_locs)
self%vars = oops_variables(c_vars)
self%fmpi = fckit_mpi_comm(c_comm)

! Call Fortran
call aq_geovals_setup(self,locs%nlocs())

end subroutine aq_geovals_setup_c
! ------------------------------------------------------------------------------
!> Create GeoVals and do nothing
subroutine aq_geovals_create_c(c_key_self,c_vars,c_comm) bind(c,name='aq_geovals_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GeoVals
type(c_ptr),value,intent(in) :: c_vars     !< Variables
type(c_ptr),intent(in),value :: c_comm

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%init()
call aq_geovals_registry%add(c_key_self)
call aq_geovals_registry%get(c_key_self,self)
self%vars = oops_variables(c_vars)
self%fmpi = fckit_mpi_comm(c_comm)

end subroutine aq_geovals_create_c
! ------------------------------------------------------------------------------
!> Delete GeoVals
subroutine aq_geovals_delete_c(c_key_self) bind(c,name='aq_geovals_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GeoVals

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_delete(self)

! Clear interface
call aq_geovals_registry%remove(c_key_self)

end subroutine aq_geovals_delete_c
! ------------------------------------------------------------------------------
!> Copy GeoVals
subroutine aq_geovals_copy_c(c_key_self,c_key_other) bind(c,name='aq_geovals_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GeoVals
integer(c_int),intent(in) :: c_key_other   !< Other GeoVals

! Local variables
type(aq_geovals),pointer :: self
type(aq_geovals),pointer :: other

! Interface
call aq_geovals_registry%get(c_key_self,self)
call aq_geovals_registry%get(c_key_other,other)

! Call Fortran
call aq_geovals_copy(self,other)

end subroutine aq_geovals_copy_c
! ------------------------------------------------------------------------------
subroutine aq_geovals_fill_c(c_key, c_nloc, c_indx, c_nval, c_vals) bind(c, name="aq_geovals_fill_f90")
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nval
real(c_double), intent(in) :: c_vals(c_nval)

type(aq_geovals), pointer :: self

call aq_geovals_registry%get(c_key,self)

call aq_geovals_fill(self, c_nloc, c_indx, c_nval, c_vals)

end subroutine aq_geovals_fill_c
! ------------------------------------------------------------------------------
subroutine aq_geovals_fillad_c(c_key, c_nloc, c_indx, c_nval, c_vals) bind(c, name="aq_geovals_fillad_f90")
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nval
real(c_double), intent(inout) :: c_vals(c_nval)

type(aq_geovals),pointer :: self

call aq_geovals_registry%get(c_key, self)

call aq_geovals_fillad(self, c_nloc, c_indx, c_nval, c_vals)

end subroutine aq_geovals_fillad_c
! ------------------------------------------------------------------------------
!> Set GeoVals to zero
subroutine aq_geovals_zero_c(c_key_self) bind(c,name='aq_geovals_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_zero(self)

end subroutine aq_geovals_zero_c
! ------------------------------------------------------------------------------
!> Get GeoVals absolute value
subroutine aq_geovals_abs_c(c_key_self) bind(c,name='aq_geovals_abs_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_abs(self)

end subroutine aq_geovals_abs_c
! ------------------------------------------------------------------------------
!> Generate random GeoVals
subroutine aq_geovals_random_c(c_key_self) bind(c,name='aq_geovals_random_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_random(self)

end subroutine aq_geovals_random_c
! ------------------------------------------------------------------------------
!> Multiply GeoVals with a scalar
subroutine aq_geovals_mult_c(c_key_self,c_zz) bind(c,name='aq_geovals_mult_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
real(c_double),intent(in) :: c_zz       !< Multiplier

! Local variables
type(aq_geovals),pointer :: self
integer :: jo,jv

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_mult(self,c_zz)

end subroutine aq_geovals_mult_c
! ------------------------------------------------------------------------------
!> Add GeoVals
subroutine aq_geovals_add_c(c_key_self,c_key_other) bind(c,name='aq_geovals_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GeoVals
integer(c_int),intent(in) :: c_key_other !< Other GeoVals

! Local variables
type(aq_geovals),pointer :: self
type(aq_geovals),pointer :: other

! Interface
call aq_geovals_registry%get(c_key_self,self)
call aq_geovals_registry%get(c_key_other,other)

! Call Fortran
call aq_geovals_add(self,other)

end subroutine aq_geovals_add_c
! ------------------------------------------------------------------------------
!> Subtract GeoVals
subroutine aq_geovals_diff_c(c_key_self,c_key_other) bind(c,name='aq_geovals_diff_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GeoVals
integer(c_int),intent(in) :: c_key_other !< Other GeoVals

! Local variables
type(aq_geovals),pointer :: self
type(aq_geovals),pointer :: other

! Interface
call aq_geovals_registry%get(c_key_self,self)
call aq_geovals_registry%get(c_key_other,other)

! Call Fortran
call aq_geovals_diff(self,other)

end subroutine aq_geovals_diff_c
! ------------------------------------------------------------------------------
!> Schur product for GeoVals
subroutine aq_geovals_schurmult_c(c_key_self,c_key_other) bind(c,name='aq_geovals_schurmult_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GeoVals
integer(c_int),intent(in) :: c_key_other !< Other GeoVals

! Local variables
type(aq_geovals),pointer :: self
type(aq_geovals),pointer :: other

! Interface
call aq_geovals_registry%get(c_key_self,self)
call aq_geovals_registry%get(c_key_other,other)

! Call Fortran
call aq_geovals_schurmult(self,other)

end subroutine aq_geovals_schurmult_c
! ------------------------------------------------------------------------------
!> Schur division for GeoVals
subroutine aq_geovals_divide_c(c_key_self,c_key_other) bind(c,name='aq_geovals_divide_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GeoVals
integer(c_int),intent(in) :: c_key_other !< Other GeoVals

! Local variables
type(aq_geovals),pointer :: self
type(aq_geovals),pointer :: other

! Interface
call aq_geovals_registry%get(c_key_self,self)
call aq_geovals_registry%get(c_key_other,other)

! Call Fortran
call aq_geovals_divide(self,other)

end subroutine aq_geovals_divide_c
! ------------------------------------------------------------------------------
!> Compute GeoVals RMS
subroutine aq_geovals_rms_c(c_key_self,c_rms) bind(c,name='aq_geovals_rms_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
real(c_double),intent(inout) :: c_rms   !< RMS

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_rms(self,c_rms)

end subroutine aq_geovals_rms_c
! ------------------------------------------------------------------------------
!> GeoVals dot product
subroutine aq_geovals_dotprod_c(c_key_geovals1,c_key_geovals2,c_prod) bind(c,name='aq_geovals_dotprod_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_geovals1 !< GeoVals 1
integer(c_int),intent(in) :: c_key_geovals2 !< GeoVals 2
real(c_double),intent(inout) :: c_prod  !< Dot product

! Local variables
type(aq_geovals),pointer :: geovals1,geovals2

! Interface
call aq_geovals_registry%get(c_key_geovals1,geovals1)
call aq_geovals_registry%get(c_key_geovals2,geovals2)

! Call Fortran
call aq_geovals_dotprod(geovals1,geovals2,c_prod)

end subroutine aq_geovals_dotprod_c
! ------------------------------------------------------------------------------
!> Compute GeoVals statistics
subroutine aq_geovals_stats_c(c_key_self,c_kobs,c_pmin,c_pmax,c_pave,c_pstd) bind(c,name='aq_geovals_stats_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
integer(c_int),intent(inout) :: c_kobs  !< Number of observations
real(c_double),intent(inout) :: c_pmin  !< Minimum value
real(c_double),intent(inout) :: c_pmax  !< Maximum value
real(c_double),intent(inout) :: c_pave  !< Mean
real(c_double),intent(inout) :: c_pstd  !< StdDev

! Local variables
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_stats(self,c_kobs,c_pmin,c_pmax,c_pave,c_pstd)

end subroutine aq_geovals_stats_c
! ------------------------------------------------------------------------------
!> Find and locate GeoVals max. value
subroutine aq_geovals_maxloc_c(c_key_self,c_mxval,c_mxloc,c_mxvar) bind(c,name='aq_geovals_maxloc_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
real(c_double),intent(inout) :: c_mxval !< Maximum value
integer(c_int),intent(inout) :: c_mxloc !< Location of maximum value
type(c_ptr),value,intent(in) :: c_mxvar !< Variable of maximum value

! Local variables
type(aq_geovals),pointer :: self
type(oops_variables) :: mxvar

! Interface
call aq_geovals_registry%get(c_key_self,self)
mxvar = oops_variables(c_mxvar)

! Call Fortran
call aq_geovals_maxloc(self,c_mxval,c_mxloc,mxvar)

end subroutine aq_geovals_maxloc_c
! ------------------------------------------------------------------------------
!> Read GeoVals from file
subroutine aq_geovals_read_file_c(c_key_self,c_conf) bind(c,name='aq_geovals_read_file_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GeoVals
type(c_ptr),value,intent(in) :: c_conf     !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_geovals),pointer :: self

! Interface
call aq_geovals_registry%get(c_key_self,self)
f_conf = fckit_configuration(c_conf)

! Call Fortran
call aq_geovals_read_file(self,f_conf)

end subroutine aq_geovals_read_file_c
! ------------------------------------------------------------------------------
!> Write GeoVals to file
subroutine aq_geovals_write_file_c(c_key_self,c_conf) bind(c,name='aq_geovals_write_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_geovals),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call aq_geovals_registry%get(c_key_self,self)

! Call Fortran
call aq_geovals_write_file(self,f_conf)

end subroutine aq_geovals_write_file_c
! ------------------------------------------------------------------------------
!> GeoVals analytic initialization
subroutine aq_geovals_analytic_init_c(c_key_self,c_locs,c_conf) bind(c,name='aq_geovals_analytic_init_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GeoVals
type(c_ptr),value,intent(in) :: c_locs  !< Locations
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(aq_geovals),pointer :: self
type(aq_locs) :: locs

! Interface
f_conf = fckit_configuration(c_conf)
call aq_geovals_registry%get(c_key_self,self)
locs = aq_locs(c_locs)

! Call Fortran
call aq_geovals_analytic_init(self,locs,f_conf)

end subroutine aq_geovals_analytic_init_c
! ------------------------------------------------------------------------------
end module aq_geovals_interface
