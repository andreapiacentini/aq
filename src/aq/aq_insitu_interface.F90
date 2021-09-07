! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_insitu_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use aq_geovals_mod
use aq_obsvec_mod
use aq_insitu_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Get equivalent
subroutine aq_insitu_equiv_c(c_key_geovals,c_key_insitu,c_bias) bind(c,name='aq_insitu_equiv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_geovals !< GeoVals
integer(c_int),intent(in) :: c_key_insitu    !< Observation vector
real(c_double),intent(in) :: c_bias        !< Bias

! Local variables
type(aq_geovals),pointer :: geovals
type(aq_obsvec),pointer :: insitu

! Interface
call aq_geovals_registry%get(c_key_geovals,geovals) 
call aq_obsvec_registry%get(c_key_insitu,insitu)

! Call Fortran
call aq_insitu_equiv(geovals,insitu,c_bias)

end subroutine aq_insitu_equiv_c
! ------------------------------------------------------------------------------
!> Get equivalent - tangent linear
subroutine aq_insitu_equiv_tl_c(c_key_geovals,c_key_insitu,c_bias) bind(c,name='aq_insitu_equiv_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_geovals  !< GeoVals
integer(c_int),intent(in) :: c_key_insitu !< Observation vector
real(c_double),intent(in) :: c_bias     !< Bias

! Local variables
type(aq_geovals),pointer  :: geovals
type(aq_obsvec),pointer :: insitu

! Interface
call aq_geovals_registry%get(c_key_geovals,geovals) 
call aq_obsvec_registry%get(c_key_insitu,insitu)

! Call Fortran
call aq_insitu_equiv(geovals,insitu,c_bias)

end subroutine aq_insitu_equiv_tl_c
! ------------------------------------------------------------------------------
!> Get equivalent - adjoint
subroutine aq_insitu_equiv_ad_c(c_key_geovals,c_key_insitu,c_bias) bind(c,name='aq_insitu_equiv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_geovals  !< GeoVals
integer(c_int),intent(in) :: c_key_insitu !< Observation vector
real(c_double),intent(inout) :: c_bias  !< Bias

! Local variables
type(aq_geovals),pointer  :: geovals
type(aq_obsvec),pointer :: insitu

! Interface
call aq_geovals_registry%get(c_key_geovals,geovals)
call aq_obsvec_registry%get(c_key_insitu,insitu)

! Call Fortran
call aq_insitu_equiv_ad(geovals,insitu,c_bias)

end subroutine aq_insitu_equiv_ad_c
! ------------------------------------------------------------------------------
end module aq_insitu_interface
