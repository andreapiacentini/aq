! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_stream_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use aq_gom_mod
use aq_obsvec_mod
use aq_stream_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction
subroutine aq_stream_equiv_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='aq_stream_equiv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(in) :: c_bias     !< Bias

! Local variables
type(aq_gom),pointer  :: gom
type(aq_obsvec),pointer :: hofx

! Interface
call aq_gom_registry%get(c_key_gom,gom) 
call aq_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call aq_stream_equiv(gom,hofx,c_bias)

end subroutine aq_stream_equiv_c
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction - tangent linear
subroutine aq_stream_equiv_tl_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='aq_stream_equiv_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(in) :: c_bias     !< Bias

! Local variables
type(aq_gom),pointer  :: gom
type(aq_obsvec),pointer :: hofx

! Interface
call aq_gom_registry%get(c_key_gom,gom) 
call aq_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call aq_stream_equiv(gom,hofx,c_bias)

end subroutine aq_stream_equiv_tl_c
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction - adjoint
subroutine aq_stream_equiv_ad_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='aq_stream_equiv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(inout) :: c_bias  !< Bias

! Local variables
type(aq_gom),pointer  :: gom
type(aq_obsvec),pointer :: hofx

! Interface
call aq_gom_registry%get(c_key_gom,gom)
call aq_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call aq_stream_equiv_ad(gom,hofx,c_bias)

end subroutine aq_stream_equiv_ad_c
! ------------------------------------------------------------------------------
end module aq_stream_interface
