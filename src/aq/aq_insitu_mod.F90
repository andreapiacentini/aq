! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_insitu_mod

use kinds
use aq_geovals_mod
use aq_obsvec_mod

implicit none

private
public :: aq_insitu_equiv,aq_insitu_equiv_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent (TL calls this subroutine too)
subroutine aq_insitu_equiv(geovals,insitu,bias)

implicit none

! Passed variables
type(aq_geovals),intent(in) :: geovals        !< GeoVals
type(aq_obsvec),intent(inout) :: insitu !< Observation vector
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,geovals%nobs
  insitu%values(1,iobs) = geovals%x(iobs)+bias
enddo

end subroutine aq_insitu_equiv
! ------------------------------------------------------------------------------
!> Get equivalent - adjoint
subroutine aq_insitu_equiv_ad(geovals,insitu,bias)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: geovals     !< GeoVals
type(aq_obsvec),intent(in) :: insitu    !< Observation vector
real(kind_real),intent(inout) :: bias !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,geovals%nobs
  geovals%x(iobs) = insitu%values(1,iobs)
  bias = bias+insitu%values(1,iobs)
enddo

end subroutine aq_insitu_equiv_ad
! ------------------------------------------------------------------------------
end module aq_insitu_mod
