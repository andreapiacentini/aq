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
  if (insitu%values(1,iobs) /= insitu%missing ) then
    insitu%values(1,iobs) = geovals%x(iobs)+bias
  endif
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
  if (insitu%values(1,iobs) /= insitu%missing ) then
    bias = bias+insitu%values(1,iobs)
  endif
enddo

end subroutine aq_insitu_equiv_ad
! ------------------------------------------------------------------------------
end module aq_insitu_mod
