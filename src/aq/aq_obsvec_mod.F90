! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_obsvec_mod

use iso_c_binding
use kinds
use random_mod
use fckit_module
use aq_constants_mod

implicit none

private
public :: aq_obsvec
public :: aq_obsvec_registry
public :: aq_obsvec_setup,aq_obsvec_clone,aq_obsvec_delete,aq_obsvec_copy,aq_obsvec_zero, &
        & aq_obsvec_settomissing_ith,aq_obsvec_ones,aq_obsvec_mask,aq_obsvec_mask_with_missing, aq_obsvec_threshold_check, &
        & aq_obsvec_mul_scal,aq_obsvec_add,aq_obsvec_sub,aq_obsvec_mul,aq_obsvec_div, &
        & aq_obsvec_axpy,aq_obsvec_invert,aq_obsvec_random,aq_obsvec_dotprod,aq_obsvec_stats, &
        & aq_obsvec_size,aq_obsvec_nobs,aq_obsvec_nobs_withmask,aq_obsvec_get_withmask
! ------------------------------------------------------------------------------
interface
  subroutine aq_obsvec_random_i(odb,nn,zz) bind(c,name='aq_obsvec_random_f')
  use iso_c_binding
  implicit none
  type(c_ptr),intent(in) :: odb
  integer(c_int),intent(in) :: nn
  real(c_double),intent(inout) :: zz
  end subroutine aq_obsvec_random_i
end interface
! ------------------------------------------------------------------------------
type aq_obsvec
  integer :: nobs = 0                        !< Number of observations
  integer :: nlev = 0                        !< Number of levels
  real(kind_real),allocatable :: values(:,:) !< Values
  real(kind_real) :: missing                 !< Missing value
end type aq_obsvec

#define LISTED_TYPE aq_obsvec

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_obsvec_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup observation vector
subroutine aq_obsvec_setup(self,nlev,nobs)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
integer,intent(in) :: nlev            !< Number of levels
integer,intent(in) :: nobs            !< Number of observations

! Set sizes
self%nlev = nlev
self%nobs = nobs

! Release memory
if (allocated(self%values)) deallocate(self%values)

! Allocation
allocate(self%values(self%nlev,self%nobs))

! Initialization
self%values = 0.0_kind_real
self%missing = aq_missing_value

end subroutine aq_obsvec_setup
! ------------------------------------------------------------------------------
!> Clone observation vector
subroutine aq_obsvec_clone(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

! Set sizes
self%nlev = other%nlev
self%nobs = other%nobs
self%missing = other%missing

! Allocation
allocate(self%values(self%nlev,self%nobs))

end subroutine aq_obsvec_clone
! ------------------------------------------------------------------------------
!> Delete observation vector
subroutine aq_obsvec_delete(self)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector

! Release memory
deallocate(self%values)

end subroutine aq_obsvec_delete
! ------------------------------------------------------------------------------
!> Copy observation vector
subroutine aq_obsvec_copy(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((other%nlev/=self%nlev).or.(other%nobs/=self%nobs)) then
  ! Release memory
  deallocate(self%values)

  ! Set sizes
  self%nlev = other%nlev
  self%nobs = other%nobs
  self%missing = other%missing

  ! Allocation
  allocate(self%values(self%nlev,self%nobs))
endif

! Copy data
self%values = other%values

end subroutine aq_obsvec_copy
! ------------------------------------------------------------------------------
!> Set observation vector to zero
subroutine aq_obsvec_zero(self)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector

! Set observation vector to zero
self%values = 0.0

end subroutine aq_obsvec_zero
! ------------------------------------------------------------------------------
!> Set i-th value of observation vector to missing value
subroutine aq_obsvec_settomissing_ith(self, i)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
integer, intent(in) :: i

! Set observation vector to zero
self%values(:,i) = self%missing

end subroutine aq_obsvec_settomissing_ith
! ------------------------------------------------------------------------------
!> Set observation vector to ones
subroutine aq_obsvec_ones(self)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector

! Set observation vector to ones
self%values = 1.0

end subroutine aq_obsvec_ones
! ------------------------------------------------------------------------------
!> Mask observation vector (set values to missing values where mask == 1)
subroutine aq_obsvec_mask(self,mask)
implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: mask    !< mask

if ((self%nobs/=mask%nobs).or.(self%nlev/=mask%nlev)) call abor1_ftn('aq_obsvec_mask: inconsistent sizes')

where(mask%values == 1) self%values = self%missing

end subroutine aq_obsvec_mask
! ------------------------------------------------------------------------------
!> Mask observation vector (set values to missing values where mask == missing value)
subroutine aq_obsvec_mask_with_missing(self,mask)
implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: mask    !< mask

if ((self%nobs/=mask%nobs).or.(self%nlev/=mask%nlev)) call abor1_ftn('aq_obsvec_mask_with_missing: inconsistent sizes')

where(mask%values == mask%missing) self%values = self%missing

end subroutine aq_obsvec_mask_with_missing
! ------------------------------------------------------------------------------
!> Mask observations and set qcflag to one where observation vectors differ by more than lambda
subroutine aq_obsvec_threshold_check(self,other,qcflag,config)
implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self          !< Simulated obs vector
type(aq_obsvec),intent(inout) :: other         !< Obs vector
type(aq_obsvec),intent(inout) :: qcflag        !< Mask
type(fckit_Configuration),intent(in) :: config !< Filter config

real(aq_real) :: threshold

call config%get_or_die("threshold",threshold)
if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_mask: inconsistent sizes')

where(abs(self%values - other%values) > threshold) 
  self%values = self%missing
  qcflag%values = 1
end where

end subroutine aq_obsvec_threshold_check
! ------------------------------------------------------------------------------
!> Multiply observation vector with a scalar
subroutine aq_obsvec_mul_scal(self,zz)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier

! Multiply observation vector with a scalar
where(self%values /= self%missing) self%values = zz*self%values

end subroutine aq_obsvec_mul_scal
! ------------------------------------------------------------------------------
!> Add observation vector
subroutine aq_obsvec_add(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_add: inconsistent sizes')

! Add observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values+other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine aq_obsvec_add
! ------------------------------------------------------------------------------
!> Subtract observation vector
subroutine aq_obsvec_sub(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_sub: inconsistent sizes')

! Subtract observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values-other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine aq_obsvec_sub
! ------------------------------------------------------------------------------
!> Multiply observation vector
subroutine aq_obsvec_mul(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_mul: inconsistent sizes')

! Multiply observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values*other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine aq_obsvec_mul
! ------------------------------------------------------------------------------
!> Divide observation vector
subroutine aq_obsvec_div(self,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_div: inconsistent sizes')

! Divide observation vector
where(self%values /= self%missing .and. other%values /= other%missing .and. abs(other%values)>1.e-30)
  self%values = self%values/other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine aq_obsvec_div
! ------------------------------------------------------------------------------
!> Apply axpy on observation vector
subroutine aq_obsvec_axpy(self,zz,other)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier
type(aq_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('aq_obsvec_axpy: inconsistent sizes')

! Apply axpy on observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values+zz*other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine aq_obsvec_axpy
! ------------------------------------------------------------------------------
!> Invert observation vector
subroutine aq_obsvec_invert(self)

implicit none

! Passed variables
type(aq_obsvec),intent(inout) :: self !< Observation vector

! Invert observation vector
where(self%values /= self%missing .and. abs(self%values)>1.e-30) self%values = 1.0/self%values

end subroutine aq_obsvec_invert
! ------------------------------------------------------------------------------
!> Generate random observation vector
subroutine aq_obsvec_random(c_odb,self)

implicit none

! Passed variables
type(c_ptr),intent(in) :: c_odb       !< Observation data base
type(aq_obsvec),intent(inout) :: self !< Observation vector

! Local variables
integer :: nval

! Compute total size
nval = self%nobs*self%nlev

! Get random values
call aq_obsvec_random_i(c_odb,nval,self%values(1,1))

end subroutine aq_obsvec_random
! ------------------------------------------------------------------------------
!> Compute dot product between observation vectors
subroutine aq_obsvec_dotprod(obsvec1,obsvec2,zz)

implicit none

! Passed variables
type(aq_obsvec),intent(in) :: obsvec1 !< Observation vector 1
type(aq_obsvec),intent(in) :: obsvec2 !< Observation vector 2
real(kind_real),intent(inout) :: zz   !< Dot product

! Local variables
integer :: jlev,jobs

! Check sizes
if ((obsvec1%nobs/=obsvec2%nobs).or.(obsvec1%nlev/=obsvec2%nlev)) call abor1_ftn('aq_obsvec_dotprod: inconsistent sizes')

! Initalization
zz = 0.0

! Loop over values
do jobs=1,obsvec1%nobs
  do jlev=1,obsvec1%nlev
    if (obsvec1%values(jlev, jobs) /= obsvec1%missing .and. &
        obsvec2%values(jlev, jobs) /= obsvec2%missing)      &
      zz = zz+obsvec1%values(jlev,jobs)*obsvec2%values(jlev,jobs)
  enddo
enddo

end subroutine aq_obsvec_dotprod
! ------------------------------------------------------------------------------
!> Compute observation vector statistics
subroutine aq_obsvec_stats(self,zmin,zmax,zavg)

implicit none

! Passed variables
type(aq_obsvec),intent(in):: self        !< Observation vector
real(kind_real),intent(inout) :: zmin    !< Minimum
real(kind_real),intent(inout) :: zmax    !< Maximum
real(kind_real),intent(inout) :: zavg    !< Average

if (self%nobs*self%nlev>0) then
  ! Compute statistics
  if (.not.allocated(self%values)) call abor1_ftn('aq_obsvec_stats: obs vector not allocated')
  zmin = minval(self%values, mask = (self%values /= self%missing))
  zmax = maxval(self%values, mask = (self%values /= self%missing))
  zavg = sum(self%values, mask = (self%values /= self%missing)) /    &
         count(mask = (self%values /= self%missing))
else
  ! Empty observation vector
  zmin = 0.0
  zmax = 0.0
  zavg = 0.0
endif

end subroutine aq_obsvec_stats

! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine aq_obsvec_size(self,kobs)
implicit none
type(aq_obsvec),intent(in) :: self !< Observation vector
integer,intent(inout) :: kobs      !< Observation vector size
kobs = size(self%values) + 2
end subroutine aq_obsvec_size

! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine aq_obsvec_nobs(self,kobs)

implicit none

! Passed variables
type(aq_obsvec),intent(in) :: self !< Observation vector
integer,intent(inout) :: kobs      !< Observation vector size

! Get observation vector size
kobs = count(mask = (self%values /= self%missing))

end subroutine aq_obsvec_nobs

! ------------------------------------------------------------------------------
!> Get observation vector size (only non-masked observations)
subroutine aq_obsvec_nobs_withmask(self,obsmask,kobs)

implicit none

! Passed variables
type(aq_obsvec),intent(in) :: self    !< Observation vector
type(aq_obsvec),intent(in) :: obsmask !< mask
integer,intent(inout) :: kobs         !< Observation vector size

! Get observation vector size
kobs = count(mask = (self%values /= self%missing) .and.     &
                    (obsmask%values /= obsmask%missing))

end subroutine aq_obsvec_nobs_withmask

! ------------------------------------------------------------------------------
!> Get non-missing values from observation vector into vals array
subroutine aq_obsvec_get_withmask(self,obsmask,vals,nvals)

implicit none

! Passed variables
type(aq_obsvec),intent(in) :: self !< Observation vector
type(aq_obsvec),intent(in) :: obsmask !< mask
integer,intent(in) :: nvals        !< Number of non-missing values
real(kind_real), dimension(nvals), intent(out) :: vals!< returned value

integer :: jobs, jlev, jval

jval = 1
! Loop over values
do jobs=1,self%nobs
  do jlev=1,self%nlev
    if ((self%values(jlev, jobs) /= self%missing) .and.           &
        (obsmask%values(jlev, jobs) /= obsmask%missing)) then
      if (jval > nvals) call abor1_ftn('aq_obsvec_get: inconsistent vector size')
      vals(jval) = self%values(jlev, jobs)
      jval = jval + 1
    endif
  enddo
enddo

end subroutine aq_obsvec_get_withmask

! ------------------------------------------------------------------------------
end module aq_obsvec_mod
