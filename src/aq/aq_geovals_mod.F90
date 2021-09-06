! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_geovals_mod

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use netcdf
use oops_variables_mod
use aq_geom_mod
use aq_locs_mod
use aq_tools_mod, only: ncerr
use random_mod

implicit none
private
public :: aq_geovals
public :: aq_geovals_registry
public :: aq_geovals_setup,aq_geovals_delete,aq_geovals_copy,aq_geovals_zero,aq_geovals_abs,aq_geovals_random,aq_geovals_mult, &
        & aq_geovals_add,aq_geovals_diff,aq_geovals_schurmult,aq_geovals_divide,aq_geovals_rms,aq_geovals_dotprod,aq_geovals_stats,aq_geovals_maxloc, &
        & aq_geovals_read_file, aq_geovals_write_file,aq_geovals_analytic_init
! ------------------------------------------------------------------------------
type :: aq_geovals
  integer :: nobs                      !< Number of observations
  real(kind_real), allocatable :: x(:) !< Chemical observations values
  logical :: lalloc = .false.          !< Allocation flag
  type(oops_variables) :: vars         !< Variables
end type aq_geovals

#define LISTED_TYPE aq_geovals

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_geovals_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup GeoVals
subroutine aq_geovals_setup(self,nobs)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self      !< GeoVals
integer, intent(in) :: nobs             !< Number of observations

! Set attributes
self%nobs = nobs

! Allocation
allocate(self%x(self%nobs))

self%lalloc = .true.

end subroutine aq_geovals_setup
! ------------------------------------------------------------------------------
!> Delete GeoVals
subroutine aq_geovals_delete(self)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals

! Release memory
if (allocated(self%x)) deallocate(self%x)

self%lalloc = .false.

end subroutine aq_geovals_delete
! ------------------------------------------------------------------------------
!> Copy GeoVals
subroutine aq_geovals_copy(self,other)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self            !< GeoVals
type(aq_geovals),intent(in) :: other              !< Other GeoVals

! Copy attributes
self%nobs = other%nobs

! Allocation
if (.not.self%lalloc) then
  allocate(self%x(self%nobs))
  self%lalloc = .true.
endif

! Copy
!AQ Temporary stub waiting for geovals files I/O
if (allocated(other%x)) then
self%x = other%x
end if
!AQ till here

end subroutine aq_geovals_copy
! ------------------------------------------------------------------------------
!> Set GeoVals to zero
subroutine aq_geovals_zero(self)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals

! Set to zero
self%x = 0.0

end subroutine aq_geovals_zero
! ------------------------------------------------------------------------------
!> Get GeoVals absolute value
subroutine aq_geovals_abs(self)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals

! Get absolute value
self%x = abs(self%x)

end subroutine aq_geovals_abs
! ------------------------------------------------------------------------------
!> Generate random GeoVals values
subroutine aq_geovals_random(self)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals

! Local variables
integer :: nv
real(kind_real),allocatable :: values(:,:)

! TODO(Benjamin): change that in a following PR
nv = 0
nv = nv+1
allocate(values(nv,self%nobs))

! Generate random GeoVals values
call normal_distribution(values,0.0_kind_real,1.0_kind_real,seed=14871,reset=.true.)

! Split random values
nv = 0
nv = nv+1
self%x = values(nv,:)

end subroutine aq_geovals_random
! ------------------------------------------------------------------------------
!> Multiply GeoVals with a scalar
subroutine aq_geovals_mult(self,zz)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
real(kind_real),intent(in) :: zz   !< Multiplier

! Multiply GeoVals with a scalar
self%x = zz*self%x

end subroutine aq_geovals_mult
! ------------------------------------------------------------------------------
!> Add GeoVals
subroutine aq_geovals_add(self,other)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
type(aq_geovals),intent(in) :: other   !< Other GeoVals

! Add GeoVals
self%x = self%x+other%x

end subroutine aq_geovals_add
! ------------------------------------------------------------------------------
!> Subtract GeoVals
subroutine aq_geovals_diff(self,other)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
type(aq_geovals),intent(in) :: other   !< Other GeoVals

! Subtract GeoVals
self%x = self%x-other%x

end subroutine aq_geovals_diff
! ------------------------------------------------------------------------------
!> Schur product for GeoVals
subroutine aq_geovals_schurmult(self,other)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
type(aq_geovals),intent(in) :: other   !< Other GeoVals

! Multiply GeoVals
self%x = self%x*other%x

end subroutine aq_geovals_schurmult
! ------------------------------------------------------------------------------
!> Schur division for GeoVals
subroutine aq_geovals_divide(self,other)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
type(aq_geovals),intent(in) :: other   !< Other GeoVals

! Local variables
real(kind_real) :: tol
integer :: jloc

! Set tolerance
tol = epsilon(tol)

! Conditional division
do jloc=1,self%nobs
  if (abs(other%x(jloc))>tol) then
    self%x(jloc) = self%x(jloc)/other%x(jloc)
  else
    self%x(jloc) = 0.0
  endif
enddo

end subroutine aq_geovals_divide
! ------------------------------------------------------------------------------
!> Compute GeoVals RMS
subroutine aq_geovals_rms(self,rms)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self   !< GeoVals
real(kind_real),intent(inout) :: rms !< RMS

! Local variables
integer :: nv

! Initialization
rms = 0.0
nv = 0

! Loop over values
rms = rms+sum(self%x**2)
nv = nv+1

! Normalize and take square-root
rms = sqrt(rms/real(self%nobs*nv,kind_real))

end subroutine aq_geovals_rms
! ------------------------------------------------------------------------------
!> GeoVals dot product
subroutine aq_geovals_dotprod(geovals1,geovals2,prod)

implicit none

! Passed variables
type(aq_geovals),intent(in) :: geovals1       !< GeoVals 1
type(aq_geovals),intent(in) :: geovals2       !< GeoVals 2
real(kind_real),intent(inout) :: prod !< Dot product

! Local variables
integer :: jo,jv

! Check
if (geovals1%nobs/=geovals2%nobs) call abor1_ftn('aq_geovals_dotprod: inconsistent GeoVals sizes')

! Initialization
prod = 0.0

! Dot product
!AQ Temporary stub waiting for geovals files I/O
if (allocated(geovals1%x).and.allocated(geovals2%x)) then
prod = prod+sum(geovals1%x*geovals2%x)
!AQ here too
else
  prod = 1.0 ! Not zero for tests using ratios
end if
!AQ till here

end subroutine aq_geovals_dotprod
! ------------------------------------------------------------------------------
!> Compute GeoVals stats
subroutine aq_geovals_stats(self,kobs,pmin,pmax,prms)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self       !< GeoVals
integer,intent(inout) :: kobs            !< Number of observations
real(kind_real),intent(inout) :: pmin    !< Minimum value
real(kind_real),intent(inout) :: pmax    !< Maximum value
real(kind_real),intent(inout) :: prms    !< RMS

! Local variables
integer :: nv

! Compute GeoVals stats
kobs = self%nobs
if (self%nobs>0) then
  pmin = huge(1.0)
  pmax = -huge(1.0)
  prms = 0.0
  nv = 0
  pmin = min(pmin,minval(self%x))
  pmax = max(pmax,maxval(self%x))
  prms = prms+sum(self%x**2)
  nv = nv+1
  prms = sqrt(prms/real(self%nobs*nv,kind_real))
else
  pmin = 0.0
  pmax = 0.0
  prms = 0.0
end if

end subroutine aq_geovals_stats
! ------------------------------------------------------------------------------
!> Find and locate GeoVals max. value
subroutine aq_geovals_maxloc(self,mxval,mxloc,mxvar)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self          !< GeoVals
real(kind_real),intent(inout) :: mxval      !< Maximum value
integer,intent(inout) :: mxloc              !< Location of maximum value
type(oops_variables),intent(inout) :: mxvar !< Variable of maximum value

! Local variables
integer :: mxloc_arr(1),mxval_tmp
character(len=1) :: var

! Initialization
mxval = -huge(1.0)

!AQ Temporary stub waiting for geovals files I/O
if (allocated(self%x)) then
! Find GeoVals max. value
mxval_tmp = maxval(self%x)
if (mxval_tmp>mxval) then
  mxval = mxval
  mxloc_arr = maxloc(self%x)
  var = 'x'
endif

! Locate GeoVals max. value
mxloc = mxloc_arr(1)
!AQ here too
else
   mxval = 9999.0
   mxloc = 9999.0
end if
!AQ till here
! Set GeoVals max. variable
call mxvar%push_back(var)

end subroutine aq_geovals_maxloc
! ------------------------------------------------------------------------------
!> Read GeoVals from file
subroutine aq_geovals_read_file(self,f_conf)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self             !< GeoVals
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,nobs,x_id,q_id,u_id,v_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Get filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('aq_geovals_read_file: reading '//trim(filename))

! Open NetCDF file
call ncerr(nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))

! Get dimension id
call ncerr(nf90_inq_dimid(ncid,'nobs',nobs_id))

! Get dimension
call ncerr(nf90_inquire_dimension(ncid,nobs_id,len=nobs))

! GeoVals setup
call aq_geovals_setup(self,nobs)

! Get variables ids
call ncerr(nf90_inq_varid(ncid,'x',x_id))

! Get variables
call ncerr(nf90_get_var(ncid,x_id,self%x))

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine aq_geovals_read_file
! ------------------------------------------------------------------------------
!> Write GeoVals to file
subroutine aq_geovals_write_file(self,f_conf)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self !< GeoVals
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,x_id,q_id,u_id,v_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Check allocation
if (.not.self%lalloc) call abor1_ftn('aq_geovals_write_file: geovals not allocated')

! Set filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('aq_geovals_write_file: writing '//trim(filename))
#ifdef AQ_FIX_GeoVals_MPI
! Create NetCDF file
call ncerr(nf90_create(trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nobs',self%nobs,nobs_id))

! Define variables
call ncerr(nf90_def_var(ncid,'x',nf90_double,(/nobs_id/),x_id))

! End definitions
call ncerr(nf90_enddef(ncid))

! Put variables
call ncerr(nf90_put_var(ncid,x_id,self%x))

! Close NetCDF file
call ncerr(nf90_close(ncid))
#endif
end subroutine aq_geovals_write_file
! ------------------------------------------------------------------------------
!> GeoVals analytic initialization
subroutine aq_geovals_analytic_init(self,locs,f_conf)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self             !< GeoVals
type(aq_locs),intent(inout) :: locs            !< Locations
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: iloc
real(kind_real) :: x,y
character(len=30) :: ic
character(len=:),allocatable :: str
real(kind_real), pointer :: lonlat(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

z_field = locs%altitude()
call z_field%data(z)

! Check allocation
if (.not.self%lalloc) call abor1_ftn('aq_geovals_analytic init: geovals not allocated')

! Get analytic configuration
call f_conf%get_or_die("analytic_init",str)
ic = str
call fckit_log%info('aq_geovals_analytic_init: ic = '//trim(ic))
! do iloc=1,locs%nlocs()
!   select case (trim(ic))
!   case ('baroclinic-instability')
!     ! Go to cartesian coordinates
!     call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

!     ! Compute values for baroclinic instability
!     if (self%vars%has('x')) call baroclinic_instability(x,y,z(iloc),'x',self%x(iloc))
!     if (self%vars%has('q')) call baroclinic_instability(x,y,z(iloc),'q',self%q(iloc))
!     if (self%vars%has('u')) call baroclinic_instability(x,y,z(iloc),'u',self%u(iloc))
!     if (self%vars%has('v')) call baroclinic_instability(x,y,z(iloc),'v',self%v(iloc))
!   case ('large-vortices')
!     ! Go to cartesian coordinates
!     call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

!     ! Compute values for large vortices
!     if (self%vars%has('x')) call large_vortices(x,y,z(iloc),'x',self%x(iloc))
!     if (self%vars%has('q')) call large_vortices(x,y,z(iloc),'q',self%q(iloc))
!     if (self%vars%has('u')) call large_vortices(x,y,z(iloc),'u',self%u(iloc))
!     if (self%vars%has('v')) call large_vortices(x,y,z(iloc),'v',self%v(iloc))
!   case default
!     call abor1_ftn('aq_geovals_analytic_init: unknown initialization')
!   endselect
! enddo

call lonlat_field%final()
call z_field%final()

end subroutine aq_geovals_analytic_init
! ------------------------------------------------------------------------------
end module aq_geovals_mod
