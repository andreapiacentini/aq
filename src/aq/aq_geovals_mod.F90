! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2021-2022 CERFACS.
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
use fckit_mpi_module
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
        & aq_geovals_add,aq_geovals_diff,aq_geovals_schurmult,aq_geovals_divide,aq_geovals_rms,aq_geovals_dotprod, &
        & aq_geovals_stats,aq_geovals_maxloc,aq_geovals_read_file, aq_geovals_write_file,aq_geovals_analytic_init, &
        & aq_geovals_fill, aq_geovals_fillad
! ------------------------------------------------------------------------------
type :: aq_geovals
  integer :: nobs                      !< Number of observations
  real(kind_real), allocatable :: x(:) !< Chemical observations values
  logical :: lalloc = .false.          !< Allocation flag
  type(oops_variables) :: vars         !< Variables
  type(fckit_mpi_comm) :: fmpi         !< Communicator
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
subroutine aq_geovals_fill(self, c_nloc, c_indx, c_nval, c_vals)
implicit none
type(aq_geovals), intent(inout) :: self
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nval
real(c_double), intent(in) :: c_vals(c_nval)

integer :: jvar, jloc, iloc, ii

if (.not.self%lalloc) call abor1_ftn('aq_geovals_fill: gom not allocated')

ii = 0
do jvar=1,self%vars%nvars()
  do jloc=1,c_nloc
    iloc = c_indx(jloc)
    ii = ii + 1
    self%x(iloc) = c_vals(ii)
  enddo
enddo
! if (ii /= c_nval) call abor1_ftn('aq_geovals_fill: error size')

end subroutine aq_geovals_fill
! ------------------------------------------------------------------------------
subroutine aq_geovals_fillad(self, c_nloc, c_indx, c_nval, c_vals)
implicit none
type(aq_geovals), intent(in) :: self
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nval
real(c_double), intent(inout) :: c_vals(c_nval)

integer :: jvar, jloc, iloc, ii

if (.not.self%lalloc) call abor1_ftn('aq_geovals_fillad: gom not allocated')

ii = 0
do jvar=1,self%vars%nvars()
  do jloc=1,c_nloc
    iloc = c_indx(jloc)
    ii = ii + 1
    c_vals(ii) = self%x(iloc)
  enddo
enddo
! if (ii /= c_nval) call abor1_ftn('aq_geovals_fillad: error size')

end subroutine aq_geovals_fillad
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
real(kind_real) :: norm

! Initialization
rms = 0.0
nv = 0

if (self%nobs>0) then
   ! Loop over values
   rms = rms+sum(self%x**2)
   nv = nv+1

   ! Total number of values
   norm = real(self%nobs*nv,kind_real)
else
   ! No value
   rms = 0.0
   norm = 0.0
end if

! Allreduce
call self%fmpi%allreduce(rms,fckit_mpi_sum())
call self%fmpi%allreduce(norm,fckit_mpi_sum())

! Normalize and take square-root
rms = sqrt(rms/norm)

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

! Dot product
if ((geovals1%nobs>0).and.(geovals2%nobs>0)) then
   prod = sum(geovals1%x*geovals2%x)
else
   prod = 0.0
end if
call geovals1%fmpi%allreduce(prod,fckit_mpi_sum())

end subroutine aq_geovals_dotprod
! ------------------------------------------------------------------------------
!> Compute GeoVals stats
subroutine aq_geovals_stats(self,kobs,pmin,pmax,pave,pstd)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self       !< GeoVals
integer,intent(inout) :: kobs            !< Number of observations
real(kind_real),intent(inout) :: pmin    !< Minimum value
real(kind_real),intent(inout) :: pmax    !< Maximum value
real(kind_real),intent(inout) :: pave    !< Mean
real(kind_real),intent(inout) :: pstd    !< StdDev

! Compute GeoVals stats
kobs = self%nobs
call self%fmpi%allreduce(kobs,fckit_mpi_sum())
if (kobs>0) then
  pmin = huge(1.0)
  if (self%nobs>0) pmin = min(pmin,minval(self%x))
  call self%fmpi%allreduce(pmin,fckit_mpi_min())
  pmax = -huge(1.0)
  if (self%nobs>0) pmax = max(pmax,maxval(self%x))
  call self%fmpi%allreduce(pmax,fckit_mpi_max())
  pave = 0.0
  if (self%nobs>0) pave = sum(self%x)
  call self%fmpi%allreduce(pave,fckit_mpi_sum())
  pave = pave/real(kobs,kind_real)
  pstd = 0.0
  if (self%nobs>0) pstd = sum((self%x-pave)**2)
  call self%fmpi%allreduce(pstd,fckit_mpi_sum())
  pstd = sqrt(pstd/real(kobs,kind_real))
else
  pmin = 0.0
  pmax = 0.0
  pave = 0.0
  pstd = 0.0
end if

end subroutine aq_geovals_stats
! ------------------------------------------------------------------------------
!> Find and locate GeoVals max. value
subroutine aq_geovals_maxloc(self,mxval,mxloc,mxvar)

implicit none

! Passed variables
type(aq_geovals),intent(inout) :: self      !< GeoVals
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

! Input on master proc only
if (self%fmpi%rank() == 0) then
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
else
   nobs = 0
end if

! GeoVals setup
call aq_geovals_setup(self,nobs)

if (self%fmpi%rank() == 0) then
   ! Get variables ids
   call ncerr(nf90_inq_varid(ncid,'x',x_id))

   ! Get variables
   call ncerr(nf90_get_var(ncid,x_id,self%x))

   ! Close NetCDF file
   call ncerr(nf90_close(ncid))
end if

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

! Output from master proc only
if ( self%fmpi%rank() == 0) then
   ! Set filename
   call f_conf%get_or_die("filename",str)
   filename = str
   call fckit_log%info('aq_geovals_write_file: writing '//trim(filename))

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
end if

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
real(kind_real), pointer :: xy(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

real(kind_real), parameter :: dp_pi=3.14159265359
real(kind_real), parameter :: dLon0 = 6.3
real(kind_real), parameter :: dLat0 = 0.8
real(kind_real), parameter :: dR0   = 3.0
real(kind_real), parameter :: dD    = 10.0
real(kind_real), parameter :: dT    = 6.0
real(kind_real) :: dp_conv

real(kind_real) :: dSinC, dCosC, dCosT, dSinT
real(kind_real) :: dTrm, dX, dY, dZ
real(kind_real) :: dlon, dlat
real(kind_real) :: dRho, dVt, dOmega

dp_conv = dp_pi/180._kind_real
dSinC = sin( dLat0 )
dCosC = cos( dLat0 )

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(xy)

z_field = locs%altitude()
call z_field%data(z)

! Check allocation
if (.not.self%lalloc) call abor1_ftn('aq_geovals_analytic init: geovals not allocated')

do iloc=1,locs%nlocs()
   ! Find the rotated longitude and latitude of a point on a sphere
   !    with pole at (dLon0, dLat0).
   dCosT = cos( xy(2,iloc)*dp_conv )
   dSinT = sin( xy(2,iloc)*dp_conv )

   dTrm = dCosT * cos( xy(1,iloc)*dp_conv - dLon0 )
   dX   = dSinC * dTrm - dCosC * dSinT
   dY   = dCosT * sin( xy(1,iloc)*dp_conv - dLon0 )
   dZ   = dSinC * dSinT + dCosC * dTrm

   dlon = atan2( dY, dX )
   if( dlon < 0.0_kind_real ) dlon = dlon + 2.0_kind_real * dp_pi
   dlat = asin( dZ )

   dRho = dR0 * cos(dlat)
   dVt = 3.0_kind_real * sqrt(3.0_kind_real)/2.0_kind_real/cosh(dRho)/cosh(dRho)*tanh(dRHo)
   if (dRho == 0.0_kind_real) then
      dOmega = 0.0_kind_real
   else
      dOmega = dVt / dRho
   end if

   self%x(iloc) = 2.0_kind_real * &
      & ( 1.0_kind_real + tanh( dRho / dD * sin( dLon - dOmega * dT ) ) )
end do

call lonlat_field%final()
call z_field%final()

end subroutine aq_geovals_analytic_init
! ------------------------------------------------------------------------------
end module aq_geovals_mod
