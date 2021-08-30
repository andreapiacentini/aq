! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_tools_mod

use netcdf

implicit none

private
public :: ncerr
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Check NetCDF status
subroutine ncerr(info)

implicit none

! Passed variables
integer,intent(in) :: info !< Info index

! Check status
if (info/=nf90_noerr) call abor1_ftn(trim(nf90_strerror(info)))

end subroutine ncerr

end module aq_tools_mod
