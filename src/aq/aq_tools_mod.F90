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
