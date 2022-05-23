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

module aq_constants_mod
   use kinds
   use mpi

   implicit none

   private
   public :: aq_strlen, aq_varlen
   public :: oops_int, oops_long, oops_single, oops_real
   public :: aq_int, aq_long, aq_single, aq_real, aq_double
   public :: pi, deg_to_rad, rad_to_deg, req
   public :: aq_missing_value

   integer, parameter :: aq_strlen = 256
   integer, parameter :: aq_varlen = 100 ! as MAXVARLEN in oops_variables_mod

   integer, parameter :: oops_int    = kind_int
   integer, parameter :: oops_long   = kind_long
   integer, parameter :: oops_single = kind_single
   integer, parameter :: oops_real   = kind_real

   integer, parameter :: aq_int = kind_int
   integer, parameter :: aq_long = kind_long
   integer, parameter :: aq_single = kind_single
   integer, parameter :: aq_real = kind_real
   integer, parameter :: aq_double = kind_real

   ! Mathematical parameters
   real(aq_real),parameter :: pi = acos(-1.0)                    !< Pi
   real(aq_real),parameter :: deg_to_rad = pi/180.0_aq_real      !< Degrees to radians
   real(aq_real),parameter :: rad_to_deg = 180.0_aq_real/pi      !< Radians to degrees

   ! Geophysical parameters
   real(aq_real),parameter :: req = 6371229.0_aq_real            !< Earth radius at equator (m)

   ! Missing values
   real(kind_real),parameter :: aq_missing_value = -1.0d+38

end module aq_constants_mod
