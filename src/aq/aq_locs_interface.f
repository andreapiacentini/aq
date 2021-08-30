! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  This Define interfaces for accessing C++ LocationsAQ objects from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
function aq_locs_nlocs_c(locs) bind(C,name="aq_locs_nlocs_f90")
    use iso_c_binding, only: c_ptr, c_int
    integer(c_int) :: aq_locs_nlocs_c
    type(c_ptr), value :: locs
end function

function aq_locs_lonlat_c(locs) bind(C,name="aq_locs_lonlat_f90")
    use iso_c_binding, only: c_ptr
    type(c_ptr) :: aq_locs_lonlat_c
    type(c_ptr), value :: locs
end function

function aq_locs_altitude_c(locs) bind(C,name="aq_locs_altitude_f90")
    use iso_c_binding, only: c_ptr
    type(c_ptr) :: aq_locs_altitude_c
    type(c_ptr), value :: locs
end function

function aq_locs_times_c(locs, idx) bind(C,name="aq_locs_times_f90")
    use iso_c_binding, only: c_ptr, c_size_t
    type(c_ptr) :: aq_locs_times_c
    type(c_ptr), value :: locs
    integer(c_size_t) idx
end function

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
