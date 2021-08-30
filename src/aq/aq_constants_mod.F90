module aq_constants_mod
   use kinds
   use mpi
   
   implicit none

   private
   public :: aq_strlen, aq_varlen
   public :: oops_int, oops_long, oops_single, oops_real
   public :: aq_int, aq_long, aq_single, aq_real, aq_double
   public :: pi, deg_to_rad, rad_to_deg, req

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

end module aq_constants_mod

