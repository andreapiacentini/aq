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

module aq_blas_mod
   use aq_constants_mod

   implicit none

   private
   public :: aq_dot_product, aq_scal, aq_copy, aq_axpy

   interface aq_dot_product
      module procedure aq_ddot_1d
      module procedure aq_ddot_2d
      module procedure aq_ddot_3d
      module procedure aq_ddot_4d
      module procedure aq_sdot_1d
      module procedure aq_sdot_2d
      module procedure aq_sdot_3d
      module procedure aq_sdot_4d
   end interface aq_dot_product

   interface aq_scal
      module procedure aq_sscal_1d
      module procedure aq_dscal_1d
      module procedure aq_sscal_2d
      module procedure aq_dscal_2d
      module procedure aq_sscal_3d
      module procedure aq_dscal_3d
   end interface aq_scal

   interface aq_copy
      module procedure aq_scopy_1d
      module procedure aq_dcopy_1d
      module procedure aq_scopy_2d
      module procedure aq_dcopy_2d
      module procedure aq_scopy_3d
      module procedure aq_dcopy_3d
      module procedure aq_scopy_serialize
      module procedure aq_dcopy_serialize
      module procedure aq_scopy_deserialize
      module procedure aq_dcopy_deserialize
   end interface aq_copy

   interface aq_axpy
      module procedure aq_saxpy_1d
      module procedure aq_daxpy_1d
      module procedure aq_saxpy_2d
      module procedure aq_daxpy_2d
      module procedure aq_saxpy_3d
      module procedure aq_daxpy_3d
   end interface aq_axpy

contains

   real(kind=oops_real) function aq_ddot_1d(n, x1, x2) result(zprod)
      integer           , intent(in) :: n
      real(kind=aq_real), intent(in) :: x1(:)
      real(kind=aq_real), intent(in) :: x2(:)

      real(kind=aq_real) :: ddot
      zprod = ddot(n, x1, 1, x2, 1)

   end function aq_ddot_1d

   real(kind=oops_real) function aq_sdot_1d(n, x1, x2) result(zprod)
      integer             , intent(in) :: n
      real(kind=aq_single), intent(in) :: x1(:)
      real(kind=aq_single), intent(in) :: x2(:)

      real(kind=aq_single) :: sdot
      zprod = real(sdot(n, x1, 1, x2, 1), kind=oops_real)

   end function aq_sdot_1d

   real(kind=oops_real) function aq_ddot_2d(n, x1, x2) result(zprod)
      integer           , intent(in) :: n
      real(kind=aq_real), intent(in) :: x1(:,:)
      real(kind=aq_real), intent(in) :: x2(:,:)

      real(kind=aq_real) :: ddot
      zprod = ddot(n, x1, 1, x2, 1)

   end function aq_ddot_2d

   real(kind=oops_real) function aq_sdot_2d(n, x1, x2) result(zprod)
      integer             , intent(in) :: n
      real(kind=aq_single), intent(in) :: x1(:,:)
      real(kind=aq_single), intent(in) :: x2(:,:)

      real(kind=aq_single) :: sdot
      zprod = real(sdot(n, x1, 1, x2, 1), kind=oops_real)

   end function aq_sdot_2d

   real(kind=oops_real) function aq_ddot_3d(n, x1, x2) result(zprod)
      integer           , intent(in) :: n
      real(kind=aq_real), intent(in) :: x1(:,:,:)
      real(kind=aq_real), intent(in) :: x2(:,:,:)

      real(kind=aq_real) :: ddot
      zprod = ddot(n, x1, 1, x2, 1)

   end function aq_ddot_3d

   real(kind=oops_real) function aq_sdot_3d(n, x1, x2) result(zprod)
      integer             , intent(in) :: n
      real(kind=aq_single), intent(in) :: x1(:,:,:)
      real(kind=aq_single), intent(in) :: x2(:,:,:)

      real(kind=aq_single) :: sdot
      zprod = real(sdot(n, x1, 1, x2, 1), kind=oops_real)

   end function aq_sdot_3d

   real(kind=oops_real) function aq_ddot_4d(n, x1, x2) result(zprod)
      integer           , intent(in) :: n
      real(kind=aq_real), intent(in) :: x1(:,:,:,:)
      real(kind=aq_real), intent(in) :: x2(:,:,:,:)

      real(kind=aq_real) :: ddot
      zprod = ddot(n, x1, 1, x2, 1)

   end function aq_ddot_4d

   real(kind=oops_real) function aq_sdot_4d(n, x1, x2) result(zprod)
      integer             , intent(in) :: n
      real(kind=aq_single), intent(in) :: x1(:,:,:,:)
      real(kind=aq_single), intent(in) :: x2(:,:,:,:)

      real(kind=aq_single) :: sdot
      zprod = real(sdot(n, x1, 1, x2, 1), kind=oops_real)

   end function aq_sdot_4d

   subroutine aq_sscal_1d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(inout) :: x1(:)
      call sscal(n, real(coeff,kind=aq_single), x1, 1)
   end subroutine aq_sscal_1d

   subroutine aq_dscal_1d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(inout) :: x1(:)
      call dscal(n, real(coeff,kind=aq_double), x1, 1)
   end subroutine aq_dscal_1d

   subroutine aq_sscal_2d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(inout) :: x1(:,:)
      call sscal(n, real(coeff,kind=aq_single), x1, 1)
   end subroutine aq_sscal_2d

   subroutine aq_dscal_2d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(inout) :: x1(:,:)
      call dscal(n, real(coeff,kind=aq_double), x1, 1)
   end subroutine aq_dscal_2d

   subroutine aq_sscal_3d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(inout) :: x1(:,:,:)
      call sscal(n, real(coeff,kind=aq_single), x1, 1)
   end subroutine aq_sscal_3d

   subroutine aq_dscal_3d(n, coeff, x1)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(inout) :: x1(:,:,:)
      call dscal(n, real(coeff,kind=aq_double), x1, 1)
   end subroutine aq_dscal_3d

   subroutine aq_scopy_1d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_single), intent(in)    :: x(:)
      real(kind=aq_single), intent(inout) :: y(:)
      call scopy(n, x, 1, y, 1)
   end subroutine aq_scopy_1d

   subroutine aq_dcopy_1d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_double), intent(in)    :: x(:)
      real(kind=aq_double), intent(inout) :: y(:)
      call dcopy(n, x, 1, y, 1)
   end subroutine aq_dcopy_1d

   subroutine aq_scopy_2d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_single), intent(in)    :: x(:,:)
      real(kind=aq_single), intent(inout) :: y(:,:)
      call scopy(n, x, 1, y, 1)
   end subroutine aq_scopy_2d

   subroutine aq_dcopy_2d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_double), intent(in)    :: x(:,:)
      real(kind=aq_double), intent(inout) :: y(:,:)
      call dcopy(n, x, 1, y, 1)
   end subroutine aq_dcopy_2d

   subroutine aq_scopy_3d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_single), intent(in)    :: x(:,:,:)
      real(kind=aq_single), intent(inout) :: y(:,:,:)
      call scopy(n, x, 1, y, 1)
   end subroutine aq_scopy_3d

   subroutine aq_dcopy_3d(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_double), intent(in)    :: x(:,:,:)
      real(kind=aq_double), intent(inout) :: y(:,:,:)
      call dcopy(n, x, 1, y, 1)
   end subroutine aq_dcopy_3d

   subroutine aq_scopy_serialize(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_single), intent(in)    :: x(:,:)
      real(kind=aq_single), intent(inout) :: y(:)
      call scopy(n, x, 1, y, 1)
   end subroutine aq_scopy_serialize

   subroutine aq_dcopy_serialize(n, x, y)
      integer             , intent(in)    :: n
      real(kind=aq_real),   intent(in)    :: x(:,:)
      real(kind=oops_real), intent(inout) :: y(:)
      call dcopy(n, x, 1, y, 1)
   end subroutine aq_dcopy_serialize

   subroutine aq_scopy_deserialize(n, x, y)
      integer             , intent(in)  :: n
      real(kind=aq_single), intent(in)  :: x(:)
      real(kind=aq_single), intent(out) :: y(:,:)
      call scopy(n, x, 1, y, 1)
   end subroutine aq_scopy_deserialize

   subroutine aq_dcopy_deserialize(n, x, y)
      integer             , intent(in)  :: n
      real(kind=oops_real), intent(in)  :: x(:)
      real(kind=aq_real),   intent(out) :: y(:,:)
      call dcopy(n, x, 1, y, 1)
   end subroutine aq_dcopy_deserialize

   subroutine aq_saxpy_1d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(in)    :: x(:)
      real(kind=aq_single), intent(inout) :: y(:)
      call saxpy(n, real(coeff,kind=aq_single), x, 1, y, 1)
   end subroutine aq_saxpy_1d

   subroutine aq_daxpy_1d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(in)    :: x(:)
      real(kind=aq_double), intent(inout) :: y(:)
      call daxpy(n, real(coeff,kind=aq_double), x, 1, y, 1)
   end subroutine aq_daxpy_1d

   subroutine aq_saxpy_2d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(in)    :: x(:,:)
      real(kind=aq_single), intent(inout) :: y(:,:)
      call saxpy(n, real(coeff,kind=aq_single), x, 1, y, 1)
   end subroutine aq_saxpy_2d

   subroutine aq_daxpy_2d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(in)    :: x(:,:)
      real(kind=aq_double), intent(inout) :: y(:,:)
      call daxpy(n, real(coeff,kind=aq_double), x, 1, y, 1)
   end subroutine aq_daxpy_2d

   subroutine aq_saxpy_3d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_single), intent(in)    :: x(:,:,:)
      real(kind=aq_single), intent(inout) :: y(:,:,:)
      call saxpy(n, real(coeff,kind=aq_single), x, 1, y, 1)
   end subroutine aq_saxpy_3d

   subroutine aq_daxpy_3d(n, coeff, x, y)
      integer             , intent(in)    :: n
      real(kind=oops_real), intent(in)    :: coeff
      real(kind=aq_double), intent(in)    :: x(:,:,:)
      real(kind=aq_double), intent(inout) :: y(:,:,:)
      call daxpy(n, real(coeff,kind=aq_double), x, 1, y, 1)
   end subroutine aq_daxpy_3d

end module aq_blas_mod
