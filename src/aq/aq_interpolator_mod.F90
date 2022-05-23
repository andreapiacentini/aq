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

module aq_interpolator_mod

use aq_constants_mod
use aq_geom_mod
use aq_fields_mod

use fckit_module
use fckit_mpi_module
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use oops_variables_mod

use interp_matrix_structure_mod, only : csr_format
use space_time_operator_mod, only : observ_operator
use matrix_manipulations, only : multiply_matrix_csr_vector, addmult_matrixt_csr_vector, deallocate_operator

implicit none

private
public :: aq_interpolator

type :: aq_interpolator
   integer :: loc_nlocs
   type(csr_format) :: Hmat
contains
   procedure, public :: create  => aq_interpolator_create
   procedure, public :: final   => aq_interpolator_delete
   procedure, public :: delete  => aq_interpolator_delete
   procedure, public :: apply   => aq_interpolator_apply
   procedure, public :: applyAD => aq_interpolator_applyAD
   !
   final :: aq_interpolator_final_auto
   !
end type aq_interpolator
! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------
!> Create and store an interpolator
subroutine aq_interpolator_create(self, geom, loc_nlocs, lats, lons)
   !
   implicit none
   !
   class(aq_interpolator), intent(inout) :: self
   type(aq_geom), intent(in)             :: geom
   integer, intent(in)                   :: loc_nlocs
   real(kind_real), intent(in)           :: lats(loc_nlocs)
   real(kind_real), intent(in)           :: lons(loc_nlocs)

   ! Local variables
   character(len=aq_strlen) :: msg
   integer :: ib
   real(kind_real), dimension(1,1) :: dummylev
   real(kind_real), dimension(1) :: dummycoord
   real(kind_real), dimension(:), allocatable :: dummytime

   ! Interpolator needed only where there are observations
   self%loc_nlocs = loc_nlocs
   if ( loc_nlocs == 0 ) return
   !
   allocate(dummytime(loc_nlocs))
   dummytime(:) = 0_kind_real

   call observ_operator ( &
      &   geom%ny, &
      &   geom%nx, &
      &   1, & ! Only on input level in the gathered surface field
      &   1, & ! Only one exact time (no time interpolation)
      &   geom%lats, &
      &   geom%lons, &
      &   dummycoord, & ! Vert coord not relevat
      &   dummycoord, & ! Obs time not relevant
      &   loc_nlocs, &
      &   1, & ! Levels in and out are 1
      &   1, &
      &   lats, &
      &   lons, &
      &   .false., & ! aq grid not considered as lon periodic
      &   dummytime, &
      &   dummylev, &
      &   2, & ! It is the ground interpolator option
      &   .false., & ! no input scaling
      &   .false., & ! no averaging kernel
      &   .false., & ! no time interpolation
      &   self%Hmat)

   write(msg,'(A)') 'Built interpolator'
   call fckit_log%debug(msg)

   deallocate(dummytime)

end subroutine aq_interpolator_create
! ------------------------------------------------------------------------------
!> Delete interpolators
subroutine aq_interpolator_final_auto(self)
   !
   type(aq_interpolator), intent(inout) :: self
   !
   call self%delete()
   !
end subroutine aq_interpolator_final_auto

subroutine aq_interpolator_delete(this)
   !
   class(aq_interpolator), intent(inout) :: this
   !
   if (this%loc_nlocs > 0) call deallocate_operator(this%Hmat)
   this%loc_nlocs = 0
   !
end subroutine aq_interpolator_delete
! ------------------------------------------------------------------------------
!> Apply the direct interpolator (acting both as linear and nonlinear)
subroutine aq_interpolator_apply(self, field, vars, mask, vals)
   implicit none
   class(aq_interpolator), intent(in) :: self
   class(aq_fields), intent(in)       :: field
   type(oops_variables), intent(in)   :: vars
   integer(c_int), intent(in)         :: mask(:)
   real(c_double), intent(inout)      :: vals(:)

   ! Local variables
   integer       :: nlev = 1
   integer       :: glo_nlocs, offset, jvar, ib
   real(aq_real), allocatable, dimension(:,:) :: surf_fld
   character(len=aq_strlen) :: fname
   !
   if (trim(field%geom%orientation) == 'down') nlev = field%geom%levels
   call field%fmpi%allreduce(self%loc_nlocs,glo_nlocs,fckit_mpi_sum())
   if ( glo_nlocs == 0 ) return
   allocate(surf_fld(field%geom%grid%nx(1),field%geom%grid%ny()))

   offset = 0
   do jvar=1,vars%nvars()
      fname = vars%variable(jvar)

      ! Gather field on master
      call field%gather_var_at_lev(trim(fname), nlev, surf_fld, 0)
      ! If the geometry as a halo, perform the interpolation in parallel
      ! otherwise assume that by construction loc_nlocs is >0 only on the master
      if ( field%geom%halo >= 1 ) call field%fmpi%broadcast(surf_fld, root=0)

      ! Local Interpolation
      if ( self%loc_nlocs > 0 ) then
         call multiply_matrix_csr_vector( &
            &   self%Hmat, &
            &   pack(surf_fld,.true.), &
            &   1, &
            &   self%loc_nlocs, &
            &   vals(offset+1:offset+self%loc_nlocs))
         do ib = 1, self%loc_nlocs
            if (mask(ib)==0) vals(offset+ib) = aq_missing_value
         end do
         ! Update offset
         offset = offset+self%loc_nlocs
      end if
   enddo

   if (size(vals) /= offset) call abor1_ftn('aq_interpolator_apply: error size')

   ! Release memory
   deallocate(surf_fld)

end subroutine aq_interpolator_apply
! ------------------------------------------------------------------------------
!> Apply the adjoint interpolator
subroutine aq_interpolator_applyAD(self, field, vars, mask, vals)
   implicit none
   class(aq_interpolator), intent(in) :: self
   class(aq_fields), intent(inout)    :: field
   type(oops_variables), intent(in)   :: vars
   integer(c_int), intent(in)         :: mask(:)
   real(c_double), intent(in)         :: vals(:)

   ! Local variables
   integer       :: nlev = 1
   integer       :: loc_nlocs, glo_nlocs, offset, jvar, ib
   real(aq_real), allocatable, dimension(:) :: surf_1d(:)
   real(aq_real), allocatable, dimension(:,:) :: surf_fld
   real(aq_real), allocatable, dimension(:) :: masked_vals(:)
   character(len=aq_strlen) :: fname
   real(aq_real) :: filter_val

   if (trim(field%geom%orientation) == 'down') nlev = field%geom%levels

   call field%fmpi%allreduce(self%loc_nlocs,glo_nlocs,fckit_mpi_sum())
   if ( glo_nlocs == 0 ) return
   allocate(surf_fld(field%geom%grid%nx(1),field%geom%grid%ny()))

   offset = 0
   if ( self%loc_nlocs > 0 ) then
      allocate(masked_vals(self%loc_nlocs))
      allocate(surf_1d(field%geom%grid%nx(1)*field%geom%grid%ny()))
   end if
   do jvar=1,vars%nvars()
      fname = vars%variable(jvar)

      surf_fld(:,:) = 0_kind_real
      if ( self%loc_nlocs > 0 ) then
         do ib = 1, self%loc_nlocs
            if (mask(ib)>0) then
               masked_vals(ib) = vals(offset+ib)
            else
               masked_vals(ib) = aq_missing_value
            end if
         end do
         surf_1d = 0_kind_real
         call addmult_matrixt_csr_vector( &
            &   self%Hmat, &
            &   masked_vals, &
            &   1, &
            &   self%loc_nlocs, &
            &   surf_1d, &
            &   aq_missing_value)

         surf_fld = unpack(surf_1d,surf_fld==0_kind_real,surf_fld)
      end if

      if ( field%geom%halo >= 1 ) call field%fmpi%allreduce(surf_fld,fckit_mpi_sum())
      call field%scatteradd_var_at_lev(trim(fname), nlev, surf_fld, 0)

      ! Update offset
      offset = offset+self%loc_nlocs
   enddo

   if (size(vals) /= offset) call abor1_ftn('aq_interpolator_applyAD: error size')

   ! Release memory
   deallocate(surf_fld)
   if ( self%loc_nlocs > 0 ) then
      deallocate(masked_vals)
      deallocate(surf_1d)
   end if

end subroutine aq_interpolator_applyAD
!
end module aq_interpolator_mod
