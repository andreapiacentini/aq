! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_transform_mod

   use atlas_module
   use aq_constants_mod
   use fckit_configuration_module, only: fckit_configuration, fckit_YAMLConfiguration
   use iso_c_binding

   implicit none

   private
   public :: aq_transform

   type aq_transform_params
      real(aq_real), dimension(:), allocatable :: p
   end type aq_transform_params

   type aq_transform
      integer(aq_int) :: nb_vars = 0
      character(len=aq_varlen), dimension(:), allocatable :: var_names
      character(len=aq_strlen), dimension(:), allocatable :: transform
      type(aq_transform_params), dimension(:), allocatable :: params

   contains
      procedure, public :: setup => aq_transform_setup
      procedure, public :: delete => aq_transform_delete
      procedure, public :: copy => aq_transform_copy
      procedure, public :: apply => aq_transform_direct
      procedure, public :: apply_inverse => aq_transform_inverse

   end type aq_transform

contains

   subroutine aq_transform_setup(self, vars, config)
      class(aq_transform), intent(inout)    :: self
      character(len=*),       intent(in)    :: vars(:)
      type(fckit_Configuration), intent(in) :: config

      integer :: ib, ib_t, idx
      character(len=aq_strlen) :: key
      character(len=:), allocatable :: transform
      type(fckit_Configuration) :: lconf
      character(kind=c_char,len=:), allocatable :: json_str
      character(len=aq_varlen) :: var, rvar

      call self%delete()

      ! Specific treatment for strings with dots as in PM2.5
      json_str = config%json()
      do ib = 1, size(vars)
         var = trim(vars(ib))
         rvar = var
         do
            idx = index(rvar, '.') ; if (idx == 0) exit
            rvar(idx:idx) = '_'
         end do
         if (trim(var) /= trim(rvar)) then
            do
               idx = index(json_str, trim(var)) ; if (idx == 0) exit
               json_str(idx:idx+len(trim(var))-1) = trim(rvar)
            end do
         end if
      end do
      lconf = fckit_YAMLConfiguration(json_str)

      self%nb_vars = 0
      do ib = 1, size(vars)
         var = trim(vars(ib))
         rvar = var
         do
            idx = index(rvar, '.') ; if (idx == 0) exit
            rvar(idx:idx) = '_'
         end do
         if (lconf%has(trim(rvar))) self%nb_vars = self%nb_vars + 1
      end do

      allocate(self%var_names(self%nb_vars))
      allocate(self%transform(self%nb_vars))
      allocate(self%params(self%nb_vars))

      ib_t = 0
      do ib = 1, size(vars)
         var = trim(vars(ib))
         rvar = var
         do
            idx = index(rvar, '.') ; if (idx == 0) exit
            rvar(idx:idx) = '_'
         end do
         if (lconf%has(trim(rvar))) then
            ib_t = ib_t + 1
            self%var_names(ib_t) = trim(vars(ib))
            key = trim(rvar)//'.method'
            call lconf%get_or_die(trim(key),transform)
            self%transform(ib_t) = trim(transform)
            key = trim(rvar)//'.parameters'
            if (lconf%has(trim(key))) then
               call lconf%get_or_die(trim(key),self%params(ib_t)%p)
            end if
         end if
      end do

   end subroutine aq_transform_setup

   subroutine aq_transform_delete(self)
      class(aq_transform), intent(inout) :: self

      integer :: ib

      if (allocated(self%var_names)) then
         deallocate(self%var_names)
         deallocate(self%transform)
         do ib = 1, size(self%params)
            deallocate(self%params(ib)%p)
         end do
         deallocate(self%params)
      end if
      self%nb_vars = 0

   end subroutine aq_transform_delete

   subroutine aq_transform_copy(self, other)
      class(aq_transform), intent(inout) :: self
      class(aq_transform), intent(in)    :: other

      integer :: ib_t, ib

      call self%delete()

      if (other%nb_vars == 0) return

      self%nb_vars = other%nb_vars

      allocate(self%var_names(self%nb_vars))
      allocate(self%transform(self%nb_vars))
      allocate(self%params(self%nb_vars))

      self%var_names(:) = other%var_names(:)
      self%transform(:) = other%transform(:)
      do ib_t = 1, self%nb_vars
         allocate(self%params(ib_t)%p(size(other%params(ib_t)%p)))
         self%params(ib_t)%p(:) = other%params(ib_t)%p(:)
      end do

   end subroutine aq_transform_copy

   subroutine aq_transform_direct(self,afset,vars)
      class(aq_transform), intent(in) :: self
      class(atlas_FieldSet),  intent(inout) :: afset
      character(len=*),       intent(in)    :: vars(:)

      type(atlas_Field) :: afld
      real(aq_single), pointer :: flds(:,:)
      real(aq_real), pointer :: fldd(:,:)

      integer(atlas_kind_idx) :: ib_t, ib_i, ib_j, ib_k, ib_var, il_lev
      logical :: ll_sgl

      if (self%nb_vars == 0) return

      afld = afset%field(trim(vars(1)))
      ll_sgl = afld%kind() == aq_single
      do ib_var = 1, size(vars)
         do ib_t = 1, self%nb_vars
            if (trim(vars(ib_var)) == trim(self%var_names(ib_t))) then
               afld = afset%field(trim(vars(ib_var)))
               select case(trim(self%transform(ib_t)))
               case('log_threshold')
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = log10(flds(:,:)+real(self%params(ib_t)%p(1),kind=aq_single))
                  else
                     call afld%data(fldd)
                     fldd(:,:) = log10(fldd(:,:)+self%params(ib_t)%p(1))
                  end if
               case('log_bounded')
                  !AQ CHECK FORMULAE
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = log10(flds(:,:)/ &
                        & (real(self%params(ib_t)%p(1),kind=aq_single)-flds(:,:)))
                  else
                     call afld%data(fldd)
                     fldd(:,:) = log10(fldd(:,:)/(self%params(ib_t)%p(1)-fldd(:,:)))
                  end if
               case('scale')
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = flds(:,:)*real(self%params(ib_t)%p(1),kind=aq_single)
                  else
                     call afld%data(fldd)
                     fldd(:,:) = fldd(:,:)*self%params(ib_t)%p(1)
                  end if
               case default
                  call abor1_ftn('Transformation '//trim(self%transform(ib_t))//&
                     &' for variable '//trim(vars(ib_var))//' not recognised')
               end select
               call afld%final()
               exit
            end if
         end do
      end do

   end subroutine aq_transform_direct

   subroutine aq_transform_inverse(self,afset,vars)
      class(aq_transform), intent(in) :: self
      class(atlas_FieldSet),  intent(inout) :: afset
      character(len=*),       intent(in)    :: vars(:)

      type(atlas_Field) :: afld
      real(aq_single), pointer :: flds(:,:)
      real(aq_real), pointer :: fldd(:,:)

      integer(atlas_kind_idx) :: ib_t, ib_i, ib_j, ib_k, ib_var, il_lev
      logical :: ll_sgl

      if (self%nb_vars == 0) return

      afld = afset%field(trim(vars(1)))
      ll_sgl = afld%kind() == aq_single
      do ib_var = 1, size(vars)
         do ib_t = 1, self%nb_vars
            if (trim(vars(ib_var)) == trim(self%var_names(ib_t))) then
               afld = afset%field(trim(vars(ib_var)))
               select case(trim(self%transform(ib_t)))
               case('log_threshold')
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = 10.0**flds(:,:)-real(self%params(ib_t)%p(1),kind=aq_single)
                  else
                     call afld%data(fldd)
                     fldd(:,:) = 10.0**fldd(:,:)-self%params(ib_t)%p(1)
                  end if
               case('log_bounded')
                  !AQ CHECK FORMULAE
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = (10.0**flds(:,:)*real(self%params(ib_t)%p(1),kind=aq_single))/&
                        & (1.0+10.0**flds(:,:))
                  else
                     call afld%data(fldd)
                     fldd(:,:) = (10.0**fldd(:,:)*self%params(ib_t)%p(1))/(1.0+10.0**fldd(:,:))
                  end if
               case('scale')
                  if (ll_sgl) then
                     call afld%data(flds)
                     flds(:,:) = flds(:,:)/real(self%params(ib_t)%p(1),kind=aq_single)
                  else
                     call afld%data(fldd)
                     fldd(:,:) = fldd(:,:)/self%params(ib_t)%p(1)
                  end if
               case default
                  call abor1_ftn('Transformation '//trim(self%transform(ib_t))//&
                     &' for variable '//trim(vars(ib_var))//' not recognised')
               end select
               call afld%final()
               exit
            end if
         end do
      end do

   end subroutine aq_transform_inverse

end module aq_transform_mod
