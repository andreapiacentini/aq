! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_fields_mod

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_module
use fckit_mpi_module
use mpi
use datetime_mod
use duration_mod
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use missing_values_mod
!$ use omp_lib
use oops_variables_mod
use aq_geom_mod
use aq_constants_mod
use aq_blas_mod
use aq_field_io_mod
use aq_transform_mod
!AP use aq_interp_mod
!AP use aq_locs_mod
use random_mod

implicit none

private
public :: aq_fields
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 7 !< Random seed (for reproducibility)

type aq_fld_ptr_d
  real(kind=aq_real), pointer :: fld(:,:)
end type aq_fld_ptr_d

type aq_fld_ptr_s
  real(kind=aq_single), pointer :: fld(:,:)
end type aq_fld_ptr_s

type, extends(atlas_FieldSet) :: aq_fields
  integer(atlas_kind_idx)       :: n_vars
  character(len=:), allocatable :: fs_name
  type(datetime)                :: date
  character(len=:), allocatable :: var_name(:)
  type(aq_transform)            :: var_transf
  type(aq_geom)                 :: geom
  type(fckit_mpi_comm)          :: fmpi
  integer                       :: prec
  integer                       :: mpi_kind
  type(aq_fld_ptr_d), allocatable :: fldsd(:)
  type(aq_fld_ptr_s), allocatable :: fldss(:)
  integer                       :: locsize
contains
  procedure, public :: create                => aq_field_create
  procedure, public :: create_from           => aq_field_create_from_other
  procedure, public :: final                 => aq_field_delete
  procedure, public :: delete                => aq_field_delete
  procedure, public :: zero                  => aq_field_zero
  procedure, public :: ones                  => aq_field_ones
  procedure, public :: dirac                 => aq_field_dirac
  procedure, public :: random                => aq_field_random
  procedure, public :: copy                  => aq_field_copy
  procedure, public :: self_add              => aq_field_self_add
  procedure, public :: self_sub              => aq_field_self_sub
  procedure, public :: self_mul              => aq_field_self_mul
  procedure, public :: add_incr              => aq_field_add_incr
  procedure, public :: diff_incr             => aq_field_diff_incr
  procedure, public :: axpy                  => aq_field_self_axpy
  procedure, public :: schur                 => aq_field_self_schur
  procedure, public :: dot_prod_with         => aq_field_dot_prod_with
  procedure, public :: norm                  => aq_field_norm2
  generic,   public :: stats                 => aq_field_stats_tot, &
     &                                          aq_field_stats_per_var, &
     &                                          aq_field_stats_per_var_lev
  procedure, public :: info                  => aq_field_info
  procedure, public :: write                 => aq_field_write
  procedure, public :: read                  => aq_field_read
  procedure, public :: analytic_IC           => aq_field_ana_IC
  procedure, public :: var                   => aq_field_var
  procedure, public :: gather_var_at_lev     => aq_field_gather_var_at_lev
  procedure, public :: scatteradd_var_at_lev => aq_field_scatteradd_var_at_lev
  procedure, public :: serialsize            => aq_field_serial_size
  procedure, public :: serialize             => aq_field_serialize
  procedure, public :: deserialize           => aq_field_deserialize
  generic,   public :: serialize_prec        => aq_field_serialize_single, &
     &                                          aq_field_serialize_real
  generic,   public :: deserialize_prec      => aq_field_deserialize_single, &
     &                                          aq_field_deserialize_real
  procedure, public :: set_atlas             => aq_field_set_atlas
  procedure, public :: to_atlas              => aq_field_to_atlas
  !
  procedure, public :: kind                  => aq_field_kind
  procedure, public :: name                  => aq_field_name
  procedure, public :: rename                => aq_field_rename
  !
  procedure, private :: idx_var              => aq_field_idx_var
  procedure, private :: aq_field_stats_tot, &
     &                  aq_field_stats_per_var, &
     &                  aq_field_stats_per_var_lev
  procedure, private :: aq_field_serialize_single, &
     &                  aq_field_serialize_real, &
     &                  aq_field_deserialize_single, &
     &                  aq_field_deserialize_real
  !
  final :: aq_field_final_auto
  !
end type aq_fields

interface aq_field_stats
  module procedure aq_field_stats_tot
  module procedure aq_field_stats_per_var
  module procedure aq_field_stats_per_var_lev
end interface aq_field_stats
! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------
!> Create fields from geometry and variables
subroutine aq_field_create(self, geom, vars, name, date, kind)
  !
  class(aq_fields),           intent(inout) :: self
  type(aq_geom),              intent(in)    :: geom
  type(oops_variables),       intent(in)    :: vars
  character(len=*), optional, intent(in)    :: name
  type(datetime), optional,   intent(in)    :: date
  integer, optional,          intent(in)    :: kind
  !
  integer(atlas_kind_idx) :: ib_var
  type(atlas_Field) :: afld
  !
  if (present(kind)) then
     self%prec = kind
     if (kind == aq_single) then
        self%mpi_kind = MPI_FLOAT
     else
        self%mpi_kind = MPI_DOUBLE
     end if
  else
     self%prec = aq_real
     self%mpi_kind = MPI_DOUBLE
  end if
  !
  if (present(name)) then
     self%fs_name = trim(name)
  else
     self%fs_name = "aq_fieldset"
  end if
  !
  if (present(date)) then
     self%date = date
  else
     call datetime_create("1970-01-01T00:00:00Z", self%date)
  end if
  !
  if ( self%is_null() ) &
     & self = atlas_FieldSet(self%fs_name)
  !
  call self%geom%clone(geom)
  self%fmpi = geom%fmpi
  self%locsize = self%geom%fs%size()*self%geom%levels
  !
  self%n_vars = vars%nvars()
  if (allocated(self%var_name) .and. self%fmpi%rank() == 0) &
     & print '(A)','Warning: creating an aq field on top a previous non finalized one'
  self%var_name = vars%varlist()
  !
  select case(self%prec)
  case(aq_single)
     allocate(self%fldss(self%n_vars))
     do ib_var = 1, self%n_vars
        afld = self%geom%fs%create_field(name=trim(self%var_name(ib_var)),  &
           &                                     kind=atlas_real(aq_single))
        call self%add(afld)
        call afld%data(self%fldss(ib_var)%fld)
        call afld%final()
     end do
  case(aq_real)
     allocate(self%fldsd(self%n_vars))
     do ib_var = 1, self%n_vars
        afld = self%geom%fs%create_field(name=trim(self%var_name(ib_var)),  &
           &                                     kind=atlas_real(aq_real))
        call self%add(afld)
        call afld%data(self%fldsd(ib_var)%fld)
        call afld%final()
     end do
  end select
  !
  call self%zero()
  !
end subroutine aq_field_create
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine aq_field_create_from_other(self, other, name, date, kind)
  class(aq_fields),           intent(inout) :: self
  type(aq_fields),            intent(in)    :: other
  character(len=*), optional, intent(in)    :: name
  type(datetime), optional,   intent(in)    :: date
  integer, optional,          intent(in)    :: kind
  !
  integer(atlas_kind_idx) :: ib_var
  type(atlas_Field) :: afld
  !
  if (present(name)) then
     self%fs_name = trim(name)
  else
     self%fs_name = other%name()
  end if
  !
  if (present(date)) then
     self%date = date
  else
     self%date = other%date
  end if
  !
  if ( self%is_null() ) &
     & self = atlas_FieldSet(self%fs_name)
  !
  call self%geom%clone(other%geom)
  self%fmpi = self%geom%fmpi
  self%locsize = self%geom%fs%size()*self%geom%levels
  !
  self%n_vars = other%n_vars
  if (allocated(self%var_name)) then
     if (self%fmpi%rank() == 0) &
        & print '(A)','Warning: create_from_other''ing an aq field on top a previous non finalized one'
     deallocate(self%var_name)
  end if
  allocate(character(aq_varlen)::self%var_name(self%n_vars))
  self%var_name(:) = other%var_name(:)
  !
  call self%var_transf%copy(other%var_transf)
  !
  if (present(kind)) then
     self%prec = kind
     if (kind == aq_single) then
        self%mpi_kind = MPI_FLOAT
     else
        self%mpi_kind = MPI_DOUBLE
     end if
  else
     self%prec = other%prec
     self%mpi_kind = other%mpi_kind
  end if
  select case(self%prec)
  case(aq_single)
     if (allocated(self%fldss)) deallocate(self%fldss)
     allocate(self%fldss(self%n_vars))
     do ib_var = 1, self%n_vars
        afld = self%geom%fs%create_field(name=trim(self%var_name(ib_var)),  &
           &                                     kind=atlas_real(aq_single))
        call self%add(afld)
        call afld%data(self%fldss(ib_var)%fld)
        call afld%final()
     end do
  case(aq_real)
     if (allocated(self%fldsd)) deallocate(self%fldsd)
     allocate(self%fldsd(self%n_vars))
     do ib_var = 1, self%n_vars
        afld = self%geom%fs%create_field(name=trim(self%var_name(ib_var)),  &
           &                                     kind=atlas_real(aq_real))
        call self%add(afld)
        call afld%data(self%fldsd(ib_var)%fld)
        call afld%final()
     end do
  end select
  !
  call self%zero()
  !
end subroutine aq_field_create_from_other
! ------------------------------------------------------------------------------
!> Delete fields
subroutine aq_field_final_auto(self)
  type(aq_fields), intent(inout) :: self
  !
  call self%delete()
  !
end subroutine aq_field_final_auto

subroutine aq_field_delete(this)
  class(aq_fields), intent(inout) :: this
  !
  integer(atlas_kind_idx) :: ib_var
  !
  call this%var_transf%delete()
  if (allocated(this%var_name)) deallocate(this%var_name)
  if (allocated(this%fldsd)) then
     do ib_var = 1, this%n_vars
        nullify(this%fldsd(ib_var)%fld)
     end do
     deallocate(this%fldsd)
  end if
  if (allocated(this%fldss)) then
     do ib_var = 1, this%n_vars
        nullify(this%fldss(ib_var)%fld)
     end do
     deallocate(this%fldss)
  end if
  call this%atlas_FieldSet%final()
  !
end subroutine aq_field_delete

function aq_field_kind(this) result(kind)
  class(aq_fields), intent(in) :: this
  integer                      :: kind
  !
  kind = this%prec
  !
end function aq_field_kind

function aq_field_name(this) result(name)
  class(aq_fields), intent(in)  :: this
  character(len=:), allocatable :: name
  !
  name = trim(this%fs_name)
  !
end function aq_field_name

subroutine aq_field_rename(self, name)
  class(aq_fields), intent(inout) :: self
  character(len=*), intent(in)    :: name
  !
  self%fs_name = trim(name)
  !
end subroutine aq_field_rename

subroutine aq_field_zero(self)
  class(aq_fields), intent(inout) :: self
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec == aq_single) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        self%fldss(ib_var)%fld(:,:) = 0.0_aq_single
     end do
!$omp end parallel do
  else
!$omp parallel do
     do ib_var = 1, self%n_vars
        self%fldsd(ib_var)%fld(:,:) = 0.0_aq_real
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_zero

subroutine aq_field_ones(self)
  class(aq_fields), intent(inout) :: self
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec == aq_single) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        self%fldss(ib_var)%fld(:,:) = 1.0_aq_single
     end do
!$omp end parallel do
  else
!$omp parallel do
     do ib_var = 1, self%n_vars
        self%fldsd(ib_var)%fld(:,:) = 1.0_aq_real
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_ones

subroutine aq_field_dirac(self, config)
  class(aq_fields),          intent(inout) :: self
  type(fckit_Configuration), intent(in)    :: config
  !
  integer :: ndir
  integer, allocatable :: ixdir(:)
  integer, allocatable :: iydir(:)
  integer, allocatable :: ildir(:)
  real(aq_real), allocatable :: coorddir(:)
  character(len=:), allocatable :: ifdir(:)
  !
  integer :: ib_dir
  !
  if (config%has("londir")) then
     ndir = config%get_size("londir")
     call config%get_or_die("londir",coorddir)
     allocate(ixdir(ndir))
     do ib_dir = 1, ndir
        ixdir(ib_dir) = 1+nint((coorddir(ib_dir) - self%geom%xmin)/self%geom%deltax)
     enddo
     deallocate(coorddir)
  else if (config%has("ixdir")) then
     ndir = config%get_size("ixdir")
     call config%get_or_die("ixdir", ixdir)
  else
     call abor1_ftn("Missing lon or ix information for Dirac")
  end if
  !
  if (config%has("latdir")) then
     if (ndir /= config%get_size("latdir")) &
        & call abor1_ftn("Incoherent nr of coordinates for Dirac")
     call config%get_or_die("latdir",coorddir)
     allocate(iydir(ndir))
     do ib_dir = 1, ndir
        iydir(ib_dir) = 1+nint((coorddir(ib_dir) - self%geom%ymin)/self%geom%deltay)
     enddo
     deallocate(coorddir)
  else if (config%has("iydir")) then
     if (ndir /= config%get_size("iydir")) &
        & call abor1_ftn("Incoherent nr of coordinates for Dirac")
     call config%get_or_die("iydir", iydir)
  else
     call abor1_ftn("Missing lat or iy information for Dirac")
  end if
  !
  if (config%has("levdir")) then
     if (ndir /= config%get_size("levdir")) &
        & call abor1_ftn("Incoherent nr of coordinates for Dirac")
     call config%get_or_die("levdir", ildir)
     if ( self%geom%orientation == "up" ) then
        do ib_dir = 1,ndir
           ildir(ib_dir) = self%geom%mod_levels-ildir(ib_dir)+1
        end do
     else
        do ib_dir = 1,ndir
           ildir(ib_dir) = self%geom%mod_levels-self%geom%levels+ildir(ib_dir)
        end do
     end if
  else if (config%has("ildir")) then
     if (ndir /= config%get_size("ildir")) &
        & call abor1_ftn("Incoherent nr of coordinates for Dirac")
     call config%get_or_die("ildir", ildir)
  else
     call abor1_ftn("Missing level information for Dirac")
  end if
  !
  call config%get_or_die("var",   ifdir)
  !
  call self%zero()
  !
  do ib_dir = 1, ndir
     if (ixdir(ib_dir) > self%geom%grid%nx(1) .or. ixdir(ib_dir) < 1) &
        & call abor1_ftn('X index for dirac out of bounds')
     if (iydir(ib_dir) > self%geom%grid%ny() .or. iydir(ib_dir) < 1) &
        & call abor1_ftn('Y index for dirac out of bounds')
     if (ildir(ib_dir) > self%geom%levels .or. ildir(ib_dir) < 1) &
        & call abor1_ftn('Level for dirac out of bounds')
     if ( iydir(ib_dir) >= self%geom%fs%j_begin() .and. &
        & iydir(ib_dir) <= self%geom%fs%j_end()) then
        if ( ixdir(ib_dir) >= self%geom%fs%i_begin(iydir(ib_dir)) .and. &
           & ixdir(ib_dir) <= self%geom%fs%i_end(iydir(ib_dir))) then
           if (self%prec == aq_single) then
              self%fldss(self%idx_var(trim(ifdir(ib_dir))))%fld(ildir(ib_dir),&
                 &      self%geom%fs%index(ixdir(ib_dir),iydir(ib_dir))) = 1.0_aq_single
           else
              self%fldsd(self%idx_var(trim(ifdir(ib_dir))))%fld(ildir(ib_dir),&
                 &      self%geom%fs%index(ixdir(ib_dir),iydir(ib_dir))) = 1.0_aq_real
           end if
        end if
     end if
  end do
  !
  call self%halo_exchange()
  !
  deallocate(ixdir)
  deallocate(iydir)
  deallocate(ildir)
  deallocate(ifdir)
  !
end subroutine aq_field_dirac

subroutine aq_field_random(self)
  class(aq_fields), intent(inout) :: self
  !
  integer, parameter :: rseed = 7
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec == aq_single) then
     do ib_var = 1, self%n_vars
        call normal_distribution(self%fldss(ib_var)%fld(:,:), 0.0_aq_single, 1.0_aq_single, rseed)
     end do
  else
     do ib_var = 1, self%n_vars
        call normal_distribution(self%fldsd(ib_var)%fld(:,:), 0.0_aq_real, 1.0_aq_real, rseed)
     end do
  end if
  !
end subroutine aq_field_random

subroutine aq_field_copy(self, other)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: other
  !
  integer(atlas_kind_idx) :: ib_var
  integer(atlas_kind_idx) :: ib_pos
  logical :: ll_copy
  !
  ll_copy = self%n_vars == other%n_vars
  if (ll_copy) ll_copy = all(self%var_name == other%var_name)
  if (ll_copy) then
     call self%var_transf%copy(other%var_transf)
     select case(self%prec)
     case(aq_single)
        select case(other%prec)
        case(aq_single)
!$omp parallel do
           do ib_var = 1, self%n_vars
              call aq_copy(self%locsize, other%fldss(ib_var)%fld, self%fldss(ib_var)%fld)
           end do
!$omp end parallel do
        case default
           ! Demote
!$omp parallel do
           do ib_var = 1, self%n_vars
              self%fldss(ib_var)%fld(:,:) = real(other%fldsd(ib_var)%fld(:,:),kind=aq_single)
           end do
!$omp end parallel do
        end select
     case default
        select case(other%prec)
        case(aq_single)
           ! Promote
!$omp parallel do
           do ib_var = 1, self%n_vars
              self%fldsd(ib_var)%fld(:,:) = real(other%fldss(ib_var)%fld(:,:),kind=aq_real)
           end do
!$omp end parallel do
        case default
!$omp parallel do
           do ib_var = 1, self%n_vars
              call aq_copy(self%locsize, other%fldsd(ib_var)%fld, self%fldsd(ib_var)%fld)
           end do
!$omp end parallel do
        end select
     end select
  else
     if (self%n_vars <= other%n_vars) then
        ! Extraction of subset or variable reshuffling
        select case(self%prec)
        case(aq_single)
           select case(other%prec)
           case(aq_single)
!$omp parallel do private(ib_pos)
              do ib_var = 1, self%n_vars
                 ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
                 call aq_copy(self%locsize, other%fldss(ib_pos)%fld, self%fldss(ib_var)%fld)
              end do
!$omp end parallel do
           case default
              ! Demote
!$omp parallel do private(ib_pos)
              do ib_var = 1, self%n_vars
                 ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
                 self%fldss(ib_var)%fld(:,:) = real(other%fldsd(ib_pos)%fld(:,:),kind=aq_single)
              end do
!$omp end parallel do
           end select
        case default
           select case(other%prec)
           case(aq_single)
              ! Promote
!$omp parallel do private(ib_pos)
              do ib_var = 1, self%n_vars
                 ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
                 self%fldsd(ib_var)%fld(:,:) = real(other%fldss(ib_pos)%fld(:,:),kind=aq_real)
              end do
!$omp end parallel do
           case default
!$omp parallel do private(ib_pos)
              do ib_var = 1, self%n_vars
                 ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
                 call aq_copy(self%locsize, other%fldsd(ib_pos)%fld, self%fldsd(ib_var)%fld)
              end do
!$omp end parallel do
           end select
        end select
     else
        ! Injection into larger set
        select case(self%prec)
        case(aq_single)
           select case(other%prec)
           case(aq_single)
!$omp parallel do private(ib_pos)
              do ib_var = 1, other%n_vars
                 ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: injection source variable ' &
                    & //trim(other%var_name(ib_var))//' not in target variables')
                 call aq_copy(self%locsize, other%fldss(ib_var)%fld, self%fldss(ib_pos)%fld)
              end do
!$omp end parallel do
           case default
              ! Demote
!$omp parallel do private(ib_pos)
              do ib_var = 1, other%n_vars
                 ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: injection source variable ' &
                    & //trim(other%var_name(ib_var))//' not in target variables')
                 self%fldss(ib_pos)%fld(:,:) = real(other%fldsd(ib_var)%fld(:,:),kind=aq_single)
              end do
!$omp end parallel do
           end select
        case default
           select case(other%prec)
           case(aq_single)
              ! Promote
!$omp parallel do private(ib_pos)
              do ib_var = 1, other%n_vars
                 ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: injection source variable ' &
                    & //trim(other%var_name(ib_var))//' not in target variables')
                 self%fldsd(ib_pos)%fld(:,:) = real(other%fldss(ib_var)%fld(:,:),kind=aq_real)
              end do
!$omp end parallel do
           case default
!$omp parallel do private(ib_pos)
              do ib_var = 1, other%n_vars
                 ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
                 if (ib_pos < 1) &
                    & call abor1_ftn('aq_field_copy: injection source variable ' &
                    & //trim(other%var_name(ib_var))//' not in target variables')
                 call aq_copy(self%locsize, other%fldsd(ib_var)%fld, self%fldsd(ib_pos)%fld)
              end do
!$omp end parallel do
           end select
        end select
     end if
  end if
  !
end subroutine aq_field_copy

subroutine aq_field_self_add(self, other)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: other
  !
  integer(atlas_kind_idx) :: ib_var
  integer(atlas_kind_idx) :: ib_pos
  logical :: ll_direct
  !
  if (self%prec /= other%prec) &
     & call abor1_ftn('aq_field_self_add only allowed for equal precisions')
  ll_direct = self%n_vars == other%n_vars
  if (ll_direct) ll_direct = all(self%var_name == other%var_name)
  if (ll_direct) then
     if (self%prec == aq_single) then
!$omp parallel do
        do ib_var = 1, self%n_vars
           call aq_axpy(self%locsize, 1.0_oops_real, other%fldss(ib_var)%fld, self%fldss(ib_var)%fld)
        end do
!$omp end parallel do
     else
!$omp parallel do
        do ib_var = 1, self%n_vars
           call aq_axpy(self%locsize, 1.0_oops_real, other%fldsd(ib_var)%fld, self%fldsd(ib_var)%fld)
        end do
!$omp end parallel do
     end if
  else
     if (self%n_vars <= other%n_vars) then
        ! Extraction of subset or variable reshuffling
        if (self%prec == aq_single) then
!$omp parallel do private(ib_pos)
           do ib_var = 1, self%n_vars
              ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
              if (ib_pos < 1) &
                 & call abor1_ftn('aq_field_self_add: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
              call aq_axpy(self%locsize, 1.0_oops_real, other%fldss(ib_pos)%fld, self%fldss(ib_var)%fld)
           end do
!$omp end parallel do
        else
!$omp parallel do private(ib_pos)
           do ib_var = 1, self%n_vars
              ib_pos = findloc(other%var_name, self%var_name(ib_var), dim=1)
              if (ib_pos < 1) &
                 & call abor1_ftn('aq_field_self_add: extraction target variable ' &
                    & //trim(self%var_name(ib_var))//' not in source variables')
              call aq_axpy(self%locsize, 1.0_oops_real, other%fldsd(ib_pos)%fld, self%fldsd(ib_var)%fld)
           end do
!$omp end parallel do
        end if
     else
        ! Injection into larger set
        if (self%prec == aq_single) then
!$omp parallel do private(ib_pos)
           do ib_var = 1, other%n_vars
              ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
              if (ib_pos < 1) &
                 & call abor1_ftn('aq_field_self_add: injection source variable ' &
                    & //trim(self%var_name(ib_var))//' not in target variables')
              call aq_axpy(self%locsize, 1.0_oops_real, other%fldss(ib_var)%fld, self%fldss(ib_pos)%fld)
           end do
!$omp end parallel do
        else
!$omp parallel do private(ib_pos)
           do ib_var = 1, other%n_vars
              ib_pos = findloc(self%var_name, other%var_name(ib_var), dim=1)
              if (ib_pos < 1) &
                 & call abor1_ftn('aq_field_self_add: injection source variable ' &
                    & //trim(self%var_name(ib_var))//' not in target variables')
              call aq_axpy(self%locsize, 1.0_oops_real, other%fldsd(ib_var)%fld, self%fldsd(ib_pos)%fld)
           end do
!$omp end parallel do
        end if
     end if
  end if
  !
end subroutine aq_field_self_add

subroutine aq_field_self_sub(self, other)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: other
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec /= other%prec) &
     & call abor1_ftn('aq_field_self_sub only allowed for equal precisions')
  if (self%prec == aq_single) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_axpy(self%locsize, -1.0_oops_real, other%fldss(ib_var)%fld, self%fldss(ib_var)%fld)
     end do
!$omp end parallel do
  else
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_axpy(self%locsize, -1.0_oops_real, other%fldsd(ib_var)%fld, self%fldsd(ib_var)%fld)
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_self_sub

subroutine aq_field_self_mul(self, coeff)
  class(aq_fields), intent(inout) :: self
  real(oops_real),    intent(in)    :: coeff
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec == aq_single) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_scal(self%locsize, coeff, self%fldss(ib_var)%fld)
     end do
!$omp end parallel do
  else
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_scal(self%locsize, coeff, self%fldsd(ib_var)%fld)
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_self_mul

subroutine aq_field_add_incr(self, incr)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: incr
  !
  integer(atlas_kind_idx) :: ib_var, ib_pos
  !
  if (self%prec == aq_single .or. incr%prec == aq_single) &
     & call abor1_ftn('aq_field_add_incr not implemented for single precision')
!$omp parallel do private(ib_pos)
  do ib_var = 1, incr%n_vars
     ib_pos = findloc(self%var_name, incr%var_name(ib_var), dim=1)
     if (ib_pos < 1) &
        & call abor1_ftn('aq_field_add_incr: incr variable '//trim(incr%var_name(ib_var))//' not in state variables')
     call aq_axpy(self%locsize, 1.0_oops_real, incr%fldsd(ib_var)%fld, self%fldsd(ib_pos)%fld)
  end do
!$omp end parallel do
  !
end subroutine aq_field_add_incr

subroutine aq_field_diff_incr(self, fld1, fld2)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: fld1
  class(aq_fields), intent(in)    :: fld2
  !
  integer(atlas_kind_idx) :: ib_var, ib_pos
  !
  if (self%prec == aq_single .or. fld1%prec == aq_single .or. fld2%prec == aq_single) &
     & call abor1_ftn('aq_field_diff_incr not implemented for single precision')
  !
!$omp parallel do private(ib_pos)
  do ib_var = 1, self%n_vars
     ib_pos = findloc(fld1%var_name, self%var_name(ib_var), dim=1)
     if (ib_pos < 1) &
        & call abor1_ftn('aq_field_diff_incr: incr variable '//trim(self%var_name(ib_var))//' not in state variables')
     call aq_copy(self%locsize, fld1%fldsd(ib_pos)%fld, self%fldsd(ib_var)%fld)
     call aq_axpy(self%locsize, -1.0_oops_real, fld2%fldsd(ib_pos)%fld, self%fldsd(ib_var)%fld)
  end do
!$omp end parallel do
  !
end subroutine aq_field_diff_incr

subroutine aq_field_self_axpy(self, coeff, other)
  class(aq_fields), intent(inout) :: self
  real(oops_real),  intent(in)    :: coeff
  class(aq_fields), intent(in)    :: other
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec /= other%prec) &
     & call abor1_ftn('aq_field_self_axpy only allowed for equal precisions')
  if (self%prec == aq_single) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_axpy(self%locsize, coeff, other%fldss(ib_var)%fld, self%fldss(ib_var)%fld)
     end do
!$omp end parallel do
  else
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_axpy(self%locsize, coeff, other%fldsd(ib_var)%fld, self%fldsd(ib_var)%fld)
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_self_axpy

subroutine aq_field_self_schur(self, other)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(in)    :: other
  !
  integer(atlas_kind_idx) :: ib_var
  !
  if (self%prec == aq_single) then
     if (other%prec == aq_single) then
!$omp parallel do
        do ib_var = 1, self%n_vars
           self%fldss(ib_var)%fld(:,:) = self%fldss(ib_var)%fld(:,:) * other%fldss(ib_var)%fld(:,:)
        end do
!$omp end parallel do
     else
!$omp parallel do
        do ib_var = 1, self%n_vars
           self%fldss(ib_var)%fld(:,:) = self%fldss(ib_var)%fld(:,:) * other%fldsd(ib_var)%fld(:,:)
        end do
!$omp end parallel do
     end if
  else
     if (other%prec == aq_single) then
!$omp parallel do
        do ib_var = 1, self%n_vars
           self%fldsd(ib_var)%fld(:,:) = self%fldsd(ib_var)%fld(:,:) * other%fldss(ib_var)%fld(:,:)
        end do
!$omp end parallel do
     else
!$omp parallel do
        do ib_var = 1, self%n_vars
           self%fldsd(ib_var)%fld(:,:) = self%fldsd(ib_var)%fld(:,:) * other%fldsd(ib_var)%fld(:,:)
        end do
!$omp end parallel do
     end if
  end if
  !
end subroutine aq_field_self_schur

subroutine aq_field_dot_prod_with(self, other, zprod)
  class(aq_fields), intent(inout) :: self
  class(aq_fields), intent(inout) :: other
  real(oops_real),    intent(out) :: zprod
  !
  integer(atlas_kind_idx) :: ib_i, ib_j, ib_k, ib_var
  !
  if (self%prec /= other%prec) &
     & call abor1_ftn('aq_field_dot_prod_with only allowed for equal precisions')
  !
  call self%halo_exchange()
  call other%halo_exchange()
  !
  if (self%geom%halo > 0) then
     zprod = 0.0_oops_real
     if (self%prec == aq_single) then
!$omp parallel do reduction(+:zprod) private(ib_j, ib_i, ib_k)
        do ib_var = 1, self%n_vars
           do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
              do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
                 do ib_k = 1, self%geom%levels
                    zprod = zprod + &
                       &    self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j)) * &
                       &    other%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 end do
              end do
           end do
        end do
!$omp end parallel do
     else
!$omp parallel do reduction(+:zprod) private(ib_j, ib_i, ib_k)
        do ib_var = 1, self%n_vars
           do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
              do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
                 do ib_k = 1, self%geom%levels
                    zprod = zprod + &
                       &    self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j)) * &
                       &    other%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 end do
              end do
           end do
        end do
!$omp end parallel do
     end if
  else
     zprod = 0.0_oops_real
     if (self%prec == aq_single) then
!!$omp parallel do reduction(+:zprod)
        do ib_var = 1, self%n_vars
           zprod = zprod + aq_dot_product(self%locsize, self%fldss(ib_var)%fld, other%fldss(ib_var)%fld)
        end do
!!$omp end parallel do
     else
!!$omp parallel do reduction(+:zprod)
        do ib_var = 1, self%n_vars
           zprod = zprod + aq_dot_product(self%locsize, self%fldsd(ib_var)%fld, other%fldsd(ib_var)%fld)
        end do
!!$omp end parallel do
     end if
  end if
  call self%fmpi%allreduce(zprod,fckit_mpi_sum())
  !
end subroutine aq_field_dot_prod_with

subroutine aq_field_norm2(self, norm)
  class(aq_fields), intent(inout) :: self
  real(oops_real),  intent(out)   :: norm
  !
  call self%dot_prod_with(self, norm)
  norm = sqrt(norm)
  !
end subroutine aq_field_norm2

subroutine aq_field_stats_tot(self, valmin, valmax, mean, stddev, divnm1)
  class(aq_fields),  intent(inout) :: self
  real(aq_real),     intent(out)   :: valmin
  real(aq_real),     intent(out)   :: valmax
  real(aq_real),     intent(out)   :: mean
  real(aq_real),     intent(out)   :: stddev
  logical, optional, intent(in)    :: divnm1
  !
  integer(atlas_kind_idx) :: ib_i, ib_j, ib_k, ib_var
  real(aq_real) :: den
  logical :: sgl
  !
  sgl = self%prec == aq_single
  !
  call self%halo_exchange()
  !
  if (sgl) then
     valmin = minval(self%fldss(1)%fld)
     valmax = maxval(self%fldss(1)%fld)
     do ib_var = 2, self%n_vars
        valmin = min(valmin,minval(self%fldss(ib_var)%fld))
        valmax = max(valmax,maxval(self%fldss(ib_var)%fld))
     end do
  else
     valmin = minval(self%fldsd(1)%fld)
     valmax = maxval(self%fldsd(1)%fld)
     do ib_var = 2, self%n_vars
        valmin = min(valmin,minval(self%fldsd(ib_var)%fld))
        valmax = max(valmax,maxval(self%fldsd(ib_var)%fld))
     end do
  end if
  call self%fmpi%allreduce(valmin,fckit_mpi_min())
  call self%fmpi%allreduce(valmax,fckit_mpi_max())
  !
  mean = 0.0_aq_real
  stddev = 0.0_aq_real
  den = real(self%geom%grid%size()*self%geom%levels*self%n_vars,kind=aq_real)
  !
  if ( self%geom%halo > 0) then
     do ib_var = 1, self%n_vars
        do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
           do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
              do ib_k = 1, self%geom%levels
                 if (sgl) then
                    mean = mean + self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 else
                    mean = mean + self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_var = 1, self%n_vars
           mean = mean + sum(self%fldss(ib_var)%fld)
        end do
     else
        do ib_var = 1, self%n_vars
           mean = mean + sum(self%fldsd(ib_var)%fld)
        end do
     end if
  end if
  mean = mean / den
  call self%fmpi%allreduce(mean,fckit_mpi_sum())
  !
  if ( self%geom%halo > 0) then
     do ib_var = 1, self%n_vars
        do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
           do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
              do ib_k = 1, self%geom%levels
                 if (sgl) then
                    stddev = stddev + &
                       &     (self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean)**2
                 else
                    stddev = stddev + &
                       &     (self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean)**2
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_var = 1, self%n_vars
           stddev = stddev + sum((self%fldss(ib_var)%fld-mean)**2)
        end do
     else
        do ib_var = 1, self%n_vars
           stddev = stddev + sum((self%fldsd(ib_var)%fld-mean)**2)
        end do
     end if
  end if
  if (present(divnm1)) then
     if (divnm1) den = den - 1.0_aq_real
  end if
  stddev = sqrt(stddev / den)
  call self%fmpi%allreduce(stddev,fckit_mpi_sum())
  !
end subroutine aq_field_stats_tot

subroutine aq_field_stats_per_var(self, valmin, valmax, mean, stddev, divnm1)
  class(aq_fields),  intent(inout) :: self
  real(aq_real),     intent(out)   :: valmin(:)
  real(aq_real),     intent(out)   :: valmax(:)
  real(aq_real),     intent(out)   :: mean(:)
  real(aq_real),     intent(out)   :: stddev(:)
  logical, optional, intent(in)    :: divnm1
  !
  integer(atlas_kind_idx) :: ib_i, ib_j, ib_k, ib_var
  real(aq_real) :: den
  logical :: sgl
  !
  sgl = self%prec == aq_single
  !
  if (size(valmin) /= self%n_vars) call abor1_ftn('valmin per var wrong argument size')
  if (size(valmax) /= self%n_vars) call abor1_ftn('valmax per var wrong argument size')
  if (size(mean) /= self%n_vars) call abor1_ftn('mean per var wrong argument size')
  if (size(stddev) /= self%n_vars) call abor1_ftn('stddev per var wrong argument size')
  !
  call self%halo_exchange()
  !
  if (sgl) then
     do ib_var = 1, self%n_vars
        valmin(ib_var) = minval(self%fldss(ib_var)%fld)
        valmax(ib_var) = maxval(self%fldss(ib_var)%fld)
     end do
  else
     do ib_var = 1, self%n_vars
        valmin(ib_var) = minval(self%fldsd(ib_var)%fld)
        valmax(ib_var) = maxval(self%fldsd(ib_var)%fld)
     end do
  end if
  call self%fmpi%allreduce(valmin,fckit_mpi_min())
  call self%fmpi%allreduce(valmax,fckit_mpi_max())
  !
  mean(:) = 0.0_aq_real
  stddev(:) = 0.0_aq_real
  den = real(self%geom%grid%size()*self%geom%levels,kind=aq_real)
  !
  if ( self%geom%halo > 0) then
     do ib_var = 1, self%n_vars
        do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
           do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
              do ib_k = 1, self%geom%levels
                 if (sgl) then
                    mean(ib_var) = mean(ib_var) + &
                       &           self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 else
                    mean(ib_var) = mean(ib_var) + &
                       &           self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_var = 1, self%n_vars
           mean(ib_var) = sum(self%fldss(ib_var)%fld)
        end do
     else
        do ib_var = 1, self%n_vars
           mean(ib_var) = sum(self%fldsd(ib_var)%fld)
        end do
     end if
  end if
  mean(:) = mean(:) / den
  call self%fmpi%allreduce(mean,fckit_mpi_sum())
  !
  if ( self%geom%halo > 0) then
     do ib_var = 1, self%n_vars
        do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
           do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
              do ib_k = 1, self%geom%levels
                 if (sgl) then
                    stddev(ib_var) = stddev(ib_var) + &
                       &   (self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean(ib_var))**2
                 else
                    stddev(ib_var) = stddev(ib_var) + &
                       &   (self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean(ib_var))**2
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_var = 1, self%n_vars
           stddev(ib_var) = sum((self%fldss(ib_var)%fld-mean(ib_var))**2)
        end do
     else
        do ib_var = 1, self%n_vars
           stddev(ib_var) = sum((self%fldsd(ib_var)%fld-mean(ib_var))**2)
        end do
     end if
  end if
  if (present(divnm1)) then
     if (divnm1) den = den - 1.0
  end if
  stddev(:) = sqrt(stddev(:) / den)
  call self%fmpi%allreduce(stddev,fckit_mpi_sum())
  !
end subroutine aq_field_stats_per_var

subroutine aq_field_stats_per_var_lev(self, valmin, valmax, mean, stddev, divnm1)
  class(aq_fields),  intent(inout) :: self
  real(aq_real),     intent(out)   :: valmin(:,:)
  real(aq_real),     intent(out)   :: valmax(:,:)
  real(aq_real),     intent(out)   :: mean(:,:)
  real(aq_real),     intent(out)   :: stddev(:,:)
  logical, optional, intent(in)    :: divnm1
  !
  integer(atlas_kind_idx) :: ib_i, ib_j, ib_k, ib_var
  real(aq_real) :: den
  logical :: sgl
  !
  sgl = self%prec == aq_single
  !
  if (any(shape(valmin) /= [self%n_vars, self%geom%levels])) &
     & call abor1_ftn('valmin per var and lev wrong argument size')
  if (any(shape(valmax) /= [self%n_vars, self%geom%levels])) &
     & call abor1_ftn('valmax per var and lev wrong argument size')
  if (any(shape(mean) /= [self%n_vars, self%geom%levels])) &
     & call abor1_ftn('mean per var and lev wrong argument size')
  if (any(shape(stddev) /= [self%n_vars, self%geom%levels])) &
     & call abor1_ftn('stddev per var and lev wrong argument size')
  !
  call self%halo_exchange()
  !
  if (sgl) then
     do ib_k = 1, self%geom%levels
        do ib_var = 1, self%n_vars
           valmin(ib_var, :) = minval(self%fldss(ib_var)%fld,dim=2)
           valmax(ib_var, :) = maxval(self%fldss(ib_var)%fld,dim=2)
        end do
     end do
  else
     do ib_k = 1, self%geom%levels
        do ib_var = 1, self%n_vars
           valmin(ib_var, :) = minval(self%fldsd(ib_var)%fld,dim=2)
           valmax(ib_var, :) = maxval(self%fldsd(ib_var)%fld,dim=2)
        end do
     end do
  end if
  call self%fmpi%allreduce(valmin,fckit_mpi_min())
  call self%fmpi%allreduce(valmax,fckit_mpi_max())
  !
  mean(:,:) = 0.0_aq_real
  stddev(:,:) = 0.0_aq_real
  den = real(self%geom%grid%size(),kind=aq_real)
  !
  if ( self%geom%halo > 0) then
     do ib_k = 1, self%geom%levels
        do ib_var = 1, self%n_vars
           do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
              do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
                 if (sgl) then
                    mean(ib_var,ib_k) = mean(ib_var,ib_k) + &
                       &     self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 else
                    mean(ib_var,ib_k) = mean(ib_var,ib_k) + &
                       &     self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_k = 1, self%geom%levels
           do ib_var = 1, self%n_vars
              mean(ib_var,ib_k) = sum(self%fldss(ib_var)%fld(ib_k,:))
           end do
        end do
     else
        do ib_k = 1, self%geom%levels
           do ib_var = 1, self%n_vars
              mean(ib_var,ib_k) = sum(self%fldsd(ib_var)%fld(ib_k,:))
           end do
        end do
     end if
  end if
  mean(:,:) = mean(:,:) / den
  call self%fmpi%allreduce(mean,fckit_mpi_sum())
  !
  if ( self%geom%halo > 0) then
     do ib_k = 1, self%geom%levels
        do ib_var = 1, self%n_vars
           do ib_j = self%geom%fs%j_begin(), self%geom%fs%j_end()
              do ib_i = self%geom%fs%i_begin(ib_j), self%geom%fs%i_end(ib_j)
                 if (sgl) then
                    stddev(ib_var,ib_k) = stddev(ib_var,ib_k) + &
                       &   (self%fldss(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean(ib_var,ib_k))**2
                 else
                    stddev(ib_var,ib_k) = stddev(ib_var,ib_k) + &
                       &   (self%fldsd(ib_var)%fld(ib_k,self%geom%fs%index(ib_i,ib_j))-mean(ib_var,ib_k))**2
                 end if
              end do
           end do
        end do
     end do
  else
     if (sgl) then
        do ib_k = 1, self%geom%levels
           do ib_var = 1, self%n_vars
              stddev(ib_var,ib_k) = sum((self%fldss(ib_var)%fld(ib_k,:)-mean(ib_var,ib_k))**2)
           end do
        end do
     else
        do ib_k = 1, self%geom%levels
           do ib_var = 1, self%n_vars
              stddev(ib_var,ib_k) = sum((self%fldsd(ib_var)%fld(ib_k,:)-mean(ib_var,ib_k))**2)
           end do
        end do
     end if
  end if
  if (present(divnm1)) then
     if (divnm1) den = den - 1.0_aq_real
  end if
  stddev(:,:) = sqrt(stddev(:,:) / den)
  call self%fmpi%allreduce(stddev,fckit_mpi_sum())
  !
end subroutine aq_field_stats_per_var_lev

subroutine aq_field_info(self, config)
  class(aq_fields),          intent(inout) :: self
  type(fckit_Configuration), intent(inout) :: config

  type(fckit_Configuration) :: geom_config
  type(fckit_Configuration) :: stat_config
  type(fckit_Configuration) :: species_config

  integer :: year, month, day, hour, minute, second
  integer :: ib_var, nx, ny, nz, halo, mod_levels
  real(aq_real) :: deltax, deltay, xmin, ymin
  character(len=:), allocatable :: domname
  character(len=:), allocatable :: orientation
  character(len=:), allocatable :: model
  real(aq_real), allocatable, dimension(:) :: minvar, maxvar, meanvar, stddevvar

  call config%set('name',self%fs_name)
  call config%set('state variables',self%var_name)
  call config%set('n_vars',self%n_vars)
  call datetime_to_yyyymmddhhmmss(self%date, year, month, day, hour, minute, second)
  call config%set('date',[year,month,day,hour,minute,second])
  call config%set('precision',self%prec)

  geom_config = fckit_configuration()
  call self%geom%info(nx, ny, nz, deltax, deltay, xmin, ymin, mod_levels, halo, domname, orientation, model)
  call geom_config%set('nx',nx)
  call geom_config%set('ny',ny)
  call geom_config%set('nz',nz)
  call geom_config%set('deltax',deltax)
  call geom_config%set('deltay',deltay)
  call geom_config%set('xmin',xmin)
  call geom_config%set('ymin',ymin)
  call geom_config%set('modlev', mod_levels)
  call geom_config%set('halo',halo)
  call geom_config%set('domname',domname)
  call geom_config%set('orientation',orientation)
  call geom_config%set('model',model)
  call config%set('geometry',geom_config)
  call geom_config%final()

  stat_config = fckit_configuration()
  allocate(minvar(self%n_vars))
  allocate(maxvar(self%n_vars))
  allocate(meanvar(self%n_vars))
  allocate(stddevvar(self%n_vars))
  call self%stats(minvar, maxvar, meanvar, stddevvar)

  do ib_var = 1, self%n_vars
     species_config = fckit_configuration()
     call species_config%set('min',minvar(ib_var))
     call species_config%set('max',maxvar(ib_var))
     call species_config%set('mean',meanvar(ib_var))
     call species_config%set('stddev',stddevvar(ib_var))
     call stat_config%set(trim(self%var_name(ib_var)),species_config)
     call species_config%final()
  end do

  deallocate(minvar)
  deallocate(maxvar)
  deallocate(meanvar)
  deallocate(stddevvar)
  call config%set('statistics',stat_config)
  call stat_config%final()

end subroutine aq_field_info

subroutine aq_field_read(self, config, date)
  class(aq_fields),          intent(inout) :: self
  type(fckit_Configuration), intent(in)    :: config
  type(datetime), optional,  intent(in)    :: date

  character(len=:), allocatable :: file
  type(datetime) :: vdate
  integer :: strlen, dotpos
  type(fckit_Configuration) :: transform_config

  vdate = self%date
  if (present(date)) vdate = date

  call config%get_or_die("filename", file)
  strlen = len_trim(file)
  dotpos = index(trim(file),'.',back=.true.)

  if (file(dotpos+1:strlen) == 'nc') then
     if (trim(self%geom%model) == "MOCAGE") then
        call aq_read_mocage_nc(self, self%var_name, self%geom, file, vdate)
     else
        call abor1_ftn('NetCDF input only coded for MOCAGE')
     end if
  else
     call abor1_ftn('Input format '//file(dotpos+1:strlen)//' not recognised')
  end if

  if (config%has("transfvar")) then
     call config%get_or_die("transfvar", transform_config)
     call self%var_transf%setup(self%var_name, transform_config)
     call self%var_transf%apply(self, self%var_name)
  end if

end subroutine aq_field_read

subroutine aq_field_write(self, config, date)
  class(aq_fields),          intent(inout) :: self
  type(fckit_Configuration), intent(in)    :: config
  type(datetime), optional,  intent(in)    :: date

  type(aq_fields) :: output

  if (self%var_transf%nb_vars > 0) then
     call output%create_from(self)
     call output%copy(self)
     call output%var_transf%apply_inverse(output, output%var_name)
     call aq_field_do_write(output, config, date)
     call output%delete()
  else
     call aq_field_do_write(self, config, date)
  end if

end subroutine aq_field_write

subroutine aq_field_do_write(self, config, date)
  class(aq_fields),          intent(inout) :: self
  type(fckit_Configuration), intent(in)    :: config
  type(datetime), optional,  intent(in)    :: date

  character(len=:), allocatable :: fileprefix
  character(len=:), allocatable :: datadir
  character(len=:), allocatable :: filetype
  character(len=10) :: filedate
  character(len=:), allocatable :: file
  type(datetime) :: vdate
  integer :: year, month, day, hour, minute, second
  integer :: strlen, dotpos

  vdate = self%date
  if (present(date)) vdate = date

  if (config%has("filename")) then
     call config%get_or_die("filename", file)
  else
     call datetime_to_yyyymmddhhmmss(vdate, year, month, day, hour, &
        & minute, second)
     write(filedate,'(I4.4,3I2.2)') year, month, day, hour
     call config%get_or_die("fileprefix", fileprefix)
     if (config%has("filetype")) then
        call config%get_or_die("filetype", filetype)
     else
        filetype="nc"
     end if
     if (config%has("datadir")) then
        call config%get_or_die("datadir", datadir)
        file=trim(datadir)//'/'//trim(fileprefix)//'+'//trim(filedate)//&
           & '.'//trim(filetype)
     else
        file=trim(fileprefix)//'+'//trim(filedate)//'.'//trim(filetype)
     end if
  end if
  strlen = len(trim(file))
  dotpos = index(trim(file),'.',back=.true.)

  if (file(dotpos+1:strlen) == 'nc') then
     if (trim(self%geom%model) == "MOCAGE") then
!AP MISSING atlas_FieldSet %name() interfaces
!AP        call aq_write_mocage_nc(self, self%var_name, self%geom, file, vdate)
        call aq_write_mocage_nc(self, self%fs_name, self%var_name, self%geom, file, vdate)
!AP END
     else
        call abor1_ftn('NetCDF output only coded for MOCAGE')
     end if
  else if (file(dotpos+1:strlen) == 'msh') then
     call aq_write_field_gmsh(self, self%geom, file)
  else
     call abor1_ftn('Output format '//file(dotpos+1:strlen)//' not recognised')
  end if

end subroutine aq_field_do_write

subroutine aq_field_ana_IC(self, config)
  class(aq_fields), intent(inout)       :: self
  type(fckit_Configuration), intent(in) :: config

  type(atlas_Field) :: field_xy
  real(atlas_kind_real64), pointer :: xy(:,:)
  integer(atlas_kind_idx) :: ib
  type(fckit_Configuration) :: transform_config

  real(aq_real), parameter :: dp_pi=3.14159265359
  real(aq_real), parameter :: dLon0 = 6.3
  real(aq_real), parameter :: dLat0 = 0.8
  real(aq_real), parameter :: dR0   = 3.0
  real(aq_real), parameter :: dD    = 10.0
  real(aq_real), parameter :: dT    = 6.0
  real(aq_real) :: dp_conv

  real(aq_real) :: dSinC, dCosC, dCosT, dSinT
  real(aq_real) :: dTrm, dX, dY, dZ
  real(aq_real) :: dlon, dlat
  real(aq_real) :: dRho, dVt, dOmega

  !
  ! Fill level 1 of variable 1

  dp_conv = dp_pi/180._aq_real
  dSinC = sin( dLat0 )
  dCosC = cos( dLat0 )

  field_xy = self%geom%fs%xy()
  call field_xy%data(xy)

  do ib = 1, self%geom%fs%size()
     ! Find the rotated longitude and latitude of a point on a sphere
     !    with pole at (dLon0, dLat0).
     dCosT = cos( xy(2,ib)*dp_conv )
     dSinT = sin( xy(2,ib)*dp_conv )

     dTrm = dCosT * cos( xy(1,ib)*dp_conv - dLon0 )
     dX   = dSinC * dTrm - dCosC * dSinT
     dY   = dCosT * sin( xy(1,ib)*dp_conv - dLon0 )
     dZ   = dSinC * dSinT + dCosC * dTrm

     dlon = atan2( dY, dX )
     if( dlon < 0.0_aq_real ) dlon = dlon + 2.0_aq_real * dp_pi
     dlat = asin( dZ )

     dRho = dR0 * cos(dlat)
     dVt = 3.0_aq_real * sqrt(3.0_aq_real)/2.0_aq_real/cosh(dRho)/cosh(dRho)*tanh(dRHo)
     if (dRho == 0.0_aq_real) then
        dOmega = 0.0_aq_real
     else
        dOmega = dVt / dRho
     end if

     if (self%prec == aq_single) then
        self%fldss(1)%fld(1,ib) = real(2.0_aq_real * &
           & ( 1.0_aq_real + tanh( dRho / dD * sin( dLon - dOmega * dT ) ) ), kind=aq_single)
     else
        self%fldsd(1)%fld(1,ib) = 2.0_aq_real * &
           & ( 1.0_aq_real + tanh( dRho / dD * sin( dLon - dOmega * dT ) ) )
     end if
  end do

  call field_xy%final()

  ! Propagate on the vertical
  if (self%prec == aq_single) then
     do ib = 2, self%geom%levels
        self%fldss(1)%fld(ib,:) = 1.0_aq_single/real(ib,kind=aq_single) * self%fldss(1)%fld(1,:)
     end do
  else
     do ib = 2, self%geom%levels
        self%fldsd(1)%fld(ib,:) = 1.0_aq_real/real(ib,kind=aq_real) * self%fldsd(1)%fld(1,:)
     end do
  end if

  ! Propagate to other species
  if (self%prec == aq_single) then
     do ib = 2, self%n_vars
        self%fldss(ib)%fld(:,:) = self%fldss(ib-1)%fld(:,:)*1.2_aq_single
     end do
  else
     do ib = 2, self%n_vars
        self%fldsd(ib)%fld(:,:) = self%fldsd(ib-1)%fld(:,:)*1.2_aq_real
     end do
  end if

  call self%halo_exchange()
  !
  if (config%has("transfvar")) then
     call config%get_or_die("transfvar", transform_config)
     call self%var_transf%setup(self%var_name, transform_config)
     call self%var_transf%apply(self, self%var_name)
  end if
  !
end subroutine aq_field_ana_IC

integer(atlas_kind_idx) function aq_field_idx_var(self, var) result(ib_var)

  class(aq_fields), intent(in)  :: self
  character(len=*), intent(in)  :: var

  if (.not. self%has_field(trim(var))) &
     & call abor1_ftn('Idx var: variable '//trim(var)// ' not found')

  do ib_var = 1, self%n_vars
     if ( trim(self%var_name(ib_var)) == trim(var)) exit
  end do

end function aq_field_idx_var

subroutine aq_field_var(self, var, fld_var)
  class(aq_fields), intent(in)  :: self
  character(len=*), intent(in)  :: var
  real(aq_real),    intent(out) :: fld_var(:,:)

  integer :: il_var
  !
  il_var = self%idx_var(var)

  if ( any(shape(fld_var) /= [self%geom%levels,self%geom%fs%size()])) &
     & call abor1_ftn('Field get var: shape of output field not compliant')

  if (self%prec == aq_single) then
     ! Promote
     fld_var(:,:) = real(self%fldss(il_var)%fld(:,:),kind=aq_real)
  else
     fld_var(:,:) = self%fldsd(il_var)%fld(:,:)
  end if
  !
end subroutine aq_field_var

subroutine aq_field_gather_var_at_lev(self, var, lev, fld_2d, owner)
  class(aq_fields),                  intent(in)  :: self
  character(len=*),                  intent(in)  :: var
  integer(atlas_kind_idx),           intent(in)  :: lev
  integer(atlas_kind_idx), optional, intent(in)  :: owner
  real(aq_real),                     intent(out) :: fld_2d(:,:)

  integer(atlas_kind_idx) :: il_owner

  type(atlas_Field)  :: aloc_2d
  type(atlas_Field)  :: aglo_2d
  real(aq_single), pointer  :: loc_2ds(:)
  real(aq_single), pointer  :: glo_2ds(:,:)
  real(aq_real), pointer  :: loc_2d(:)
  real(aq_real), pointer  :: glo_2d(:,:)
  integer(atlas_kind_idx) :: il_var

  !
  il_owner = 0
  if (present(owner)) il_owner = owner

  if ( self%fmpi%rank() == il_owner ) then
     if ( any(shape(fld_2d) /= [self%geom%grid%nx(1),self%geom%grid%ny()])) &
        & call abor1_ftn('Gather var at lev: shape of global output field not compliant')
  end if

  if (lev > self%geom%levels .or. lev < 1) &
     & call abor1_ftn('Gather var at lev: level out of bounds')

  if (.not.self%has_field(trim(var))) &
     & call abor1_ftn('Var: variable '//trim(var)// ' not found')

  il_var = self%idx_var(trim(var))

  aloc_2d = self%geom%fs%create_field(name=trim(var),  &
     &                                   kind=atlas_real(self%prec), &
     &                                   levels=0)
  !
  if (self%prec == aq_single) then
     call aloc_2d%data(loc_2ds)
     loc_2ds(:) = self%fldss(il_var)%fld(lev,:)
  else
     call aloc_2d%data(loc_2d)
     loc_2d(:) = self%fldsd(il_var)%fld(lev,:)
  end if

  aglo_2d = self%geom%fs%create_field(name=trim(var),  &
     &                                   kind=atlas_real(self%prec), &
     &                                   global = .true., &
     &                                   owner = il_owner, &
     &                                   levels=0)
  call self%geom%fs%gather(aloc_2d, aglo_2d)
  if ( self%fmpi%rank() == il_owner ) then
     if (self%prec == aq_single) then
        call aglo_2d%data(glo_2ds, shape=[self%geom%grid%nx(1),self%geom%grid%ny()])
        ! Promote
        fld_2d(:,:) = real(glo_2ds(:,:),kind = aq_real)
     else
        call aglo_2d%data(glo_2d, shape=[self%geom%grid%nx(1),self%geom%grid%ny()])
        fld_2d(:,:) = glo_2d(:,:)
     end if
  end if

  call aglo_2d%final()
  call aloc_2d%final()
  !
end subroutine aq_field_gather_var_at_lev

subroutine aq_field_scatteradd_var_at_lev(self, var, lev, fld_2d, owner)
  class(aq_fields),                  intent(inout) :: self
  character(len=*),                  intent(in)    :: var
  integer(atlas_kind_idx),           intent(in)    :: lev
  integer(atlas_kind_idx), optional, intent(in)    :: owner
  real(aq_real),                     intent(in)    :: fld_2d(:,:)

  integer(atlas_kind_idx) :: il_owner

  type(atlas_Field)  :: aloc_2d
  type(atlas_Field)  :: aglo_2d
  real(aq_single), pointer  :: loc_2ds(:)
  real(aq_single), pointer  :: glo_2ds(:,:)
  real(aq_real), pointer  :: loc_2d(:)
  real(aq_real), pointer  :: glo_2d(:,:)
  integer(atlas_kind_idx) :: il_var

  !
  il_owner = 0
  if (present(owner)) il_owner = owner

  if ( self%fmpi%rank() == il_owner ) then
     if ( any(shape(fld_2d) /= [self%geom%grid%nx(1),self%geom%grid%ny()])) &
        & call abor1_ftn('Gather var at lev: shape of global output field not compliant')
  end if

  if (lev > self%geom%levels .or. lev < 1) &
     & call abor1_ftn('Gather var at lev: level out of bounds')

  il_var = self%idx_var(trim(var))

  aglo_2d = self%geom%fs%create_field(name=trim(var),  &
     &                                   kind=atlas_real(self%prec), &
     &                                   global = .true., &
     &                                   owner = il_owner, &
     &                                   levels=0)
  if ( self%fmpi%rank() == il_owner ) then
     if (self%prec == aq_single) then
        call aglo_2d%data(glo_2ds, shape=[self%geom%grid%nx(1),self%geom%grid%ny()])
        ! Demote
        glo_2ds(:,:) = real(fld_2d(:,:), kind=aq_single)
     else
        call aglo_2d%data(glo_2d, shape=[self%geom%grid%nx(1),self%geom%grid%ny()])
        glo_2d(:,:) = fld_2d(:,:)
     end if
  end if

  aloc_2d = self%geom%fs%create_field(name=trim(var),  &
     &                                   kind=atlas_real(self%prec), &
     &                                   levels=0)
  !
  call self%geom%fs%scatter(aglo_2d, aloc_2d)
  !
  if (self%prec == aq_single) then
     call aloc_2d%data(loc_2ds)
     self%fldss(il_var)%fld(lev,:) = self%fldss(il_var)%fld(lev,:) + loc_2ds(:)
  else
     call aloc_2d%data(loc_2d)
     self%fldsd(il_var)%fld(lev,:) = self%fldsd(il_var)%fld(lev,:) + loc_2d(:)
  end if
  !
  call aglo_2d%final()
  call aloc_2d%final()
  !
  call self%halo_exchange()
  !
end subroutine aq_field_scatteradd_var_at_lev

integer(oops_int) function aq_field_serial_size(self) result(ser_size)

  class(aq_fields), intent(in)  :: self

  ser_size = self%locsize * self%n_vars

end function aq_field_serial_size

subroutine aq_field_serialize(self, buff)
  class(aq_fields), intent(in)  :: self
  real(oops_real),  intent(out) :: buff(:)

  integer :: ib_var, ib_lev, ib_node, icounter
  !
  if (self%prec == oops_real) then
!$omp parallel do
     do ib_var = 1, self%n_vars
        call aq_copy(self%locsize, self%fldsd(ib_var)%fld, &
           & buff((ib_var-1)*self%locsize+1:ib_var*self%locsize))
     end do
!$omp end parallel do
  else
!$omp parallel do private(ib_lev, ib_node, icounter)
     do ib_var = 1, self%n_vars
        icounter = (ib_var-1)*self%locsize+1
        do ib_node = 1, self%geom%fs%size()
           do ib_lev = 1, self%geom%levels
              buff(icounter) = real(self%fldss(ib_var)%fld(ib_lev, ib_node), kind=oops_real)
              icounter = icounter + 1
           end do
        end do
     end do
!$omp end parallel do
  end if
  !
end subroutine aq_field_serialize

subroutine aq_field_deserialize(self, buff, offset)
   class(aq_fields), intent(inout) :: self
   real(oops_real),  intent(in)    :: buff(:)
   integer,          intent(inout) :: offset

   integer :: ib_var, ib_lev, ib_node, icounter

   !
   if (self%prec == oops_real) then
!$omp parallel do
      do ib_var = 1, self%n_vars
         call aq_copy(self%locsize, buff(offset+(ib_var-1)*self%locsize+1:offset+ib_var*self%locsize), &
            & self%fldsd(ib_var)%fld)
      end do
!$omp end parallel do
   else
!$omp parallel do private(ib_lev, ib_node, icounter)
      do ib_var = 1, self%n_vars
         icounter = offset + (ib_var-1)*self%locsize+1
         do ib_node = 1, self%geom%fs%size()
            do ib_lev = 1, self%geom%levels
               self%fldss(ib_var)%fld(ib_lev, ib_node) = real(buff(icounter), kind=aq_single)
               icounter = icounter + 1
            end do
         end do
      end do
!$omp end parallel do
   end if
   !
   offset = offset + self%n_vars * self%locsize
   !
end subroutine aq_field_deserialize

subroutine aq_field_serialize_single(self, buff)
  class(aq_fields), intent(in)  :: self
  real(aq_single),  intent(out) :: buff(:)

  integer :: ib_var
  !
!$omp parallel do
  do ib_var = 1, self%n_vars
     call aq_copy(self%locsize, self%fldss(ib_var)%fld, &
        & buff((ib_var-1)*self%locsize+1:ib_var*self%locsize))
  end do
!$omp end parallel do
  !
end subroutine aq_field_serialize_single

subroutine aq_field_serialize_real(self, buff)
  class(aq_fields), intent(in)  :: self
  real(aq_real),    intent(out) :: buff(:)

  integer :: ib_var
  !
!$omp parallel do
  do ib_var = 1, self%n_vars
     call aq_copy(self%locsize, self%fldsd(ib_var)%fld, &
        & buff((ib_var-1)*self%locsize+1:ib_var*self%locsize))
  end do
!$omp end parallel do
  !
end subroutine aq_field_serialize_real

subroutine aq_field_deserialize_single(self, buff, offset)
   class(aq_fields), intent(inout) :: self
   real(aq_single),  intent(in)    :: buff(:)
   integer,          intent(inout) :: offset

   integer :: ib_var
   !
!$omp parallel do
   do ib_var = 1, self%n_vars
      call aq_copy(self%locsize, buff(offset+(ib_var-1)*self%locsize+1:offset+ib_var*self%locsize), &
         & self%fldss(ib_var)%fld)
   end do
!$omp end parallel do
   !
   offset = offset + self%n_vars * self%locsize
   !
end subroutine aq_field_deserialize_single
  !
subroutine aq_field_deserialize_real(self, buff, offset)
   class(aq_fields), intent(inout) :: self
   real(aq_real),    intent(in)    :: buff(:)
   integer,          intent(inout) :: offset

   integer :: ib_var
   !
!$omp parallel do
   do ib_var = 1, self%n_vars
      call aq_copy(self%locsize, buff(offset+(ib_var-1)*self%locsize+1:offset+ib_var*self%locsize), &
         & self%fldsd(ib_var)%fld)
   end do
!$omp end parallel do
   !
   offset = offset + self%n_vars * self%locsize
   !
end subroutine aq_field_deserialize_real

subroutine aq_field_set_atlas(self, vars, fieldset)
   class(aq_fields),     intent(in)    :: self
   type(oops_variables), intent(in)    :: vars
   type(atlas_fieldset), intent(inout) :: fieldset
   !
   integer(atlas_kind_idx) :: ib_var
   character(len=aq_varlen) :: fieldname
   type(atlas_Field) :: afld
   !
   do ib_var = 1, vars%nvars()
      fieldname = vars%variable(ib_var)
      if (self%has_field(trim(fieldname))) then
         afld = self%field(trim(fieldname))
         call fieldset%add(afld)
         call afld%final()
      else
         call abor1_ftn('Variable '//trim(fieldname)//' not in field')
      end if
   end do
   !
end subroutine aq_field_set_atlas

subroutine aq_field_to_atlas(self, vars, fieldset)
   class(aq_fields),     intent(inout) :: self
   type(oops_variables), intent(in)    :: vars
   type(atlas_fieldset), intent(inout) :: fieldset
   !
   integer(atlas_kind_idx) :: ib_var
   character(len=aq_varlen) :: fieldname
   type(atlas_Field) :: afld_s, afld_t
   real(kind=aq_single), pointer :: flds(:,:)
   real(kind=aq_real), pointer :: fldd(:,:)
   !
   do ib_var = 1, vars%nvars()
      fieldname = vars%variable(ib_var)
      if (self%has_field(trim(fieldname))) then
         if (fieldset%has_field(trim(fieldname))) then
            afld_s = self%field(trim(fieldname))
            afld_t = fieldset%field(trim(fieldname))
            if (afld_t /= afld_s) then
               if (self%prec == aq_single) then
                  call afld_t%data(flds)
                  flds(:,:) = self%fldss(self%idx_var(trim(fieldname)))%fld(:,:)
               else
                  call afld_t%data(fldd)
                  fldd(:,:) = self%fldsd(self%idx_var(trim(fieldname)))%fld(:,:)
               end if
!AQ This output should only be on debug channel
!AQ            else
!AQ               if (self%fmpi%rank() == 0) &
!AQ                  & print '(3A)', 'Skipping copy of ',trim(fieldname),&
!AQ                  &' because src and tgt share the same pointer'
            end if
            call afld_s%final()
            call afld_t%final()
         else
            call abor1_ftn('Variable '//trim(fieldname)//' not in destination fieldset')
         end if
      else
         call abor1_ftn('Variable '//trim(fieldname)//' not in source field')
      end if
   end do
   !
end subroutine aq_field_to_atlas

! ------------------------------------------------------------------------------
end module aq_fields_mod
