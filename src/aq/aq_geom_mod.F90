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

module aq_geom_mod
   use aq_constants_mod
   use fckit_module
   use atlas_module

   implicit none

   private
   public :: aq_geom

   type aq_geom
      type(atlas_StructuredGrid)                  :: grid
      type(fckit_mpi_comm)                        :: fmpi
      character(len=:), allocatable               :: domname
      character(len=:), allocatable               :: model
      character(len=:), allocatable               :: orientation
      type(atlas_functionspace_StructuredColumns) :: fs
      integer(atlas_kind_idx)                     :: levels
      integer(atlas_kind_idx)                     :: mod_levels
      integer(atlas_kind_idx)                     :: halo
      integer(atlas_kind_idx)                     :: bbox_imin, bbox_imax, bbox_isz
      integer(atlas_kind_idx)                     :: bbox_jmin, bbox_jmax, bbox_jsz
      real(aq_real), allocatable                  :: halo_mask(:)
      real(aq_real)                               :: deltax, deltay
      real(aq_real)                               :: xmin, ymin
      integer(atlas_kind_idx)                     :: nx, ny, nz
      real(aq_real), allocatable                  :: lons(:)
      real(aq_real), allocatable                  :: lats(:)
   contains
      procedure, public :: create => aq_geom_create
      procedure, public :: fill_atlas_fieldset  => aq_geom_fill_atlas_fieldset
      procedure, public :: clone  => aq_geom_clone
      procedure, public :: info   => aq_geom_info
      procedure, public :: delete => aq_geom_delete
      procedure, public :: final  => aq_geom_delete
   end type aq_geom

contains

   subroutine aq_geom_create(self, config, fckit_mpi)
      class(aq_geom), intent(inout)    :: self
      type(fckit_Configuration)        :: config
      type(fckit_mpi_comm), intent(in) :: fckit_mpi
      !
      character(len=:), allocatable :: str
      !
      integer, allocatable :: ila_work(:)
      integer(atlas_kind_idx) :: ib_i, ib_j
      !
      self%fmpi = fckit_mpi
      call config%get_or_die("nx",self%nx)
      call config%get_or_die("ny",self%ny)
      call config%get_or_die("dx",self%deltax)
      call config%get_or_die("dy",self%deltay)
      call config%get_or_die("xmin",self%xmin)
      call config%get_or_die("ymin",self%ymin)
      call config%get_or_die("domname",self%domname)
      !
      self%halo = 0
      if (config%has("halo")) call config%get_or_die("halo",self%halo)
      !
      self%model = "MOCAGE"
      if (config%has("model")) call config%get_or_die("model",self%model)
      !
      self%orientation = "up"
      if (config%has("orientation")) then
         call config%get_or_die("orientation",str)
         select case(str)
            case("up")
               self%orientation = str
            case("UP")
               self%orientation = "up"
            case("down")
               self%orientation = str
            case("DOWN")
               self%orientation = "down"
            case default
               call abor1_ftn("Geom orientation can only be 'up' or 'down'")
            end select
      end if
      !
      call config%get_or_die("levels",self%levels)
      !
      self%mod_levels = self%levels
      if (config%has("model levels")) &
         & call config%get_or_die("model levels",self%mod_levels)
      if (self%mod_levels < 0) self%mod_levels = self%levels
      !
      self%nz = self%levels
      !
      ! Bounding box global indexes. Notice that C++ atlas has already a similar function
      self%bbox_jmin = self%fs%j_begin()
      self%bbox_jmax = self%fs%j_end()
      allocate(ila_work(self%bbox_jmin:self%bbox_jmax))
      do ib_j = self%bbox_jmin, self%bbox_jmax
         ila_work(ib_j) = self%fs%i_begin(ib_j)
      end do
      self%bbox_imin = minval(ila_work)
      do ib_j = self%bbox_jmin, self%bbox_jmax
         ila_work(ib_j) = self%fs%i_end(ib_j)
      end do
      self%bbox_imax = maxval(ila_work)
      deallocate(ila_work)
      !
      allocate(self%lons(self%nx))
      allocate(self%lats(self%ny))
      !
      do ib_i = 1, self%nx
         self%lons(ib_i) = self%grid%x(ib_i,1)
      end do
      do ib_j = 1, self%ny
         self%lats(ib_j) = self%grid%y(ib_j)
      end do
      !
      if ( self%halo>0 ) then
         allocate( self%halo_mask(self%fs%size()) )
         self%halo_mask(:) = 0.0
         do ib_j = self%fs%j_begin(), self%fs%j_end()
            do ib_i = self%fs%i_begin(ib_j), self%fs%i_end(ib_j)
               self%halo_mask(self%fs%index(ib_i,ib_j)) = 1.0
            end do
         end do
      end if
      !
   end subroutine aq_geom_create

   subroutine aq_geom_fill_atlas_fieldset(self,afieldset)
      class(aq_geom),intent(inout)       :: self
      type(atlas_fieldset),intent(inout) :: afieldset
      !
      integer :: ix,iy,iz,inode
      real(aq_real) :: lonlat(2), dx, dy
      real(aq_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
      type(atlas_field) :: afield
      !
      afield = self%fs%create_field(name='area',kind=atlas_real(aq_real),levels=0)
      call afield%data(real_ptr_1)
      dy = self%deltay*deg_to_rad*req
      inode = 0
      do iy=self%fs%j_begin(),self%fs%j_end()
         do ix=self%fs%i_begin(iy),self%fs%i_end(iy)
            inode = inode+1
            lonlat = self%grid%lonlat(ix,iy)
            dx = self%deltax*deg_to_rad*req*cos(lonlat(2)*deg_to_rad)
            real_ptr_1(inode) = dx*dy
         end do
      end do
      call afieldset%add(afield)
      call afield%final()
      !
      afield = self%fs%create_field(name='vunit',kind=atlas_real(aq_real))
      call afield%data(real_ptr_2)
      do iz=1,self%nz
        real_ptr_2(iz,:) = real(iz,aq_real)
      end do
      call afieldset%add(afield)
      call afield%final()
      !
   end subroutine aq_geom_fill_atlas_fieldset

   subroutine aq_geom_clone(self, other)
      class(aq_geom), intent(inout) :: self
      class(aq_geom), intent(in)    :: other
      !
      self%grid = atlas_structuredgrid(other%grid%c_ptr())
      self%fs = atlas_functionspace_structuredcolumns(other%fs%c_ptr())
      !
      self%nx = other%nx
      self%ny = other%ny
      self%nz = other%nz
      self%deltax = other%deltax
      self%deltay = other%deltay
      self%xmin = other%xmin
      self%ymin = other%ymin
      self%fmpi = other%fmpi
      self%domname = other%domname
      self%model = other%model
      self%orientation = other%orientation
      self%levels = other%levels
      self%mod_levels = other%mod_levels
      self%halo = other%halo
      self%bbox_imin = other%bbox_imin
      self%bbox_imax = other%bbox_imax
      self%bbox_isz  = other%bbox_isz
      self%bbox_jmin = other%bbox_jmin
      self%bbox_jmax = other%bbox_jmax
      self%bbox_jsz  = other%bbox_jsz
      allocate(self%lons(self%nx))
      allocate(self%lats(self%ny))
      self%lons(:) = other%lons(:)
      self%lats(:) = other%lats(:)
      if ( self%halo>0 ) then
         allocate( self%halo_mask(size(other%halo_mask)) )
         self%halo_mask(:) = other%halo_mask(:)
      end if
      !
   end subroutine aq_geom_clone

   subroutine aq_geom_delete(self)
      class(aq_geom), intent(inout)  :: self
      if (allocated(self%lons)) deallocate(self%lons)
      if (allocated(self%lats)) deallocate(self%lats)
      if (allocated(self%halo_mask)) deallocate(self%halo_mask)
   end subroutine aq_geom_delete

   subroutine aq_geom_info(self, nx, ny, nz, deltax, deltay, xmin, ymin, mod_levels, halo, domname, orientation, model)
      class(aq_geom), intent(in)                  :: self
      integer(aq_int), intent(out)                :: nx, ny, nz
      real(aq_real), intent(out)                  :: deltax, deltay
      real(aq_real), intent(out)                  :: xmin, ymin
      integer(aq_int), intent(out)                :: mod_levels
      integer(aq_int), intent(out)                :: halo
      character(len=:), allocatable , intent(out) :: domname
      character(len=:), allocatable , intent(out) :: orientation
      character(len=:), allocatable , intent(out) :: model
      !
      nx = self%nx
      ny = self%ny
      nz = self%nz
      xmin = self%xmin
      ymin = self%ymin
      deltax = self%deltax
      deltay = self%deltay
      mod_levels = self%mod_levels
      halo = self%halo
      domname = self%domname
      orientation = self%orientation
      model = self%model
      !
   end subroutine aq_geom_info

end module aq_geom_mod
