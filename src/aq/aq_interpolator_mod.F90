module aq_interpolator_mod

use aq_constants_mod
use aq_geom_mod
use aq_fields_mod

use atlas_module
use fckit_module
use fckit_configuration_module, only: fckit_configuration
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
   logical :: use_atlas
   type(csr_format) :: Hmat
   type(atlas_Interpolation) :: interp
   type(atlas_functionspace_PointCloud) :: locs_ptcloud
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
subroutine aq_interpolator_create(self, config, geom, loc_nlocs, lats, lons)
   !
   implicit none
   !
   class(aq_interpolator), intent(inout) :: self
   type(fckit_Configuration), intent(in) :: config
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
   type(atlas_field) :: afield
   real(kind_real), dimension(:,:), pointer :: data_ptr
   type(atlas_Config) :: int_config

   ! Number of local obs locations
   self%loc_nlocs = loc_nlocs

   ! Select the interpolator
   if (config%has("use atlas")) then
      call config%get_or_die("use atlas",self%use_atlas)
   else
      self%use_atlas = .false.
   end if
   !
   if (self%use_atlas) then
      ! Observation locations described as an Atlas PointCloud functionspace
      afield = atlas_field(name='lonlat',kind=atlas_real(kind_real),&
         &                 shape=(/2,loc_nlocs/))
      call afield%data(data_ptr)
      do ib=1, loc_nlocs
         data_ptr(1,ib) = lons(ib)
         data_ptr(2,ib) = lats(ib)
      end do
      self%locs_ptcloud = atlas_functionspace_pointcloud(afield)

      ! Create the interpolator
      int_config = atlas_Config()
      call int_config%set("type", "structured-linear2D")
      call int_config%set("adjoint", .true.)

      !AQ WAITING FOR ATLAS FIX
      !AQ self%interp = atlas_Interpolation(int_config, &
      !AQ               & geom%fs_surf, self%locs_ptcloud)

      call afield%final()
      !
   else
      ! Interpolator needed only where there are observations
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
   end if

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
   if (this%use_atlas) then
      call this%interp%final()
      call this%locs_ptcloud%final()
   else
      if (this%loc_nlocs > 0) call deallocate_operator(this%Hmat)
   end if
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
   integer       :: glo_nlocs, offset, jvar, ib, ib_i, ib_j
   real(aq_real), allocatable, dimension(:,:) :: surf_fld
   character(len=aq_strlen) :: fname
   type(atlas_Field)  :: aloc_2d
   type(atlas_Field)  :: avals
   real(aq_single), pointer  :: loc_2ds(:)
   real(aq_real), pointer  :: loc_2d(:)
   real(aq_real), pointer  :: pvals(:)
   !
   if (trim(field%geom%orientation) == 'down') nlev = field%geom%levels
   call field%fmpi%allreduce(self%loc_nlocs,glo_nlocs,fckit_mpi_sum())
   if ( glo_nlocs == 0 ) return

   if (self%use_atlas) then
      offset = 0
      do jvar=1,vars%nvars()
         fname = vars%variable(jvar)
         aloc_2d = field%geom%fs_surf%create_field(name=fname,  &
            &                                     kind=atlas_real(field%prec), &
            &                                     levels=0)
         !
         ! Inject fields from the 3d representation into the 2D with halo
         !
         if (field%prec == aq_single) then
            call aloc_2d%data(loc_2ds)
            do ib_j = field%geom%fs%j_begin(), field%geom%fs%j_end()
               do ib_i =  field%geom%fs%i_begin(ib_j), field%geom%fs%i_end(ib_j)
                  loc_2d(field%geom%fs_surf%index(ib_i,ib_j)) = &
                     & field%fldss(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j))
               end do
            end do
         else
            call aloc_2d%data(loc_2d)
            do ib_j = field%geom%fs%j_begin(), field%geom%fs%j_end()
               do ib_i =  field%geom%fs%i_begin(ib_j), field%geom%fs%i_end(ib_j)
                  loc_2d(field%geom%fs_surf%index(ib_i,ib_j)) = &
                     & field%fldsd(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j))
               end do
            end do
         end if
         call field%geom%fs_surf%halo_exchange(aloc_2d)
         !
         ! Prepare the vals storage
         !
         avals = self%locs_ptcloud%create_field(name=fname,&
            &                                   kind=aq_real, &
            &                                   levels=0)

         !AQ WAITING FOR ATLAS FIX
         !AQ      call self%interp%execute(aloc_2d,avals)
         call avals%data(pvals)
         if ( self%loc_nlocs > 0 ) then
            vals(offset+1:offset+self%loc_nlocs) = pvals(:)
            do ib = 1, self%loc_nlocs
               if (mask(ib)==0) vals(offset+ib) = aq_missing_value
            end do
            ! Update offset
            offset = offset+self%loc_nlocs
         end if

         call aloc_2d%final()
         call avals%final()

      end do

      if (size(vals) /= offset) call abor1_ftn('aq_field_getvals: error size')

   else

      allocate(surf_fld(field%geom%grid%nx(1),field%geom%grid%ny()))

      offset = 0
      do jvar=1,vars%nvars()
         fname = vars%variable(jvar)

         ! Gather surface field on master
         call field%gather_var_at_lev(trim(fname), nlev, surf_fld, 0)
         ! Brodcast it to all processes
         call field%fmpi%broadcast(surf_fld, root=0)

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

   end if

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
   integer       :: loc_nlocs, glo_nlocs, offset, jvar, ib, ib_i, ib_j
   real(aq_real), allocatable, dimension(:) :: surf_1d(:)
   real(aq_real), allocatable, dimension(:,:) :: surf_fld
   real(aq_real), allocatable, dimension(:) :: masked_vals(:)
   character(len=aq_strlen) :: fname
   type(atlas_Field)  :: aloc_2d
   type(atlas_Field)  :: avals
   real(aq_single), pointer  :: loc_2ds(:)
   real(aq_real), pointer  :: loc_2d(:)
   real(aq_real), pointer  :: pvals(:)

   if (trim(field%geom%orientation) == 'down') nlev = field%geom%levels

   call field%fmpi%allreduce(self%loc_nlocs,glo_nlocs,fckit_mpi_sum())
   if ( glo_nlocs == 0 ) return

   if ( self%use_atlas ) then
      offset = 0
      do jvar=1,vars%nvars()
         fname = vars%variable(jvar)
         !
         ! Prepare the vals storage
         !
         avals = self%locs_ptcloud%create_field(name=fname,&
            &                                   kind=aq_real, &
            &                                   levels=0)

         call avals%data(pvals)
         if ( self%loc_nlocs > 0 ) then
            do ib = 1, self%loc_nlocs
               if ( vals(offset + ib) /= aq_missing_value .and. &
                  & mask(ib) /= 0 ) then
                  pvals(ib) = vals(offset + ib)
               else
                  pvals(ib) = 0.
               end if
            end do
            ! Update offset
            offset = offset+self%loc_nlocs
         end if

         aloc_2d = field%geom%fs_surf%create_field(name=fname,  &
            &                                      kind=atlas_real(field%prec), &
            &                                      levels=0)

         !AQ WAITING FOR ATLAS FIX
         !AQ      call self%interp%execute_adjoint(aloc_2d,avals)
         !
         ! Inject fields from the 3d representation into the 2D with halo
         !
         if (field%prec == aq_single) then
            call aloc_2d%data(loc_2ds)
            do ib_j = field%geom%fs%j_begin(), field%geom%fs%j_end()
               do ib_i =  field%geom%fs%i_begin(ib_j), field%geom%fs%i_end(ib_j)
                  field%fldss(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j)) = &
                     & field%fldss(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j)) + &
                     & loc_2d(field%geom%fs_surf%index(ib_i,ib_j))
               end do
            end do
         else
            call aloc_2d%data(loc_2d)
            do ib_j = field%geom%fs%j_begin(), field%geom%fs%j_end()
               do ib_i =  field%geom%fs%i_begin(ib_j), field%geom%fs%i_end(ib_j)
                  field%fldsd(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j)) = &
                     & field%fldsd(jvar)%fld(nlev,field%geom%fs%index(ib_i,ib_j)) + &
                     & loc_2d(field%geom%fs_surf%index(ib_i,ib_j))
               end do
            end do
         end if

         call aloc_2d%final()
         call avals%final()

      end do

      if (size(vals) /= offset) call abor1_ftn('aq_field_getvalsad: error size')

   else

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

         ! Collect contributions from all the processes
         call field%fmpi%allreduce(surf_fld,fckit_mpi_sum())
         ! Add the result to the adjoint field on all processors
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

   end if

end subroutine aq_interpolator_applyAD
!
end module aq_interpolator_mod
