! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module aq_obsdb_mod

use atlas_module
use datetime_mod
use duration_mod
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
! use netcdf
use aq_locs_mod
use aq_obsvec_mod
use random_mod
use string_f_c_mod
use H5_UTILS_MOD, ONLY : ip_hdf_namelen, ig_hdfverb, ip_hid_t, h5state_t, &
  open_h5group, close_h5group, close_h5space, check_h5file, Get_h5dset_size, &
  create_h5file, create_h5group_instrdom, create_h5group_instr, create_h5group, &
  create_attrib_string
use H5_READ_MOD, only : open_h5file_rdonly, close_h5file
USE HDF5

implicit none

private
public :: aq_obsdb
public :: aq_obsdb_registry
public :: aq_obsdb_setup,aq_obsdb_delete,aq_obsdb_read,aq_obsdb_get,aq_obsdb_put,aq_obsdb_locations,aq_obsdb_generate,aq_obsdb_nobs
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 1 !< Random seed (for reproducibility)
type(datetime), save :: obs_ref_time 

type column_data
  character(len=50) :: colname                !< Column name
  type(column_data),pointer :: next => null() !< Next column
  integer :: nlev                             !< Number of levels
  real(kind_real),allocatable :: values(:,:)  !< Values
end type column_data

type group_data
  character(len=ip_hdf_namelen) :: grpname       !< Group name
  type(group_data),pointer :: next => null()     !< Next group
  integer :: nobs                                !< Number of observations
  type(datetime),allocatable :: times(:)         !< Time-slots
  type(column_data),pointer :: colhead => null() !< Head column
end type group_data

type aq_obsdb
  integer :: ngrp = 1                           !< Number of groups (1 obs space == 1 instrument == 1 chemical spc in aq)
  character(len=:),allocatable :: instrname     !< Instrument name 
  character(len=:),allocatable :: spcname       !< Name of the chemical species
  character(len=1024) :: filein                 !< Input filename
  character(len=1024) :: fileout                !< Output filename
  type(datetime) :: winbgn            !< Start of window
  type(datetime) :: winend            !< End of window
  integer(kind=ip_hid_t) :: il_hdat_id = 0
  integer(kind=ip_hid_t) :: il_hstat_id = 0
  type(h5state_t) :: h5statein
  type(h5state_t) :: h5stateout
  type(group_data),pointer :: grphead => null() !< Head group
end type aq_obsdb

#define LISTED_TYPE aq_obsdb

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: aq_obsdb_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup observation data
subroutine aq_obsdb_setup(self,f_conf,winbgn,winend,openfile,other)
use string_utils

implicit none

! Passed variables
type(aq_obsdb),intent(inout) :: self           !< Observation data
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(in) :: winbgn            !< Start of window
type(datetime),intent(in) :: winend            !< End of window
logical,intent(in) :: openfile
type(aq_obsdb),intent(in) :: other             !< Used only for file id

! Local variables
character(len=1024) :: fin,fout
character(len=:),allocatable :: str
integer :: il_err
logical :: ll_exists
character(len=12) :: cl_obsgrp = 'OBSERVATIONS'

! Reference time for hdf5 I/O files
call datetime_create("1970-01-01T00:00:00Z",obs_ref_time)
ig_hdfverb = 1
! Input file
if (f_conf%has("obsdatain")) then
  call f_conf%get_or_die("obsdatain.obsfile",str)
  fin = str
else
  fin = ''
endif
call fckit_log%info('aq_obsdb_setup: file in = '//trim(fin))

! Output file
if (f_conf%has("obsdataout")) then
  call f_conf%get_or_die("obsdataout.obsfile",str)
  call swap_name_member(f_conf, str)

  fout = str
  call fckit_log%info('aq_obsdb_setup: file out = '//trim(fout))
else
  fout = ''
endif
! Set attributes
self%filein = fin
self%fileout = fout
self%winbgn = winbgn
self%winend = winend
call f_conf%get_or_die("instr name",self%instrname)
call f_conf%get_or_die("obs type",self%spcname)

if (openfile) then
  ! Init the hdf5 fortran library 
  ! should go in our HDF5 library in a better place,
  ! cause it was hidden in CREATE_H5FILE, which makes no sense
  call H5open_f (il_err)
  ! Open input hdf5 file
  call open_h5file_rdonly(self%filein, self%il_hdat_id)
  call check_h5file(self%il_hdat_id, trim(self%instrname), "Surface")
  call fckit_log%info('aq_obsdb_setup: file in opened correctly')

  ! Create HSTAT.h5 output file
  call create_h5file(self%fileout, self%il_hstat_id)
else
  self%il_hdat_id = other%il_hdat_id
  self%il_hstat_id = other%il_hstat_id
endif

call open_h5group(self%il_hdat_id, '/'//trim(self%instrname), self%h5statein%instr_id)

call H5Lexists_f(self%il_hstat_id, trim(self%instrname), ll_exists, il_err)

if (.not.ll_exists) then
  call H5Gcreate_f(self%il_hstat_id, trim(self%instrname), &
         & self%h5stateout%instr_id, il_err)
  call create_attrib_string(self%h5stateout%instr_id, 'MeasurementType', "Surface")

  ! Create instrument/domain  sub-groups
  call create_h5group(self%h5stateout%instr_id, 'GEOLOCALIZATION')
  call create_h5group(self%h5stateout%instr_id, 'OBSERVATIONS')
  call create_h5group(self%h5stateout%instr_id,trim(cl_obsgrp)//'/'//trim(self%spcname))
endif

end subroutine aq_obsdb_setup
! ------------------------------------------------------------------------------
!> Delete observation data
subroutine aq_obsdb_delete(self,closefile)

implicit none

! Passed variables
type(aq_obsdb),intent(inout) :: self !< Observation data
logical,intent(in) :: closefile

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs

! Write observation data
if (self%fileout/='') call aq_obsdb_write(self)

if (closefile) then
  call close_h5group(self%h5statein%instr_id)
  call close_h5group(self%h5stateout%instr_id)
  call close_h5file(self%il_hdat_id, .false.)
  call close_h5file(self%il_hstat_id, .true.)
endif

! Release memory
do while (associated(self%grphead))
  jgrp => self%grphead
  self%grphead => jgrp%next
  do jobs=1,jgrp%nobs
    call datetime_delete(jgrp%times(jobs))
  enddo
  deallocate(jgrp%times)
  do while (associated(jgrp%colhead))
    jcol => jgrp%colhead
    jgrp%colhead => jcol%next
    deallocate(jcol%values)
    deallocate(jcol)
  enddo
  deallocate(jgrp)
enddo

end subroutine aq_obsdb_delete
! ------------------------------------------------------------------------------
!> Get observation data
subroutine aq_obsdb_get(self,grp,col,ovec)

implicit none

! Passed variables
type(aq_obsdb),intent(in) :: self     !< Observation data
character(len=*),intent(in) :: grp    !< Group
character(len=*),intent(in) :: col    !< Column
type(aq_obsvec),intent(inout) :: ovec !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Print Group and column
call fckit_log%debug('aq_obsdb_get: grp = '//trim(grp))
call fckit_log%debug('aq_obsdb_get: col = '//trim(col))

! Find observation group
call aq_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('aq_obsdb_get: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('aq_obsdb_get: cannot find '//trim(grp))
  call abor1_ftn('aq_obsdb_get: obs group not found')
endif

! Find observation column
call aq_obsdb_find_column(jgrp,col,jcol)

if (.not.associated(jcol)) then
  call fckit_log%error('aq_obsdb_get: cannot find '//trim(col))
  call abor1_ftn('aq_obsdb_get: obs column not found')
endif

! Get observation data
if (allocated(ovec%values)) deallocate(ovec%values)
ovec%nlev = jcol%nlev
! get all the obs
ovec%nobs = jgrp%nobs
allocate(ovec%values(ovec%nlev,ovec%nobs))
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    ovec%values(jlev,jobs) = jcol%values(jlev,jobs)
  enddo
enddo

end subroutine aq_obsdb_get
! ------------------------------------------------------------------------------
!> Put observations data
subroutine aq_obsdb_put(self,grp,col,ovec)

implicit none

! Passed variables
type(aq_obsdb),intent(inout) :: self !< Observation data
character(len=*),intent(in) :: grp   !< Group
character(len=*),intent(in) :: col   !< Column
type(aq_obsvec),intent(in) :: ovec   !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Find observation group
call aq_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('aq_obsdb_put: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('aq_obsdb_put: cannot find '//trim(grp))
  call abor1_ftn('aq_obsdb_put: obs group not found')
endif

! Find observation column (and add it if not there)
call aq_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) then
  if (.not.associated(jgrp%colhead)) call abor1_ftn('aq_obsdb_put: no locations')
  jcol => jgrp%colhead
  do while (associated(jcol%next))
    jcol => jcol%next
  enddo
  allocate(jcol%next)
  jcol => jcol%next
  jcol%colname = col
  jcol%nlev = ovec%nlev
  allocate(jcol%values(jcol%nlev,jgrp%nobs))
endif

! Put observation data
if (ovec%nobs/=jgrp%nobs) call abor1_ftn('aq_obsdb_put: error obs number')
if (ovec%nlev/=jcol%nlev) call abor1_ftn('aq_obsdb_put: error col number')
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    jcol%values(jlev,jobs) = ovec%values(jlev,jobs)
  enddo
enddo

end subroutine aq_obsdb_put
! ------------------------------------------------------------------------------
!> Get locations from observation data
subroutine aq_obsdb_locations(self,grp,fields,c_times)

implicit none

! Passed variables
type(aq_obsdb),intent(in) :: self   !< Observation data
character(len=*),intent(in) :: grp  !< Group
type(atlas_fieldset), intent(inout) :: fields !< Locations FieldSet
type(c_ptr), intent(in), value :: c_times !< pointer to times array in C++

! Local variables
integer :: nlocs, jo
character(len=8),parameter :: col = 'Location'
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
type(atlas_field) :: field_z, field_lonlat
real(kind_real), pointer :: z(:), lonlat(:,:)

! Find observation group
call aq_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  call fckit_log%error('aq_obsdb_get: cannot find '//trim(grp))
  call abor1_ftn('aq_obsdb_locations: obs group not found')
endif
nlocs = jgrp%nobs

! Find observation column
call aq_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) call abor1_ftn('aq_obsdb_locations: obs column not found')

! Set number of observations

field_lonlat = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=[2,nlocs])
field_z = atlas_field(name="altitude", kind=atlas_real(kind_real), shape=[nlocs])

call field_lonlat%data(lonlat)
call field_z%data(z)

! Copy coordinates
do jo = 1, nlocs
  lonlat(1,jo) = jcol%values(1,jo)
  lonlat(2,jo) = jcol%values(2,jo)
  z(jo) = jcol%values(3,jo)
  call f_c_push_to_datetime_vector(c_times, jgrp%times(jo))
enddo

call fields%add(field_lonlat)
call fields%add(field_z)

! release pointers
call field_lonlat%final()
call field_z%final()

end subroutine aq_obsdb_locations
! ------------------------------------------------------------------------------
!> Generate observation data
subroutine aq_obsdb_generate(self,grp,f_conf,bgn,step,ktimes,kobs)

implicit none

! Passed variables
type(aq_obsdb),intent(inout) :: self           !< Observation data
character(len=*),intent(in) :: grp             !< Group
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(in) :: bgn               !< Start time
type(duration),intent(in) :: step              !< Time-step
integer,intent(in) :: ktimes                   !< Number of time-slots
integer,intent(inout) :: kobs                  !< Number of observations

! Local variables
integer :: nlev,nlocs
real(kind_real) :: err
type(datetime),allocatable :: times(:)
type(aq_obsvec) :: obsloc,obserr

! Get number of observations
call f_conf%get_or_die("obs_density",nlocs)
kobs = nlocs*ktimes

! Allocation
allocate(times(kobs))

! Generate locations
call aq_obsdb_generate_locations(nlocs,ktimes,bgn,step,times,obsloc)

! Create observations data
call aq_obsdb_create(self,trim(grp),times,obsloc)

! Create observation error
call f_conf%get_or_die("obs_error",err)
call f_conf%get_or_die("nval",nlev)
call aq_obsvec_setup(obserr,nlev,kobs)
obserr%values(:,:) = err
call aq_obsdb_put(self,trim(grp),'ObsError',obserr)

! Release memory
deallocate(times)
deallocate(obsloc%values)
deallocate(obserr%values)

end subroutine aq_obsdb_generate
! ------------------------------------------------------------------------------
!> Get observation data size
subroutine aq_obsdb_nobs(self,grp,kobs)

implicit none

! Passed variables
type(aq_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
integer,intent(inout) :: kobs      !< Number of observations

! Local variables
type(group_data),pointer :: jgrp

! Find group
call aq_obsdb_find_group(self,grp,jgrp)

! Get observation data size
if (associated(jgrp)) then
  kobs = jgrp%nobs
else
  kobs = 0
endif

end subroutine aq_obsdb_nobs
! ------------------------------------------------------------------------------
!> Read observation data
subroutine aq_obsdb_read(self)

use H5_READ_MOD, only : readslice_h5dset
use H5_SELECTION_MOD, only : Get_number_selected_timeelts

implicit none

! Passed variables
type(aq_obsdb),intent(inout) :: self !< Observation data

! Local variables
integer :: igrp,icol,iobs,ncol,nobs,jobs
integer :: ncid,grpname_id,ngrp_id,nobs_id,ncol_id,times_id,nlev_id,colname_id,values_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=6) :: igrpchar
character(len=50) :: stime
character(len=1024) :: record
logical, allocatable :: inwindow(:)
type(datetime) :: tobs
type(duration) :: dtwinbgn, dtwinend, dt
character(len=1024) :: timestr1, timestr2
real(kind_real),allocatable :: readbuf(:,:)
integer(kind=ip_hid_t) :: il_instr_idin
integer :: id_tmin, id_tmax
integer :: il_err
! integer :: ncols = 2
! character(len=50), dimension(ncols) :: colnames = ['Location','ObsValue']
! integer, dimension(ncols) :: coldims = [ 2, 1 ]
real(kind=4), dimension(:), allocatable :: rla_lats, rla_lons, rla_obs
integer(kind=4), dimension(:), allocatable :: ila_times
type(datetime),allocatable :: times(:)
type(aq_obsvec) :: obsloc,obsval,obserr

! Get the window begin and window end time in seconds using the obs_ref_time
call datetime_to_string(self%winbgn,timestr1)
call datetime_to_string(self%winend,timestr2)
call datetime_diff(self%winbgn,obs_ref_time,dtwinbgn)
call datetime_diff(self%winend,obs_ref_time,dtwinend)
! Count the obs for the given instrument in the hdf5 file
! Loss of precision here cause obs were initially based on mocage datetime, 
! this will stop working sometime during this century
id_tmin = int(duration_seconds(dtwinbgn),kind=4)
id_tmax = int(duration_seconds(dtwinend),kind=4)

call fckit_log%info('aq_obsdb_read: reading = '//trim(self%instrname))

nobs = Get_number_selected_timeelts(self%h5statein,id_tmin,id_tmax)
! do not read 0 size obs !!!
if (nobs == 0) then
  call fckit_log%info('aq_obsdb_read: no obs found between '//trim(timestr1)//' and '//trim(timestr2))
else
  ! Read the data from the hdf5
  allocate(ila_times(nobs),rla_lats(nobs),rla_lons(nobs),rla_obs(nobs))

  call readslice_h5dset(self%h5statein, 'GEOLOCALIZATION/Timestamp', ila_times)
  call readslice_h5dset(self%h5statein, 'GEOLOCALIZATION/Latitude', rla_lats)
  call readslice_h5dset(self%h5statein, 'GEOLOCALIZATION/Longitude', rla_lons)
  call readslice_h5dset(self%h5statein, 'OBSERVATIONS/'//trim(self%spcname)//'/Y', rla_obs)
  ! Setup observation vector for the locations
  call aq_obsvec_setup(obsloc,3,nobs)

  ! Setup observation vector for the observations
  call aq_obsvec_setup(obsval,1,nobs)
  call aq_obsvec_setup(obserr,1,nobs)
  allocate(times(nobs))

  ! Fill the arrays
  do iobs=1,nobs
    tobs = obs_ref_time
    dt = int(ila_times(iobs))
    call datetime_update(tobs,dt)
    times(iobs) = tobs
    obsloc%values(:,iobs) = (/real(rla_lons(iobs),kind=kind_real),real(rla_lats(iobs),kind=kind_real),0.0d0/)
    obsval%values(:,iobs) = rla_obs(iobs)
    obserr%values(:,iobs) = rla_obs(iobs)
  enddo

  ! Store observations data in the obsdb structure
  call aq_obsdb_create(self,trim(self%spcname),times,obsloc)
  call aq_obsdb_put(self,trim(self%spcname),'ObsValue',obsval)
  call aq_obsdb_put(self,trim(self%spcname),'ObsError',obserr) ! This should not be mandatory but it is asked by InSitu!!!!

  deallocate(ila_times,rla_lats,rla_lons,rla_obs,times)
endif

call close_h5space(self%h5statein%memspace_id)
call close_h5space(self%h5statein%dataspace_id) ! Verify if these calls are really needed here

end subroutine aq_obsdb_read
! ------------------------------------------------------------------------------
!> Write observation data
subroutine aq_obsdb_write(self)

use H5_WRITE_MOD, only : writeslice_h5dset_scalar

implicit none

! Passed variables
type(aq_obsdb),intent(in) :: self !< Observation data

! Local variables
integer :: igrp,icol,iobs,ncol,nlevmax
integer :: ncid,nstrmax_id,grpname_id,ngrp_id,nobs_id,ncol_id,nlevmax_id,times_id,nlev_id,colname_id,values_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=6) :: igrpchar
character(len=50) :: stime
integer(kind=ip_hid_t) :: il_hstat_id, il_instr_id
character(len=15) :: cl_geogrp = 'GEOLOCALIZATION'
character(len=12) :: cl_obsgrp = 'OBSERVATIONS'
character(len=19), dimension(:), allocatable :: cla_timehuman
integer :: il
integer :: il_err
integer, dimension(:), allocatable :: timestamp                 
type(duration) :: dtdiff                   

CALL H5open_f(il_err)

jgrp => self%grphead
if (jgrp%nobs > 0) then
  allocate(cla_timehuman(jgrp%nobs),timestamp(jgrp%nobs))
  do iobs=1,jgrp%nobs
    call datetime_to_string(jgrp%times(iobs),stime)
    cla_timehuman(iobs) = stime(:19)
    call datetime_diff(jgrp%times(iobs), obs_ref_time, dtdiff)
    timestamp(iobs) = int(duration_seconds(dtdiff),kind=4)
  end do
  call writeslice_h5dset_scalar(self%h5stateout, trim(cl_geogrp)//'/TimeHuman', cla_timehuman)
  call writeslice_h5dset_scalar(self%h5stateout, trim(cl_geogrp)//'/Timestamp', timestamp)
  deallocate(cla_timehuman,timestamp)

  ! call writeslice_h5dset_scalar(il_instr_id, trim(cl_geogrp)//'/Timestamp', IDINT(self%time))

  ! Loop over columns
  icol = 0
  jcol => jgrp%colhead
  do while (associated(jcol))
    icol = icol+1
    select case(trim(jcol%colname))
    case ('Location')
      call writeslice_h5dset_scalar(self%h5stateout, trim(cl_geogrp)//'/Longitude', jcol%values(1,:))
      call writeslice_h5dset_scalar(self%h5stateout, trim(cl_geogrp)//'/Latitude', jcol%values(2,:))
    case ('hofx')
      call writeslice_h5dset_scalar(self%h5stateout, trim(cl_obsgrp)//'/'//trim(self%spcname)//'/Hx', jcol%values(1,:))
    case ('ObsValue')
      call writeslice_h5dset_scalar(self%h5stateout, trim(cl_obsgrp)//'/'//trim(self%spcname)//'/Y', jcol%values(1,:))
    case ('ObsError')
       call writeslice_h5dset_scalar(self%h5stateout, trim(cl_obsgrp)//'/'//trim(self%spcname)//'/Covariance', jcol%values(1,:))
    case ('EffectiveQC')
       call writeslice_h5dset_scalar(self%h5stateout, trim(cl_obsgrp)//'/'//trim(self%spcname)//'/EffectiveQC', jcol%values(1,:))
    case ('EffectiveError')
       call writeslice_h5dset_scalar(self%h5stateout, trim(cl_obsgrp)//'/'//trim(self%spcname)//'/EffectiveError', jcol%values(1,:))
    case ('ObsBias')
       call fckit_log%info('Warning: ObsBias not implemented in aq_obsdb_write')
    case default    
       call fckit_log%info('Warning: '//trim(jcol%colname)//' not known in aq_obsdb_write')
    end select
    ! Update
    jcol => jcol%next
  enddo
endif

end subroutine aq_obsdb_write
! ------------------------------------------------------------------------------
!  Private
! ------------------------------------------------------------------------------
!> Find observation data group
subroutine aq_obsdb_find_group(self,grp,find)

implicit none

! Passed variables
type(aq_obsdb),intent(in) :: self              !< Observation data
character(len=*),intent(in) :: grp             !< Group
type(group_data),pointer,intent(inout) :: find !< Result

! Initialization
find => self%grphead

! Loop
do while (associated(find))
  if (find%grpname==grp) exit
  find => find%next
enddo

end subroutine aq_obsdb_find_group
! ------------------------------------------------------------------------------
!> Find observation data column
subroutine aq_obsdb_find_column(grp,col,find)

implicit none

! Passed variables
type(group_data),intent(in) :: grp              !< Observation data
character(len=*),intent(in) :: col              !< Column
type(column_data),pointer,intent(inout) :: find !< Result

! Initialization
find=>grp%colhead

! Loop
do while (associated(find))
  if (find%colname==col) exit
  find => find%next
enddo

end subroutine aq_obsdb_find_column
! ------------------------------------------------------------------------------
!> Generate random locations
subroutine aq_obsdb_generate_locations(nlocs,ntimes,bgn,step,times,obsloc)

implicit none

! Passed variables
integer,intent(in) :: nlocs                         !< Number of locations
integer,intent(in) :: ntimes                        !< Number of time-slots
type(datetime),intent(in) :: bgn                    !< Start time
type(duration),intent(in) :: step                   !< Time-step
type(datetime),intent(inout) :: times(nlocs*ntimes) !< Time-slots
type(aq_obsvec),intent(inout) :: obsloc             !< Observation locations

! Local variables
integer :: jobs,iobs,jstep
real(kind_real) :: x(nlocs),y(nlocs),z(nlocs),lon(nlocs),lat(nlocs)
type(datetime) :: now

! Generate random locations
! call uniform_distribution(x,0.0_kind_real,domain_zonal,rseed)
! call uniform_distribution(y,0.0_kind_real,domain_meridional,rseed)
! call uniform_distribution(z,0.0_kind_real,domain_depth,rseed)

! Convert to lon/lat
do jobs=1,nlocs
  ! call xy_to_lonlat(x(jobs),y(jobs),lon(jobs),lat(jobs))
enddo

! Setup observation vector
call aq_obsvec_setup(obsloc,3,nlocs*ntimes)

! Set observation locations
now = bgn
iobs=0
do jstep=1,ntimes
  do jobs=1,nlocs
    iobs = iobs+1
    times(iobs) = now
    obsloc%values(:,iobs) = (/lon(jobs),lat(jobs),z(jobs)/)
  enddo
  call datetime_update(now,step)
enddo

! Release memory
call datetime_delete(now)

end subroutine aq_obsdb_generate_locations
! ------------------------------------------------------------------------------
!> Create observation data
subroutine aq_obsdb_create(self,grp,times,locs)

implicit none

! Passed varaibles
type(aq_obsdb),intent(inout) :: self  !< Observation data
character(len=*),intent(in) :: grp    !< Group
type(datetime),intent(in) :: times(:) !< Time-slots
type(aq_obsvec),intent(in) :: locs    !< Locations

! Local variables
type(group_data),pointer :: igrp
integer :: jobs,jlev

! Find observation group
call aq_obsdb_find_group(self,grp,igrp)
if (associated(igrp)) call abor1_ftn('aq_obsdb_create: obs group already exists')
if (associated(self%grphead)) then
  igrp => self%grphead
  do while (associated(igrp%next))
    igrp => igrp%next
  enddo
  allocate(igrp%next)
  igrp => igrp%next
else
  allocate(self%grphead)
  igrp => self%grphead
endif

! Create observation data
igrp%grpname = grp
igrp%nobs = size(times)
allocate(igrp%times(igrp%nobs))
igrp%times(:) = times(:)
allocate(igrp%colhead)
igrp%colhead%colname = 'Location'
igrp%colhead%nlev = 3
allocate(igrp%colhead%values(3,igrp%nobs))
if (locs%nlev/=3) call abor1_ftn('aq_obsdb_create: error locations not 3D')
if (locs%nobs/=igrp%nobs) call abor1_ftn('aq_obsdb_create: error locations number')
do jobs=1,igrp%nobs
  do jlev=1,3
    igrp%colhead%values(jlev,jobs) = locs%values(jlev,jobs)
  enddo
enddo
self%ngrp = self%ngrp+1

end subroutine aq_obsdb_create
! ------------------------------------------------------------------------------
end module aq_obsdb_mod
