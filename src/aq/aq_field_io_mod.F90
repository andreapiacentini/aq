module aq_field_io_mod
   use aq_constants_mod
   use aq_geom_mod
   use datetime_mod
   use duration_mod
   use atlas_module
   use fckit_module
   use netcdf
   use, intrinsic :: iso_c_binding

   implicit none

   private
   public :: aq_write_field_gmsh, &
      &      aq_read_mocage_nc, &
      &      aq_write_mocage_nc

contains

   subroutine aq_nferr(info)

      ! Passed variables
      integer,intent(in) :: info !< Info index

      ! Check status
      if (info/=nf90_noerr) call abor1_ftn(trim(nf90_strerror(info)))

   end subroutine aq_nferr

   subroutine aq_write_field_gmsh(afset, geom, file)

      class(atlas_FieldSet),  intent(in) :: afset
      type(aq_geom),          intent(in) :: geom
      character(len=*),       intent(in) :: file

      type(atlas_Output) :: gmsh
      type(atlas_Mesh) :: mesh
      type(atlas_Partitioner) :: partitioner
      type(atlas_GridDistribution) :: griddistribution
      type(atlas_MeshGenerator) :: meshgenerator

      partitioner = atlas_Partitioner(type="checkerboard")
      griddistribution = partitioner%partition(geom%grid)
      call partitioner%final()
      meshgenerator = atlas_MeshGenerator()
      mesh = meshgenerator%generate(geom%grid, griddistribution)
      call meshgenerator%final()
      call griddistribution%final()
      gmsh = atlas_output_Gmsh(trim(file))
      call gmsh%write(mesh)
      call mesh%final()

      call gmsh%write(afset)

      call gmsh%final()

   end subroutine aq_write_field_gmsh

   subroutine aq_read_mocage_nc(afset, vars, geom, file, date)

      class(atlas_FieldSet),  intent(inout) :: afset
      character(len=*),       intent(in)    :: vars(:)
      type(aq_geom),          intent(in)    :: geom
      character(len=*),       intent(in)    :: file
      type(datetime),         intent(in)    :: date

      integer :: ncid, varid, dimid, nblevels
      integer, allocatable :: mod_levs(:)
      type(atlas_Field) :: afld
      real(aq_single), pointer :: flds(:,:)
      real(aq_real), pointer :: fldd(:,:)
      real(aq_single), allocatable :: ncbuff3ds(:,:,:)
      real(aq_real), allocatable :: ncbuff3dd(:,:,:)
      integer(atlas_kind_idx) :: ib_i, ib_j, ib_k, ib_var, il_lev
      logical :: ll_sgl
      character(len=:), allocatable :: cl_pos
      integer :: poslen
      logical :: ll_up
      integer :: mod_levels

      if (geom%fmpi%size() == 1) then
         call aq_nferr(nf90_open(trim(file), NF90_NOWRITE, ncid))
      else
         call aq_nferr(nf90_open_par(trim(file), &
         &                           ior(NF90_NOWRITE, NF90_MPIIO), &
         &                           comm = geom%fmpi%communicator(), &
         &                           info = geom%fmpi%info_null(), ncid=ncid))
      end if

      afld = afset%field(trim(vars(1)))
      ll_sgl = afld%kind() == aq_single
      call afld%final()
      allocate(ncbuff3ds(geom%bbox_imin:geom%bbox_imax, &
         &               geom%bbox_jmin:geom%bbox_jmax, &
         &               geom%levels))
      if (.not.ll_sgl) then
         allocate(ncbuff3dd(geom%bbox_imin:geom%bbox_imax, &
            &               geom%bbox_jmin:geom%bbox_jmax, &
            &               geom%levels))
      end if
      call aq_nferr(nf90_inq_dimid(ncid, "lev", dimid))
      call aq_nferr(nf90_inquire_dimension(ncid, dimid, len = nblevels))
      allocate(mod_levs(nblevels))
      call aq_nferr(nf90_inq_varid(ncid, "lev", varid))
      call aq_nferr(nf90_get_var(ncid, varid, mod_levs))
      mod_levels = maxval(mod_levs)
      if (mod_levels /= geom%mod_levels) &
         & call abor1_ftn("Number of model levels in file inconsistent with 'model levels' geometry entry")
      deallocate(mod_levs)
      if (nf90_inquire_attribute(ncid, varid, "positive", len=poslen) == nf90_noerr) then
         allocate(character(len=poslen)::cl_pos)
         call aq_nferr(nf90_get_att(ncid, varid, "positive", cl_pos))
      else
         cl_pos = "down"
      end if
      ll_up = trim(geom%orientation) /= trim(cl_pos)
      do ib_var = 1, size(vars)
         afld = afset%field(trim(vars(ib_var)))
         if (ll_sgl) then
            call afld%data(flds)
         else
            call afld%data(fldd)
         end if
         call aq_nferr(nf90_inq_varid(ncid, trim(vars(ib_var)), varid))
         call aq_nferr(nf90_get_var(ncid, varid, ncbuff3ds, &
            &                       start = [geom%bbox_imin, geom%bbox_jmin, nblevels-geom%levels+1], &
            &                       count = shape(ncbuff3ds)))
         if (.not.ll_sgl) ncbuff3dd = real(ncbuff3ds,kind=aq_real)
         do ib_j = geom%fs%j_begin(), geom%fs%j_end()
            do ib_i = geom%fs%i_begin(ib_j), geom%fs%i_end(ib_j)
               do ib_k = 1,geom%levels
                  if ( ll_up ) then
                     il_lev = geom%levels-ib_k+1
                  else
                     il_lev = ib_k
                  end if
                  if (ll_sgl) then
                     flds(il_lev,geom%fs%index(ib_i, ib_j)) = &
                        & ncbuff3ds(ib_i, ib_j, ib_k)
                  else
                     fldd(il_lev,geom%fs%index(ib_i, ib_j)) = &
                        & ncbuff3dd(ib_i, ib_j, ib_k)
                  end if
               end do
            end do
         end do
         call afld%final()
      end do
      deallocate(ncbuff3ds)
      if (.not.ll_sgl) deallocate(ncbuff3dd)
      call aq_nferr(nf90_close(ncid))
      call afset%halo_exchange()

   end subroutine aq_read_mocage_nc

   subroutine aq_write_mocage_nc(afset, vars, geom, file, date)

      class(atlas_FieldSet),  intent(in) :: afset
      character(len=*),       intent(in) :: vars(:)
      type(aq_geom),          intent(in) :: geom
      character(len=*),       intent(in) :: file
      type(datetime),         intent(in) :: date

      integer :: ncid, il_create_mode
      integer :: il_dlat, il_dlon, il_dlev, il_dtime
      integer :: il_vlat, il_vlon, il_vlev, il_vtime
      integer, allocatable :: ila_varid(:)
      integer :: ib, ib_var, ib_i, ib_j, ib_k
      type(atlas_Field) :: aloc_3d, aglo_3d
      real(aq_single), pointer :: loc_3ds(:,:), glo_3ds(:,:,:), glo_3dn(:,:,:)
      real(aq_real), pointer :: loc_3dd(:,:), glo_3dd(:,:,:)
      character(len=20) :: cl_date
      type(datetime) :: ref_date
      type(duration) :: dt
      integer(c_int64_t) :: timesec
      logical :: ll_sgl

      aloc_3d = afset%field(trim(vars(1)))
      ll_sgl = aloc_3d%kind() == aq_single
      call aloc_3d%final()

      if( geom%fmpi%rank() == 0 ) then
         il_create_mode =      NF90_NETCDF4
         il_create_mode = ior( NF90_CLOBBER       , il_create_mode)
         call aq_nferr(nf90_create(trim(file), il_create_mode, ncid))
         call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.0'))
         call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'source','AQ-4DEnVAR'))
         call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'field',trim(afset%name())))
         call datetime_create("1970-01-01T00:00:00Z",ref_date)
         call datetime_diff(date, ref_date, dt)
         timesec = duration_seconds(dt)
         call datetime_to_string(date, cl_date)
         if (cl_date /= "1970-01-01T00:00:00Z") &
            & call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'date',cl_date))
         call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'model','MOCAGE'))
         call aq_nferr(nf90_put_att(ncid, NF90_GLOBAL, 'institution','CERFACS - SEEDS'))

         call aq_nferr(nf90_def_dim(ncid,'lat',geom%grid%ny(), il_dlat))
         call aq_nferr(nf90_def_dim(ncid,'lon',geom%grid%nx(1), il_dlon))
         call aq_nferr(nf90_def_dim(ncid,'lev',geom%levels, il_dlev))
#ifdef NDEBUG
         call aq_nferr(nf90_def_dim(ncid,'time', NF90_UNLIMITED, il_dtime))
#endif
         !
         call aq_nferr(nf90_def_var(ncid,'lat',NF90_FLOAT,[il_dlat],il_vlat))
         call aq_nferr(nf90_put_att(ncid,il_vlat,'standard_name','latitude'))
         call aq_nferr(nf90_put_att(ncid,il_vlat,'units','degrees_north'))
         call aq_nferr(nf90_def_var(ncid,'lon',NF90_FLOAT,[il_dlon],il_vlon))
         call aq_nferr(nf90_put_att(ncid,il_vlon,'standard_name','longitude'))
         call aq_nferr(nf90_put_att(ncid,il_vlon,'units','degrees_east'))
         call aq_nferr(nf90_def_var(ncid,'lev',NF90_FLOAT,[il_dlev],il_vlev))
         call aq_nferr(nf90_put_att(ncid,il_vlev,'axis','z'))
         call aq_nferr(nf90_put_att(ncid,il_vlev, 'units', 'levels'))
         if (trim(geom%orientation) == "up" ) then
            call aq_nferr(nf90_put_att(ncid,il_vlev, 'positive', 'up'))
         else
            call aq_nferr(nf90_put_att(ncid,il_vlev, 'positive', 'down'))
            call aq_nferr(nf90_put_att(ncid,il_vlev,'standard_name','depth'))
            call aq_nferr(nf90_put_att(ncid,il_vlev,'formula_terms',&
               & 'ap: a_hybr_coord b: b_hybr_coord ps: air_pressure_at_surface'))
         end if
#ifdef NDEBUG
         call aq_nferr(nf90_def_var(ncid,'time',NF90_INT64,[il_dtime],il_vtime))
         call aq_nferr(nf90_put_att(ncid,il_vtime,'axis','t'))
         call aq_nferr(nf90_put_att(ncid,il_vtime,'standard_name','Time axis'))
         call aq_nferr(nf90_put_att(ncid,il_vtime,'units',&
            &'seconds since 1970-01-01 00:00:00'))
         call aq_nferr(nf90_put_att(ncid,il_vtime,'time_origin','1970-01-01 00:00:00'))
         call aq_nferr(nf90_put_att(ncid,il_vtime,'calendar','standard'))
#endif
         allocate(ila_varid(size(vars)))

         do ib_var = 1, size(vars)
#ifdef NDEBUG
            call aq_nferr(nf90_def_var(ncid,trim(vars(ib_var)),NF90_FLOAT,&
               &                   [il_dlon,il_dlat,il_dlev,il_dtime],ila_varid(ib_var)))
#else
            call aq_nferr(nf90_def_var(ncid,trim(vars(ib_var)),NF90_FLOAT,&
               &                   [il_dlon,il_dlat,il_dlev],ila_varid(ib_var)))
#endif
            call aq_nferr(nf90_put_att(ncid,ila_varid(ib_var),'units','ppb'))
         end do

         call aq_nferr(nf90_enddef(ncid))

         call aq_nferr(nf90_put_var(ncid,il_vlat, &
            & real([(geom%grid%y(ib),ib=1,geom%grid%ny())],kind=aq_single)))

         call aq_nferr(nf90_put_var(ncid,il_vlon, &
            & real([(geom%grid%x(ib,1),ib=1,geom%grid%nx(1))],kind=aq_single)))

         if ( geom%orientation == "up" ) then
            call aq_nferr(nf90_put_var(ncid,il_vlev, &
               & real([(geom%mod_levels-ib+1,ib=1,geom%levels)],kind=aq_single)))
         else
            call aq_nferr(nf90_put_var(ncid,il_vlev, &
               & real([(geom%mod_levels-geom%levels+ib,ib=1,geom%levels)],kind=aq_single)))
         end if
#ifdef NDEBUG
         call aq_nferr(nf90_put_var(ncid,il_vtime,timesec,start=[1]))
#endif
      end if

      if( geom%fmpi%rank() == 0 ) &
         & allocate(glo_3dn(geom%grid%nx(1),geom%grid%ny(),geom%levels))

      if (ll_sgl) then
         aglo_3d = geom%fs%create_field(name='globuff',    &
            &                              kind=atlas_real(aq_single), &
            &                              global = .true.)
         if( geom%fmpi%rank() == 0 ) &
            & call aglo_3d%data(glo_3ds,shape=[geom%levels,geom%grid%nx(1),geom%grid%ny()])

         do ib_var = 1, size(vars)
            aloc_3d = afset%field(trim(vars(ib_var)))
            call aloc_3d%data(loc_3ds)
            call geom%fs%gather(aloc_3d, aglo_3d)
            call aloc_3d%final()
            if( geom%fmpi%rank() == 0 ) then
               do ib_j = 1,geom%grid%ny()
                  do ib_i = 1,geom%grid%nx(1)
                     do ib_k = 1,geom%levels
                        glo_3dn(ib_i,ib_j,ib_k) = glo_3ds(ib_k,ib_i,ib_j)
                     end do
                  end do
               end do
               call aq_nferr(nf90_put_var(ncid,ila_varid(ib_var),glo_3dn))
            end if
         end do

      else
         aglo_3d = geom%fs%create_field(name='globuff',    &
            &                              kind=atlas_real(aq_real), &
            &                              global = .true.)
         if( geom%fmpi%rank() == 0 ) &
            & call aglo_3d%data(glo_3dd,shape=[geom%levels,geom%grid%nx(1),geom%grid%ny()])

         do ib_var = 1, size(vars)
            aloc_3d = afset%field(trim(vars(ib_var)))
            call aloc_3d%data(loc_3dd)
            call geom%fs%gather(aloc_3d, aglo_3d)
            call aloc_3d%final()
            if( geom%fmpi%rank() == 0 ) then
               do ib_j = 1,geom%grid%ny()
                  do ib_i = 1,geom%grid%nx(1)
                     do ib_k = 1,geom%levels
                        glo_3dn(ib_i,ib_j,ib_k) = real(glo_3dd(ib_k,ib_i,ib_j))
                     end do
                  end do
               end do
               call aq_nferr(nf90_put_var(ncid,ila_varid(ib_var),glo_3dn))
            end if
         end do

      end if

      if( geom%fmpi%rank() == 0 ) deallocate(glo_3dn)

      call aglo_3d%final()

      if( geom%fmpi%rank() == 0 ) then
         call aq_nferr(nf90_close(ncid))
      end if

      call geom%fmpi%barrier()

   end subroutine aq_write_mocage_nc

end module aq_field_io_mod
