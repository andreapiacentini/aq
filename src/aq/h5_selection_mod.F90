! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!! category: observation operator
!! summary:  Prepare selection of elements in HDF5 dataset
!! author:   CERFACS and CNRM (G. Jonville)
!!
MODULE H5_SELECTION_MOD
   !! category: observation operator
   !! author:   CERFACS and CNRM (G. Jonville)
   !!
   !! ----------------------------------------------------------------------
   !!
   !! Functions to prepare selection of data from a time or geographical restriction
   !!
   !! ----------------------------------------------------------------------
   !!
   USE H5_READ_MOD

   IMPLICIT NONE

   INTEGER(hsize_t), DIMENSION(1) :: ila_start

CONTAINS
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Get_number_selected_timeelts(h5state, id_tmin, id_tmax)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Function to get number of selected elements in time,
      !! and to initialize data and memory spaces HDF5 identifiers
      !!
      !! ----------------------------------------------------------------------
      !!
      USE H5LT

      TYPE(h5state_t), INTENT(INOUT) :: h5state
      INTEGER,         INTENT(IN) :: id_tmin
      INTEGER,         INTENT(IN) :: id_tmax

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_memspace_id
      ! Dimensions
      INTEGER(hsize_t), DIMENSION(1) :: ila_dims
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      ! Inquiries
      INTEGER :: il_err
      ! Temporary buffers
      CHARACTER(LEN=ip_hdf_namelen) :: cl_dsname
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size
      INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ila_time
      INTEGER, DIMENSION(1) :: ila_minidx
      INTEGER, DIMENSION(1) :: ila_maxidx
      INTEGER :: il_tselsize

      ! Read time dataset in the full
      cl_dsname = 'GEOLOCALIZATION/Timestamp'
      CALL H5LTget_dataset_info_f(h5state%instr_id, cl_dsname, &
         & ila_dims, il_type_class, il_type_size, il_err)

      CALL H5Dopen_f(h5state%instr_id, cl_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,*) 'Dataset opening error ', il_err

      CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
      IF (il_err /= 0) WRITE(*,*) 'Dataspace get error ', il_err

      IF (ila_dims(1).GT.0) THEN
         ALLOCATE(ila_time(ila_dims(1)))
         CALL H5LTread_dataset_f(h5state%instr_id, cl_dsname, H5T_NATIVE_INTEGER, &
            & ila_time, ila_dims, il_err)
         CALL H5Dclose_f(il_dataset_id, il_err)
         ! Get min and max index of time selection:
         ila_minidx = MINLOC(ila_time,ila_time>=id_tmin)
         ila_maxidx = MINLOC(ila_time,ila_time>=id_tmax)   ! pb si on utilisait MAXLOC car si tmax est répété dans ila_time
         ! alors MAXLOC retourne le premier indice de la répétition,
         ! or on veux le dernier indice de la répétition.
         IF ( id_tmin.LE.ila_time(ila_dims(1)) .AND. id_tmax.GT.ila_time(ila_dims(1)) ) ila_maxidx = ila_dims(1)+1   ! sinon dans ce cas ila_maxidx=0 ce qui ne faut pas
         DEALLOCATE(ila_time)

         ! Number of selected elements in time
         il_tselsize = ila_maxidx(1) - ila_minidx(1)
         IF (ig_hdfverb.GE.1) &
            & WRITE(*,'(1X,A,I6,A)') '> [I/Oa] Selected ',il_tselsize,' elements in time'
         Get_number_selected_timeelts = il_tselsize

      ELSE

         il_tselsize = 0
         Get_number_selected_timeelts = 0

      END IF

      IF (il_tselsize .NE. 0) THEN
         ! Initialize the dataspace for the selected elements in time.
         ! The H5S API for defining dataset dataspace.
         ! Selects a hyperslab region to add to the current selected region.
         ila_start(1) = ila_minidx(1)-1
         ila_count(1) = il_tselsize
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_start, ila_count, il_err)
         IF (il_err /= 0) WRITE(*,*) 'Select hyperslab error ', il_err

         ! Creates a new simple dataspace and opens it for access
         il_memrank = 1
         ila_memdims(1) = il_tselsize   ! number of selected elements in time
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         IF (il_err /= 0) WRITE(*,*) 'Create mem dataspace error ', il_err

         h5state%memspace_id = il_memspace_id
      END IF

      h5state%dataspace_id = il_dataspace_id
      h5state%totnobs = ila_dims

   END FUNCTION Get_number_selected_timeelts
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Get_number_selected_geoelts(h5state, id_tselsize, rd_latmin, rd_latmax, rd_lonmin, rd_lonmax)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Function to get number of geographically selected elements in latitude-longitude,
      !! and to change memory space HDF5 identifiers
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(INOUT) :: h5state
      INTEGER,         INTENT(IN) :: id_tselsize
      REAL,            INTENT(IN) :: rd_latmin
      REAL,            INTENT(IN) :: rd_latmax
      REAL,            INTENT(IN) :: rd_lonmin
      REAL,            INTENT(IN) :: rd_lonmax

      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      ! Dimensions
      INTEGER(hsize_t) :: il_hselsize
      INTEGER(hsize_t), DIMENSION(:,:), ALLOCATABLE :: ila_hselcoord
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      ! Inquiries
      INTEGER :: il_err
      ! Temporary buffers
      CHARACTER(LEN=ip_hdf_namelen) :: cl_dsname
      REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rla_lat
      REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rla_lon
      INTEGER :: ib
      INTEGER :: il_hselcount

      IF (ig_hdfverb.GE.1) THEN
         WRITE(*,'(1X,A,I6,A)') '> [I/Oa] Restrain obs selection to geographical boundaries with:'
         WRITE(*,'(1X,2(A,F8.1))') '                         Latmin = ',rd_latmin, '   LatMax = ',rd_latmax
         WRITE(*,'(1X,2(A,F8.1))') '                         Lonmin = ',rd_lonmin, '   LonMax = ',rd_lonmax
      END IF

      ALLOCATE(rla_lat(id_tselsize))
      ALLOCATE(rla_lon(id_tselsize))

      ! Read latitude and longitude datasets in the time selection
      cl_dsname = 'GEOLOCALIZATION/Latitude'
      CALL READSLICE_H5DSET_R1(h5state, cl_dsname, rla_lat)
      cl_dsname = 'GEOLOCALIZATION/Longitude'
      CALL READSLICE_H5DSET_R1(h5state, cl_dsname, rla_lon)

      ! Restrain the selection on the horizontal dimension (geographical):
      ! 1-Search the size of the horizontal selection (geographical)
      il_hselsize = 0
      DO ib = 1, id_tselsize
         IF ( REAL(rla_lat(ib),KIND(rd_latmin)) .GE. rd_latmin .AND. &
            & REAL(rla_lat(ib),KIND(rd_latmax)) .LE. rd_latmax .AND. &
            & REAL(rla_lon(ib),KIND(rd_lonmin)) .GE. rd_lonmin .AND. &
            & REAL(rla_lon(ib),KIND(rd_lonmax)) .LE. rd_lonmax ) il_hselsize = il_hselsize + 1
      END DO

      ! 2-Search the coord of the horizontal selection (geographical)
      ALLOCATE(ila_hselcoord(1,il_hselsize))
      il_hselcount = 0
      DO ib = 1, id_tselsize
         IF ( REAL(rla_lat(ib),KIND(rd_latmin)) .GE. rd_latmin .AND. &
            & REAL(rla_lat(ib),KIND(rd_latmax)) .LE. rd_latmax .AND. &
            & REAL(rla_lon(ib),KIND(rd_lonmin)) .GE. rd_lonmin .AND. &
            & REAL(rla_lon(ib),KIND(rd_lonmax)) .LE. rd_lonmax ) THEN

            !WRITE(*,*) ' Lat, lon ',rla_lat(ib), rla_lon(ib)
            il_hselcount = il_hselcount + 1
            ila_hselcoord(1,il_hselcount) = ila_start(1)+ib
         END IF
      END DO

      DEALLOCATE(rla_lat)
      DEALLOCATE(rla_lon)

      IF (ig_hdfverb.GE.1) &
         & WRITE(*,'(1X,A,I6,A)') '> [I/Oa] Refined selection to ', il_hselsize, ' elements in time and space'
      !WRITE(*,*) ila_hselcoord(1,:)
      Get_number_selected_geoelts = il_hselsize

      IF (il_hselsize .NE. 0) THEN
         ! Select array elements to be included in the selection of the dataspace
         il_memrank = 1
         CALL H5Sselect_elements_f(h5state%dataspace_id, H5S_SELECT_SET_F, il_memrank, il_hselsize, ila_hselcoord, il_err)
         IF (il_err /= 0) WRITE(*,*) 'Select elements error ', il_err
         DEALLOCATE(ila_hselcoord)

         ! Creates a new simple dataspace and opens it for access
         il_memspace_id = h5state%memspace_id
         CALL H5Sclose_f(il_memspace_id, il_err)
         ila_memdims(1) = il_hselsize
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         IF (il_err /= 0) WRITE(*,*) 'Create mem dataspace error ', il_err

         h5state%memspace_id = il_memspace_id
      END IF

   END FUNCTION Get_number_selected_geoelts
   !---------------------------------------------------------------------------
END MODULE H5_SELECTION_MOD
