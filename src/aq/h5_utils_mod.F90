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

module H5_UTILS_MOD
   !! category: observation operator
   !! author:   CERFACS and CNRM (G. Jonville)
   !!
   !! ----------------------------------------------------------------------
   !!
   !! Tools to manipulate and check HDF5 files, groups, spaces and attributes
   !!
   !! ----------------------------------------------------------------------
   !!
   USE HDF5
   USE ISO_C_BINDING

   IMPLICIT NONE

   INTEGER, PARAMETER :: ip_hid_t = hid_t !! kind for integer types in datdat

   ! This type contains the identifiers to give access to an hdf5 file by hdf5 functions
   TYPE h5state_t
      INTEGER(hid_t) :: instr_id
      INTEGER(hid_t) :: dataspace_id
      INTEGER(hid_t) :: memspace_id
      INTEGER(hsize_t), DIMENSION(1) :: totnobs
   END TYPE h5state_t

   INTEGER, PARAMETER :: ip_hdf_namelen=128
   INTEGER, PARAMETER :: ip_missing_obs = -2                   ! The missing observation mask value
   REAL(KIND=8), PARAMETER :: rp_hdf5_missing_value = -1.d+38  ! The missing real value in hdf5 files
   REAL(KIND=8), PARAMETER :: rp_roundup_factor = 0.99d0       ! The factor to set an upper bound of rp_hdf5_missing_value
   INTEGER, PARAMETER      :: ip_hdf5_missing_value = -huge(1) ! The missing integer value in hdf5 files
   INTEGER :: ig_hdfverb = 0
   INTEGER :: ig_recursion
   CHARACTER(LEN=ip_hdf_namelen), DIMENSION(:), ALLOCATABLE :: cga_names_in_level

CONTAINS
   !---------------------------------------------------------------------------
   ! HDF5 File
   !---------------------------------------------------------------------------
   SUBROUTINE CREATE_H5FILE(cd_fname, id_file_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Create hdf5 file
      !!
      !! ----------------------------------------------------------------------
      !!
      CHARACTER(LEN=*), INTENT(IN) :: cd_fname
      INTEGER(hid_t), INTENT(OUT) :: id_file_id

      INTEGER :: il_err

      IF (ig_hdfverb.GE.2) THEN
         WRITE(*,'(A)') '>> [I/Oa] "INPUT / OUTPUT for assimilation"'
         ! Activate HDF5 fortran
         WRITE(*,'(A)') ' > [I/Oa] Initializing HDF5 library'
      END IF
      CALL H5open_f (il_err)

      ! Create hdf5 file
      IF (ig_hdfverb.GE.1) &
         & WRITE(*,'(2A)') ' > [I/Oa] Creating HDF5 output file ',TRIM(cd_fname)
      CALL H5Fcreate_f(cd_fname, H5F_ACC_TRUNC_F, &
         & id_file_id, il_err)
      IF (il_err /=0) WRITE(*,'(A,I3)') '!!!!! Create HDF5 file error ',il_err

   END SUBROUTINE CREATE_H5FILE
   !---------------------------------------------------------------------------
   SUBROUTINE OPEN_H5FILE_RDONLY(cd_fname, id_file_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Open hdf5 file read only
      !!
      !! ----------------------------------------------------------------------
      !!
      CHARACTER(LEN=*), INTENT(IN) :: cd_fname
      INTEGER(hid_t), INTENT(OUT) :: id_file_id

      INTEGER :: il_err

      ! Open HDF5 file
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(2A)') ' > [I/Oa] Opening HDF5 file ',TRIM(cd_fname)
      CALL H5Fopen_f(cd_fname, H5F_ACC_RDONLY_F, &
         & id_file_id, il_err)
      IF (il_err /=0) THEN
         WRITE(*,'(A,I3)') '!!!!! Open HDF5 file error ',il_err
         ! CALL MOCAGE_STOP
      END IF
   END SUBROUTINE OPEN_H5FILE_RDONLY
   !---------------------------------------------------------------------------
   SUBROUTINE CLOSE_H5FILE(id_file_id, ld_closeh5)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Close hdf5 file
      !!
      !! ----------------------------------------------------------------------
      !!
      LOGICAL, INTENT(IN) :: ld_closeh5
      INTEGER(hid_t), INTENT(IN) :: id_file_id

      CHARACTER(LEN=ip_hdf_namelen) :: cl_fname
      INTEGER(size_t) :: il_size
      INTEGER :: il_err

      CALL H5Fget_name_f(id_file_id, cl_fname, il_size, il_err)
      ! Close file
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,2A,/)') '> [I/Oa] Closing HDF5 file ',TRIM(cl_fname)
      CALL H5Fclose_f(id_file_id, il_err)

      ! Shutdown HDF5
      IF (ld_closeh5) THEN
         IF (ig_hdfverb.GE.2) WRITE(*,'(A)') '>> [I/Oa] Closing HDF5'
         CALL H5close_f(il_err)
      END IF

   END SUBROUTINE CLOSE_H5FILE
   !---------------------------------------------------------------------------
   SUBROUTINE CHECK_H5FILE(id_file_id, cd_instrname, cd_instrtype)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Check contents of hdf5 file
      !!
      !! ----------------------------------------------------------------------
      !!
      ! USE SYSTEM_UTILS_MOD, ONLY : DAIMON_Abort

      INTEGER(hid_t), INTENT(IN) :: id_file_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_instrname
      CHARACTER(LEN=*), INTENT(IN) :: cd_instrtype

      ! Identifiers
      INTEGER(hid_t) :: il_group_id
      ! Inquiries
      INTEGER :: il_err
      ! Temporary buffers
      CHARACTER(LEN=ip_hdf_namelen) :: cl_dsname
      CHARACTER(LEN=ip_hdf_namelen) :: cl_attribut

      cl_dsname='/'//TRIM(cd_instrname)
      CALL OPEN_H5GROUP(id_file_id, cl_dsname, il_group_id)

      ! Read measurement type attribute and check it with namelist
      CALL GET_INSTR_ATTR(il_group_id, '.', 'MeasurementType', cl_attribut)
      IF ( TRIM(cl_attribut) /= TRIM(cd_instrtype) ) THEN
         WRITE(*,'(1X,7A)') '/!\ For instrument ',TRIM(cd_instrname),' measurement type in the namelist (',TRIM(cd_instrtype), &
            &') is not the same than in the hdf5 file (',TRIM(cl_attribut),')'
         WRITE(*,*) '/!\ Check your namelist and/or your HDAT*.h5 file'
         WRITE(*,*) 'ABORT'
         ! CALL DAIMON_Abort
      END IF

      ! Check if groups exist
      CALL CHECK_H5GROUP(il_group_id, 'GEOLOCALIZATION')
      CALL CHECK_H5GROUP(il_group_id, 'GEOLOCALIZATION/Latitude')
      CALL CHECK_H5GROUP(il_group_id, 'GEOLOCALIZATION/Longitude')
      CALL CHECK_H5GROUP(il_group_id, 'GEOLOCALIZATION/Timestamp')
      CALL CHECK_H5GROUP(il_group_id, 'OBSERVATIONS')

      CALL CLOSE_H5GROUP(il_group_id)

   END SUBROUTINE CHECK_H5FILE
   !---------------------------------------------------------------------------
   ! HDF5 Group
   !---------------------------------------------------------------------------
   SUBROUTINE OPEN_H5GROUP(id_file_id, cd_grpname, id_grp_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Open hdf5 group
      !!
      !! ----------------------------------------------------------------------
      !!
      ! USE SYSTEM_UTILS_MOD, ONLY : DAIMON_Abort

      INTEGER(hid_t), INTENT(IN)      :: id_file_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_grpname
      INTEGER(hid_t), INTENT(OUT)  :: id_grp_id

      INTEGER :: il_err

      CALL H5Gopen_f(id_file_id, cd_grpname, &
         & id_grp_id, il_err)
      IF (il_err /= 0) THEN
         WRITE(*,'(/,1X,3A,I3,A)') '/!\ Group ',TRIM(cd_grpname),' opening error ',il_err,' in hdf5 file'
         WRITE(*,*) '/!\ Check your namelist and/or your HDAT*.h5 file'
         WRITE(*,*) 'ABORT'
         ! CALL DAIMON_Abort
      END IF

   END SUBROUTINE OPEN_H5GROUP
   !---------------------------------------------------------------------------
   SUBROUTINE CLOSE_H5GROUP(id_grp_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Close hdf5 group
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN)  :: id_grp_id

      INTEGER :: il_err

      CALL H5Gclose_f(id_grp_id, il_err)
      IF (il_err /= 0) WRITE(*,*) 'Group close error ',il_err

   END SUBROUTINE CLOSE_H5GROUP
   !---------------------------------------------------------------------------
   SUBROUTINE CHECK_H5GROUP(id_file_id, cd_grpname)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Check hdf5 group
      !!
      !! ----------------------------------------------------------------------
      !!
      ! USE SYSTEM_UTILS_MOD, ONLY : DAIMON_Abort

      INTEGER(hid_t), INTENT(IN)      :: id_file_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_grpname

      LOGICAL :: ll_exists
      INTEGER :: il_err

      CALL H5Lexists_f(id_file_id, cd_grpname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         WRITE(*,'(1X,4A)') '/!\ For this instrument, group ',TRIM(cd_grpname),' does not exist!'
         WRITE(*,*) '/!\ Check your HDAT*.h5 file'
         WRITE(*,*) 'ABORT'
         ! CALL DAIMON_Abort
      END IF

   END SUBROUTINE CHECK_H5GROUP
   !---------------------------------------------------------------------------
   SUBROUTINE CREATE_H5GROUP(id_loc_id, cd_grpname)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Create hdf5 group
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN) :: id_loc_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_grpname

      INTEGER(hid_t) :: il_group_id
      LOGICAL :: ll_exists
      INTEGER :: il_err

      ! Check if group already exists
      CALL H5Lexists_f(id_loc_id, cd_grpname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create group
         IF (ig_hdfverb.GE.2) &
            & WRITE(*,'(2A)') ' > [I/Oa] Creating HDF5 group ',TRIM(cd_grpname)
         CALL H5Gcreate_f(id_loc_id, cd_grpname, &
            & il_group_id, il_err)
         IF (il_err /=0) WRITE(*,'(/,1X,3A,I3,A)') '/!\ Create group ',TRIM(cd_grpname),' error ',il_err,' in hdf5 file'
         CALL H5Gclose_f(il_group_id, il_err)
      END IF

   END SUBROUTINE CREATE_H5GROUP
   !---------------------------------------------------------------------------
   SUBROUTINE CREATE_H5GROUP_INSTR(id_file_id, cd_grpname, cd_attrib_meastype)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Create hdf5 instrument group
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN) :: id_file_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_grpname
      CHARACTER(LEN=*), INTENT(IN) :: cd_attrib_meastype

      INTEGER(hid_t) :: il_instr_id
      INTEGER :: il_err

      ! Create group
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(2A)') ' > [I/Oa] Creating HDF5 group ',TRIM(cd_grpname)
      CALL H5Gcreate_f(id_file_id, cd_grpname, &
         & il_instr_id, il_err)
      IF (il_err /=0) WRITE(*,*) 'Create group error ',il_err

      ! Write instrument attribut
      CALL CREATE_ATTRIB_STRING(il_instr_id, 'MeasurementType', cd_attrib_meastype)

   END SUBROUTINE CREATE_H5GROUP_INSTR
   !---------------------------------------------------------------------------
   SUBROUTINE CREATE_H5GROUP_INSTRDOM(id_file_id, cd_grpname, id_instr_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Create hdf5 domain group and subgroups
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN) :: id_file_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_grpname
      INTEGER(hid_t), INTENT(OUT) :: id_instr_id

      INTEGER :: il_err

      ! Create group
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(2A)') ' > [I/Oa] Creating HDF5 group ',TRIM(cd_grpname)
      CALL H5Gcreate_f(id_file_id, cd_grpname, &
         & id_instr_id, il_err)
      IF (il_err /=0) WRITE(*,*) 'Create group error ',il_err

      ! Create instrument/domain  sub-groups
      CALL CREATE_H5GROUP(id_instr_id, 'GEOLOCALIZATION')
      CALL CREATE_H5GROUP(id_instr_id, 'OBSERVATIONS')

   END SUBROUTINE CREATE_H5GROUP_INSTRDOM
   !---------------------------------------------------------------------------
   ! HDF5 Space
   !---------------------------------------------------------------------------
   SUBROUTINE CLOSE_H5SPACE(id_space_id)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Close hdf5 space
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN)  :: id_space_id

      INTEGER :: il_err

      CALL H5Sclose_f(id_space_id, il_err)
      IF (il_err /= 0) WRITE(*,*) 'Space close error ',il_err

   END SUBROUTINE CLOSE_H5SPACE
   !---------------------------------------------------------------------------
   ! HDF5 Dataset
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Get_h5dset_size(h5state, cd_dsname)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Get the size of a dataset
      !!
      !! ----------------------------------------------------------------------
      !!
      USE H5LT
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER(hsize_t), DIMENSION(1) :: ila_dims
      ! Inquiries
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
         & ila_dims, il_type_class, il_type_size, il_err)

      Get_h5dset_size = ila_dims(1)

   END FUNCTION Get_h5dset_size
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Get_h5data_size(h5state, cd_dsname)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Get the size of a data
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE :: ila_shape
      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_rank
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
      ! Get the rank of the datatype
      CALL H5Tget_array_ndims_f(il_datatype_id, il_rank, il_err)
      ! Get the dimensions of the datatype
      ALLOCATE(ila_shape(il_rank))
      CALL H5Tget_array_dims_f(il_datatype_id, ila_shape, il_err)
      Get_h5data_size = PRODUCT(ila_shape)

      DEALLOCATE(ila_shape)
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

   END FUNCTION Get_h5data_size
   !---------------------------------------------------------------------------
   ! HDF5 Attribut
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Get_int_attribute(id_instr_id, cd_objname, cd_attrname)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read an integer attribute
      !!
      !! ----------------------------------------------------------------------
      !!
      USE H5LT

      INTEGER(HID_T) :: id_instr_id
      CHARACTER(LEN=*) :: cd_objname
      CHARACTER(LEN=*) :: cd_attrname

      INTEGER, DIMENSION(1) :: ila_int_attr
      INTEGER :: il_err

      CALL H5LTget_attribute_int_f(id_instr_id, cd_objname, cd_attrname, ila_int_attr, il_err)
      IF (il_err /= 0) WRITE(*,*) 'Get attribute integer error ',il_err

      Get_int_attribute = ila_int_attr(1)

   END FUNCTION Get_int_attribute
   !---------------------------------------------------------------------------
   SUBROUTINE GET_INSTR_ATTR(id_instr_id, cd_objname, cd_attrname, cd_string)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read a string attribute.
      !! Equivalent of H5LTget_attribute_string_f to avoid extra char reading.
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t) :: id_instr_id
      CHARACTER(LEN=*) :: cd_objname
      CHARACTER(LEN=*) :: cd_attrname
      CHARACTER(LEN=1), DIMENSION(ip_hdf_namelen), TARGET :: cd_string
      INTEGER(hid_t) :: il_instr_attr_id
      INTEGER(hid_t) :: il_attr_datatype
      TYPE(C_PTR) :: ptr
      INTEGER :: il_err

      CALL H5Aopen_by_name_f(id_instr_id, TRIM(cd_objname), TRIM(cd_attrname), &
         & il_instr_attr_id, il_err)
      CALL H5Aget_type_f(il_instr_attr_id, &
         & il_attr_datatype, il_err)
      cd_string(:)=''
      ptr=C_LOC(cd_string)
      CALL H5Aread_f(il_instr_attr_id, il_attr_datatype, &
         & ptr, il_err)
      CALL H5Aclose_f(il_instr_attr_id, il_err)

   END SUBROUTINE GET_INSTR_ATTR
   !---------------------------------------------------------------------------
   SUBROUTINE CREATE_ATTRIB_STRING(id_loc_id, cd_attribname, cd_attribdata)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Create string attribute
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), INTENT(IN) :: id_loc_id
      CHARACTER(LEN=*), INTENT(IN) :: cd_attribname
      CHARACTER(LEN=*), INTENT(IN) :: cd_attribdata

      INTEGER(size_t) :: il_attriblen
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      INTEGER(hid_t) :: il_attrib_id
      INTEGER(hsize_t), DIMENSION(1) :: ila_bufdims
      INTEGER :: il_err

      il_attriblen = LEN(TRIM(cd_attribdata))

      ! Create scalar data space (H5S_SCALAR_F=single element) for the attribute
      ! Create array data space (H5S_SIMPLE_F=regular array of elements) for the attribute
      CALL H5Screate_f(H5S_SCALAR_F, il_dataspace_id, il_err)
      !autre forme :  CALL h5screate_simple_f(1, (/1/), aspace_id, hdferr)
      ! Create datatype for the attribute
      CALL H5Tcopy_f(H5T_NATIVE_CHARACTER, il_datatype_id, il_err)
      ! set the size of the datatype to be the size of the attribut string
      CALL H5Tset_size_f(il_datatype_id, il_attriblen, il_err)
      ! set string as NULLTERM to imitate C-style strings
      CALL H5Tset_strpad_f(il_datatype_id, H5T_STR_NULLTERM_F, il_err)
      !Create attribute
      CALL H5Acreate_f(id_loc_id, TRIM(cd_attribname), il_datatype_id, il_dataspace_id, &
         & il_attrib_id, il_err)
      ! Write the attribute data
      ila_bufdims(1) = 1      ! DIMENSION for scalars
      CALL H5Awrite_f(il_attrib_id, il_datatype_id, TRIM(cd_attribdata), ila_bufdims, il_err)
      ! Close the dataspace
      CALL H5Sclose_f(il_dataspace_id, il_err)
      ! Close the attribute
      CALL H5Aclose_f(il_attrib_id, il_err)

   END SUBROUTINE CREATE_ATTRIB_STRING
   !---------------------------------------------------------------------------
   INTEGER FUNCTION Op_func(loc_id, name, info, cptr) BIND(C)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Operator function
      !!
      !! ----------------------------------------------------------------------
      !!
      INTEGER(hid_t), VALUE :: loc_id
      CHARACTER(LEN=1), DIMENSION(1:ip_hdf_namelen) :: name ! We must have LEN=1 for bind(C) strings
      ! in order to be standard compliant
      TYPE(C_PTR) :: info
      CHARACTER(LEN=ip_hdf_namelen) :: name_string = ' '
      TYPE(C_PTR) :: cptr
      INTEGER      :: i
      TYPE(h5o_info_t), TARGET :: infobuf
      INTEGER :: status

      name_string = " "
      DO i = 1, ip_hdf_namelen
         IF (name(i)(1:1) .EQ. C_NULL_CHAR) EXIT ! Read up to the C NULL termination
         name_string(i:i) = name(i)(1:1)
      END DO

      CALL H5Oget_info_by_name_f(loc_id, name_string, infobuf, status)

      IF (name(1)(1:1) .EQ. '.') THEN            ! Root group, do not print '.'
         WRITE(*,"('  (Group)')")
      ELSE
         IF (infobuf%TYPE .EQ. H5O_type_GROUP_F) THEN
            WRITE(*,'(" *)",A,"  (Group)")') TRIM(name_string)
         ELSE IF (infobuf%TYPE .EQ. H5O_type_DATASET_F) THEN
            WRITE(*,'(" *)",A,"  (Dataset)")') TRIM(name_string)
         ELSE IF (infobuf%TYPE .EQ. H5O_type_NAMED_DATAtype_F) THEN
            WRITE(*,'(" *)",A,"  (Datatype)")') TRIM(name_string)
         ELSE
            WRITE(*,'(" *)",A,"  (Unknown)")') TRIM(name_string)
         END IF
      END IF

      ig_recursion = ig_recursion+1
      cga_names_in_level(ig_recursion) = TRIM(name_string)

      Op_func = 0 ! return successful

   END FUNCTION Op_func
   !---------------------------------------------------------------------------
END MODULE H5_UTILS_MOD
