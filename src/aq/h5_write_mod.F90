! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!! category: observation operator
!! summary:  Write a slice of HDF5 dataset
!! author:   CERFACS and CNRM (G. Jonville)
!!
MODULE H5_WRITE_MOD
   !! category: observation operator
   !! author:   CERFACS and CNRM (G. Jonville)
   !!
   !! ----------------------------------------------------------------------
   !!
   !! Interface to write a slice of HDF5 dataset
   !!
   !! ----------------------------------------------------------------------
   !!
   USE H5_UTILS_MOD
   USE H5LT

   IMPLICIT NONE

   INTERFACE WRITESLICE_H5DSET_SCALAR
      MODULE PROCEDURE WRITESLICE_H5DSET_C0
      MODULE PROCEDURE WRITESLICE_H5DSET_I0
      MODULE PROCEDURE WRITESLICE_H5DSET_R0
   END INTERFACE WRITESLICE_H5DSET_SCALAR

   INTERFACE WRITESLICE_H5DSET
      MODULE PROCEDURE WRITESLICE_H5DSET_C1
      MODULE PROCEDURE WRITESLICE_H5DSET_I1
      MODULE PROCEDURE WRITESLICE_H5DSET_I2
      MODULE PROCEDURE WRITESLICE_H5DSET_R1
      MODULE PROCEDURE WRITESLICE_H5DSET_R2
      MODULE PROCEDURE WRITESLICE_H5DSET_R3
   END INTERFACE WRITESLICE_H5DSET

   INTERFACE WRITESLICE_H5DSET_CONSTANT
      MODULE PROCEDURE WRITESLICE_H5DSET_CR2
   END INTERFACE WRITESLICE_H5DSET_CONSTANT

CONTAINS
   !---------------------------------------------------------------------------
   ! HDF5 Dataset
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_C0(h5state, cd_dsname, cda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of string data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: cda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size
      INTEGER(size_t) :: il_slength

      il_memrank = 1
      ila_memdims(1) = SIZE(cda_dset, dim=1)
      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! set size of fortran sring
      il_slength = LEN(TRIM(cda_dset(1)))
      CALL H5Tcreate_f(H5T_STRING_F, il_slength, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! string scalar datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, cda_dset, ila_memdims, il_err, & ! string scalar datatype
         & il_memspace_id, il_dataspace_id)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,3A,I5,A,T92,5A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',TRIM(cda_dset(1)),'  Max= ',TRIM(cda_dset(SIZE(cda_dset))),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_C0
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_I0(h5state, cd_dsname, ida_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of integer data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER, DIMENSION(:), INTENT(IN) :: ida_dset

      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(ida_dset, dim=1)
      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, H5T_STD_I32LE, il_memspace_id, & ! integer scalar datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, H5T_NATIVE_INTEGER, ida_dset, ila_memdims, il_err, & ! integer scalar datatype
         & il_memspace_id, il_dataspace_id)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,I12,A,I12,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(ida_dset),'  Max= ',MAXVAL(ida_dset),')'

      ! Close objects
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_I0
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_R0(h5state, cd_dsname, rda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of real data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(rda_dset, dim=1)

      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, H5T_IEEE_F32LE, il_memspace_id, & ! float scalar datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, H5T_NATIVE_REAL, REAL(rda_dset,kind=4), ila_memdims, il_err, & ! float scalar datatype
         & il_memspace_id, il_dataspace_id)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

      ! Close objects
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_R0
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_C1(h5state, cd_dsname, cda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of string array(1) data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: cda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size
      INTEGER(size_t) :: il_slength

      il_memrank = 1
      ila_memdims(1) = SIZE(cda_dset, dim=1)
      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype of rank 1 and size 1 for dataset
      il_datarank = 1
      ila_datadims(1) = 1
      CALL H5Tarray_create_f(H5T_FORTRAN_S1, il_datarank, ila_datadims, il_datatype_id, il_err)
      ! set size of fortran sring
      il_slength = LEN(TRIM(cda_dset(1)))
      CALL H5Tset_size_f(il_datatype_id, il_slength, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! string array(1) datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, cda_dset, ila_memdims, il_err, & ! string array(1) datatype
         & il_memspace_id, il_dataspace_id)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,3A,I5,A,T92,5A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',TRIM(cda_dset(1)),'  Max= ',TRIM(cda_dset(SIZE(cda_dset))),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_C1
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_I1(h5state, cd_dsname, ida_dset, ida_mask)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of integer array(1) data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER, DIMENSION(:), INTENT(IN) :: ida_dset
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ida_mask

      INTEGER, DIMENSION(:), ALLOCATABLE :: ila_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(ida_dset, dim=1)
      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype of rank 1 and size 1 for dataset
      il_datarank = 1
      ila_datadims(1) = 1
      CALL H5Tarray_create_f(H5T_STD_I32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! integer array(1) datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ALLOCATE(ila_dset, SOURCE=ida_dset) ! allocate to the same dimensions and copy elements of source to target
      ! Rewrite the missing value in output file
      IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) ila_dset=ip_hdf5_missing_value
      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, ila_dset, ila_memdims, il_err, & ! integer array(1) datatype
         & il_memspace_id, il_dataspace_id)
      DEALLOCATE(ila_dset)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,I12,A,I12,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(ida_dset),'  Max= ',MAXVAL(ida_dset),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_I1
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_I2(h5state, cd_dsname, ida_dset, ida_mask)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of integer array(n) data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ida_dset
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ida_mask

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ila_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(ida_dset, dim=2)
      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype for dataset
      il_datarank = 1
      ila_datadims(1) = SIZE(ida_dset, dim=1)
      CALL H5Tarray_create_f(H5T_STD_I32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, &
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ALLOCATE(ila_dset, SOURCE=ida_dset) ! allocate to the same dimensions and copy elements of source to target
      ! Rewrite the missing value in output file
      IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) ila_dset=ip_hdf5_missing_value
      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, ila_dset, ila_memdims, il_err, &
         & il_memspace_id, il_dataspace_id)
      DEALLOCATE(ila_dset)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,I12,A,I12,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(ida_dset),'  Max= ',MAXVAL(ida_dset),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_I2
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_R1(h5state, cd_dsname, rda_dset, ida_mask)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of real array(1) data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rda_dset
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ida_mask

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rla_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(rda_dset, dim=1)

      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype of rank 1 and size 1 for dataset
      il_datarank = 1
      ila_datadims(1) = 1
      CALL H5Tarray_create_f(H5T_IEEE_F32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! float array(1) datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ALLOCATE(rla_dset, SOURCE=rda_dset) ! allocate to the same dimensions and copy elements of source to target
      ! Rewrite the missing value in output file
      IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) rla_dset=rp_hdf5_missing_value
      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, REAL(rla_dset, KIND=4), ila_memdims, il_err, & ! float array(1) datatype
         & il_memspace_id, il_dataspace_id)
      DEALLOCATE(rla_dset)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_R1
   !---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_R2(h5state, cd_dsname, rda_dset, ida_mask)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of real array(n) data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rda_dset
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ida_mask

      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rla_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(rda_dset, dim=2)

      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype for dataset
      il_datarank = 1
      ila_datadims(1) = SIZE(rda_dset, dim=1)
      CALL H5Tarray_create_f(H5T_IEEE_F32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! float array(n) datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
            & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ALLOCATE(rla_dset, SOURCE=rda_dset) ! allocate to the same dimensions and copy elements of source to target
      ! Rewrite the missing value in output file
      IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) rla_dset=rp_hdf5_missing_value
      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, REAL(rla_dset, KIND=4), ila_memdims, il_err, & ! float array(n) datatype
         & il_memspace_id, il_dataspace_id)
      DEALLOCATE(rla_dset)
      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,A,I5,A,T92,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
         & ' of ',ila_memdims,' elements in HSTAT*.h5', &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_R2
!---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_R3(h5state, cd_dsname, rda_dset, ida_mask)
   !! category: observation operator
   !! author:   CERFACS and CNRM (E. Emili)
   !!
   !! ----------------------------------------------------------------------
   !!
   !! Write HDF5 dataset slice of real 3d-array data type
   !!
   !! ----------------------------------------------------------------------
   !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: rda_dset
      INTEGER, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: ida_mask

      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: rla_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(2) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(rda_dset, dim=3)

      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
                       & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype for dataset
      il_datarank = 2
      ila_datadims(1) = SIZE(rda_dset, dim=1)
      ila_datadims(2) = SIZE(rda_dset, dim=2)
      CALL H5Tarray_create_f(H5T_IEEE_F32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
                                 & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, &
                          & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

      ELSE
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         ! Get current dimensions of dataset (rank is 1 in the hdf5 file)
         CALL H5LTget_dataset_info_f(h5state%instr_id, cd_dsname, &
                                     & ila_datasetdims, il_type_class, il_type_size, il_err)

         ila_begin(1) = ila_datasetdims(1)
         ila_count(1) = ila_memdims(1)
         ! Extend the dataset
         CALL H5Dset_extent_f(il_dataset_id, ila_begin(1:1)+ila_count(1:1), il_err)
         ! Create simple dataspace
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
                                 & il_memspace_id, il_err)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)
         ! Selects a hyperslab region to add to the current selected region
         CALL H5Sselect_hyperslab_f(il_dataspace_id, H5S_SELECT_SET_F, ila_begin, ila_count, il_err)

      END IF

      ALLOCATE(rla_dset, SOURCE=rda_dset) ! allocate to the same dimensions and copy elements of source to target
      ! Rewrite the missing value in output file
      IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) rla_dset=rp_hdf5_missing_value
      ! Write data to dataset / Write extended part to dataset
      CALL H5Dwrite_f(il_dataset_id, il_datatype_id, REAL(rla_dset, KIND=4), ila_memdims, il_err, &
                      & il_memspace_id, il_dataspace_id)
      DEALLOCATE(rla_dset)
      IF (ig_hdfverb.GE.2) &
      & WRITE(*,'(1X,A,A,A,I5,A,T92,A,E12.4,A,E12.4,A)') &
            & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
            & ' of ',ila_memdims,' elements in HSTAT*.h5', &
            & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

      ! Close objects
      CALL H5Tclose_f(il_datatype_id, il_err)
      CALL H5Sclose_f(il_dataspace_id, il_err)
      CALL H5Sclose_f(il_memspace_id, il_err)
      CALL H5Dclose_f(il_dataset_id, il_err)
      CALL H5Pclose_f(il_prop_id, il_err)

   END SUBROUTINE WRITESLICE_H5DSET_R3
!---------------------------------------------------------------------------
   SUBROUTINE WRITESLICE_H5DSET_CR2(h5state, cd_dsname, rda_dset, ida_mask)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Write HDF5 dataset slice of real array(n) data type only once on multiple calls
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rda_dset
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ida_mask

      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rla_dset
      ! Identifiers
      INTEGER(hid_t) :: il_memspace_id
      INTEGER(hid_t) :: il_prop_id
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_dataspace_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_memrank
      INTEGER(hsize_t), DIMENSION(1) :: ila_memdims
      INTEGER(hsize_t), DIMENSION(1) :: ila_begin
      INTEGER(hsize_t), DIMENSION(1) :: ila_count
      INTEGER :: il_datarank
      INTEGER(hsize_t), DIMENSION(1) :: ila_datadims
      INTEGER(hsize_t), DIMENSION(1) :: ila_datasetdims
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err
      ! Temporary buffers
      INTEGER(hsize_t), DIMENSION(1) :: ila_maxdims
      INTEGER :: il_type_class
      INTEGER(size_t) :: il_type_size

      il_memrank = 1
      ila_memdims(1) = SIZE(rda_dset, dim=2)

      ! Modify dataset creation properties, i.e. enable chunking
      CALL H5Pcreate_f(H5P_DATASET_CREATE_F, &
         & il_prop_id, il_err)
      CALL H5Pset_chunk_f(il_prop_id, il_memrank, ila_memdims, il_err)

      ! Create an array datatype for dataset
      il_datarank = 1
      ila_datadims(1) = SIZE(rda_dset, dim=1)
      CALL H5Tarray_create_f(H5T_IEEE_F32LE, il_datarank, ila_datadims, il_datatype_id, il_err)

      ! Check if dataset already exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (.NOT.ll_exists) THEN
         ! Create the data space with unlimited dimensions
         ila_maxdims = (/H5S_UNLIMITED_f/)
         CALL H5Screate_simple_f(il_memrank, ila_memdims, &
            & il_memspace_id, il_err, ila_maxdims)
         ! Create the dataset using prop_id creation properties
         CALL H5Dcreate_f(h5state%instr_id, cd_dsname, il_datatype_id, il_memspace_id, & ! float array(n) datatype
            & il_dataset_id, il_err, il_prop_id)
         ! Get dataset space
         CALL H5Dget_space_f(il_dataset_id, il_dataspace_id, il_err)

         ALLOCATE(rla_dset, SOURCE=rda_dset) ! allocate to the same dimensions and copy elements of source to target
         ! Rewrite the missing value in output file
         IF (PRESENT(ida_mask)) WHERE (ida_mask.EQ.ip_missing_obs) rla_dset=rp_hdf5_missing_value
         ! Write data to dataset / Write extended part to dataset
         CALL H5Dwrite_f(il_dataset_id, il_datatype_id, REAL(rla_dset, KIND=4), ila_memdims, il_err, & ! float array(n) datatype
            & il_memspace_id, il_dataspace_id)
         DEALLOCATE(rla_dset)
         IF (ig_hdfverb.GE.2) &
            & WRITE(*,'(1X,A,A,A,I5,A,T92,A,E12.4,A,E12.4,A)') &
            & '> [I/Oa] Writing/Extending dataset ',TRIM(cd_dsname), &
            & ' of ',ila_memdims,' elements in HSTAT*.h5', &
            & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

         ! Close objects
         CALL H5Tclose_f(il_datatype_id, il_err)
         CALL H5Sclose_f(il_dataspace_id, il_err)
         CALL H5Sclose_f(il_memspace_id, il_err)
         CALL H5Dclose_f(il_dataset_id, il_err)
         CALL H5Pclose_f(il_prop_id, il_err)
      END IF

   END SUBROUTINE WRITESLICE_H5DSET_CR2
!---------------------------------------------------------------------------
END MODULE H5_WRITE_MOD
