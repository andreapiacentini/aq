!! category: observation operator
!! summary:  Read a slice of HDF5 dataset
!! author:   CERFACS and CNRM (G. Jonville)
!!
MODULE H5_READ_MOD
   !! category: observation operator
   !! author:   CERFACS and CNRM (G. Jonville)
   !!
   !! ----------------------------------------------------------------------
   !!
   !! Interface to read a slice of HDF5 dataset
   !!
   !! ----------------------------------------------------------------------
   !!
   USE H5_UTILS_MOD

   IMPLICIT NONE

   INTERFACE READSIMPLE_H5DSET
      MODULE PROCEDURE READSIMPLE_H5DSET_R1
      MODULE PROCEDURE READSIMPLE_H5DSET_I1
   END INTERFACE READSIMPLE_H5DSET

   INTERFACE READSLICE_H5DSET
      MODULE PROCEDURE READSLICE_H5DSET_C1
      MODULE PROCEDURE READSLICE_H5DSET_I1
      MODULE PROCEDURE READSLICE_H5DSET_R1
      MODULE PROCEDURE READSLICE_H5DSET_R2
      MODULE PROCEDURE READSLICE_H5DSET_R3
   END INTERFACE READSLICE_H5DSET

CONTAINS
   !---------------------------------------------------------------------------
   ! HDF5 Dataset
   !---------------------------------------------------------------------------
   SUBROUTINE READSIMPLE_H5DSET_R1(h5state, cd_dsname, rda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read the first data in the dataset (used to read constant pressure/wavelength array)
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: rda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER(hsize_t), DIMENSION(1) :: ila_dims
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
      ! Read first data in the dataset
      ila_dims(1)=1
      CALL H5Dread_f(il_dataset_id, il_datatype_id, rda_dset, ila_dims, &
         & il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Reading one data of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

   END SUBROUTINE READSIMPLE_H5DSET_R1
   !---------------------------------------------------------------------------
   SUBROUTINE READSIMPLE_H5DSET_I1(h5state, cd_dsname, ida_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read the first data in the dataset (used to read constant satellite channels array)
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER(KIND=4), DIMENSION(:), INTENT(OUT) :: ida_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER(hsize_t), DIMENSION(1) :: ila_dims
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
      ! Read first data in the dataset
      ila_dims(1)=1
      CALL H5Dread_f(il_dataset_id, il_datatype_id, ida_dset, ila_dims, &
         & il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Reading one data of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(ida_dset),'  Max= ',MAXVAL(ida_dset),')'

   END SUBROUTINE READSIMPLE_H5DSET_I1
   !---------------------------------------------------------------------------
   SUBROUTINE READSLICE_H5DSET_C1(h5state, cd_dsname, cda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read HDF5 dataset slice of string 1d-array data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: cda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Inquiries
      LOGICAL :: ll_exists
      INTEGER :: il_err

      ! Check if dataset exists
      CALL H5Lexists_f(h5state%instr_id, cd_dsname, ll_exists, il_err)

      IF (ll_exists) THEN
         ! Open dataset
         CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
         IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
         ! Get dataset type
         CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
         ! Read dataset
         CALL H5Dread_f(il_dataset_id, il_datatype_id, cda_dset, h5state%totnobs, &
            & il_err, &
            & h5state%memspace_id, h5state%dataspace_id)
         IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
         ! Close dataset
         CALL H5Dclose_f(il_dataset_id, il_err)

         IF (ig_hdfverb.GE.2) &
            & WRITE(*,'(1X,A,A,T67,A,A,A,A,A)') &
            & '> [I/Oa] Reading slice of dataset ',TRIM(cd_dsname), &
            & '(Min= ',TRIM(cda_dset(1)),'  Max= ',TRIM(cda_dset(SIZE(cda_dset))),')'
      END IF

   END SUBROUTINE READSLICE_H5DSET_C1
   !---------------------------------------------------------------------------
   SUBROUTINE READSLICE_H5DSET_I1(h5state, cd_dsname, ida_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read HDF5 dataset slice of integer 1d-array data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      INTEGER(KIND=4), DIMENSION(:), INTENT(OUT) :: ida_dset   ! With INTEGER of kind=8 (int64),
      ! there is no matching specific SUBROUTINE for this generic SUBROUTINE call.   [H5DREAD_F]
      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
      ! Read dataset
      CALL H5Dread_f(il_dataset_id, il_datatype_id, ida_dset, h5state%totnobs, &
         & il_err, &
         & h5state%memspace_id, h5state%dataspace_id)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,I12,A,I12,A)') &
         & '> [I/Oa] Reading slice of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(ida_dset),'  Max= ',MAXVAL(ida_dset),')'

   END SUBROUTINE READSLICE_H5DSET_I1
   !---------------------------------------------------------------------------
   SUBROUTINE READSLICE_H5DSET_R1(h5state, cd_dsname, rda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read HDF5 dataset slice of real 1d-array data type
      !!
      !! ----------------------------------------------------------------------
      !!
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: rda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)
      ! Read dataset
      CALL H5Dread_f(il_dataset_id, il_datatype_id, rda_dset, h5state%totnobs, &
         & il_err, &
         & h5state%memspace_id, h5state%dataspace_id)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Reading slice of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

   END SUBROUTINE READSLICE_H5DSET_R1
   !---------------------------------------------------------------------------
   SUBROUTINE READSLICE_H5DSET_R2(h5state, cd_dsname, rda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read HDF5 dataset slice of real 2d-array data type
      !!
      !! ----------------------------------------------------------------------
      !!
      ! USE SYSTEM_UTILS_MOD, ONLY : DAIMON_Abort
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=4), DIMENSION(:,:), INTENT(OUT) :: rda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_rank
      INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE :: ila_shape
      INTEGER :: il_ni
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)

      ! Get the shape of the data to check dimensions
      CALL H5Tget_array_ndims_f(il_datatype_id, il_rank, il_err)
      IF (ig_hdfverb.GE.2) WRITE(*,*) 'Rank of the datatype ', il_rank
      ALLOCATE(ila_shape(il_rank))
      CALL H5Tget_array_dims_f(il_datatype_id, ila_shape, il_err)
      IF (ig_hdfverb.GE.2) WRITE(*,*) 'Dimensions of the datatype ', ila_shape
      ! Check array first dimension versus data dimension
      il_ni = SIZE(rda_dset,DIM=1)
      IF (il_ni.NE.ila_shape(1)) THEN
         WRITE(*,'(/,1X,A,I3,3A,I3)') '/!\ The array first dimension (',il_ni, &
            & ' from namelist) is NOT equal to the data dimension in ',TRIM(cd_dsname), &
            & ' dataset (',ila_shape(1),' in HDAT)'
         WRITE(*,*) '/!\ Check your namelist and/or your HDATxxx.h5 file'
         ! CALL DAIMON_Abort
      END IF
      DEALLOCATE(ila_shape)

      ! Read dataset
      CALL H5Dread_f(il_dataset_id, il_datatype_id, rda_dset, h5state%totnobs, &
         & il_err, &
         & h5state%memspace_id, h5state%dataspace_id)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Reading slice of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

   END SUBROUTINE READSLICE_H5DSET_R2
   !---------------------------------------------------------------------------
   SUBROUTINE READSLICE_H5DSET_R3(h5state, cd_dsname, rda_dset)
      !! category: observation operator
      !! author:   CERFACS and CNRM (G. Jonville)
      !!
      !! ----------------------------------------------------------------------
      !!
      !! Read HDF5 dataset slice of real 3d-array data type
      !!
      !! ----------------------------------------------------------------------
      !!
      ! USE SYSTEM_UTILS_MOD, ONLY : DAIMON_Abort
      TYPE(h5state_t), INTENT(IN) :: h5state
      CHARACTER(LEN=*), INTENT(IN) :: cd_dsname
      REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT) :: rda_dset

      ! Identifiers
      INTEGER(hid_t) :: il_dataset_id
      INTEGER(hid_t) :: il_datatype_id
      ! Dimensions
      INTEGER :: il_rank
      INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE :: ila_shape
      INTEGER :: il_ni, il_nj
      ! Inquiries
      INTEGER :: il_err

      ! Open dataset
      CALL H5Dopen_f(h5state%instr_id, cd_dsname, il_dataset_id, il_err)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' opening error ', il_err
      ! Get dataset type
      CALL H5Dget_type_f(il_dataset_id, il_datatype_id, il_err)

      ! Get the shape of the data to check dimensions
      CALL H5Tget_array_ndims_f(il_datatype_id, il_rank, il_err)
      IF (ig_hdfverb.GE.2) WRITE(*,*) 'Rank of the datatype ', il_rank
      ALLOCATE(ila_shape(il_rank))
      CALL H5Tget_array_dims_f(il_datatype_id, ila_shape, il_err)
      IF (ig_hdfverb.GE.2) WRITE(*,*) 'Dimensions of the datatype ', ila_shape
      ! Check array first and second dimension versus data dimensions
      il_ni = SIZE(rda_dset,DIM=1)
      il_nj = SIZE(rda_dset,DIM=2)
      IF (il_ni.NE.ila_shape(1) .OR. il_nj.NE.ila_shape(2)) THEN
         WRITE(*,'(/,1X,2(A,I3),3A,I3,A,I3)') '/!\ The array first or second dimension (',il_ni,' x ',il_nj, &
            & ' from namelist) is NOT equal to the data dimension in ',TRIM(cd_dsname), &
            & ' dataset (',ila_shape(1),' x ',ila_shape(2),' in HDAT)'
         WRITE(*,*) '/!\ Check your namelist and/or your HDATxxx.h5 file'
         ! CALL DAIMON_Abort
      END IF
      DEALLOCATE(ila_shape)

      ! Read dataset
      CALL H5Dread_f(il_dataset_id, il_datatype_id, rda_dset, h5state%totnobs, &
         & il_err, &
         & h5state%memspace_id, h5state%dataspace_id)
      IF (il_err /= 0) WRITE(*,'(1X,3A,I3)') 'Dataset ',TRIM(cd_dsname),' read error ',il_err
      ! Close dataset
      CALL H5Dclose_f(il_dataset_id, il_err)

      IF (ig_hdfverb.GE.2) &
         & WRITE(*,'(1X,A,A,T67,A,E12.4,A,E12.4,A)') &
         & '> [I/Oa] Reading slice of dataset ',TRIM(cd_dsname), &
         & '(Min= ',MINVAL(rda_dset),'  Max= ',MAXVAL(rda_dset),')'

   END SUBROUTINE READSLICE_H5DSET_R3
   !---------------------------------------------------------------------------
END MODULE H5_READ_MOD
