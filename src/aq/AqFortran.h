/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_AQFORTRAN_H_
#define AQ_AQFORTRAN_H_

#define AQ_STRLEN 256

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace aq {
  class LocationsAQ;
  class ObsSpaceAQ;

// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Gom key type
typedef int F90gom;
// Fields key type
typedef int F90flds;
// Observation vector key type
typedef int F90ovec;
// Observation data base type
typedef int F90odb;
// Obs op matrix
typedef int F90hmat;

/// Interface to Fortran AQ model
/*!
 * The core of the AQ model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void aq_fields_create_f90(F90flds &, const F90geom &, const oops::Variables &);
  void aq_fields_create_from_other_f90(F90flds &, const F90flds &);
  void aq_fields_delete_f90(F90flds &);
  void aq_fields_zero_f90(const F90flds &);
  void aq_fields_ones_f90(const F90flds &);
  void aq_fields_dirac_f90(const F90flds &, const eckit::Configuration &);
  void aq_fields_random_f90(const F90flds &, const oops::Variables &);
  void aq_fields_copy_f90(const F90flds &, const F90flds &);
  void aq_fields_self_add_f90(const F90flds &, const F90flds &);
  void aq_fields_self_sub_f90(const F90flds &, const F90flds &);
  void aq_fields_self_mul_f90(const F90flds &, const double &);
  void aq_fields_axpy_f90(const F90flds &, const double &, const F90flds &);
  void aq_fields_self_schur_f90(const F90flds &, const F90flds &);
  void aq_fields_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void aq_fields_add_incr_f90(const F90flds &, const F90flds &);
  void aq_fields_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);
  // void aq_fields_change_resol_f90(const F90flds &, const F90flds &);
  void aq_fields_info_f90(const F90flds &, const eckit::Configuration &);
  void aq_fields_read_file_f90(const F90flds &, const eckit::Configuration &,
                               util::DateTime &);
  void aq_fields_write_file_f90(const F90flds &, const eckit::Configuration &,
                                const util::DateTime &);
  void aq_fields_analytic_init_f90(const F90flds &, const eckit::Configuration &);
  void aq_fields_gpnorm_f90(const F90flds &, int[], double[], double[], double[]);
  void aq_fields_rms_f90(const F90flds &, double &);
  void aq_fields_sizes_f90(const F90flds &, int &, int &, int &);
  void aq_fields_lbc_f90(const F90flds &, int &);
  void aq_fields_set_atlas_f90(const F90flds &, const oops::Variables &,
                               atlas::field::FieldSetImpl *);
  void aq_fields_to_atlas_f90(const F90flds &, const oops::Variables &,
                              atlas::field::FieldSetImpl *);
// AQ  void aq_fields_from_atlas_f90(const F90flds &, const oops::Variables &,
// AQ                              atlas::field::FieldSetImpl *);
  void aq_fields_getpoint_f90(const F90flds&, const F90iter&, const int &, double &);
  void aq_fields_setpoint_f90(const F90flds&, const F90iter&, const int &, const double &);
  void aq_fields_serialsize_f90(const F90flds &, int &);
  void aq_fields_serialize_f90(const F90flds &, const int &, double[]);
  void aq_fields_deserialize_f90(const F90flds &, const int &, const double[], int &);

// -----------------------------------------------------------------------------
//  GetValues
// -----------------------------------------------------------------------------
  void aq_getvalues_setup_f90(const F90hmat &);
  void aq_getvalues_delete_f90(const F90hmat &);
  void aq_getvalues_interp_f90(const LocationsAQ &, const F90flds &,
                               const util::DateTime &,
                               const util::DateTime &, const F90gom &);
  void aq_getvalues_build_f90(const LocationsAQ &, const F90flds &,
                              const util::DateTime &,
                              const util::DateTime &, const F90hmat &);
  void aq_getvalues_interp_tl_f90(const LocationsAQ &, const F90flds &,
                                  const util::DateTime &,
                                  const util::DateTime &, const F90hmat &, const F90gom &);
  void aq_getvalues_interp_ad_f90(const LocationsAQ &, const F90flds &,
                                  const util::DateTime &,
                                  const util::DateTime &, const F90hmat &, const F90gom &);

// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void aq_geom_setup_f90(F90geom &, const eckit::Configuration &, const eckit::mpi::Comm *,
                         atlas::grid::GridImpl *,
                         atlas::functionspace::FunctionSpaceImpl *);
  void aq_geom_fill_atlas_fieldset_f90(const F90geom &, atlas::field::FieldSetImpl *);
  void aq_geom_clone_f90(F90geom &, const F90geom &);
  void aq_geom_info_f90(const F90geom &, int &, int &, int &,
                        double &, double &, int &, char *, char *, char *);
  void aq_geom_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Geometry iterator
// -----------------------------------------------------------------------------
  void aq_geom_iter_setup_f90(F90iter &, const F90geom &, const int &);
  void aq_geom_iter_clone_f90(F90iter &, const F90iter &);
  void aq_geom_iter_delete_f90(F90iter &);
  void aq_geom_iter_equals_f90(const F90iter &, const F90iter&, int &);
  void aq_geom_iter_current_f90(const F90iter &, double &, double &);
  void aq_geom_iter_next_f90(const F90iter &);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void aq_gom_setup_f90(F90gom &, const LocationsAQ &, const oops::Variables &);
  void aq_gom_create_f90(F90gom &, const oops::Variables &);
  void aq_gom_delete_f90(F90gom &);
  void aq_gom_copy_f90(const F90gom &, const F90gom &);
  void aq_gom_zero_f90(const F90gom &);
  void aq_gom_abs_f90(const F90gom &);
  void aq_gom_random_f90(const F90gom &);
  void aq_gom_mult_f90(const F90gom &, const double &);
  void aq_gom_add_f90(const F90gom &, const F90gom &);
  void aq_gom_diff_f90(const F90gom &, const F90gom &);
  void aq_gom_schurmult_f90(const F90gom &, const F90gom &);
  void aq_gom_divide_f90(const F90gom &, const F90gom &);
  void aq_gom_rms_f90(const F90gom &, double &);
  void aq_gom_dotprod_f90(const F90gom &, const F90gom &, double &);
  void aq_gom_stats_f90(const F90gom &, int &, double &, double &, double &);
  void aq_gom_maxloc_f90(const F90gom &, double &, int &, const oops::Variables &);
  void aq_gom_read_file_f90(const F90gom &, const eckit::Configuration &);
  void aq_gom_write_file_f90(const F90gom &, const eckit::Configuration &);
  void aq_gom_analytic_init_f90(const F90gom &, const LocationsAQ &,
                                const eckit::Configuration &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void aq_model_setup_f90(F90model &, const eckit::Configuration &);
  void aq_model_delete_f90(F90model &);
  void aq_model_propagate_f90(const F90model &, const F90flds &);
  void aq_model_propagate_tl_f90(const F90model &, const F90flds &, const F90flds &);
  void aq_model_propagate_ad_f90(const F90model &, const F90flds &, const F90flds &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void aq_obsdb_setup_f90(F90odb &, const eckit::Configuration &,
                          const util::DateTime &, const util::DateTime &,
                          const eckit::mpi::Comm *);
  void aq_obsdb_delete_f90(F90odb &);
  void aq_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void aq_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void aq_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              atlas::field::FieldSetImpl *, std::vector<util::DateTime> &);
  void aq_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration &, const util::DateTime &,
                             const util::Duration &, const int &, int &);
  void aq_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);

// -----------------------------------------------------------------------------
//  Observation vector
// -----------------------------------------------------------------------------
  void aq_obsvec_setup_f90(F90ovec &, const int &, const int &);
  void aq_obsvec_clone_f90(F90ovec &, const F90ovec &);
  void aq_obsvec_delete_f90(F90ovec &);
  void aq_obsvec_copy_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_zero_f90(const F90ovec &);
  void aq_obsvec_settomissing_ith_f90(const F90ovec &, const int &);
  void aq_obsvec_ones_f90(const F90ovec &);
  /// set ObsVector (with key \p obsvector_key) values to missing values where
  /// mask ObsVector (with key \p mask_key) values are set to 1
  void aq_obsvec_mask_f90(const F90ovec & obsvector_key, const F90ovec & mask_key);
  /// set ObsVector (with key \p obsvector_key) values to missing values where
  /// mask ObsVector (with key \p mask_key) values are set to missing value
  void aq_obsvec_mask_with_missing_f90(const F90ovec & obsvector_key,
                                       const F90ovec & mask_key);
  void aq_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void aq_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void aq_obsvec_invert_f90(const F90ovec &);
  void aq_obsvec_random_f90(const ObsSpaceAQ &, const F90ovec &);
  void aq_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void aq_obsvec_stats_f90(const F90ovec &, double &, double &, double &);
  void aq_obsvec_nobs_f90(const F90ovec &, int &);
  void aq_obsvec_size_f90(const F90ovec &, int &);
  /// fill \p data (size \p nobs) with all non-masked out (non-missing) values
  void aq_obsvec_get_withmask_f90(const F90ovec &, const F90ovec & mask_key,
                                  double * data, const int & nobs);
  void aq_obsvec_nobs_withmask_f90(const F90ovec &, const F90ovec & mask_key, int &);


// -----------------------------------------------------------------------------
//  Streamfunction observations
// -----------------------------------------------------------------------------
  void aq_stream_equiv_f90(const F90gom &, const F90ovec &, const double &);
  void aq_stream_equiv_tl_f90(const F90gom &, const F90ovec &, const double &);
  void aq_stream_equiv_ad_f90(const F90gom &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind observations
// -----------------------------------------------------------------------------
  void aq_wind_equiv_f90(const F90gom &, const F90ovec &, const double &);
  void aq_wind_equiv_tl_f90(const F90gom &, const F90ovec &, const double &);
  void aq_wind_equiv_ad_f90(const F90gom &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind speed observations
// -----------------------------------------------------------------------------
  void aq_wspeed_equiv_f90(const F90gom &, const F90ovec &, const double &);
  void aq_wspeed_equiv_tl_f90(const F90gom &, const F90ovec &, const F90gom &, const double &);
  void aq_wspeed_equiv_ad_f90(const F90gom &, const F90ovec &, const F90gom &, double &);
  void aq_wspeed_gettraj_f90(const int &, const oops::Variables &, const F90gom &);
  void aq_wspeed_settraj_f90(const F90gom &, const F90gom &);

}

}  // namespace aq
#endif  // AQ_AQFORTRAN_H_
