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

#ifndef AQ_AQ_OBSVEC_INTERFACE_H_
#define AQ_AQ_OBSVEC_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

// Forward declarations
namespace aq {
  class ObsSpace;

extern "C" {
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
  void aq_obsvec_threshold_check_f90(const F90ovec & obsvector_key, const F90ovec & othervector_key,
   const F90ovec & mask_key, const eckit::Configuration &);
  void aq_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void aq_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void aq_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void aq_obsvec_invert_f90(const F90ovec &);
  void aq_obsvec_random_f90(const ObsSpace &, const F90ovec &);
  void aq_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void aq_obsvec_stats_f90(const F90ovec &, double &, double &, double &);
  void aq_obsvec_nobs_f90(const F90ovec &, int &);
  void aq_obsvec_size_f90(const F90ovec &, int &);
  /// fill \p data (size \p nobs) with all non-masked out (non-missing) values
  void aq_obsvec_get_withmask_f90(const F90ovec &, const F90ovec & mask_key,
                                  double * data, const int & nobs);
  void aq_obsvec_nobs_withmask_f90(const F90ovec &, const F90ovec & mask_key, int &);
}
}  // namespace aq
#endif  // AQ_AQ_OBSVEC_INTERFACE_H_
