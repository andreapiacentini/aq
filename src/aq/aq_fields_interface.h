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

#ifndef AQ_AQ_FIELDS_INTERFACE_H_
#define AQ_AQ_FIELDS_INTERFACE_H_

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
}

namespace aq {
extern "C" {
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
}
}  // namespace aq
#endif  // AQ_AQ_FIELDS_INTERFACE_H_
